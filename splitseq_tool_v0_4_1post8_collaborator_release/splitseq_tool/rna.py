from __future__ import annotations

import re
import gzip
import shlex
from pathlib import Path

from .fastq import iter_fastq, extract_splitseq_cell_umi_from_header
from .utils import ensure_dir, run_cmd, cleanup_intermediates, require_output_files


def _base_qname(header: str) -> str:
    """Return a whitespace-free base read name without prior SPLiT-seq suffixes."""
    token = header.strip().split()[0].lstrip("@")
    # Remove a previous rebuilt-tool suffix if one is already present.
    token = re.sub(r"_[ACGTNacgtn]{24}_[ACGTNacgtn]+$", "", token)
    token = re.sub(r":CB_[ACGTNacgtn]{24}_UMI_[ACGTNacgtn]+$", "", token)
    return token


def prepare_umi_tools_fastq(input_fastq: Path, output_fastq: Path) -> dict:
    """Rewrite demultiplexed FASTQ headers into UMI-tools-compatible form.

    UMI-tools' standard single-cell extraction convention appends the cell
    barcode and then UMI to the read name, separated by underscores.  We use:

        @original_read_id_<24nt_cell_barcode>_<UMI>

    Critically, the SPLiT-seq metadata is placed in the first whitespace-
    delimited FASTQ header token so STAR preserves it in BAM QNAMEs.
    """
    stats = {"records": 0, "rewritten": 0, "missing_cell_or_umi": 0}
    with open(output_fastq, "w") as out:
        for rec in iter_fastq(input_fastq):
            stats["records"] += 1
            cell, umi = extract_splitseq_cell_umi_from_header(rec.name)
            if not cell or not umi:
                stats["missing_cell_or_umi"] += 1
                raise ValueError(
                    f"Could not recover cell barcode and UMI from demultiplexed FASTQ header: {rec.name}"
                )
            name = f"@{_base_qname(rec.name)}_{cell}_{umi}"
            out.write(f"{name}\n{rec.seq}\n+\n{rec.qual}\n")
            stats["rewritten"] += 1
    return stats


def run_star_pipeline(args) -> None:
    """Run STAR -> featureCounts -> UMI-tools after SPLiT-seq demux.

    The input MergedCells_passing.fastq may use either the legacy header form
    or the rebuilt-tool CB_/UMI_ form.  Before STAR, headers are normalized to
    UMI-tools' standard ``readid_CELL_UMI`` convention so STAR preserves the
    cell barcode and UMI in BAM QNAMEs.
    """
    out_dir = ensure_dir(args.output_dir)
    work = ensure_dir(out_dir / args.sample_id / "rna_star")
    log = work / "rna_star_pipeline.log"
    merged = Path(args.merged_fastq)
    if not merged.exists():
        raise FileNotFoundError(merged)

    star_input = work / "MergedCells_passing.umi_tools_headers.fastq"
    prep = prepare_umi_tools_fastq(merged, star_input)
    with open(log, "a") as lh:
        lh.write(
            f"Prepared STAR/UMI-tools FASTQ: records={prep['records']} "
            f"rewritten={prep['rewritten']} missing={prep['missing_cell_or_umi']}\n"
        )

    run_cmd([
        "STAR", "--runThreadN", str(args.threads), "--readFilesIn", str(star_input),
        "--outFilterMismatchNoverLmax", str(args.star_mismatch_nover_lmax),
        "--genomeDir", args.star_index,
        "--alignIntronMax", str(args.align_intron_max),
        "--outSAMtype", "BAM", "SortedByCoordinate",
    ], cwd=work, log_path=log)

    bam = work / "Aligned.sortedByCoord.out.bam"
    derived_intermediates = run_featurecounts_and_umi(args, bam, work, log)

    # STAR alignment BAM and the rewritten FASTQ are large staging products for
    # the final gene-by-cell UMI matrix. Delete them only after the RNA outputs
    # have been successfully produced and validated.
    cleanup_intermediates(
        [star_input, bam] + derived_intermediates,
        work / "cleanup_summary.tsv",
        keep_all=bool(getattr(args, "keep_all_intermediates", False)),
        reason="RNA intermediate after successful STAR/featureCounts/UMI-tools completion",
    )


def resume_featurecounts(args) -> None:
    work = ensure_dir(args.output_dir)
    log = work / "resume_featurecounts.log"
    bam = Path(args.aligned_bam)
    if not bam.exists():
        raise FileNotFoundError(bam)
    derived_intermediates = run_featurecounts_and_umi(args, bam, work, log)

    # Never delete --aligned-bam here: it is user-supplied input to this resume
    # command. Only delete intermediates created by this command itself.
    cleanup_intermediates(
        derived_intermediates,
        work / "cleanup_summary.tsv",
        keep_all=bool(getattr(args, "keep_all_intermediates", False)),
        reason="featureCounts/UMI-tools intermediate after successful resume completion",
    )


def run_featurecounts_and_umi(args, bam: Path, work: Path, log: Path) -> list[Path]:
    run_cmd([
        "featureCounts", "-F", "SAF", "-a", args.saf_file,
        "-o", "gene_assigned", "-R", "BAM", str(bam),
        "-T", str(args.threads), "-M", "-O",
    ], cwd=work, log_path=log)

    fc_bam = work / (bam.name + ".featureCounts.bam")
    assigned_sorted = work / "assigned_sorted.bam"
    run_cmd(["samtools", "sort", str(fc_bam), "-o", str(assigned_sorted)], cwd=work, log_path=log)
    run_cmd(["samtools", "index", str(assigned_sorted)], cwd=work, log_path=log)

    assigned_gene = work / "assigned_gene_only.bam"
    assigned_sorted_q = shlex.quote(str(assigned_sorted))
    assigned_gene_q = shlex.quote(str(assigned_gene))
    shell_cmd = (
        f"samtools view -h {assigned_sorted_q} | "
        "awk 'BEGIN{OFS=\"\t\"} /^@/ || /XS:Z:Assigned/' | "
        f"samtools view -b -o {assigned_gene_q} -"
    )
    run_cmd(["bash", "-lc", shell_cmd], cwd=work, log_path=log)
    run_cmd(["samtools", "index", str(assigned_gene)], cwd=work, log_path=log)

    # Sanity-check that BAM QNAMEs retain the cell barcode + UMI suffix before
    # invoking UMI-tools. This provides a useful failure message instead of a
    # cryptic downstream parsing exception.
    assigned_gene_check_q = shlex.quote(str(assigned_gene))
    qname_check = (
        f"samtools view {assigned_gene_check_q} | head -n 100 | "
        "awk '$1 ~ /_[ACGTNacgtn]{24}_[ACGTNacgtn]+$/ {n++} END {if (n==0) exit 1}'"
    )
    try:
        run_cmd(["bash", "-lc", qname_check], cwd=work, log_path=log)
    except Exception as exc:
        raise RuntimeError(
            "Assigned BAM read names do not contain the expected _CELL_UMI suffix. "
            "Cell/UMI metadata was likely lost before or during alignment."
        ) from exc

    counts = work / "counts.tsv.gz"
    run_cmd([
        "umi_tools", "count", "--wide-format-cell-counts", "--per-gene",
        "--gene-tag=XT", "--assigned-status-tag=XS", "--per-cell",
        "--extract-umi-method=read_id", "--umi-separator=_",
        "-I", str(assigned_gene), "-S", str(counts),
    ], cwd=work, log_path=log)

    # Generate cellIDs_keep.txt directly from the wide count-table header.
    # This avoids shell quoting problems and does not depend on an external zcat.
    cell_keep = work / "cellIDs_keep.txt"
    opener = gzip.open if str(counts).endswith(".gz") else open
    with opener(counts, "rt") as inp:
        header = inp.readline().rstrip("\r\n").split("\t")
    if len(header) < 2:
        raise RuntimeError(f"UMI-tools count table has no cell columns: {counts}")
    with open(cell_keep, "w") as out:
        for cell in header[1:]:
            if cell:
                out.write(cell + "\n")

    require_output_files([counts], allow_empty=False)
    require_output_files([cell_keep], allow_empty=True)

    return [
        fc_bam,
        assigned_sorted,
        Path(str(assigned_sorted) + ".bai"),
        assigned_gene,
        Path(str(assigned_gene) + ".bai"),
    ]
