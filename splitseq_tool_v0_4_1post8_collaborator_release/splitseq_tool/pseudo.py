from __future__ import annotations

import gzip
import shlex
from collections import Counter
from pathlib import Path

from .fastq import iter_fastq, extract_splitseq_cell_umi_from_header
from .utils import ensure_dir, open_text, run_cmd, write_kv, cleanup_intermediates, require_output_files


def _qual(n: int) -> str:
    return "I" * n


def prepare_kb_fastqs_from_merged(args) -> None:
    """Create kb/kallisto-bustools compatible paired FASTQs from MergedCells_passing.fastq.

    Read 1 of the synthetic pair contains corrected cell barcode + UMI.
    Read 2 contains the biological cDNA sequence from the demultiplexed FASTQ.

    This lets kb/kallisto parse an already-corrected SPLiT-seq cell barcode/UMI layout
    without relying on its raw-read barcode correction behavior.
    """
    out_dir = ensure_dir(args.output_dir)
    merged = Path(args.merged_fastq)
    if not merged.exists():
        raise FileNotFoundError(merged)

    prefix = args.sample_id or merged.stem
    barcode_fastq = out_dir / f"{prefix}.kb_barcode_umi_R1.fastq"
    cdna_fastq = out_dir / f"{prefix}.kb_cdna_R2.fastq"
    whitelist = out_dir / f"{prefix}.kb_cell_whitelist.txt"
    summary = Counter()
    cell_counts = Counter()
    umi_missing = 0

    requested_cells = None
    if getattr(args, "cell_whitelist", None):
        requested_cells = set()
        with open_text(args.cell_whitelist, "rt") as inp:
            for line_number, line in enumerate(inp, start=1):
                cell = line.strip().split()[0] if line.strip() else ""
                if not cell:
                    continue
                if len(cell) != int(args.cell_barcode_len) or not all(c in "ACGTNacgtn" for c in cell):
                    raise ValueError(
                        f"Invalid cell barcode in --cell-whitelist at line {line_number}: {cell!r}; "
                        f"expected {int(args.cell_barcode_len)} nucleotide characters"
                    )
                requested_cells.add(cell.upper())
        if not requested_cells:
            raise ValueError("--cell-whitelist did not contain any valid cell barcodes")

    # Preserve gzip output if requested.
    def _writer(path: Path):
        if str(path).endswith(".gz"):
            return gzip.open(path, "wt")
        return open(path, "w")

    if getattr(args, "gzip_output", False):
        barcode_fastq = Path(str(barcode_fastq) + ".gz")
        cdna_fastq = Path(str(cdna_fastq) + ".gz")

    with _writer(barcode_fastq) as bout, _writer(cdna_fastq) as rout:
        for rec in iter_fastq(merged):
            summary["records_scanned"] += 1
            cell, umi = extract_splitseq_cell_umi_from_header(rec.name)
            if not cell:
                summary["missing_cell_id"] += 1
                continue
            if not umi:
                umi_missing += 1
                raise ValueError(
                    f"Missing UMI in demultiplexed FASTQ header: {rec.name}. "
                    "Pseudoalignment requires a real UMI for every record."
                )
            cell = cell.upper()
            umi = umi.upper()
            if len(cell) != int(args.cell_barcode_len):
                raise ValueError(
                    f"Unexpected cell-barcode length in demultiplexed FASTQ header: {cell!r}; "
                    f"expected {int(args.cell_barcode_len)}"
                )
            if len(umi) != int(args.umi_len):
                raise ValueError(
                    f"Unexpected UMI length in demultiplexed FASTQ header: {umi!r}; "
                    f"expected {int(args.umi_len)}. Set --umi-len to the demultiplexed UMI length."
                )
            if requested_cells is not None and cell not in requested_cells:
                summary["records_not_in_cell_whitelist"] += 1
                continue

            cell_counts[cell] += 1
            bc_umi = cell + umi
            name = rec.name.split()[0]
            bout.write(f"{name}\n{bc_umi}\n+\n{_qual(len(bc_umi))}\n")
            rout.write(f"{name}\n{rec.seq}\n+\n{rec.qual}\n")
            summary["records_written"] += 1

    if summary["records_written"] == 0:
        raise ValueError(
            "No records were written to pseudoaligner FASTQs. Check the demultiplexed input, "
            "cell-barcode/UMI lengths, and any supplied --cell-whitelist."
        )

    # If an explicit cell whitelist was supplied, reads were filtered to it above.
    # Write only cells that are both requested and actually observed, preventing kb
    # barcode correction from reassigning excluded cells to nearby on-list barcodes.
    with open(whitelist, "w") as out:
        for cell in sorted(cell_counts):
            out.write(cell + "\n")

    summary["unique_cells_written"] = len(cell_counts)
    summary["records_missing_umi"] = umi_missing
    summary["cell_barcode_len"] = int(args.cell_barcode_len)
    summary["umi_len"] = int(args.umi_len)
    write_kv(out_dir / "kb_prepare_summary.tsv", list(summary.items()))

    with open(out_dir / "kb_inputs.tsv", "w") as out:
        out.write("key\tpath\n")
        out.write(f"barcode_umi_fastq\t{barcode_fastq}\n")
        out.write(f"cdna_fastq\t{cdna_fastq}\n")
        out.write(f"cell_whitelist\t{whitelist}\n")
        cb_len = int(args.cell_barcode_len)
        umi_len = int(args.umi_len)
        tech = f"0,0,{cb_len}:0,{cb_len},{cb_len + umi_len}:1,0,0"
        out.write(f"technology_string\t{tech}\n")

    print("Prepared kb-compatible FASTQs:")
    print(f"  barcode/UMI FASTQ: {barcode_fastq}")
    print(f"  cDNA FASTQ:        {cdna_fastq}")
    print(f"  whitelist:         {whitelist}")
    print(f"  technology string: {tech}")


def run_kb_pipeline(args) -> None:
    """Run kb-python/kallisto-bustools from demultiplexed SPLiT-seq FASTQ.

    This command uses the safe two-step route:
      1. Convert MergedCells_passing.fastq into synthetic kb-compatible FASTQs.
      2. Run kb count using the synthetic barcode/UMI read and the cDNA read.
    """
    out_dir = ensure_dir(args.output_dir)
    prepare_dir = ensure_dir(out_dir / "kb_prepared_fastqs")

    # Reuse preparation code with a shallow args-like object.
    class PrepArgs:
        pass

    p = PrepArgs()
    p.output_dir = str(prepare_dir)
    p.merged_fastq = args.merged_fastq
    p.sample_id = args.sample_id
    p.cell_whitelist = getattr(args, "cell_whitelist", None)
    p.cell_barcode_len = getattr(args, "cell_barcode_len", 24)
    p.umi_len = getattr(args, "umi_len", 10)
    p.gzip_output = getattr(args, "gzip_fastqs", False)
    prepare_kb_fastqs_from_merged(p)

    inputs = {}
    with open(prepare_dir / "kb_inputs.tsv") as f:
        next(f)
        for line in f:
            k, v = line.rstrip("\n").split("\t", 1)
            inputs[k] = v

    kb_out = ensure_dir(out_dir / "kb_count")
    log = kb_out / "kb_count.log"
    tech = getattr(args, "technology_string", None) or inputs["technology_string"]

    cmd = [
        "kb", "count",
        "-i", args.kallisto_index,
        "-g", args.t2g_file,
        "-x", tech,
        "-w", inputs["cell_whitelist"],
        "-o", str(kb_out),
        "-t", str(args.threads),
    ]

    workflow = getattr(args, "workflow", None)
    if workflow:
        cmd += ["--workflow", workflow]
    strand = getattr(args, "strand", None)
    if strand:
        cmd += ["--strand", strand]
    if getattr(args, "h5ad", False):
        cmd += ["--h5ad"]
    if getattr(args, "keep_tmp", False) or getattr(args, "keep_all_intermediates", False):
        cmd += ["--keep-tmp"]
    extra = getattr(args, "kb_extra_args", None)
    if extra:
        cmd += shlex.split(str(extra))

    # kb custom technology uses first FASTQ as barcode/UMI read and second as cDNA read.
    cmd += [inputs["barcode_umi_fastq"], inputs["cdna_fastq"]]

    run_cmd(cmd, cwd=kb_out, log_path=log)

    with open(out_dir / "kb_pipeline_outputs.tsv", "w") as out:
        out.write("key\tpath_or_value\n")
        out.write(f"prepared_fastq_dir\t{prepare_dir}\n")
        out.write(f"kb_count_dir\t{kb_out}\n")
        out.write(f"technology_string\t{tech}\n")
        out.write(f"matrix_dir\t{kb_out}\n")

    outputs_manifest = out_dir / "kb_pipeline_outputs.tsv"
    require_output_files([outputs_manifest], allow_empty=False)
    meaningful_kb_outputs = [p for p in kb_out.iterdir() if p.name != log.name]
    if not meaningful_kb_outputs:
        raise RuntimeError(
            "Refusing intermediate-file cleanup because kb count completed without producing any "
            f"non-log outputs in {kb_out}"
        )

    # The two synthetic FASTQs can be very large; they are only adapters for
    # kb/kallisto-bustools. Keep their small whitelist/summary metadata for
    # reproducibility, but remove the synthetic reads after a successful run.
    cleanup_intermediates(
        [inputs["barcode_umi_fastq"], inputs["cdna_fastq"]],
        out_dir / "cleanup_summary.tsv",
        keep_all=bool(getattr(args, "keep_all_intermediates", False)),
        reason="synthetic kb FASTQs no longer needed after successful kb count",
    )


def write_raw_splitseq_kb_design(args) -> None:
    """Write helper files for the direct raw FASTQ kallisto-bustools route.

    This is provided as an advanced/experimental option. It creates a concatenated
    24-nt whitelist in R1+R2+R3 order and reports the raw-read technology string.
    """
    from itertools import product
    from .barcodes import load_barcode_list

    out_dir = ensure_dir(args.output_dir)
    r1 = load_barcode_list(args.round1, 8, "auto")
    r2 = load_barcode_list(args.round2, 8, "auto")
    r3 = load_barcode_list(args.round3, 8, "auto")
    whitelist = out_dir / "splitseq_R1R2R3_onlist.txt"
    with open(whitelist, "w") as out:
        for a, b, c in product(r1, r2, r3):
            out.write(a + b + c + "\n")

    # Specify barcode in R1+R2+R3 order even though the positions on the barcode read
    # are typically R3, R2, R1. Coordinates are user-configurable and zero-based.
    r1s = int(args.r1_start); r2s = int(args.r2_start); r3s = int(args.r3_start)
    us = int(args.umi_start); ul = int(args.umi_len)
    tech = f"1,{r1s},{r1s+8},1,{r2s},{r2s+8},1,{r3s},{r3s+8}:1,{us},{us+ul}:0,0,0"
    with open(out_dir / "raw_kb_design.tsv", "w") as out:
        out.write("key\tvalue\n")
        out.write(f"technology_string\t{tech}\n")
        out.write(f"whitelist\t{whitelist}\n")
        out.write("note\tExperimental: kb/bustools correction on a 24-nt concatenated barcode is not identical to per-round SPLiT-seq correction. Prefer rna-kb after splitseq-tool demux for best equivalence.\n")
    print(f"Wrote raw kb design files to {out_dir}")
