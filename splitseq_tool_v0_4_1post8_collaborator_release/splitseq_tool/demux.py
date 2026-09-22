from __future__ import annotations

from collections import Counter, defaultdict
from itertools import zip_longest
from pathlib import Path

from .barcodes import BarcodeCorrector, load_barcode_list
from .fastq import iter_fastq, extract_splitseq_cell_umi_from_header
from .utils import ensure_dir, write_kv, cleanup_intermediates, require_output_files
from .quality import QualityAccumulator, summary_rows, write_quality_qc


def _read_id(name: str) -> str:
    token = name.split()[0].lstrip("@")
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def demux_barcoded_fastq(args) -> None:
    """Streaming SPLiT-seq demultiplexing with per-round correction.

    Cell thresholding defaults to unique UMIs per cell to match the legacy ``-m``
    interpretation. ``--cell-threshold-mode reads`` is available when desired.
    """
    out_dir = ensure_dir(args.output_dir)
    r1_list = load_barcode_list(args.round1, 8, "auto")
    r2_list = load_barcode_list(args.round2, 8, "auto")
    r3_list = load_barcode_list(args.round3, 8, "auto")
    if not r1_list or not r2_list or not r3_list:
        raise ValueError(f"Empty barcode whitelist detected: R1={len(r1_list)} R2={len(r2_list)} R3={len(r3_list)}")
    c1 = BarcodeCorrector(r1_list, args.barcode_mismatches)
    c2 = BarcodeCorrector(r2_list, args.barcode_mismatches)
    c3 = BarcodeCorrector(r3_list, args.barcode_mismatches)
    cell_reads = Counter()
    cell_umis = defaultdict(set)
    stats = Counter()
    quality_threshold = int(getattr(args, "quality_threshold", 30))
    phred_offset = int(getattr(args, "phred_offset", 33))
    q_r1 = QualityAccumulator(quality_threshold, phred_offset)
    q_r2 = QualityAccumulator(quality_threshold, phred_offset)
    q_r3 = QualityAccumulator(quality_threshold, phred_offset)
    temp = out_dir / "MergedCells_1.fastq"
    passing = out_dir / "MergedCells_passing.fastq"

    with open(temp, "w") as out:
        for rec1, rec2 in zip_longest(iter_fastq(args.r1_fastq), iter_fastq(args.r2_fastq)):
            if rec1 is None or rec2 is None:
                stats["unpaired_end_of_file"] += 1
                raise ValueError("R1/R2 FASTQs contain different numbers of records")
            stats["read_pairs_scanned"] += 1
            if _read_id(rec1.name) != _read_id(rec2.name):
                stats["pair_id_mismatch"] += 1
                if not getattr(args, "allow_pair_id_mismatch", False):
                    raise ValueError(f"R1/R2 read IDs are out of sync near pair {stats['read_pairs_scanned']}: {rec1.name} vs {rec2.name}")
            seq2 = rec2.seq
            qual2 = rec2.qual
            # Sequencing-quality QC is calculated directly from the original
            # barcode-read quality string, independently of barcode matching.
            # Q30 follows the standard convention Phred >= 30.
            for region_name, start, acc in [
                ("R1", args.r1_start, q_r1),
                ("R2", args.r2_start, q_r2),
                ("R3", args.r3_start, q_r3),
            ]:
                if len(seq2) >= start + 8 and len(qual2) >= start + 8:
                    acc.add(qual2[start:start + 8])
                else:
                    stats[f"quality_too_short_for_{region_name}"] += 1
            required_len = max(
                args.umi_start + args.umi_len,
                args.r1_start + 8,
                args.r2_start + 8,
                args.r3_start + 8,
            )
            if len(seq2) < required_len:
                stats["too_short_for_barcode_layout"] += 1
                continue
            umi = seq2[args.umi_start:args.umi_start + args.umi_len]
            if len(umi) != args.umi_len:
                stats["invalid_umi_length"] += 1
                continue
            obs3 = seq2[args.r3_start:args.r3_start + 8]
            obs2 = seq2[args.r2_start:args.r2_start + 8]
            obs1 = seq2[args.r1_start:args.r1_start + 8]
            rr1 = c1.correct(obs1, args.ambiguous_policy)
            rr2 = c2.correct(obs2, args.ambiguous_policy)
            rr3 = c3.correct(obs3, args.ambiguous_policy)
            for name, res in [("R1", rr1), ("R2", rr2), ("R3", rr3)]:
                if res.ambiguous:
                    stats[f"{name}_ambiguous"] += 1
                if res.corrected and res.distance == 0:
                    stats[f"{name}_exact"] += 1
                elif res.corrected:
                    stats[f"{name}_corrected"] += 1
                else:
                    stats[f"{name}_failed"] += 1
            if not (rr1.corrected and rr2.corrected and rr3.corrected):
                stats["failed_any_round"] += 1
                continue
            cell = rr1.corrected + rr2.corrected + rr3.corrected
            cell_reads[cell] += 1
            cell_umis[cell].add(umi)
            # Keep cell barcode and UMI in the first whitespace-delimited token
            # so aligners preserve them in BAM QNAMEs. UMI-tools' single-cell
            # convention is readid_CELL_UMI.
            base_name = rec1.name.strip().split()[0]
            new_name = f"{base_name}_{cell}_{umi}"
            out.write(f"{new_name}\n{rec1.seq}\n+\n{rec1.qual}\n")
            stats["demuxed_reads"] += 1

    threshold_mode = getattr(args, "cell_threshold_mode", "umis")
    min_count = int(getattr(args, "min_cell_count", 100))
    metric = {c: (len(cell_umis[c]) if threshold_mode == "umis" else cell_reads[c]) for c in cell_reads}
    keep = {cell for cell, n in metric.items() if n >= min_count}

    with open(temp) as inp, open(passing, "w") as out:
        while True:
            h = inp.readline()
            if not h:
                break
            s = inp.readline(); plus = inp.readline(); q = inp.readline()
            cell, _umi = extract_splitseq_cell_umi_from_header(h.strip())
            if cell in keep:
                out.write(h + s + plus + q)
                stats["passing_reads_written"] += 1
            elif cell is None:
                stats["passing_header_parse_failed"] += 1

    with open(out_dir / "cell_barcode_umi_counts.tsv", "w") as out:
        out.write("cellID\tread_count\tunique_umi_count\tthreshold_value\tpasses_threshold\n")
        for cell in sorted(cell_reads):
            out.write(f"{cell}\t{cell_reads[cell]}\t{len(cell_umis[cell])}\t{metric[cell]}\t{cell in keep}\n")

    quality_rows = [
        q_r1.row("Round1", "all_barcode_reads_with_region_present", args.r1_start, 8),
        q_r2.row("Round2", "all_barcode_reads_with_region_present", args.r2_start, 8),
        q_r3.row("Round3", "all_barcode_reads_with_region_present", args.r3_start, 8),
    ]
    write_quality_qc(out_dir / "barcode_quality_qc.tsv", quality_rows)

    write_kv(out_dir / "demux_summary.tsv", list(stats.items()) + [
        ("round1_whitelist_size", len(r1_list)),
        ("round2_whitelist_size", len(r2_list)),
        ("round3_whitelist_size", len(r3_list)),
        ("unique_cellIDs_detected", len(cell_reads)),
        ("cell_threshold_mode", threshold_mode),
        ("min_cell_count", min_count),
        ("cells_passing_min", len(keep)),
        ("quality_threshold", quality_threshold),
        ("phred_offset", phred_offset),
    ] + summary_rows("round1", q_r1) + summary_rows("round2", q_r2) + summary_rows("round3", q_r3))
    with open(out_dir / "cellIDs_keep.txt", "w") as out:
        for cell in sorted(keep):
            out.write(cell + "\n")

    # MergedCells_1.fastq is only an internal staging file used to construct
    # MergedCells_passing.fastq. Remove it by default after every final demux
    # output exists; retain it only when explicitly requested for debugging.
    require_output_files([
        passing,
        out_dir / "cell_barcode_umi_counts.tsv",
        out_dir / "barcode_quality_qc.tsv",
        out_dir / "demux_summary.tsv",
        out_dir / "cellIDs_keep.txt",
    ], allow_empty=True)
    cleanup_intermediates(
        [temp],
        out_dir / "cleanup_summary.tsv",
        keep_all=bool(getattr(args, "keep_all_intermediates", False)),
        reason="demux staging FASTQ no longer needed after MergedCells_passing.fastq was created",
    )
