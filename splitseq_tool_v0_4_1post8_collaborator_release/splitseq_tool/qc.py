from __future__ import annotations

from collections import Counter
from pathlib import Path

from .barcodes import convert_cell_round1, make_round1_conversion
from .fastq import extract_splitseq_cell_umi_from_header, iter_fastq
from .utils import ensure_dir, write_kv


def count_cells_in_fastq(path: str | Path) -> Counter:
    c = Counter()
    for rec in iter_fastq(path):
        cell, _ = extract_splitseq_cell_umi_from_header(rec.name)
        if cell:
            c[cell[:24]] += 1
    return c


def compare_paired(args) -> None:
    out_dir = ensure_dir(args.output_dir)
    mrna_counts = count_cells_in_fastq(args.mrna_fastq)
    cap_raw = count_cells_in_fastq(args.capture_fastq)
    conv_map = make_round1_conversion(args.mrna_round1, args.capture_round1) if args.round1_conversion else None
    cap_conv = Counter()
    failed = 0
    for cell, n in cap_raw.items():
        if conv_map:
            new = convert_cell_round1(cell, conv_map)
            if new is None:
                failed += 1
                continue
            cap_conv[new] += n
        else:
            cap_conv[cell] += n
    mrna_set = set(mrna_counts)
    raw_set = set(cap_raw)
    conv_set = set(cap_conv)
    possible = args.possible_cell_ids
    overlap_raw = mrna_set & raw_set
    overlap_conv = mrna_set & conv_set
    expected = len(mrna_set) * len(conv_set) / possible if possible else "NA"
    ratio = len(overlap_conv) / expected if isinstance(expected, float) and expected else "NA"
    rows = [
        ("possible_cell_ids", possible),
        ("mRNA_unique_cellIDs_observed", len(mrna_set)),
        ("capture_raw_unique_cellIDs_observed", len(raw_set)),
        ("capture_converted_unique_cellIDs_observed", len(conv_set)),
        ("capture_raw_unique_failed_R1_conversion", failed),
        ("raw_overlap_unique_cellIDs", len(overlap_raw)),
        ("converted_overlap_unique_cellIDs", len(overlap_conv)),
        ("converted_overlap_pct_mRNA_observed", round(100 * len(overlap_conv) / len(mrna_set), 4) if mrna_set else 0),
        ("converted_overlap_pct_capture_converted_observed", round(100 * len(overlap_conv) / len(conv_set), 4) if conv_set else 0),
        ("expected_random_overlap", round(expected, 4) if isinstance(expected, float) else expected),
        ("observed_over_expected_overlap", round(ratio, 4) if isinstance(ratio, float) else ratio),
        ("mRNA_total_reads_in_observed_cells", sum(mrna_counts.values())),
        ("capture_total_reads_in_converted_cells", sum(cap_conv.values())),
    ]
    write_kv(out_dir / f"{args.sample_id}.paired_cellID_qc.summary.tsv", rows)
    with open(out_dir / f"{args.sample_id}.paired_cellID_qc.overlap.tsv", "w") as out:
        out.write("cellID_mRNA_space\tmRNA_read_count\tcapture_converted_read_count\tin_mRNA\tin_capture_after_R1_conversion\n")
        for cell in sorted(mrna_set | conv_set):
            out.write(f"{cell}\t{mrna_counts.get(cell,0)}\t{cap_conv.get(cell,0)}\t{cell in mrna_set}\t{cell in conv_set}\n")
