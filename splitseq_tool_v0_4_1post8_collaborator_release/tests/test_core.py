from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import pytest

from splitseq_tool.barcodes import BarcodeCorrector
from splitseq_tool.capture import (
    extract_between_flanks_with_coords,
    load_conversion_table,
    translate_dna,
)
from splitseq_tool.cli import build_parser
from splitseq_tool.config import apply_config
from splitseq_tool.fastq import iter_fastq
from splitseq_tool.cli import main as cli_main
from splitseq_tool.pseudo import prepare_kb_fastqs_from_merged


def write_fastq(path: Path, records):
    with open(path, "w") as out:
        for name, seq, qual in records:
            out.write(f"@{name}\n{seq}\n+\n{qual}\n")


def make_barcode_read(r1="AAAAAAAA", r2="GGGGGGGG", r3="TTTTTTTT", umi="ACGTACGTAA"):
    seq = ["C"] * 94
    seq[0:10] = list(umi)
    seq[10:18] = list(r3)
    seq[48:56] = list(r2)
    seq[86:94] = list(r1)
    return "".join(seq)


def test_exact_match_outranks_one_mismatch_neighbor():
    corrector = BarcodeCorrector(["AAAAAAAA", "AAAAAAAC"], 1)
    exact = corrector.correct("AAAAAAAA", "fail")
    assert exact.corrected == "AAAAAAAA"
    assert exact.distance == 0
    assert exact.ambiguous is False

    tie = corrector.correct("AAAAAAAG", "fail")
    assert tie.corrected is None
    assert tie.ambiguous is True


def test_variable_length_correction_never_crosses_lengths():
    corrector = BarcodeCorrector(["AAAA", "AAAAA"], 2)
    hit = corrector.correct("AAAT", "fail")
    assert hit.corrected == "AAAA"
    assert hit.distance == 1


def test_arbitrary_flanks_skip_internal_right_flank_when_length_known():
    # Right flank "C" also occurs inside the target. Length constraint must
    # select the later boundary that yields the intended 4-nt target.
    hit = extract_between_flanks_with_coords(
        "AAGCCGC",
        left="AA",
        right="C",
        allowed_lengths={4},
    )
    assert hit is not None
    assert hit[0] == "GCCG"


def test_conflicting_conversion_mapping_is_rejected(tmp_path):
    table = tmp_path / "map.csv"
    table.write_text("ACGT,GCTGCT\nACGT,GGCGGC\n")
    with pytest.raises(ValueError, match="Conflicting mappings"):
        load_conversion_table(table)


def test_translation_basic():
    assert translate_dna("GCTGCCGCA") == "AAA"


def test_fastq_truncation_is_error(tmp_path):
    fq = tmp_path / "bad.fastq"
    fq.write_text("@r1\nACGT\n+\n")
    with pytest.raises(ValueError, match="Truncated FASTQ"):
        list(iter_fastq(fq))


def test_fastq_sequence_quality_mismatch_is_error(tmp_path):
    fq = tmp_path / "bad.fastq"
    fq.write_text("@r1\nACGT\n+\nIII\n")
    with pytest.raises(ValueError, match="sequence/quality length mismatch"):
        list(iter_fastq(fq))


def test_demux_round_qc_and_cleanup(tmp_path):
    r1fq = tmp_path / "reads_R1.fastq"
    r2fq = tmp_path / "reads_R2.fastq"
    bseq = make_barcode_read()
    write_fastq(r1fq, [("read1/1", "ACGTACGTACGT", "I" * 12)])
    write_fastq(r2fq, [("read1/2", bseq, "I" * len(bseq))])

    wl1 = tmp_path / "r1.txt"; wl1.write_text("AAAAAAAA\n")
    wl2 = tmp_path / "r2.txt"; wl2.write_text("GGGGGGGG\n")
    wl3 = tmp_path / "r3.txt"; wl3.write_text("TTTTTTTT\n")
    out = tmp_path / "demux"

    cli_main([
        "demux",
        "--r1-fastq", str(r1fq),
        "--r2-fastq", str(r2fq),
        "--round1", str(wl1),
        "--round2", str(wl2),
        "--round3", str(wl3),
        "--min-cell-count", "1",
        "--output-dir", str(out),
    ])

    assert (out / "MergedCells_passing.fastq").stat().st_size > 0
    assert not (out / "MergedCells_1.fastq").exists()
    qc = (out / "barcode_quality_qc.tsv").read_text()
    assert "Round1" in qc and "100.0" in qc


def make_merged_fastq(path: Path, target: str):
    cell = "A" * 24
    umi = "ACGTACGTAA"
    seq = "TTAA" + target + "CCGG"
    write_fastq(path, [(f"read1_{cell}_{umi}", seq, "I" * len(seq))])
    return cell


def test_capture_multi_greenlist_and_unrestricted(tmp_path):
    merged = tmp_path / "merged.fastq"
    cell = make_merged_fastq(merged, "ACGA")
    green = tmp_path / "green.txt"; green.write_text("ACGT\n")

    out_green = tmp_path / "green_out"
    cli_main([
        "capture-multi",
        "--sample-id", "sample_01",
        "--input-fastq", str(merged),
        "--output-dir", str(out_green),
        "--cell-filter", "all",
        "--capture-mode", "greenlist",
        "--green-list", str(green),
        "--flank-left", "TTAA",
        "--flank-right", "CCGG",
        "--target-length", "4",
        "--capture-mismatches", "1",
        "--progress-every", "0",
    ])
    long_green = (out_green / "sample_01.multi_capture_counts.long.tsv").read_text()
    assert f"{cell}\tACGT\t1" in long_green

    out_free = tmp_path / "free_out"
    cli_main([
        "capture-multi",
        "--sample-id", "sample_01",
        "--input-fastq", str(merged),
        "--output-dir", str(out_free),
        "--cell-filter", "all",
        "--capture-mode", "unrestricted",
        "--flank-left", "TTAA",
        "--flank-right", "CCGG",
        "--target-length", "4",
        "--progress-every", "0",
    ])
    long_free = (out_free / "sample_01.multi_capture_counts.long.tsv").read_text()
    assert f"{cell}\tACGA\t1" in long_free


def test_capture_multi_conversion_translation(tmp_path):
    merged = tmp_path / "merged.fastq"
    cell = make_merged_fastq(merged, "ACGT")
    mapping = tmp_path / "map.csv"; mapping.write_text("barcode,mapped\nACGT,GCTGCT\n")
    out = tmp_path / "out"
    cli_main([
        "capture-multi",
        "--sample-id", "sample_01",
        "--input-fastq", str(merged),
        "--output-dir", str(out),
        "--cell-filter", "all",
        "--capture-mode", "conversion",
        "--barcode-conversion-table", str(mapping),
        "--flank-left", "TTAA",
        "--flank-right", "CCGG",
        "--target-length", "4",
        "--capture-mismatches", "0",
        "--translate-aa",
        "--progress-every", "0",
    ])
    text = (out / "sample_01.multi_capture_counts.long.tsv").read_text()
    assert f"{cell}\tAA\t1" in text


def test_boolean_cli_can_override_yaml_true(tmp_path):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(
        "defaults:\n"
        "  use_flanks: true\n"
        "capture:\n"
        "  single:\n"
        "    flank_left: TTAA\n"
        "    flank_right: CCGG\n"
    )
    parser = build_parser()
    argv = [
        "capture-single", "--config", str(cfg), "--no-use-flanks",
        "--sample-id", "sample_01", "--input-fastq", "x.fastq",
        "--output-dir", "out", "--target-sequence", "ACGT",
        "--cell-filter", "all",
    ]
    args = parser.parse_args(argv)
    args = apply_config(args, parser, argv)
    assert args.use_flanks is False


def test_prepare_kb_technology_string_tracks_custom_lengths(tmp_path):
    merged = tmp_path / "merged.fastq"
    cell = "A" * 24
    umi = "ACGTACGT"
    write_fastq(merged, [(f"read1_{cell}_{umi}", "ACGTACGT", "I" * 8)])
    out = tmp_path / "kb"
    args = Namespace(
        output_dir=str(out), merged_fastq=str(merged), sample_id="sample_01",
        cell_whitelist=None, cell_barcode_len=24, umi_len=8, gzip_output=False,
    )
    prepare_kb_fastqs_from_merged(args)
    inputs = (out / "kb_inputs.tsv").read_text()
    assert "technology_string\t0,0,24:0,24,32:1,0,0" in inputs


def test_demux_keep_all_intermediates_override(tmp_path):
    r1fq = tmp_path / "reads_R1.fastq"
    r2fq = tmp_path / "reads_R2.fastq"
    bseq = make_barcode_read()
    write_fastq(r1fq, [("read1/1", "ACGTACGT", "I" * 8)])
    write_fastq(r2fq, [("read1/2", bseq, "I" * len(bseq))])
    wl1 = tmp_path / "r1.txt"; wl1.write_text("AAAAAAAA\n")
    wl2 = tmp_path / "r2.txt"; wl2.write_text("GGGGGGGG\n")
    wl3 = tmp_path / "r3.txt"; wl3.write_text("TTTTTTTT\n")
    out = tmp_path / "demux"
    cli_main([
        "demux", "--r1-fastq", str(r1fq), "--r2-fastq", str(r2fq),
        "--round1", str(wl1), "--round2", str(wl2), "--round3", str(wl3),
        "--min-cell-count", "1", "--keep-all-intermediates",
        "--output-dir", str(out),
    ])
    assert (out / "MergedCells_1.fastq").exists()
    assert "\tkept\t" in (out / "cleanup_summary.tsv").read_text()


def test_capture_single_q30_metrics(tmp_path):
    merged = tmp_path / "merged.fastq"
    cell = "A" * 24
    umi1 = "ACGTACGTAA"
    umi2 = "ACGTACGTAC"
    seq = "TTACGTAA"
    write_fastq(merged, [
        (f"read1_{cell}_{umi1}", seq, "I" * len(seq)),
        (f"read2_{cell}_{umi2}", seq, "II" + "5" * 4 + "II"),
    ])
    out = tmp_path / "single"
    cli_main([
        "capture-single", "--sample-id", "sample_01",
        "--input-fastq", str(merged), "--output-dir", str(out),
        "--cell-filter", "all", "--target-sequence", "ACGT",
        "--quality-threshold", "30",
    ])
    qc = (out / "sample_01.capture_quality_qc.tsv").read_text().splitlines()
    assert len(qc) == 2
    fields = qc[0].split("\t")
    values = qc[1].split("\t")
    row = dict(zip(fields, values))
    assert float(row["percent_bases_Q_ge_threshold"]) == 50.0
    assert float(row["percent_observations_all_bases_Q_ge_threshold"]) == 50.0


def test_translation_requires_codon_aligned_mapped_sequences(tmp_path):
    merged = tmp_path / "merged.fastq"
    make_merged_fastq(merged, "ACGT")
    mapping = tmp_path / "map.csv"; mapping.write_text("ACGT,GCTA\n")
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="divisible by 3"):
        cli_main([
            "capture-multi", "--sample-id", "sample_01",
            "--input-fastq", str(merged), "--output-dir", str(out),
            "--cell-filter", "all", "--capture-mode", "conversion",
            "--barcode-conversion-table", str(mapping),
            "--flank-left", "TTAA", "--flank-right", "CCGG",
            "--target-length", "4", "--translate-aa",
        ])


def test_capture_all_mode_retains_zero_count_cells(tmp_path):
    cell1 = "A" * 24
    cell2 = "C" * 24
    umi1 = "ACGTACGTAA"
    umi2 = "ACGTACGTAC"
    merged = tmp_path / "merged.fastq"
    write_fastq(merged, [
        (f"r1_{cell1}_{umi1}", "TTACGTAA", "I" * 8),
        (f"r2_{cell2}_{umi2}", "TTTTTTTT", "I" * 8),
    ])
    out = tmp_path / "out"
    cli_main([
        "capture-single", "--sample-id", "sample_01",
        "--input-fastq", str(merged), "--output-dir", str(out),
        "--cell-filter", "all", "--target-sequence", "ACGT",
    ])
    lines = (out / "sample_01.single_capture_counts.wide.tsv").read_text().splitlines()
    assert lines[0] == "cellID\tACGT_count"
    rows = {line.split("\t")[0]: line.split("\t")[1] for line in lines[1:]}
    assert rows[cell1] == "1"
    assert rows[cell2] == "0"


def test_raw_kb_design_uses_configured_coordinates(tmp_path):
    wl1 = tmp_path / "r1.txt"; wl1.write_text("AAAAAAAA\n")
    wl2 = tmp_path / "r2.txt"; wl2.write_text("CCCCCCCC\n")
    wl3 = tmp_path / "r3.txt"; wl3.write_text("GGGGGGGG\n")
    out = tmp_path / "rawkb"
    cli_main([
        "raw-kb-design",
        "--round1", str(wl1), "--round2", str(wl2), "--round3", str(wl3),
        "--output-dir", str(out),
        "--umi-start", "2", "--umi-len", "8",
        "--r3-start", "20", "--r2-start", "40", "--r1-start", "60",
    ])
    text = (out / "raw_kb_design.tsv").read_text()
    assert "technology_string\t1,60,68,1,40,48,1,20,28:1,2,10:0,0,0" in text


def test_round1_conversion_rejects_duplicate_rows(tmp_path):
    from splitseq_tool.barcodes import make_round1_conversion
    rna = tmp_path / "rna.txt"; rna.write_text("AAAAAAAA\nAAAAAAAA\n")
    cap = tmp_path / "cap.txt"; cap.write_text("CCCCCCCC\nGGGGGGGG\n")
    with pytest.raises(ValueError, match="duplicate"):
        make_round1_conversion(rna, cap)


def test_prepare_kb_explicit_whitelist_filters_reads(tmp_path):
    merged = tmp_path / "merged.fastq"
    cell_keep = "A" * 24
    cell_drop = "C" * 24
    umi1 = "ACGTACGTAA"
    umi2 = "ACGTACGTAC"
    write_fastq(merged, [
        (f"r1_{cell_keep}_{umi1}", "ACGTACGT", "I" * 8),
        (f"r2_{cell_drop}_{umi2}", "TGCATGCA", "I" * 8),
    ])
    whitelist = tmp_path / "cells.txt"
    whitelist.write_text(cell_keep + "\n")
    out = tmp_path / "kb"
    args = Namespace(
        output_dir=str(out), merged_fastq=str(merged), sample_id="sample_01",
        cell_whitelist=str(whitelist), cell_barcode_len=24, umi_len=10, gzip_output=False,
    )
    prepare_kb_fastqs_from_merged(args)
    barcode_fq = out / "sample_01.kb_barcode_umi_R1.fastq"
    records = list(iter_fastq(barcode_fq))
    assert len(records) == 1
    assert records[0].seq.startswith(cell_keep)
    assert (out / "sample_01.kb_cell_whitelist.txt").read_text().strip() == cell_keep


def test_prepare_kb_rejects_umi_length_mismatch(tmp_path):
    merged = tmp_path / "merged.fastq"
    cell = "A" * 24
    umi = "ACGTACGTAA"
    write_fastq(merged, [(f"r1_{cell}_{umi}", "ACGT", "I" * 4)])
    out = tmp_path / "kb"
    args = Namespace(
        output_dir=str(out), merged_fastq=str(merged), sample_id="sample_01",
        cell_whitelist=None, cell_barcode_len=24, umi_len=8, gzip_output=False,
    )
    with pytest.raises(ValueError, match="Unexpected UMI length"):
        prepare_kb_fastqs_from_merged(args)


def test_qc_report_capture_count_threshold_uses_new_generic_flag(tmp_path):
    wide = tmp_path / "capture.tsv"
    cell1 = "A" * 24
    cell2 = "C" * 24
    wide.write_text(
        "cellID\tFEATURE_count\n"
        f"{cell1}\t3\n"
        f"{cell2}\t1\n"
    )
    out = tmp_path / "qc"
    cli_main([
        "qc-report", "--sample-id", "sample_01", "--output-dir", str(out),
        "--capture-counts-wide", str(wide), "--min-capture-count", "2",
    ])
    summary = (out / "sample_01.qc_summary.tsv").read_text()
    assert "capture_min_counts_filter\t2" in summary
    assert "capture_cells_passing_min_counts\t1" in summary
    per_cell = (out / "sample_01.per_cell_qc_metrics.tsv").read_text().splitlines()
    assert "capture_total_counts" in per_cell[0]
    assert "capture_total_reads" not in per_cell[0]


def test_legacy_output_rerun_removes_stale_cell_files(tmp_path):
    from splitseq_tool.capture import write_legacy
    out = tmp_path / "legacy"
    cell1 = "A" * 24
    cell2 = "C" * 24
    write_legacy({(cell1, "FEATURE"): 1, (cell2, "FEATURE"): 1}, out)
    assert (out / cell1).exists() and (out / cell2).exists()
    write_legacy({(cell1, "FEATURE"): 2}, out)
    assert (out / cell1).exists()
    assert not (out / cell2).exists()


def test_prepare_kb_rejects_whitelist_with_no_observed_cells(tmp_path):
    merged = tmp_path / "merged.fastq"
    observed = "A" * 24
    requested = "C" * 24
    umi = "ACGTACGTAA"
    write_fastq(merged, [(f"r1_{observed}_{umi}", "ACGT", "I" * 4)])
    whitelist = tmp_path / "cells.txt"
    whitelist.write_text(requested + "\n")
    out = tmp_path / "kb"
    args = Namespace(
        output_dir=str(out), merged_fastq=str(merged), sample_id="sample_01",
        cell_whitelist=str(whitelist), cell_barcode_len=24, umi_len=10, gzip_output=False,
    )
    with pytest.raises(ValueError, match="No records were written"):
        prepare_kb_fastqs_from_merged(args)


def test_single_target_match_prefers_exact_over_earlier_mismatch():
    from splitseq_tool.capture import find_sequence_match
    # ACGA is a 1-mismatch candidate at position 0; exact ACGT occurs later.
    hit = find_sequence_match("ACGATTTACGT", "ACGT", mismatches=1)
    assert hit == (7, 11, "forward")


def test_single_target_match_prefers_exact_reverse_complement_over_forward_mismatch():
    from splitseq_tool.capture import find_sequence_match
    # Target AGTC: forward near-match AGTA at position 0, exact reverse-complement GACT later.
    hit = find_sequence_match("AGTATTTGACT", "AGTC", mismatches=1, check_revcomp=True)
    assert hit == (7, 11, "reverse_complement")
