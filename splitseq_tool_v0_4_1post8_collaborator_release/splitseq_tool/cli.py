from __future__ import annotations

import argparse
import sys

from .capture import count_multi, count_single
from .config import apply_config, require_args, write_resolved_run_files
from .demux import demux_barcoded_fastq
from .qc import compare_paired
from .qc_report import run_qc_report
from .rna import resume_featurecounts, run_star_pipeline
from .pseudo import prepare_kb_fastqs_from_merged, run_kb_pipeline, write_raw_splitseq_kb_design


def add_config(p):
    p.add_argument(
        "--config",
        default=None,
        help=(
            "Optional YAML config. All commands can still be run with flags only. "
            "When both YAML and flags are provided, flags override YAML values."
        ),
    )


def add_cleanup_control(p):
    p.add_argument(
        "--keep-all-intermediates",
        "--keep_all_intermediates",
        dest="keep_all_intermediates",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Retain large temporary/intermediate files. By default, intermediates are "
            "removed only after the command finishes successfully and final outputs are verified."
        ),
    )


def add_common_capture(p):
    add_config(p)
    add_cleanup_control(p)
    p.add_argument("--sample-id", default=None)
    p.add_argument("--input-fastq", default=None, help="Demultiplexed MergedCells_passing.fastq")
    p.add_argument("--output-dir", default=None)
    p.add_argument("--cell-keep-list", default=None)
    p.add_argument("--cell-filter", choices=["keep-list", "all"], default="keep-list")
    p.add_argument("--round1-conversion", action=argparse.BooleanOptionalAction, default=False, help="Convert capture Round1 barcode to RNA Round1 barcode space")
    p.add_argument("--mrna-round1", default=None, help="RNA Round1 whitelist; required with --round1-conversion")
    p.add_argument("--capture-round1", default=None, help="Capture Round1 whitelist; required with --round1-conversion")
    p.add_argument("--capture-mismatches", type=int, default=1)
    p.add_argument("--count-mode", choices=["reads", "umi"], default="reads")
    p.add_argument("--legacy-output", action=argparse.BooleanOptionalAction, default=False, help="Write optional legacy per-cell compatibility files")
    p.add_argument("--include-zero-cells", action=argparse.BooleanOptionalAction, default=False, help="Include zero-count cells in legacy compatibility output")
    p.add_argument("--quality-threshold", type=int, default=30, help="Phred threshold for capture target quality QC; standard Q30 uses >=30")
    p.add_argument("--phred-offset", type=int, default=33, help="FASTQ ASCII Phred offset; modern Illumina FASTQ is normally 33")


def build_parser():
    p = argparse.ArgumentParser(prog="splitseq-tool")
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("demux", help="Streaming SPLiT-seq demultiplexing core")
    add_config(d)
    add_cleanup_control(d)
    d.add_argument("--sample-id", default=None, help="Optional sample label for reports")
    d.add_argument("--r1-fastq", default=None)
    d.add_argument("--r2-fastq", default=None)
    d.add_argument("--round1", default=None)
    d.add_argument("--round2", default=None)
    d.add_argument("--round3", default=None)
    d.add_argument("--output-dir", default=None)
    d.add_argument("--barcode-mismatches", type=int, default=1)
    d.add_argument("--min-cell-count", "--min-cell-reads", dest="min_cell_count", type=int, default=100, help="Minimum per-cell threshold; default interpretation is unique UMIs")
    d.add_argument("--cell-threshold-mode", choices=["umis", "reads"], default="umis")
    d.add_argument("--allow-pair-id-mismatch", action=argparse.BooleanOptionalAction, default=False, help="Advanced/debug only: continue when paired FASTQ read IDs differ")
    d.add_argument("--ambiguous-policy", choices=["fail", "first", "best"], default="fail")
    d.add_argument("--umi-start", type=int, default=0)
    d.add_argument("--umi-len", type=int, default=10)
    d.add_argument("--r3-start", type=int, default=10)
    d.add_argument("--r2-start", type=int, default=48)
    d.add_argument("--r1-start", type=int, default=86)
    d.add_argument("--quality-threshold", type=int, default=30, help="Phred threshold for barcode quality QC; standard Q30 uses >=30")
    d.add_argument("--phred-offset", type=int, default=33, help="FASTQ ASCII Phred offset; modern Illumina FASTQ is normally 33")
    d.set_defaults(func=demux_barcoded_fastq)

    s = sub.add_parser("capture-single", help="Count one target sequence/barcode in capture reads")
    add_common_capture(s)
    s.add_argument("--target-sequence", default=None)
    s.add_argument("--check-revcomp", action=argparse.BooleanOptionalAction, default=False, help="Also search for the reverse complement of the target")
    s.add_argument("--use-flanks", action=argparse.BooleanOptionalAction, default=False, help="Require/search target within the region delimited by flank-left/flank-right")
    s.add_argument("--flank-left", default=None, help="Optional exact left flank; required when --use-flanks is set")
    s.add_argument("--flank-right", default=None, help="Optional exact right flank; required when --use-flanks is set")
    s.set_defaults(func=count_single)

    m = sub.add_parser("capture-multi", help="Extract/count flank-delimited multiplex capture targets")
    add_common_capture(m)
    m.add_argument(
        "--capture-mode",
        "--capture_mode",
        choices=["conversion", "greenlist", "unrestricted"],
        default="conversion",
        help=(
            "Multiplex target mode. conversion (default): use a barcode-to-sequence mapping table and optional amino-acid translation; "
            "greenlist: count only pre-approved capture barcodes with no sequence mapping/translation; "
            "unrestricted: count every valid extracted capture sequence."
        ),
    )
    m.add_argument(
        "--green-list-only",
        "--green_list_only",
        dest="capture_mode",
        action="store_const",
        const="greenlist",
        help="Shortcut for --capture-mode greenlist.",
    )
    m.add_argument(
        "--green-list-free",
        "--green_list_free",
        dest="capture_mode",
        action="store_const",
        const="unrestricted",
        help="Shortcut for --capture-mode unrestricted.",
    )
    m.add_argument("--barcode-conversion-table", default=None, help="Capture-barcode -> mapped-sequence table; required in conversion mode")
    m.add_argument("--green-list", "--green_list", dest="green_list", default=None, help="Approved barcode list; required in greenlist mode")
    m.add_argument("--flank-left", default=None, help="Exact left flank; required and may be any non-zero length")
    m.add_argument("--flank-right", default=None, help="Exact right flank; required and may be any non-zero length")
    m.add_argument(
        "--target-length",
        "--target_length",
        dest="target_length",
        type=int,
        default=0,
        help="Require an exact extracted target length. 0 means use min/max length controls instead.",
    )
    m.add_argument("--min-insert-len", "--min-target-length", dest="min_insert_len", type=int, default=0, help="Minimum extracted target length; 0 disables. Ignored when --target-length > 0.")
    m.add_argument("--max-insert-len", "--max-target-length", dest="max_insert_len", type=int, default=0, help="Maximum extracted target length; 0 disables. Ignored when --target-length > 0.")
    m.add_argument("--ambiguous-policy", choices=["fail", "first", "best"], default="fail")
    m.add_argument("--translate-aa", action=argparse.BooleanOptionalAction, default=False, help="Conversion mode only: translate the mapped nucleotide sequence to amino-acid sequence")
    m.add_argument("--output-sequence", choices=["barcode", "mapped", "rna", "dna"], default="mapped", help="Conversion mode only when not translating: output the corrected capture barcode or mapped sequence. Legacy aliases rna/dna are accepted.")
    m.add_argument("--progress-every", type=int, default=1000000, help="Print progress every N reads; 0 disables progress messages")
    m.set_defaults(func=count_multi)


    pk = sub.add_parser("prepare-kb-fastqs", help="Prepare kb/kallisto-bustools compatible FASTQs from demuxed SPLiT-seq FASTQ")
    add_config(pk)
    add_cleanup_control(pk)
    pk.add_argument("--sample-id", default=None)
    pk.add_argument("--merged-fastq", default=None)
    pk.add_argument("--output-dir", default=None)
    pk.add_argument("--cell-whitelist", default=None, help="Optional cell whitelist; default derives cells from merged FASTQ")
    pk.add_argument("--cell-barcode-len", type=int, default=24)
    pk.add_argument("--umi-len", type=int, default=10)
    pk.add_argument("--gzip-output", action=argparse.BooleanOptionalAction, default=False)
    pk.set_defaults(func=prepare_kb_fastqs_from_merged)

    kb = sub.add_parser("rna-kb", help="Run pseudoalignment with kb/kallisto-bustools after SPLiT-seq demux")
    add_config(kb)
    add_cleanup_control(kb)
    kb.add_argument("--sample-id", default=None)
    kb.add_argument("--merged-fastq", default=None, help="Demultiplexed MergedCells_passing.fastq with CB/UMI in header")
    kb.add_argument("--output-dir", default=None)
    kb.add_argument("--kallisto-index", default=None)
    kb.add_argument("--t2g-file", default=None, help="Transcript-to-gene map for kb/bustools")
    kb.add_argument("--cell-whitelist", default=None, help="Optional on-list; default derives cells from merged FASTQ")
    kb.add_argument("--cell-barcode-len", type=int, default=24)
    kb.add_argument("--umi-len", type=int, default=10)
    kb.add_argument("--threads", type=int, default=10)
    kb.add_argument("--technology-string", default=None, help="Override kb technology string; default is synthetic 24nt CB + 10nt UMI")
    kb.add_argument("--workflow", default="standard")
    kb.add_argument("--strand", choices=["forward", "reverse", "unstranded"], default=None)
    kb.add_argument("--h5ad", action=argparse.BooleanOptionalAction, default=False)
    kb.add_argument("--keep-tmp", action=argparse.BooleanOptionalAction, default=False)
    kb.add_argument("--gzip-fastqs", action=argparse.BooleanOptionalAction, default=False)
    kb.add_argument("--kb-extra-args", default=None, help="Extra raw arguments appended to kb count")
    kb.set_defaults(func=run_kb_pipeline)

    rk = sub.add_parser("raw-kb-design", help="Write experimental raw FASTQ kb technology string and 24nt whitelist")
    add_config(rk)
    add_cleanup_control(rk)
    rk.add_argument("--round1", default=None)
    rk.add_argument("--round2", default=None)
    rk.add_argument("--round3", default=None)
    rk.add_argument("--output-dir", default=None)
    rk.add_argument("--umi-start", type=int, default=0)
    rk.add_argument("--umi-len", type=int, default=10)
    rk.add_argument("--r3-start", type=int, default=10)
    rk.add_argument("--r2-start", type=int, default=48)
    rk.add_argument("--r1-start", type=int, default=86)
    rk.set_defaults(func=write_raw_splitseq_kb_design)

    q = sub.add_parser("paired-qc", help="Compare mRNA and capture cell-ID spaces")
    add_config(q)
    add_cleanup_control(q)
    q.add_argument("--sample-id", default=None)
    q.add_argument("--mrna-fastq", default=None, help="Demultiplexed RNA FASTQ")
    q.add_argument("--capture-fastq", default=None, help="Demultiplexed capture FASTQ")
    q.add_argument("--output-dir", default=None)
    q.add_argument("--round1-conversion", action=argparse.BooleanOptionalAction, default=False)
    q.add_argument("--mrna-round1", default=None)
    q.add_argument("--capture-round1", default=None)
    q.add_argument("--possible-cell-ids", type=int, default=442368)
    q.set_defaults(func=compare_paired)


    qr = sub.add_parser("qc-report", help="Create QC summaries and plots from RNA/capture/paired outputs")
    add_config(qr)
    add_cleanup_control(qr)
    qr.add_argument("--sample-id", default=None)
    qr.add_argument("--output-dir", default=None)
    qr.add_argument("--demux-summary", default=None, help="Optional demux_summary.tsv")
    qr.add_argument("--rna-summary", default=None, help="Optional RNA pipeline summary TSV")
    qr.add_argument("--capture-summary", default=None, help="Optional capture summary TSV")
    qr.add_argument("--paired-summary", default=None, help="Optional paired-qc summary TSV")
    qr.add_argument("--rna-counts-matrix", default=None, help="Optional gene x cell RNA count matrix, TSV/TSV.GZ")
    qr.add_argument("--capture-counts-wide", default=None, help="Optional cell x capture-feature count table")
    qr.add_argument("--min-capture-count", "--min-capture-reads", dest="min_capture_count", type=int, default=1, help="Optional minimum capture count per cell for QC summaries; works for read- or UMI-count tables")
    qr.add_argument("--min-rna-counts", type=int, default=0, help="Optional RNA counts-per-cell filter for QC summaries")
    qr.add_argument("--min-genes", type=int, default=0, help="Optional genes-per-cell filter for QC summaries")
    qr.set_defaults(func=run_qc_report)

    r = sub.add_parser("rna-star", help="Run STAR/featureCounts/umi_tools after demux")
    add_config(r)
    add_cleanup_control(r)
    r.add_argument("--sample-id", default=None)
    r.add_argument("--merged-fastq", default=None)
    r.add_argument("--output-dir", default=None)
    r.add_argument("--star-index", default=None)
    r.add_argument("--saf-file", default=None)
    r.add_argument("--threads", type=int, default=10)
    r.add_argument("--star-mismatch-nover-lmax", type=float, default=0.05)
    r.add_argument("--align-intron-max", type=int, default=20000)
    r.set_defaults(func=run_star_pipeline)

    rf = sub.add_parser("resume-featurecounts", help="Resume RNA run from aligned BAM")
    add_config(rf)
    add_cleanup_control(rf)
    rf.add_argument("--sample-id", default=None)
    rf.add_argument("--aligned-bam", default=None)
    rf.add_argument("--output-dir", default=None)
    rf.add_argument("--saf-file", default=None)
    rf.add_argument("--threads", type=int, default=10)
    rf.set_defaults(func=resume_featurecounts)
    return p


def validate_args(args, parser):
    cmd = args.cmd
    if cmd == "demux":
        require_args(args, parser, ["r1_fastq", "r2_fastq", "round1", "round2", "round3", "output_dir"])
    elif cmd == "capture-single":
        require_args(args, parser, ["sample_id", "input_fastq", "output_dir", "target_sequence"])
        if args.cell_filter == "keep-list":
            require_args(args, parser, ["cell_keep_list"])
        if args.round1_conversion:
            require_args(args, parser, ["mrna_round1", "capture_round1"])
        if args.use_flanks:
            require_args(args, parser, ["flank_left", "flank_right"])
        if not all(c in "ACGTacgt" for c in args.target_sequence):
            parser.error("--target-sequence must contain only A/C/G/T")
        if args.use_flanks and (not all(c in "ACGTacgt" for c in args.flank_left) or not all(c in "ACGTacgt" for c in args.flank_right)):
            parser.error("--flank-left and --flank-right must contain only A/C/G/T")
    elif cmd == "capture-multi":
        require_args(args, parser, ["sample_id", "input_fastq", "output_dir"])
        if args.cell_filter == "keep-list":
            require_args(args, parser, ["cell_keep_list"])
        if args.round1_conversion:
            require_args(args, parser, ["mrna_round1", "capture_round1"])
        if not getattr(args, "flank_left", None) or not getattr(args, "flank_right", None):
            parser.error("capture-multi requires non-empty --flank-left and --flank-right; flank sequences may be any length")
        if not all(c in "ACGTacgt" for c in args.flank_left) or not all(c in "ACGTacgt" for c in args.flank_right):
            parser.error("--flank-left and --flank-right must contain only A/C/G/T")
        if getattr(args, "target_length", 0) < 0:
            parser.error("--target-length must be >= 0")
        if getattr(args, "min_insert_len", 0) < 0 or getattr(args, "max_insert_len", 0) < 0:
            parser.error("--min-target-length/--max-target-length must be >= 0")
        if getattr(args, "target_length", 0) == 0 and getattr(args, "min_insert_len", 0) and getattr(args, "max_insert_len", 0) and args.min_insert_len > args.max_insert_len:
            parser.error("--min-target-length cannot exceed --max-target-length")
        mode = getattr(args, "capture_mode", "conversion")
        if mode == "conversion":
            require_args(args, parser, ["barcode_conversion_table"])
        elif mode == "greenlist":
            require_args(args, parser, ["green_list"])
            if getattr(args, "translate_aa", False):
                parser.error("--translate-aa is only valid with --capture-mode conversion")
        elif mode == "unrestricted":
            if getattr(args, "translate_aa", False):
                parser.error("--translate-aa is only valid with --capture-mode conversion")
    elif cmd == "prepare-kb-fastqs":
        require_args(args, parser, ["merged_fastq", "output_dir"])
    elif cmd == "rna-kb":
        require_args(args, parser, ["sample_id", "merged_fastq", "output_dir", "kallisto_index", "t2g_file"])
    elif cmd == "raw-kb-design":
        require_args(args, parser, ["round1", "round2", "round3", "output_dir"])
    elif cmd == "paired-qc":
        require_args(args, parser, ["sample_id", "mrna_fastq", "capture_fastq", "output_dir"])
    elif cmd == "qc-report":
        require_args(args, parser, ["sample_id", "output_dir"])
    elif cmd == "rna-star":
        require_args(args, parser, ["sample_id", "merged_fastq", "output_dir", "star_index", "saf_file"])
    elif cmd == "resume-featurecounts":
        require_args(args, parser, ["aligned_bam", "output_dir", "saf_file"])

    if hasattr(args, "barcode_mismatches") and int(args.barcode_mismatches) < 0:
        parser.error("--barcode-mismatches must be >= 0")
    if hasattr(args, "capture_mismatches") and int(args.capture_mismatches) < 0:
        parser.error("--capture-mismatches must be >= 0")
    if hasattr(args, "min_cell_count") and int(args.min_cell_count) < 0:
        parser.error("--min-cell-count must be >= 0")
    if hasattr(args, "threads") and int(args.threads) <= 0:
        parser.error("--threads must be > 0")
    for coord_name in ("umi_start", "r1_start", "r2_start", "r3_start"):
        if hasattr(args, coord_name) and int(getattr(args, coord_name)) < 0:
            parser.error(f"--{coord_name.replace('_', '-')} must be >= 0")
    if hasattr(args, "umi_len") and int(args.umi_len) <= 0:
        parser.error("--umi-len must be > 0")
    if hasattr(args, "possible_cell_ids") and int(args.possible_cell_ids) <= 0:
        parser.error("--possible-cell-ids must be > 0")
    if hasattr(args, "progress_every") and int(args.progress_every) < 0:
        parser.error("--progress-every must be >= 0")

    if hasattr(args, "quality_threshold") and int(args.quality_threshold) < 0:
        parser.error("--quality-threshold must be >= 0")
    if hasattr(args, "phred_offset") and int(args.phred_offset) <= 0:
        parser.error("--phred-offset must be > 0")

    if getattr(args, "round1_conversion", False):
        if not getattr(args, "mrna_round1", None) or not getattr(args, "capture_round1", None):
            parser.error("--round1-conversion requires --mrna-round1 and --capture-round1, supplied either by flags or YAML")
    if getattr(args, "cell_filter", None) == "keep-list" and not getattr(args, "cell_keep_list", None):
        parser.error("--cell-filter keep-list requires --cell-keep-list; use --cell-filter all to keep all cells")


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    parser = build_parser()
    args = parser.parse_args(argv)
    args = apply_config(args, parser, argv)
    validate_args(args, parser)
    write_resolved_run_files(args, args.cmd, getattr(args, "output_dir", None))
    args.func(args)


if __name__ == "__main__":
    main()
