from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import shutil
import time
from typing import Iterable, Optional

from .barcodes import BarcodeCorrector, convert_cell_round1, load_barcode_list, make_round1_conversion
from .fastq import extract_splitseq_cell_umi_from_header, iter_fastq
from .utils import ensure_dir, hamming, revcomp, write_kv
from .quality import QualityAccumulator, summary_rows, write_quality_qc

CODON = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"STOP","TAG":"STOP","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"STOP","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def translate_dna(seq: str) -> str:
    seq = seq.upper()
    aas = []
    for i in range(0, len(seq) - 2, 3):
        aas.append(CODON.get(seq[i:i+3], "X"))
    return "".join(aas)


def load_keep_cells(path: str | Path | None) -> Optional[set[str]]:
    if not path:
        return None
    cells = set()
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s:
                cells.add(s)
    return cells


def load_conversion_table(path: str | Path) -> dict[str, str]:
    """Load capture barcode -> mapped sequence pairs.

    Duplicate barcode rows are allowed only when they map to the same sequence.
    Conflicting duplicate mappings raise an error rather than being silently
    overwritten.
    """
    d: dict[str, str] = {}
    with open(path) as f:
        for line_number, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace("\t", ",").split(",")
            if len(parts) < 2:
                continue
            barcode = parts[0].strip().upper()
            mapped = parts[1].strip().upper()
            if not barcode or not mapped:
                continue
            if not all(c in "ACGT" for c in barcode):
                continue
            if not all(c in "ACGTN" for c in mapped):
                continue
            if barcode in d and d[barcode] != mapped:
                raise ValueError(
                    f"Conflicting mappings for capture barcode {barcode!r} in {path} "
                    f"(line {line_number})"
                )
            d[barcode] = mapped
    return d


def load_green_list(path: str | Path) -> list[str]:
    """Load an approved capture-barcode list of arbitrary sequence length.

    Accepts one sequence per line or a CSV/TSV whose first field is the barcode.
    Comment/header/non-ACGT rows are ignored. Order is preserved and duplicates
    are removed.
    """
    vals: list[str] = []
    seen: set[str] = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            first_field = line.replace("\t", ",").split(",", 1)[0].strip()
            if not first_field:
                continue
            token = first_field.split()[0].upper()
            if not token or not all(c in "ACGT" for c in token):
                continue
            if token not in seen:
                seen.add(token)
                vals.append(token)
    return vals


def _validate_reference_lengths(seqs: Iterable[str], target_length: int, label: str) -> None:
    """Require reference barcodes to match --target-length when it is set."""
    if not target_length:
        return
    bad = sorted({len(s) for s in seqs if len(s) != target_length})
    if bad:
        raise ValueError(
            f"{label} contains sequence lengths incompatible with --target-length {target_length}; "
            f"other lengths observed: {bad}"
        )


def find_sequence_match(seq: str, target: str, mismatches: int = 0, check_revcomp: bool = False) -> Optional[tuple[int, int, str]]:
    """Return the best matched target span as ``(start, end, orientation)``.

    When mismatches are allowed, all candidate windows are compared and the
    lowest-Hamming-distance match is preferred. Ties are resolved by earlier
    sequence position and then forward orientation. This prevents an early
    one-mismatch window from masking a later exact target occurrence.
    """
    targets = [(target.upper(), "forward", 0)]
    if check_revcomp:
        rc = revcomp(target.upper())
        if rc != target.upper():
            targets.append((rc, "reverse_complement", 1))
    seq = seq.upper()
    best: Optional[tuple[int, int, int, int, str]] = None
    for t, orientation, orientation_rank in targets:
        L = len(t)
        if L == 0 or len(seq) < L:
            continue
        if mismatches == 0:
            pos = seq.find(t)
            if pos >= 0:
                candidate = (0, pos, orientation_rank, pos + L, orientation)
                if best is None or candidate[:3] < best[:3]:
                    best = candidate
            continue
        for i in range(0, len(seq) - L + 1):
            d = hamming(seq[i:i+L], t)
            if d <= mismatches:
                candidate = (d, i, orientation_rank, i + L, orientation)
                if best is None or candidate[:3] < best[:3]:
                    best = candidate
                    if d == 0 and i == 0 and orientation_rank == 0:
                        break
    if best is None:
        return None
    _distance, start, _orientation_rank, end, orientation = best
    return start, end, orientation


def find_sequence(seq: str, target: str, mismatches: int = 0, check_revcomp: bool = False) -> bool:
    return find_sequence_match(seq, target, mismatches, check_revcomp) is not None


def extract_between_flanks_with_coords(
    seq: str,
    left: str,
    right: str,
    allowed_lengths: Optional[set[int]] = None,
    min_len: int = 0,
    max_len: int = 0,
) -> Optional[tuple[str, int, int]]:
    """Return a flank-delimited target and coordinates.

    Flanks may be any non-empty length. When target-length constraints are
    known, all plausible left/right occurrences are searched until the first
    pair yielding an allowed target length is found. This avoids a short right
    flank that happens to occur inside the true target from prematurely ending
    extraction. With no length constraints, behavior remains the historical
    first-left / first-subsequent-right match.
    """
    seq_u = seq.upper()
    left_u = left.upper()
    right_u = right.upper()
    if not left_u or not right_u:
        return None

    p1 = seq_u.find(left_u)
    while p1 >= 0:
        start = p1 + len(left_u)
        p2 = seq_u.find(right_u, start)
        while p2 >= 0:
            length = p2 - start
            valid = length > 0
            if valid and allowed_lengths:
                valid = length in allowed_lengths
            if valid and min_len:
                valid = length >= min_len
            if valid and max_len:
                valid = length <= max_len
            if valid:
                return seq_u[start:p2], start, p2
            p2 = seq_u.find(right_u, p2 + 1)
        p1 = seq_u.find(left_u, p1 + 1)
    return None


def extract_between_flanks(seq: str, left: str, right: str) -> Optional[str]:
    hit = extract_between_flanks_with_coords(seq, left, right)
    return hit[0] if hit is not None else None


def normalize_cell_id(raw_cell: str, round1_conversion: Optional[dict[str, str]], use_conversion: bool) -> Optional[str]:
    if not raw_cell:
        return None
    cell = raw_cell[:24]
    if use_conversion and round1_conversion:
        conv = convert_cell_round1(cell, round1_conversion)
        return conv
    return cell


def write_legacy(counts: dict[tuple[str, str], int], out_dir: str | Path, include_zero_cells: bool = False, all_cells: Iterable[str] | None = None) -> None:
    # Legacy output is a complete per-run directory. Clear a previous directory
    # before rewriting so reruns cannot leave stale per-cell files behind.
    out_dir = Path(out_dir)
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir = ensure_dir(out_dir)
    by_cell: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for (cell, label), c in counts.items():
        if c > 0 or include_zero_cells:
            by_cell[cell].append((label, c))
    if include_zero_cells and all_cells:
        for cell in all_cells:
            by_cell.setdefault(cell, [])
    for cell, vals in by_cell.items():
        with open(out_dir / cell, "w") as out:
            out.write(f"AAsequence,{cell}\n")
            for label, c in sorted(vals):
                out.write(f"{label},{c}\n")


def write_long(counts: dict[tuple[str, str], int], path: str | Path) -> None:
    with open(path, "w") as out:
        out.write("cellID\tsequence\tcount\n")
        for (cell, seq), c in sorted(counts.items()):
            out.write(f"{cell}\t{seq}\t{c}\n")


def write_wide(counts: dict[tuple[str, str], int], path: str | Path, cells: Iterable[str] | None = None, labels: Iterable[str] | None = None) -> None:
    labels = sorted(set(labels or [k[1] for k in counts]))
    cells = sorted(set(cells or [k[0] for k in counts]))
    with open(path, "w") as out:
        header = ["cellID"] + [f"{l}_count" for l in labels]
        out.write("\t".join(header) + "\n")
        for cell in cells:
            out.write(cell)
            for label in labels:
                out.write(f"\t{counts.get((cell, label), 0)}")
            out.write("\n")


def count_single(args) -> None:
    out_dir = ensure_dir(args.output_dir)
    keep = load_keep_cells(args.cell_keep_list) if args.cell_filter == "keep-list" else None
    r1conv = make_round1_conversion(args.mrna_round1, args.capture_round1) if args.round1_conversion else None
    counts: dict[tuple[str, str], int] = Counter()
    umi_sets: dict[tuple[str, str], set[str]] = defaultdict(set)
    stats = Counter()
    observed_cells: set[str] = set()
    target = args.target_sequence.upper()
    quality_threshold = int(getattr(args, "quality_threshold", 30))
    phred_offset = int(getattr(args, "phred_offset", 33))
    target_quality = QualityAccumulator(quality_threshold, phred_offset)

    for rec in iter_fastq(args.input_fastq):
        stats["reads_scanned"] += 1
        raw_cell, umi = extract_splitseq_cell_umi_from_header(rec.name)
        if not raw_cell:
            stats["unparsed_header"] += 1
            continue
        cell = normalize_cell_id(raw_cell, r1conv, bool(args.round1_conversion))
        if not cell:
            stats["round1_conversion_failed"] += 1
            continue
        if keep is not None and cell not in keep:
            stats["not_in_keep_list"] += 1
            continue
        stats["reads_in_output_cells"] += 1
        observed_cells.add(cell)
        search_seq = rec.seq
        search_offset = 0
        if getattr(args, "use_flanks", False):
            extracted = extract_between_flanks_with_coords(rec.seq, args.flank_left, args.flank_right)
            if extracted is None:
                stats["missing_flanks"] += 1
                continue
            search_seq, search_offset, _search_end = extracted
        match = find_sequence_match(search_seq, target, args.capture_mismatches, args.check_revcomp)
        if match is None:
            continue
        match_start, match_end, _orientation = match
        abs_start = search_offset + match_start
        abs_end = search_offset + match_end
        if len(rec.qual) >= abs_end:
            target_quality.add(rec.qual[abs_start:abs_end])
        else:
            stats["target_quality_too_short"] += 1
        stats["target_positive_reads"] += 1
        key = (cell, target)
        if args.count_mode == "umi":
            if not umi:
                stats["missing_umi"] += 1
                continue
            umi_sets[key].add(umi)
        else:
            counts[key] += 1

    if args.count_mode == "umi":
        counts = Counter({k: len(v) for k, v in umi_sets.items()})
    all_cells = keep if keep is not None else sorted(observed_cells)
    prefix = out_dir / args.sample_id
    write_wide(counts, str(prefix) + ".single_capture_counts.wide.tsv", cells=all_cells, labels=[target])
    write_long(counts, str(prefix) + ".single_capture_counts.long.tsv")
    if args.legacy_output:
        write_legacy(counts, out_dir / f"output_PepVarCounts_{args.sample_id}_out", args.include_zero_cells, all_cells)
    write_quality_qc(
        str(prefix) + ".capture_quality_qc.tsv",
        [target_quality.row("single_capture_target", "target_positive_reads_in_output_cells", "matched_in_R1", len(target))],
    )
    write_kv(
        str(prefix) + ".single_capture_summary.tsv",
        list(stats.items()) +
        [("positive_cells", len({c for (c, _), v in counts.items() if v > 0})),
         ("quality_threshold", quality_threshold),
         ("phred_offset", phred_offset)] +
        summary_rows("target", target_quality),
    )


def count_multi(args) -> None:
    """Count flank-delimited multiplex capture targets in one of three modes.

    Modes
    -----
    conversion (default)
        Correct extracted barcode against the capture-barcode keys of a conversion
        table, then optionally output the barcode, mapped sequence, or translated amino-acid sequence.
    greenlist
        Correct extracted barcode against a user-approved list and count the
        approved barcode itself. No sequence mapping or translation is performed.
    unrestricted
        Count every valid extracted A/C/G/T barcode sequence exactly as observed.
        No whitelist correction, sequence mapping, or translation is performed.
    """
    out_dir = ensure_dir(args.output_dir)
    keep = load_keep_cells(args.cell_keep_list) if args.cell_filter == "keep-list" else None
    r1conv = make_round1_conversion(args.mrna_round1, args.capture_round1) if args.round1_conversion else None

    mode = getattr(args, "capture_mode", "conversion") or "conversion"
    target_length = int(getattr(args, "target_length", 0) or 0)

    conv: dict[str, str] = {}
    green: list[str] = []
    corrector = None

    if mode == "conversion":
        conv = load_conversion_table(args.barcode_conversion_table)
        if not conv:
            raise ValueError("No valid entries were loaded from --barcode-conversion-table")
        _validate_reference_lengths(conv.keys(), target_length, "Barcode conversion table")
        if args.translate_aa:
            bad_lengths = sorted({len(mapped) for mapped in conv.values() if len(mapped) % 3 != 0})
            if bad_lengths:
                raise ValueError(
                    "--translate-aa requires every mapped sequence length to be divisible by 3; "
                    f"incompatible lengths observed: {bad_lengths}"
                )
        corrector = BarcodeCorrector(conv.keys(), args.capture_mismatches)
        reference_n = len(conv)
        reference_lengths = {len(x) for x in conv}
        feature_type = "amino_acid" if args.translate_aa else ("capture_barcode" if args.output_sequence in {"barcode", "rna"} else "mapped_sequence")
        source_desc = f"{reference_n:,} barcode conversion entries"
    elif mode == "greenlist":
        green = load_green_list(args.green_list)
        if not green:
            raise ValueError("No valid nucleotide barcodes were loaded from --green-list")
        _validate_reference_lengths(green, target_length, "Green list")
        corrector = BarcodeCorrector(green, args.capture_mismatches)
        reference_n = len(green)
        reference_lengths = {len(x) for x in green}
        feature_type = "approved_barcode"
        source_desc = f"{reference_n:,} green-listed barcodes"
    elif mode == "unrestricted":
        reference_n = 0
        reference_lengths = set()
        feature_type = "observed_barcode"
        source_desc = "no barcode whitelist (all valid extracted targets are counted)"
    else:
        raise ValueError(f"Unknown capture mode: {mode}")

    print(
        f"capture-multi: mode={mode}; {source_desc}; "
        f"mismatches={args.capture_mismatches if mode != 'unrestricted' else 'not used'}; "
        f"count_mode={args.count_mode}; target_length={target_length or 'variable'}; "
        f"flank_lengths={len(args.flank_left)}/{len(args.flank_right)}",
        flush=True,
    )
    if keep is not None:
        print(f"capture-multi: loaded {len(keep):,} kept mRNA cell IDs", flush=True)
    if corrector is not None:
        print(
            f"capture-multi: correction lookup contains {len(corrector.lookup):,} entries",
            flush=True,
        )

    counts: dict[tuple[str, str], int] = Counter()
    umi_sets: dict[tuple[str, str], set[str]] = defaultdict(set)
    stats = Counter()
    observed_cells: set[str] = set()
    quality_threshold = int(getattr(args, "quality_threshold", 30))
    phred_offset = int(getattr(args, "phred_offset", 33))
    extracted_quality = QualityAccumulator(quality_threshold, phred_offset)
    accepted_quality = QualityAccumulator(quality_threshold, phred_offset)
    started = time.monotonic()
    raw_progress_every = getattr(args, "progress_every", 1_000_000)
    progress_every = 1_000_000 if raw_progress_every is None else int(raw_progress_every)

    for rec in iter_fastq(args.input_fastq):
        stats["reads_scanned"] += 1
        if progress_every > 0 and stats["reads_scanned"] % progress_every == 0:
            elapsed = max(time.monotonic() - started, 1e-9)
            rate = stats["reads_scanned"] / elapsed
            print(
                "capture-multi progress: "
                f"reads={stats['reads_scanned']:,}; "
                f"kept-cell reads={stats['reads_in_output_cells']:,}; "
                f"accepted={stats['accepted_reads']:,}; "
                f"missing-flanks={stats['missing_flanks']:,}; wrong-length={stats['target_wrong_length']:,}; "
                f"barcode-miss={stats['barcode_not_in_whitelist']:,}; "
                f"rate={rate:,.0f} reads/s",
                flush=True,
            )

        raw_cell, umi = extract_splitseq_cell_umi_from_header(rec.name)
        if not raw_cell:
            stats["unparsed_header"] += 1
            continue
        cell = normalize_cell_id(raw_cell, r1conv, bool(args.round1_conversion))
        if not cell:
            stats["round1_conversion_failed"] += 1
            continue
        if keep is not None and cell not in keep:
            stats["not_in_keep_list"] += 1
            continue
        stats["reads_in_output_cells"] += 1
        observed_cells.add(cell)

        allowed_lengths = {target_length} if target_length else (reference_lengths if (not args.min_insert_len and not args.max_insert_len and reference_lengths) else None)
        extracted = extract_between_flanks_with_coords(
            rec.seq,
            args.flank_left,
            args.flank_right,
            allowed_lengths=allowed_lengths,
            min_len=0 if target_length else args.min_insert_len,
            max_len=0 if target_length else args.max_insert_len,
        )
        if extracted is None:
            unconstrained = extract_between_flanks_with_coords(rec.seq, args.flank_left, args.flank_right)
            if unconstrained is None:
                stats["missing_flanks"] += 1
            else:
                stats["target_wrong_length"] += 1
            continue
        bc, bc_start, bc_end = extracted

        # Exact user-defined target length takes precedence over min/max range.
        if target_length:
            if len(bc) != target_length:
                stats["target_wrong_length"] += 1
                continue
        else:
            if args.min_insert_len and len(bc) < args.min_insert_len:
                stats["insert_too_short"] += 1
                continue
            if args.max_insert_len and len(bc) > args.max_insert_len:
                stats["insert_too_long"] += 1
                continue

        if "N" in bc:
            stats["insert_contains_N"] += 1
            continue

        target_qual = rec.qual[bc_start:bc_end] if len(rec.qual) >= bc_end else None
        if target_qual is not None:
            extracted_quality.add(target_qual)
        else:
            stats["target_quality_too_short"] += 1

        if mode == "unrestricted":
            corrected_bc = bc
        else:
            res = corrector.correct(bc, args.ambiguous_policy)
            if not res.corrected:
                stats["barcode_not_in_whitelist"] += 1
                if res.ambiguous:
                    stats["ambiguous_barcode"] += 1
                continue
            corrected_bc = res.corrected
            if res.distance:
                stats["barcode_corrected_reads"] += 1
            else:
                stats["barcode_exact_reads"] += 1

        if target_qual is not None:
            accepted_quality.add(target_qual)

        if mode == "conversion":
            dna = conv[corrected_bc]
            label = translate_dna(dna) if args.translate_aa else (corrected_bc if args.output_sequence in {"barcode", "rna"} else dna)
        else:
            # In barcode-only modes, the counted feature is the capture barcode itself.
            label = corrected_bc

        key = (cell, label)
        stats["accepted_reads"] += 1
        if args.count_mode == "umi":
            if not umi:
                stats["missing_umi"] += 1
                continue
            umi_sets[key].add(umi)
        else:
            counts[key] += 1

    elapsed = max(time.monotonic() - started, 1e-9)
    print(
        "capture-multi scan complete: "
        f"{stats['reads_scanned']:,} reads in {elapsed/60:.1f} min "
        f"({stats['reads_scanned']/elapsed:,.0f} reads/s); "
        f"accepted={stats['accepted_reads']:,}; unique_features={len({k[1] for k in counts}) if args.count_mode != 'umi' else len({k[1] for k in umi_sets}):,}",
        flush=True,
    )

    if args.count_mode == "umi":
        counts = Counter({k: len(v) for k, v in umi_sets.items()})

    all_cells = keep if keep is not None else sorted(observed_cells)
    labels = sorted({k[1] for k in counts})
    prefix = out_dir / args.sample_id
    write_long(counts, str(prefix) + ".multi_capture_counts.long.tsv")
    write_wide(counts, str(prefix) + ".multi_capture_counts.wide.tsv", cells=all_cells, labels=labels)
    if args.legacy_output:
        # Retain the historical file structure/header for downstream compatibility,
        # even in barcode-only modes where the row labels are barcodes rather than AA.
        write_legacy(counts, out_dir / f"output_PepVarCounts_{args.sample_id}_out", args.include_zero_cells, all_cells)

    expected_len = target_length or (args.min_insert_len if args.min_insert_len and args.min_insert_len == args.max_insert_len else "")
    write_quality_qc(
        str(prefix) + ".capture_quality_qc.tsv",
        [
            extracted_quality.row("multiplex_extracted_target", "valid_flank_extracted_targets_before_mode_filtering", "between_flanks", expected_len),
            accepted_quality.row("multiplex_accepted_target", "accepted_capture_reads", "between_flanks", expected_len),
        ],
    )
    write_kv(
        str(prefix) + ".multi_capture_summary.tsv",
        [
            ("capture_mode", mode),
            ("output_feature_type", feature_type),
            ("reference_barcode_count", reference_n),
            ("target_length", target_length if target_length else "variable/range"),
            ("flank_left_length", len(args.flank_left)),
            ("flank_right_length", len(args.flank_right)),
        ] +
        list(stats.items()) +
        [("positive_cells", len({c for (c, _), v in counts.items() if v > 0})),
         ("unique_counted_features", len(labels)),
         ("quality_threshold", quality_threshold),
         ("phred_offset", phred_offset)] +
        summary_rows("extracted_target", extracted_quality) +
        summary_rows("target", accepted_quality),
    )

