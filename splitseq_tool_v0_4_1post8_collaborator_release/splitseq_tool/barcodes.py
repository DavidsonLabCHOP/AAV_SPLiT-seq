from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, Optional

from .utils import hamming


def load_barcode_list(path: str | Path, length: int = 8, mode: str = "auto", deduplicate: bool = True) -> list[str]:
    """Load barcode whitelist.

    mode:
      auto: if token length is length, use token; otherwise use first length bases.
      first: always use first `length` bases.
      second: use bases length:2*length, useful for old tag+barcode files.
      raw: use whole first token.
    """
    vals: list[str] = []
    with open(path) as f:
        for line in f:
            token = line.strip().split()[0] if line.strip() else ""
            if not token:
                continue
            token = token.upper()
            if mode == "raw":
                bc = token
            elif mode == "first":
                bc = token[:length]
            elif mode == "second":
                bc = token[length:2 * length]
            else:
                if len(token) == length:
                    bc = token
                else:
                    bc = token[:length]
            if len(bc) == length and all(c in "ACGTN" for c in bc):
                vals.append(bc)
    if not deduplicate:
        return vals
    # preserve order, remove duplicates
    seen = set()
    out = []
    for v in vals:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


@dataclass
class CorrectionResult:
    corrected: Optional[str]
    distance: Optional[int]
    ambiguous: bool = False


class BarcodeCorrector:
    def __init__(self, whitelist: Iterable[str], max_mismatches: int = 1):
        self.whitelist = list(whitelist)
        self.max_mismatches = int(max_mismatches)
        self.lookup: Dict[str, CorrectionResult] = {}
        self._build()

    def _build(self):
        for bc in self.whitelist:
            self.lookup.setdefault(bc, CorrectionResult(bc, 0, False))
        alphabet = "ACGT"
        if self.max_mismatches >= 1:
            for bc in self.whitelist:
                chars = list(bc)
                for i, old in enumerate(chars):
                    for nt in alphabet:
                        if nt == old:
                            continue
                        mut = bc[:i] + nt + bc[i + 1:]
                        prev = self.lookup.get(mut)
                        if prev is None:
                            self.lookup[mut] = CorrectionResult(bc, 1, False)
                        elif prev.distance is None or 1 < prev.distance:
                            self.lookup[mut] = CorrectionResult(bc, 1, False)
                        elif prev.distance == 1 and prev.corrected != bc:
                            # Ambiguity exists only between equally good matches.
                            # An exact whitelist match (distance 0) must always
                            # outrank a one-mismatch alternative.
                            self.lookup[mut] = CorrectionResult(prev.corrected, 1, True)
                        # If prev.distance == 0, keep the exact match unchanged.
        # Note: for >1 mismatch, use scan fallback in correct().

    def correct(self, observed: str, ambiguous_policy: str = "fail") -> CorrectionResult:
        observed = observed.upper()

        # For 0/1 mismatch correction the lookup table is exhaustive:
        #   max_mismatches == 0 -> exact whitelist entries only
        #   max_mismatches == 1 -> exact entries + every 1-nt neighbor
        # Therefore a lookup miss cannot possibly be rescued within the
        # requested distance.  Returning immediately is critical for large
        # capture whitelists; scanning the whole whitelist for every miss can
        # turn a streaming job into O(reads * whitelist_size).
        if self.max_mismatches <= 1:
            res = self.lookup.get(observed)
            if res is None:
                return CorrectionResult(None, None, False)
            if res.ambiguous and ambiguous_policy == "fail":
                return CorrectionResult(None, res.distance, True)
            return res

        # Generic fallback is needed only for mismatch tolerances > 1.
        best: list[tuple[int, str]] = []
        for bc in self.whitelist:
            # Hamming-distance correction is substitution-only; sequences of
            # different lengths are never valid correction candidates.
            if len(bc) != len(observed):
                continue
            d = hamming(observed, bc)
            if d <= self.max_mismatches:
                best.append((d, bc))
        if not best:
            return CorrectionResult(None, None, False)
        best.sort()
        if len(best) > 1 and best[0][0] == best[1][0]:
            if ambiguous_policy == "fail":
                return CorrectionResult(None, best[0][0], True)
        return CorrectionResult(best[0][1], best[0][0], len(best) > 1 and best[0][0] == best[1][0])


def make_round1_conversion(mrna_r1_path: str | Path, capture_r1_path: str | Path) -> dict[str, str]:
    # Positional conversion depends on preserving every input row. Duplicate
    # Round1 entries are therefore errors rather than silently deduplicated.
    mrna = load_barcode_list(mrna_r1_path, 8, "auto", deduplicate=False)
    cap = load_barcode_list(capture_r1_path, 8, "auto", deduplicate=False)
    if len(mrna) != len(cap):
        raise ValueError(
            f"Round1 conversion requires the same number of RNA and capture barcodes; "
            f"got {len(mrna)} and {len(cap)}"
        )
    if len(set(mrna)) != len(mrna):
        raise ValueError("RNA Round1 conversion list contains duplicate barcodes")
    if len(set(cap)) != len(cap):
        raise ValueError("Capture Round1 conversion list contains duplicate barcodes")
    return dict(zip(cap, mrna))


def convert_cell_round1(cell_id: str, conversion: dict[str, str]) -> Optional[str]:
    if len(cell_id) < 24:
        return None
    r1 = cell_id[:8]
    if r1 not in conversion:
        return None
    return conversion[r1] + cell_id[8:24]


def barcode_distance_summary(barcodes: Iterable[str]) -> dict[str, int]:
    bcs = list(barcodes)
    out = {"n_barcodes": len(bcs), "min_distance": 999, "pairs_within_1": 0, "pairs_within_2": 0, "pairs_within_3": 0}
    for a, b in combinations(bcs, 2):
        d = hamming(a, b)
        out["min_distance"] = min(out["min_distance"], d)
        if d <= 1:
            out["pairs_within_1"] += 1
        if d <= 2:
            out["pairs_within_2"] += 1
        if d <= 3:
            out["pairs_within_3"] += 1
    if out["min_distance"] == 999:
        out["min_distance"] = 0
    return out
