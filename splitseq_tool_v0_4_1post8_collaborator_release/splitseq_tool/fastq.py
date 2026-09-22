from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional
import re

from .utils import open_text


@dataclass
class FastqRecord:
    name: str
    seq: str
    plus: str
    qual: str


def iter_fastq(path: str | Path) -> Iterator[FastqRecord]:
    """Yield validated FASTQ records from plain or gzipped input.

    Truncated/malformed FASTQ files raise an error instead of silently ending a
    run. Sequence and quality lengths must match for every record.
    """
    with open_text(path, "rt") as f:
        record_number = 0
        while True:
            name = f.readline()
            if not name:
                break
            record_number += 1
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            if not seq or not plus or not qual:
                raise ValueError(f"Truncated FASTQ record {record_number} in {path}")
            name = name.rstrip("\r\n")
            seq = seq.rstrip("\r\n")
            plus = plus.rstrip("\r\n")
            qual = qual.rstrip("\r\n")
            if not name.startswith("@"):
                raise ValueError(f"FASTQ record {record_number} in {path} does not start with '@'")
            if not plus.startswith("+"):
                raise ValueError(f"FASTQ record {record_number} in {path} has an invalid '+' line")
            if len(seq) != len(qual):
                raise ValueError(
                    f"FASTQ record {record_number} in {path} has sequence/quality length mismatch "
                    f"({len(seq)} != {len(qual)})"
                )
            yield FastqRecord(name, seq, plus, qual)


def extract_splitseq_cell_umi_from_header(header: str) -> tuple[Optional[str], Optional[str]]:
    """Extract corrected 24-nt cell ID and UMI from a MergedCells FASTQ header.

    Supports the rebuilt tool header form ``...:CB_<24nt>_UMI_<umi>`` anywhere
    in the header, plus the legacy field-9 layout.
    """
    # UMI-tools compatible rebuilt-tool form: ..._<24nt cell>_<UMI>
    m = re.search(r"_([ACGTNacgtn]{24})_([ACGTNacgtn]+)(?:\s|$)", header)
    if m:
        return m.group(1).upper(), m.group(2).upper()

    m = re.search(r"(?:^|[:_])CB_([ACGTNacgtn]{24})_(?:UMI|UB)_([ACGTNacgtn]+)", header)
    if m:
        return m.group(1).upper(), m.group(2).upper()
    m = re.search(r"CB_([ACGTNacgtn]{24})", header)
    if m:
        cell = m.group(1).upper()
        um = re.search(r"(?:UMI|UB)_([ACGTNacgtn]+)", header)
        return cell, um.group(1).upper() if um else None

    parts = header.split(":")
    if len(parts) > 9:
        token = parts[9].strip()
        bits = token.split("_")
        if len(bits) >= 2 and len(bits[1]) >= 24:
            cell = bits[1][:24].upper()
            umi = None
            for i, b in enumerate(bits):
                if b.upper() in {"UMI", "UB"} and i + 1 < len(bits):
                    umi = bits[i + 1].upper()
                    break
            return cell, umi

    clean = header.replace("@", " ").replace(":", " ").replace("_", " ")
    for token in clean.split():
        if len(token) == 24 and all(c in "ACGTNacgtn" for c in token):
            return token.upper(), None
    return None, None
