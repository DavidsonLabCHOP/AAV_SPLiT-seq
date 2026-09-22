#!/usr/bin/env python3
"""Standalone entry point for the SPLiT-seq QC report.

This lets users run QC without installing the full package if they launch it from
inside the repository root:
    python scripts/splitseq_qc_report.py --sample-id S1 ...
"""
import sys
from pathlib import Path

repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo))

from splitseq_tool.cli import main  # noqa: E402

if __name__ == '__main__':
    main(['qc-report'] + sys.argv[1:])
