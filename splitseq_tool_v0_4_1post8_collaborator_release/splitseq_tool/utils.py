from __future__ import annotations

import gzip
import os
import subprocess
import shutil
from pathlib import Path
from typing import Iterable, Iterator, Optional, Tuple


def open_text(path: str | Path, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def run_cmd(cmd: list[str], cwd: str | Path | None = None, log_path: str | Path | None = None) -> None:
    """Run a command, streaming stdout/stderr into optional log file."""
    msg = " ".join(cmd)
    if log_path:
        with open(log_path, "a") as log:
            log.write(f"\n$ {msg}\n")
            log.flush()
            subprocess.run(cmd, cwd=cwd, check=True, stdout=log, stderr=subprocess.STDOUT)
    else:
        print(f"$ {msg}", flush=True)
        subprocess.run(cmd, cwd=cwd, check=True)


def write_kv(path: str | Path, rows: Iterable[tuple[str, object]]) -> None:
    with open(path, "w") as out:
        for k, v in rows:
            out.write(f"{k}\t{v}\n")


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(x != y for x, y in zip(a, b))


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]



def path_size_bytes(path: str | Path) -> int:
    """Return total size of a file or directory tree, best-effort."""
    p = Path(path)
    try:
        if p.is_file() or p.is_symlink():
            return p.stat().st_size
        if p.is_dir():
            total = 0
            for child in p.rglob("*"):
                try:
                    if child.is_file() or child.is_symlink():
                        total += child.stat().st_size
                except OSError:
                    pass
            return total
    except OSError:
        pass
    return 0


def cleanup_intermediates(
    paths: Iterable[str | Path],
    summary_path: str | Path,
    keep_all: bool = False,
    reason: str = "intermediate after successful command",
) -> int:
    """Delete intermediate paths after successful completion and log every action.

    Cleanup is intentionally called only *after* a command has produced and
    validated its final outputs. If ``keep_all`` is true, paths are retained and
    recorded as such. Missing paths are also recorded. Returns bytes deleted.
    """
    rows = []
    deleted_bytes = 0
    for raw in paths:
        p = Path(raw)
        size = path_size_bytes(p) if p.exists() or p.is_symlink() else 0
        if not (p.exists() or p.is_symlink()):
            rows.append((str(p), "missing", 0, reason))
            continue
        if keep_all:
            rows.append((str(p), "kept", size, "--keep-all-intermediates"))
            continue
        try:
            if p.is_dir() and not p.is_symlink():
                shutil.rmtree(p)
            else:
                p.unlink()
            deleted_bytes += size
            rows.append((str(p), "deleted", size, reason))
        except OSError as exc:
            rows.append((str(p), "delete_failed", size, f"{reason}: {exc}"))

    sp = Path(summary_path)
    sp.parent.mkdir(parents=True, exist_ok=True)
    with open(sp, "w") as out:
        out.write("path\taction\tbytes\treason\n")
        for path, action, size, why in rows:
            out.write(f"{path}\t{action}\t{size}\t{why}\n")
        out.write(f"TOTAL_DELETED_BYTES\tdeleted_total\t{deleted_bytes}\t\n")
    return deleted_bytes


def require_output_files(paths: Iterable[str | Path], allow_empty: bool = True) -> None:
    """Verify expected final files exist before deleting intermediates."""
    missing = []
    empty = []
    for raw in paths:
        p = Path(raw)
        if not p.exists():
            missing.append(str(p))
        elif not allow_empty and p.is_file() and p.stat().st_size == 0:
            empty.append(str(p))
    if missing or empty:
        pieces = []
        if missing:
            pieces.append("missing: " + ", ".join(missing))
        if empty:
            pieces.append("empty: " + ", ".join(empty))
        raise RuntimeError(
            "Refusing intermediate-file cleanup because expected final outputs were not valid ("
            + "; ".join(pieces)
            + ")"
        )
