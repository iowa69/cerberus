"""FASTQ helpers: fast read counting, pair sync verification."""
from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

from cerberus.utils.shell import which


def count_reads(path: Path) -> int:
    """Count records in a FASTQ(.gz) file. Uses pigz/zcat | wc -l if available."""
    if not path.exists():
        return 0
    if path.stat().st_size == 0:
        return 0

    is_gz = path.suffix == ".gz"
    if is_gz and which("pigz"):
        return _line_count_pipe(["pigz", "-dc", str(path)]) // 4
    if is_gz and which("zcat"):
        return _line_count_pipe(["zcat", str(path)]) // 4
    if not is_gz and which("wc"):
        out = subprocess.check_output(["wc", "-l", str(path)]).decode().split()[0]
        return int(out) // 4
    return _count_python(path, is_gz)


def _line_count_pipe(cmd: list[str]) -> int:
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
    if p1.stdout:
        p1.stdout.close()
    out, _ = p2.communicate()
    p1.wait()
    return int(out.decode().strip())


def _count_python(path: Path, is_gz: bool) -> int:
    opener = gzip.open if is_gz else open
    n = 0
    with opener(path, "rt") as f:  # type: ignore[operator]
        for _ in f:
            n += 1
    return n // 4


def detect_paired_orphans(r1: Path, r2: Path) -> tuple[int, int]:
    """Return (r1_count, r2_count). If equal, no orphans."""
    return count_reads(r1), count_reads(r2)


def is_gzipped(path: Path) -> bool:
    if path.suffix == ".gz":
        return True
    try:
        with path.open("rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def output_name(out_dir: Path, sample_id: str, mode: str, suffix: str) -> Path:
    """Canonical output filename layout."""
    return out_dir / f"{sample_id}.{mode}.{suffix}"
