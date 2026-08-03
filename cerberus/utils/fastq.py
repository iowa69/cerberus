"""FASTQ helpers: fast read counting, pair-sync verification.

Compression is detected from the file's magic bytes rather than its
extension, so a gzipped file named ``reads.fastq`` (or plain text named
``.gz``) is counted correctly instead of silently returning nonsense.
"""
from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

from cerberus.utils.logger import get_logger
from cerberus.utils.shell import which

log = get_logger("fastq")


def is_gzipped(path: Path) -> bool:
    try:
        with path.open("rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def count_reads(path: Path) -> int:
    """Count records in a FASTQ(.gz) file.

    Uses pigz/zcat piped to ``wc -l`` when available. Returns 0 for a missing
    or empty file. A decompression failure is reported rather than silently
    yielding a plausible-looking number.
    """
    if not path.exists() or path.stat().st_size == 0:
        return 0

    lines: int | None = None
    if is_gzipped(path):
        for tool in (["pigz", "-dc"], ["zcat"], ["gzip", "-dc"]):
            if which(tool[0]):
                lines = _line_count_pipe([*tool, str(path)])
                break
    elif which("wc"):
        try:
            out = subprocess.run(
                ["wc", "-l", str(path)], capture_output=True, check=True,
            ).stdout.decode().split()[0]
            lines = int(out)
        except (subprocess.CalledProcessError, ValueError, IndexError) as e:
            log.warning("wc -l failed on %s (%s); counting in Python", path, e)

    if lines is None:
        lines = _count_python(path)

    if lines % 4:
        log.warning(
            "%s holds %d lines, which is not a multiple of 4 — the file looks truncated; "
            "reporting %d complete records", path.name, lines, lines // 4,
        )
    return lines // 4


def _line_count_pipe(cmd: list[str]) -> int | None:
    """Count lines from a decompressor. Returns None when it fails."""
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
    if p1.stdout:
        p1.stdout.close()
    out, _ = p2.communicate()
    p1.wait()
    err = p1.stderr.read().decode(errors="replace") if p1.stderr else ""
    if p1.stderr:
        p1.stderr.close()
    if p1.returncode != 0:
        log.warning("%s failed (rc=%d): %s", cmd[0], p1.returncode, err.strip())
        return None
    try:
        return int(out.decode().strip())
    except ValueError:
        return None


def _count_python(path: Path) -> int:
    opener = gzip.open if is_gzipped(path) else open
    n = 0
    with opener(path, "rb") as f:                        # type: ignore[operator]
        for _ in f:
            n += 1
    return n


def detect_paired_orphans(r1: Path, r2: Path) -> tuple[int, int]:
    """Return (r1_count, r2_count). If equal, no orphans."""
    return count_reads(r1), count_reads(r2)


def pairs_in_sync(r1: Path, r2: Path) -> bool:
    """True when R1 and R2 hold the same number of records."""
    a, b = detect_paired_orphans(r1, r2)
    return a == b


def output_name(out_dir: Path, sample_id: str, mode: str, suffix: str) -> Path:
    """Canonical output filename layout."""
    return out_dir / f"{sample_id}.{mode}.{suffix}"
