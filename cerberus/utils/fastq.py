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


# Counting is pure I/O over files that do not change once written, and the
# accounting, the pipelines and the report all ask about the same paths. Cache
# on identity (size, mtime) so a run does not decompress the same file four
# times; any rewrite invalidates the entry.
_COUNT_CACHE: dict[tuple[str, int, int], int] = {}


def clear_count_cache() -> None:
    _COUNT_CACHE.clear()


def count_reads(path: Path) -> int:
    """Count records in a FASTQ(.gz) file.

    Uses pigz/zcat piped to ``wc -l`` when available. Returns 0 for a missing
    or empty file. A decompression failure is reported rather than silently
    yielding a plausible-looking number.
    """
    try:
        st = path.stat()
    except OSError:
        return 0
    if st.st_size == 0:
        return 0

    key = (str(path), st.st_size, int(st.st_mtime_ns))
    cached = _COUNT_CACHE.get(key)
    if cached is not None:
        return cached

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
        if lines is None:
            log.warning("Could not count records in %s; reporting 0", path)
            return 0

    if lines % 4:
        log.warning(
            "%s holds %d lines, which is not a multiple of 4 — the file looks truncated; "
            "reporting %d complete records", path.name, lines, lines // 4,
        )
    result = lines // 4
    _COUNT_CACHE[key] = result
    return result


def _line_count_pipe(cmd: list[str]) -> int | None:
    """Count lines from a decompressor. Returns None when it fails."""
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
    if p1.stdout:
        p1.stdout.close()
    out, _ = p2.communicate()
    # Drain stderr BEFORE waiting. A decompressor that writes more than a pipe
    # buffer of warnings would otherwise block on write() while we block in
    # wait(), and the run would hang forever.
    err = b""
    if p1.stderr:
        try:
            err = p1.stderr.read()
        except OSError:
            pass
        finally:
            p1.stderr.close()
    p1.wait()
    if p1.returncode != 0:
        log.warning("%s failed (rc=%d): %s", cmd[0],
                    p1.returncode, err.decode(errors="replace").strip())
        return None
    try:
        return int(out.decode().strip())
    except ValueError:
        return None


def _count_python(path: Path) -> int | None:
    """Pure-Python fallback. Returns None if the file cannot be read.

    A truncated gzip raises EOFError, which is *not* an OSError — letting it
    escape would kill a finished run from inside a diagnostics helper.
    """
    opener = gzip.open if is_gzipped(path) else open
    n = 0
    try:
        with opener(path, "rb") as f:                    # type: ignore[operator]
            for _ in f:
                n += 1
    except (OSError, EOFError, gzip.BadGzipFile) as e:
        log.warning("%s could not be read (%s); it may be truncated or corrupt", path.name, e)
        return None
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
