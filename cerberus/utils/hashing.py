"""SHA256 file hashing with progress."""
from __future__ import annotations

import hashlib
from pathlib import Path

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None


_CHUNK = 1024 * 1024


def sha256_file(path: Path, *, show_progress: bool = False) -> str:
    h = hashlib.sha256()
    size = path.stat().st_size
    bar = tqdm(total=size, unit="B", unit_scale=True, desc=f"sha256 {path.name}") \
        if (show_progress and tqdm and size > 50 * 1024 * 1024) else None
    try:
        with path.open("rb") as f:
            while chunk := f.read(_CHUNK):
                h.update(chunk)
                if bar:
                    bar.update(len(chunk))
    finally:
        if bar:
            bar.close()
    return h.hexdigest()


def verify_sha256(path: Path, expected: str, *, show_progress: bool = False) -> bool:
    actual = sha256_file(path, show_progress=show_progress)
    return actual.lower() == expected.lower()
