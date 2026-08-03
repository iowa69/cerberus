"""Stream-concatenate FASTQ(.gz) files.

Concatenated gzip members form a single valid gzip stream, so joining
``.fq.gz`` files byte-wise produces a readable ``.fq.gz``. Every consumer in
the pipeline (zcat, pigz, kraken2, bbduk, samtools, Biopython) handles
multi-member gzip.
"""
from __future__ import annotations

import gzip
from pathlib import Path

from cerberus.config import CerberusConfig
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import pipe, run, which

log = get_logger("concat")


def _write_empty_gzip(path: Path) -> None:
    """An empty *valid* gzip member, not a zero-byte file.

    A zero-byte ``.gz`` fails ``gzip -t`` and makes downstream tools error
    out on what is really a legitimate "nothing survived" result.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as f:
        f.write(b"")


def concat_gz(
    cfg: CerberusConfig,
    *,
    inputs: list[Path],
    output: Path,
    log_dir: Path,
    tag: str,
) -> Path:
    inputs = [p for p in inputs if p is not None and p.exists()]
    output.parent.mkdir(parents=True, exist_ok=True)

    if not inputs:
        log.warning("No inputs to concatenate for %s; writing an empty gzip", tag)
        if not cfg.dry_run:
            _write_empty_gzip(output)
        return output

    log.info("Concatenating %d files -> %s", len(inputs), output.name)

    if cfg.dry_run:
        log_dir.mkdir(parents=True, exist_ok=True)
        (log_dir / f"{tag}.concat.log").write_text(
            "DRY-RUN concat:\n" + "\n".join(str(p) for p in inputs)
        )
        return output

    with output.open("wb") as out:
        for p in inputs:
            with p.open("rb") as fin:
                while chunk := fin.read(4 * 1024 * 1024):
                    out.write(chunk)

    log.info("Concatenated %d files (%d bytes)", len(inputs), output.stat().st_size)
    return output


def concat_mate_tagged(
    cfg: CerberusConfig,
    *,
    inputs: list[tuple[Path, str]],
    output: Path,
    log_dir: Path,
    tag: str,
) -> Path:
    """Merge FASTQ files, appending a mate suffix to each source's read names.

    ``inputs`` is a list of ``(path, suffix)``; a suffix of ``""`` leaves the
    names alone. Merging R1 and R2 into one file without this gives both mates
    of a fragment the same read ID, and any de-duplication by name then throws
    away half the data.

    Each stream is rewritten with awk and re-compressed, then appended as its
    own gzip member — concatenated members read back as one stream.
    """
    inputs = [(p, sfx) for p, sfx in inputs if p is not None and p.exists()]
    output.parent.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    if not inputs:
        log.warning("No inputs to merge for %s; writing an empty gzip", tag)
        if not cfg.dry_run:
            _write_empty_gzip(output)
        return output

    if cfg.dry_run:
        (log_dir / f"{tag}.concat.log").write_text(
            "DRY-RUN merge:\n" + "\n".join(f"{p} (+{s!r})" for p, s in inputs)
        )
        return output

    decomp = ["pigz", "-dc"] if which("pigz") else ["gzip", "-dc"]
    comp = ["pigz", "-c", "-p", str(cfg.threads)] if which("pigz") else ["gzip", "-c"]

    staged: list[Path] = []
    for i, (path, suffix) in enumerate(inputs):
        if not suffix:
            staged.append(path)
            continue
        part = output.with_name(f"{output.name}.part{i}")
        # Strip any existing comment after whitespace, then append the suffix.
        awk = r'NR%4==1 {sub(/[ \t].*/, "", $0); print $0 "' + suffix + r'"; next} {print}'
        pipe(
            [[*decomp, str(path)], ["awk", awk], comp],
            log_path=log_dir / f"{tag}.tag{i}.log",
            final_stdout=part,
        )
        staged.append(part)

    with output.open("wb") as out:
        for p in staged:
            with p.open("rb") as fin:
                while chunk := fin.read(4 * 1024 * 1024):
                    out.write(chunk)
    for p in staged:
        if p.name.startswith(output.name + ".part"):
            p.unlink(missing_ok=True)

    log.info("Merged %d files with mate tags -> %s (%d bytes)",
             len(inputs), output.name, output.stat().st_size)
    return output


def compress_to(cfg: CerberusConfig, src: Path, dst: Path, log_dir: Path, tag: str) -> Path:
    """Compress ``src`` to ``dst``.

    pigz writes to stdout with ``-c``, so the destination has to be captured
    explicitly — routing it through the step log (as an earlier version did)
    put the compressed payload in the log file and never created ``dst``.
    """
    if dst.suffix != ".gz":
        raise ValueError("compress_to: destination must end in .gz")
    dst.parent.mkdir(parents=True, exist_ok=True)
    if cfg.dry_run:
        return dst
    if which("pigz"):
        run(["pigz", "-c", "-p", str(cfg.threads), str(src)],
            log_path=log_dir / f"{tag}.pigz.log", stdout_path=dst)
    else:
        with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
            for chunk in iter(lambda: fin.read(1024 * 1024), b""):
                fout.write(chunk)
    return dst
