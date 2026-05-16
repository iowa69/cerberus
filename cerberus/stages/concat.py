"""Stream-concatenate FASTQ(.gz) files using cat (gzip files are concatenatable as gzip streams)."""
from __future__ import annotations

from pathlib import Path

from cerberus.config import CerberusConfig
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import run

log = get_logger("concat")


def concat_gz(
    cfg: CerberusConfig,
    *,
    inputs: list[Path],
    output: Path,
    log_dir: Path,
    tag: str,
) -> Path:
    inputs = [p for p in inputs if p is not None and p.exists()]
    if not inputs:
        log.warning("No inputs to concatenate for %s; writing empty output", tag)
        if not cfg.dry_run:
            output.write_bytes(b"")
        return output

    output.parent.mkdir(parents=True, exist_ok=True)
    log.info("Concatenating %d files → %s", len(inputs), output.name)

    if cfg.dry_run:
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


def compress_to(cfg: CerberusConfig, src: Path, dst: Path, log_dir: Path, tag: str) -> Path:
    from cerberus.utils.shell import which
    if dst.suffix != ".gz":
        raise ValueError("compress_to: destination must end in .gz")
    if which("pigz"):
        run(["pigz", "-c", "-p", str(cfg.threads), str(src)],
            log_path=log_dir / f"{tag}.pigz.log", dry_run=cfg.dry_run)
    else:
        import gzip
        with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
            for chunk in iter(lambda: fin.read(1024 * 1024), b""):
                fout.write(chunk)
    return dst
