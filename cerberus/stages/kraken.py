"""Kraken2 wrapper for the GDPR scrub. Extracts reads NOT classified as host.

Filename handling: kraken2's ``--unclassified-out '<dir>/foo#.fq'`` template
gets ``#`` replaced with ``_1``/``_2`` for paired data. We pass the template
unchanged and discover the actual files via glob — robust to substitution
variants between kraken2 versions.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from cerberus.config import CerberusConfig
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import require_tools, run, which

log = get_logger("kraken")

# Taxa scrubbed by the GDPR pass. Kept for reference; the DB only contains
# mammalian taxa so anything classified will be host-derived.
GDPR_DROP_TAXA = ["9605", "9527", "40674", "9606"]


@dataclass
class Kraken2Outputs:
    cleaned_r1: Path
    cleaned_r2: Path | None
    classified_r1: Path | None
    classified_r2: Path | None
    report: Path
    output: Path


def kraken2_paired(
    cfg: CerberusConfig,
    *,
    db: Path,
    r1_in: Path,
    r2_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
    keep_classified: bool = False,
) -> Kraken2Outputs:
    require_tools("kraken2")
    workdir.mkdir(parents=True, exist_ok=True)

    unclass_template = str(workdir / f"{tag}.kraken_unclass#.fq")
    report = workdir / f"{tag}.kraken2.report.txt"
    output = workdir / f"{tag}.kraken2.output.tsv"

    cmd = [
        "kraken2",
        "--db", str(db),
        "--threads", str(cfg.threads),
        "--paired",
        "--gzip-compressed",
        "--unclassified-out", unclass_template,
        "--report", str(report),
        "--output", str(output),
        "--use-names",
        "--confidence", "0.5",
    ]
    if keep_classified:
        cmd.extend([
            "--classified-out", str(workdir / f"{tag}.kraken_class#.fq"),
        ])
    cmd.extend([str(r1_in), str(r2_in)])

    run(cmd, log_path=log_dir / f"{tag}.kraken2.log", dry_run=cfg.dry_run)

    out_r1, out_r2 = _gzip_pair_outputs(workdir, f"{tag}.kraken_unclass", cfg)

    class_r1 = class_r2 = None
    if keep_classified:
        class_r1, class_r2 = _gzip_pair_outputs(workdir, f"{tag}.kraken_class", cfg)

    return Kraken2Outputs(
        cleaned_r1=out_r1,
        cleaned_r2=out_r2,
        classified_r1=class_r1,
        classified_r2=class_r2,
        report=report,
        output=output,
    )


def kraken2_single(
    cfg: CerberusConfig,
    *,
    db: Path,
    reads_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
) -> Kraken2Outputs:
    require_tools("kraken2")
    workdir.mkdir(parents=True, exist_ok=True)

    unclass = workdir / f"{tag}.kraken_unclass.fq"
    report = workdir / f"{tag}.kraken2.report.txt"
    output = workdir / f"{tag}.kraken2.output.tsv"

    cmd = [
        "kraken2",
        "--db", str(db),
        "--threads", str(cfg.threads),
        "--gzip-compressed",
        "--unclassified-out", str(unclass),
        "--report", str(report),
        "--output", str(output),
        "--use-names",
        "--confidence", "0.5",
        str(reads_in),
    ]
    run(cmd, log_path=log_dir / f"{tag}.kraken2.log", dry_run=cfg.dry_run)

    out_gz = _gzip_inplace(unclass, cfg)
    return Kraken2Outputs(
        cleaned_r1=out_gz,
        cleaned_r2=None,
        classified_r1=None,
        classified_r2=None,
        report=report,
        output=output,
    )


def _gzip_pair_outputs(
    workdir: Path, name_root: str, cfg: CerberusConfig,
) -> tuple[Path, Path]:
    """kraken2 produces ``<root>_1.fq`` / ``<root>_2.fq`` for paired data.
    Older versions used ``<root>1.fq`` / ``<root>2.fq``. Glob for either.
    Gzip both in place and return the .gz paths.
    """
    candidates_r1 = sorted(workdir.glob(f"{name_root}*1.fq")) + \
                    sorted(workdir.glob(f"{name_root}*1.fastq"))
    candidates_r2 = sorted(workdir.glob(f"{name_root}*2.fq")) + \
                    sorted(workdir.glob(f"{name_root}*2.fastq"))

    if not candidates_r1 or not candidates_r2:
        # Fall back: emit empty placeholders so downstream stages still see files.
        log.warning("kraken2 produced no R1/R2 unclassified files for %s; emitting empty placeholders", name_root)
        r1 = workdir / f"{name_root}_1.fq"
        r2 = workdir / f"{name_root}_2.fq"
        r1.touch()
        r2.touch()
    else:
        r1 = candidates_r1[0]
        r2 = candidates_r2[0]

    return _gzip_inplace(r1, cfg), _gzip_inplace(r2, cfg)


def _gzip_inplace(src: Path, cfg: CerberusConfig) -> Path:
    dst = src.with_suffix(src.suffix + ".gz")
    if cfg.dry_run:
        return dst
    if not src.exists():
        # Produce a valid empty gzip so downstream stages can read it.
        _write_empty_gzip(dst)
        return dst
    if which("pigz"):
        run(["pigz", "-f", "-p", str(cfg.threads), str(src)],
            log_path=src.with_suffix(".gzip.log"), dry_run=cfg.dry_run)
    else:
        import gzip
        with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
            for chunk in iter(lambda: fin.read(1024 * 1024), b""):
                fout.write(chunk)
        src.unlink()
    return dst


def _write_empty_gzip(path: Path) -> None:
    """Write a minimal valid empty gzip stream so bbduk/seqkit can read it."""
    import gzip
    with gzip.open(path, "wb") as f:
        f.write(b"")
