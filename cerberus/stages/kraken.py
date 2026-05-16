"""Kraken2 wrapper for the GDPR scrub. Extracts reads NOT classified as host."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from cerberus.config import CerberusConfig
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import require_tools, run

log = get_logger("kraken")

# Taxa to scrub. 9605 = Homininae; 9527 = Catarrhini (apes incl. humans);
# 40674 = Mammalia. We use Mammalia to be conservative for GDPR.
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

    unclass_r1 = workdir / f"{tag}.kraken_unclass.R1.fq"
    unclass_r2 = workdir / f"{tag}.kraken_unclass.R2.fq"
    class_r1 = workdir / f"{tag}.kraken_class.R1.fq" if keep_classified else None
    class_r2 = workdir / f"{tag}.kraken_class.R2.fq" if keep_classified else None
    report = workdir / f"{tag}.kraken2.report.txt"
    output = workdir / f"{tag}.kraken2.output.tsv"

    cmd = [
        "kraken2",
        "--db", str(db),
        "--threads", str(cfg.threads),
        "--paired",
        "--gzip-compressed",
        "--unclassified-out", str(workdir / f"{tag}.kraken_unclass.R#.fq"),
        "--report", str(report),
        "--output", str(output),
        "--use-names",
        "--confidence", "0.5",
    ]
    if keep_classified:
        cmd.extend([
            "--classified-out", str(workdir / f"{tag}.kraken_class.R#.fq"),
        ])
    cmd.extend([str(r1_in), str(r2_in)])

    run(cmd, log_path=log_dir / f"{tag}.kraken2.log", dry_run=cfg.dry_run)

    # gzip the output FASTQs to keep storage tidy
    out_r1 = unclass_r1.with_suffix(unclass_r1.suffix + ".gz")
    out_r2 = unclass_r2.with_suffix(unclass_r2.suffix + ".gz")
    _gzip_if_needed(unclass_r1, out_r1, cfg.threads, cfg.dry_run)
    _gzip_if_needed(unclass_r2, out_r2, cfg.threads, cfg.dry_run)

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

    out_gz = unclass.with_suffix(unclass.suffix + ".gz")
    _gzip_if_needed(unclass, out_gz, cfg.threads, cfg.dry_run)
    return Kraken2Outputs(
        cleaned_r1=out_gz,
        cleaned_r2=None,
        classified_r1=None,
        classified_r2=None,
        report=report,
        output=output,
    )


def _gzip_if_needed(src: Path, dst: Path, threads: int, dry_run: bool) -> None:
    if not src.exists() or dry_run:
        return
    from cerberus.utils.shell import which
    if which("pigz"):
        run(["pigz", "-f", "-p", str(threads), str(src)],
            log_path=src.with_suffix(".gzip.log"), dry_run=dry_run)
    else:
        import gzip
        with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
            for chunk in iter(lambda: fin.read(1024 * 1024), b""):
                fout.write(chunk)
        src.unlink()
