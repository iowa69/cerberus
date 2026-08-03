"""Kraken2 wrapper for the GDPR scrub. Extracts reads NOT classified as host.

Confidence
----------
Kraken2's ``--confidence`` is the fraction of a read's k-mers that must map
to the reported taxon's root-to-leaf path. The usual reason to raise it is to
suppress false *positive* identifications when the database spans all of
life. Cerberus' GDPR database is the opposite situation: it contains **only**
host taxa, so every classification is a hit on the thing we want gone, and a
high threshold merely lets host reads escape as "unclassified".

A 150 bp read has 121 31-mers; two sequencing errors already knock out ~62 of
them, so a genuine human read can score ~0.4 and slip past a 0.5 threshold.
On ONT data at 2-5% error, essentially no read reaches 0.5 at all, which made
the mechanism a complete no-op. The default is therefore low (0.05) and
exposed as ``--gdpr-confidence`` for anyone who needs to trade recall for
retained microbial signal.

Filename handling: kraken2's ``--unclassified-out '<dir>/foo#.fq'`` template
gets ``#`` replaced with ``_1``/``_2`` for paired data. We pass the template
unchanged and discover the actual files via glob — robust to substitution
variants between kraken2 versions.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

from cerberus.config import CerberusConfig
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import require_tools, run, which

log = get_logger("kraken")


class KrakenOutputError(RuntimeError):
    """Kraken2 ran but its outputs could not be located."""


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
        "--unclassified-out", unclass_template,
        "--report", str(report),
        "--output", str(output),
        "--use-names",
        "--confidence", str(cfg.gdpr_confidence),
    ]
    if _is_gzip(r1_in):
        cmd.append("--gzip-compressed")
    if keep_classified:
        cmd.extend(["--classified-out", str(workdir / f"{tag}.kraken_class#.fq")])
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
        "--unclassified-out", str(unclass),
        "--report", str(report),
        "--output", str(output),
        "--use-names",
        "--confidence", str(cfg.gdpr_confidence),
    ]
    if _is_gzip(reads_in):
        cmd.append("--gzip-compressed")
    cmd.append(str(reads_in))
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


def _is_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return str(path).endswith(".gz")


def _gzip_pair_outputs(
    workdir: Path, name_root: str, cfg: CerberusConfig,
) -> tuple[Path, Path]:
    """kraken2 produces ``<root>_1.fq`` / ``<root>_2.fq`` for paired data.

    Older versions used ``<root>1.fq`` / ``<root>2.fq``, so glob for either.
    Gzip both in place and return the .gz paths.

    A genuine "no reads survived" result and a naming mismatch used to look
    identical here — both produced empty placeholders. They are now
    distinguished: zero surviving reads is normal and yields valid empty
    gzips, but files that exist under an unrecognised name are an error,
    because silently publishing an empty GDPR release is the worst possible
    failure mode.
    """
    if cfg.dry_run:
        return (workdir / f"{name_root}_1.fq.gz", workdir / f"{name_root}_2.fq.gz")

    candidates_r1 = sorted(workdir.glob(f"{name_root}*1.fq")) + \
                    sorted(workdir.glob(f"{name_root}*1.fastq"))
    candidates_r2 = sorted(workdir.glob(f"{name_root}*2.fq")) + \
                    sorted(workdir.glob(f"{name_root}*2.fastq"))

    if not candidates_r1 or not candidates_r2:
        stray = [p.name for p in workdir.glob(f"{name_root}*")]
        if stray:
            raise KrakenOutputError(
                f"kraken2 wrote {stray} but Cerberus could not identify the R1/R2 pair for "
                f"{name_root!r}. Refusing to continue rather than emit an empty output."
            )
        log.info("kraken2 classified every read as host for %s; output is empty.", name_root)
        r1 = workdir / f"{name_root}_1.fq"
        r2 = workdir / f"{name_root}_2.fq"
        r1.write_bytes(b"")
        r2.write_bytes(b"")
    else:
        r1, r2 = candidates_r1[0], candidates_r2[0]

    return _gzip_inplace(r1, cfg), _gzip_inplace(r2, cfg)


def _gzip_inplace(src: Path, cfg: CerberusConfig) -> Path:
    dst = src.with_suffix(src.suffix + ".gz")
    if cfg.dry_run:
        return dst
    if not src.exists():
        if dst.exists():
            # Already compressed by a previous step; do not clobber real data.
            return dst
        _write_empty_gzip(dst)
        return dst
    if which("pigz"):
        run(["pigz", "-f", "-p", str(cfg.threads), str(src)],
            log_path=src.with_suffix(".gzip.log"), dry_run=cfg.dry_run)
    else:
        with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
            for chunk in iter(lambda: fin.read(1024 * 1024), b""):
                fout.write(chunk)
        src.unlink()
    return dst


def _write_empty_gzip(path: Path) -> None:
    """Write a minimal valid empty gzip stream so bbduk/seqkit can read it."""
    with gzip.open(path, "wb") as f:
        f.write(b"")
