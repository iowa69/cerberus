"""BBDuk k-mer wrapper for auxiliary host-reference matching."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import require_tools, run

log = get_logger("kmer")


@dataclass
class BBDukKmerOutputs:
    r1: Path
    r2: Path | None
    stats: Path
    contaminated: Path | None


def _empty_paired_passthrough(
    workdir: Path, tag: str, r1_in: Path, r2_in: Path,
) -> BBDukKmerOutputs:
    """If kraken2 produced no surviving reads, write empty gzipped outputs
    so the rest of the pipeline still has files to chain."""
    out_r1 = workdir / f"{tag}.kmerclean.R1.fq.gz"
    out_r2 = workdir / f"{tag}.kmerclean.R2.fq.gz"
    stats = workdir / f"{tag}.bbduk_kmer.stats.txt"
    import gzip
    for p in (out_r1, out_r2):
        with gzip.open(p, "wb") as f:
            f.write(b"")
    stats.write_text(f"# bbduk skipped: input {r1_in.name} was empty\n")
    return BBDukKmerOutputs(r1=out_r1, r2=out_r2, stats=stats, contaminated=None)


def bbduk_kmer_paired(
    cfg: CerberusConfig,
    tuned: TunedParams,
    *,
    ref: Path,
    r1_in: Path,
    r2_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
) -> BBDukKmerOutputs:
    require_tools("bbduk.sh")
    workdir.mkdir(parents=True, exist_ok=True)

    if not r1_in.exists() or r1_in.stat().st_size <= 40:
        log.warning("Skipping bbduk %s — input is empty/near-empty (%s)", tag, r1_in)
        return _empty_paired_passthrough(workdir, tag, r1_in, r2_in)

    out_r1 = workdir / f"{tag}.kmerclean.R1.fq.gz"
    out_r2 = workdir / f"{tag}.kmerclean.R2.fq.gz"
    matched_r1 = workdir / f"{tag}.host_kmer.R1.fq.gz" if cfg.keep_intermediates else None
    matched_r2 = workdir / f"{tag}.host_kmer.R2.fq.gz" if cfg.keep_intermediates else None
    stats = workdir / f"{tag}.bbduk_kmer.stats.txt"

    cmd = [
        "bbduk.sh",
        f"-Xmx{cfg.memory_gb}g",
        f"in1={r1_in}", f"in2={r2_in}",
        f"out1={out_r1}", f"out2={out_r2}",
        f"ref={ref}",
        f"k={tuned.bbduk_k}",
        f"mcf={tuned.bbduk_mcf}",
        "rcomp=t",
        f"stats={stats}",
        f"threads={cfg.threads}",
    ]
    if matched_r1:
        cmd.extend([f"outm1={matched_r1}", f"outm2={matched_r2}"])

    run(cmd, log_path=log_dir / f"{tag}.bbduk_kmer.log", dry_run=cfg.dry_run)
    return BBDukKmerOutputs(r1=out_r1, r2=out_r2, stats=stats, contaminated=matched_r1)


def bbduk_kmer_single(
    cfg: CerberusConfig,
    tuned: TunedParams,
    *,
    ref: Path,
    reads_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
) -> BBDukKmerOutputs:
    require_tools("bbduk.sh")
    workdir.mkdir(parents=True, exist_ok=True)

    if not reads_in.exists() or reads_in.stat().st_size <= 40:
        log.warning("Skipping bbduk %s — input is empty/near-empty (%s)", tag, reads_in)
        out = workdir / f"{tag}.kmerclean.fq.gz"
        stats = workdir / f"{tag}.bbduk_kmer.stats.txt"
        import gzip
        with gzip.open(out, "wb") as f:
            f.write(b"")
        stats.write_text(f"# bbduk skipped: input {reads_in.name} was empty\n")
        return BBDukKmerOutputs(r1=out, r2=None, stats=stats, contaminated=None)

    out = workdir / f"{tag}.kmerclean.fq.gz"
    matched = workdir / f"{tag}.host_kmer.fq.gz" if cfg.keep_intermediates else None
    stats = workdir / f"{tag}.bbduk_kmer.stats.txt"

    cmd = [
        "bbduk.sh",
        f"-Xmx{cfg.memory_gb}g",
        f"in={reads_in}",
        f"out={out}",
        f"ref={ref}",
        f"k={tuned.bbduk_k}",
        f"mcf={tuned.bbduk_mcf}",
        "rcomp=t",
        f"stats={stats}",
        f"threads={cfg.threads}",
    ]
    if matched:
        cmd.append(f"outm={matched}")
    run(cmd, log_path=log_dir / f"{tag}.bbduk_kmer.log", dry_run=cfg.dry_run)
    return BBDukKmerOutputs(r1=out, r2=None, stats=stats, contaminated=matched)
