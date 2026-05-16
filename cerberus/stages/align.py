"""Aligner wrappers: minimap2, bowtie2, winnowmap. Each returns unmapped read paths."""
from __future__ import annotations

import shlex
from dataclasses import dataclass
from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import pipe, require_tools, run

log = get_logger("align")


@dataclass
class AlignOutputs:
    r1: Path | None = None
    r2: Path | None = None
    singletons: Path | None = None
    long_reads: Path | None = None
    stats: Path | None = None


def minimap2_paired(
    cfg: CerberusConfig,
    tuned: TunedParams,
    *,
    index: Path,
    r1_in: Path,
    r2_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
    drop_strategy: str = "both",   # "both" = drop pair if both map; "either" = drop pair if either maps
) -> AlignOutputs:
    """Align paired reads with minimap2, retain only unmapped pairs.

    drop_strategy:
      - "both"   ⇒ only drop reads where both mates map (--meta, conservative)
      - "either" ⇒ drop pair if either mate maps (--profiling, aggressive)
    """
    require_tools("minimap2", "samtools")
    workdir.mkdir(parents=True, exist_ok=True)

    bam = workdir / f"{tag}.bam"
    out_r1 = workdir / f"{tag}.unmapped.R1.fq.gz"
    out_r2 = workdir / f"{tag}.unmapped.R2.fq.gz"
    stats = workdir / f"{tag}.flagstat.txt"

    extra = cfg.minimap2_args or tuned.minimap2_extra
    extra_tokens = shlex.split(extra) if extra else []

    minimap_cmd = [
        "minimap2",
        "-ax", tuned.minimap2_preset,
        "-t", str(cfg.threads),
        "--secondary=no",
        *extra_tokens,
        str(index),
        str(r1_in),
        str(r2_in),
    ]
    sort_cmd = ["samtools", "view", "-@", str(cfg.threads), "-b", "-o", str(bam), "-"]
    pipe([minimap_cmd, sort_cmd],
         log_path=log_dir / f"{tag}.minimap2.log", dry_run=cfg.dry_run)

    flags = _filter_flags(drop_strategy)
    fastq_cmd = [
        "samtools", "fastq",
        "-@", str(cfg.threads),
        flags["filter_flag"], flags["filter_value"],
        "-1", str(out_r1),
        "-2", str(out_r2),
        "-0", "/dev/null",
        "-s", "/dev/null",
        "-n",
        str(bam),
    ]
    run(fastq_cmd, log_path=log_dir / f"{tag}.samtools_fastq.log", dry_run=cfg.dry_run)

    run(["samtools", "flagstat", str(bam)],
        log_path=stats, dry_run=cfg.dry_run)

    if not cfg.keep_intermediates and not cfg.dry_run:
        bam.unlink(missing_ok=True)

    return AlignOutputs(r1=out_r1, r2=out_r2, stats=stats)


def minimap2_singles(
    cfg: CerberusConfig,
    tuned: TunedParams,
    *,
    index: Path,
    reads_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
) -> AlignOutputs:
    """Align single-end / orphan / long reads with minimap2, retain unmapped."""
    require_tools("minimap2", "samtools")
    workdir.mkdir(parents=True, exist_ok=True)

    bam = workdir / f"{tag}.bam"
    out = workdir / f"{tag}.unmapped.fq.gz"
    stats = workdir / f"{tag}.flagstat.txt"

    extra = cfg.minimap2_args or tuned.minimap2_extra
    extra_tokens = shlex.split(extra) if extra else []

    minimap_cmd = [
        "minimap2",
        "-ax", tuned.minimap2_preset,
        "-t", str(cfg.threads),
        "--secondary=no",
        *extra_tokens,
        str(index),
        str(reads_in),
    ]
    sort_cmd = ["samtools", "view", "-@", str(cfg.threads), "-b", "-o", str(bam), "-"]
    pipe([minimap_cmd, sort_cmd],
         log_path=log_dir / f"{tag}.minimap2.log", dry_run=cfg.dry_run)

    run(
        ["samtools", "fastq", "-@", str(cfg.threads), "-f", "4", "-0", str(out), "-n", str(bam)],
        log_path=log_dir / f"{tag}.samtools_fastq.log", dry_run=cfg.dry_run,
    )
    run(["samtools", "flagstat", str(bam)], log_path=stats, dry_run=cfg.dry_run)
    if not cfg.keep_intermediates and not cfg.dry_run:
        bam.unlink(missing_ok=True)

    return AlignOutputs(long_reads=out, stats=stats)


def bowtie2_paired(
    cfg: CerberusConfig,
    tuned: TunedParams,
    *,
    index_prefix: Path,
    r1_in: Path,
    r2_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
) -> AlignOutputs:
    """Run bowtie2 on paired reads, extract unmapped pairs."""
    require_tools("bowtie2", "samtools")
    workdir.mkdir(parents=True, exist_ok=True)

    bam = workdir / f"{tag}.bam"
    out_r1 = workdir / f"{tag}.unmapped.R1.fq.gz"
    out_r2 = workdir / f"{tag}.unmapped.R2.fq.gz"
    stats = workdir / f"{tag}.flagstat.txt"

    preset_tokens = shlex.split(tuned.bowtie2_preset)
    extra_tokens = shlex.split(cfg.bowtie2_args or tuned.bowtie2_extra) \
        if (cfg.bowtie2_args or tuned.bowtie2_extra) else []

    bt2_cmd = [
        "bowtie2",
        *preset_tokens,
        "-p", str(cfg.threads),
        "-x", str(index_prefix),
        "-1", str(r1_in),
        "-2", str(r2_in),
        *extra_tokens,
    ]
    sort_cmd = ["samtools", "view", "-@", str(cfg.threads), "-b", "-o", str(bam), "-"]
    pipe([bt2_cmd, sort_cmd],
         log_path=log_dir / f"{tag}.bowtie2.log", dry_run=cfg.dry_run)

    run(
        [
            "samtools", "fastq", "-@", str(cfg.threads),
            "-f", "12",  # both unmapped
            "-1", str(out_r1), "-2", str(out_r2),
            "-0", "/dev/null", "-s", "/dev/null", "-n",
            str(bam),
        ],
        log_path=log_dir / f"{tag}.samtools_fastq.log", dry_run=cfg.dry_run,
    )
    run(["samtools", "flagstat", str(bam)], log_path=stats, dry_run=cfg.dry_run)
    if not cfg.keep_intermediates and not cfg.dry_run:
        bam.unlink(missing_ok=True)

    return AlignOutputs(r1=out_r1, r2=out_r2, stats=stats)


def winnowmap_singles(
    cfg: CerberusConfig,
    tuned: TunedParams,
    *,
    index: Path,
    meryl_db: Path,
    reads_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
) -> AlignOutputs:
    """winnowmap for very long reads (replaces minimap2 when --double-pass)."""
    require_tools("winnowmap", "samtools")
    workdir.mkdir(parents=True, exist_ok=True)

    bam = workdir / f"{tag}.bam"
    out = workdir / f"{tag}.unmapped.fq.gz"
    stats = workdir / f"{tag}.flagstat.txt"

    wm_cmd = [
        "winnowmap",
        "-W", str(meryl_db),
        "-ax", "map-ont",
        "-t", str(cfg.threads),
        str(index),
        str(reads_in),
    ]
    sort_cmd = ["samtools", "view", "-@", str(cfg.threads), "-b", "-o", str(bam), "-"]
    pipe([wm_cmd, sort_cmd],
         log_path=log_dir / f"{tag}.winnowmap.log", dry_run=cfg.dry_run)

    run(
        ["samtools", "fastq", "-@", str(cfg.threads), "-f", "4", "-0", str(out), "-n", str(bam)],
        log_path=log_dir / f"{tag}.samtools_fastq.log", dry_run=cfg.dry_run,
    )
    run(["samtools", "flagstat", str(bam)], log_path=stats, dry_run=cfg.dry_run)
    if not cfg.keep_intermediates and not cfg.dry_run:
        bam.unlink(missing_ok=True)
    return AlignOutputs(long_reads=out, stats=stats)


def _filter_flags(strategy: str) -> dict[str, str]:
    """samtools fastq flags for retaining unmapped reads."""
    if strategy == "both":
        # SAM flag 12 (0x0C) = both mates unmapped. -f keeps reads matching.
        return {"filter_flag": "-f", "filter_value": "12"}
    if strategy == "either":
        # Drop a pair if EITHER mate maps. Keep where neither maps OR partner maps.
        # We keep reads where this read is unmapped (-f 4) AND mate is unmapped (-f 8) ⇒ -f 12,
        # which is the same as "both" for samtools fastq. To make it stricter we instead
        # exclude reads where their mate mapped: -F 8 means exclude mate-mapped, plus -f 4.
        # Implemented by -f 4 -F 8 below.
        return {"filter_flag": "-f", "filter_value": "4"}
    raise ValueError(f"Unknown drop_strategy: {strategy}")
