"""BBDuk entropy filter — strips low-complexity reads at the tail of each pipeline."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import require_tools, run

log = get_logger("entropy")


@dataclass
class EntropyOutputs:
    r1: Path
    r2: Path | None
    stats: Path


def entropy_paired(
    cfg: CerberusConfig,
    tuned: TunedParams,
    *,
    r1_in: Path,
    r2_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
) -> EntropyOutputs:
    require_tools("bbduk.sh")
    workdir.mkdir(parents=True, exist_ok=True)
    out_r1 = workdir / f"{tag}.entropy.R1.fq.gz"
    out_r2 = workdir / f"{tag}.entropy.R2.fq.gz"
    stats = workdir / f"{tag}.entropy.stats.txt"

    entropy = cfg.entropy if cfg.entropy is not None else tuned.entropy
    # The window cannot exceed the reads it is measured over.
    window = min(tuned.entropy_window, max(10, tuned.min_length))
    cmd = [
        "bbduk.sh",
        f"-Xmx{max(2, cfg.memory_gb // 2)}g",
        f"in1={r1_in}", f"in2={r2_in}",
        f"out1={out_r1}", f"out2={out_r2}",
        f"entropy={entropy}",
        f"entropywindow={window}",
        "entropyk=5",
        f"stats={stats}",
        f"threads={cfg.threads}",
    ]
    run(cmd, log_path=log_dir / f"{tag}.entropy.log", dry_run=cfg.dry_run)
    return EntropyOutputs(r1=out_r1, r2=out_r2, stats=stats)


def entropy_single(
    cfg: CerberusConfig,
    tuned: TunedParams,
    *,
    reads_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
) -> EntropyOutputs:
    require_tools("bbduk.sh")
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / f"{tag}.entropy.fq.gz"
    stats = workdir / f"{tag}.entropy.stats.txt"

    entropy = cfg.entropy if cfg.entropy is not None else tuned.entropy
    # The window cannot exceed the reads it is measured over.
    window = min(tuned.entropy_window, max(10, tuned.min_length))
    cmd = [
        "bbduk.sh",
        f"-Xmx{max(2, cfg.memory_gb // 2)}g",
        f"in={reads_in}",
        f"out={out}",
        f"entropy={entropy}",
        f"entropywindow={window}",
        "entropyk=5",
        f"stats={stats}",
        f"threads={cfg.threads}",
    ]
    run(cmd, log_path=log_dir / f"{tag}.entropy.log", dry_run=cfg.dry_run)
    return EntropyOutputs(r1=out, r2=None, stats=stats)
