"""Shared types and helpers for all pipelines."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from cerberus.config import CerberusConfig
from cerberus.utils.fastq import count_reads
from cerberus.utils.logger import get_logger

log = get_logger("pipeline")


@dataclass
class PipelineResult:
    mode: str                      # "meta" | "profiling" | "long-meta" | ...
    paired_r1: Path | None = None
    paired_r2: Path | None = None
    singletons: Path | None = None    # merged file (profiling) or orphan singles (meta)
    long_reads: Path | None = None    # long-read mode output
    stats: dict[str, int] = field(default_factory=dict)  # stage_name -> read count
    logs: list[Path] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def primary_output(self) -> Path | None:
        return self.paired_r1 or self.singletons or self.long_reads


def stage_dir(work_dir: Path, stage: str) -> Path:
    """Directory for one stage's intermediates, under the mode's work dir."""
    d = work_dir / stage
    d.mkdir(parents=True, exist_ok=True)
    return d


def publish(cfg: CerberusConfig, src: Path | None, dst: Path) -> Path | None:
    """Move a stage output to its final location.

    Nothing is unlinked until the replacement is known to exist, so a
    ``--dry-run`` (or a stage that produced nothing) can never destroy the
    results of a previous real run.
    """
    if cfg.dry_run:
        return dst
    if src is None or not src.exists():
        log.warning("Expected output %s was not produced; leaving %s untouched", src, dst.name)
        return None
    dst.parent.mkdir(parents=True, exist_ok=True)
    src.replace(dst)
    return dst


class StageTracker:
    """Collects per-stage surviving-record counts for the accounting table."""

    def __init__(self, cfg: CerberusConfig):
        self.cfg = cfg
        self.stats: dict[str, int] = {}
        self.warnings: list[str] = []

    def record(self, name: str, *paths: Path | None) -> int:
        if self.cfg.dry_run:
            return 0
        total = sum(count_reads(p) for p in paths if p is not None)
        self.stats[name] = total
        return total

    def check_not_empty(self, stage: str, before: int, after: int) -> None:
        """Flag a stage that consumed reads and emitted none."""
        if before > 0 and after == 0:
            self.warnings.append(
                f"Stage '{stage}' removed every read ({before} in, 0 out). "
                "That usually means a reference/preset mismatch rather than a clean sample."
            )
