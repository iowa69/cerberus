"""Shared types and helpers for all pipelines."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class PipelineResult:
    mode: str                      # "meta" | "profiling" | "gdpr" | "long-meta" | ...
    paired_r1: Path | None = None
    paired_r2: Path | None = None
    singletons: Path | None = None    # merged file (profiling, gdpr) or orphan singles (meta)
    long_reads: Path | None = None    # long-read mode output
    stats: dict[str, int] = field(default_factory=dict)  # stage_name → read count
    logs: list[Path] = field(default_factory=list)

    @property
    def primary_output(self) -> Path | None:
        return self.paired_r1 or self.singletons or self.long_reads

    def as_inputs_for_gdpr(self) -> dict[str, Path]:
        """Map the pipeline outputs into the input shape GDPR expects."""
        out: dict[str, Path] = {}
        if self.paired_r1 and self.paired_r2:
            out["r1"] = self.paired_r1
            out["r2"] = self.paired_r2
        if self.singletons:
            out["singletons"] = self.singletons
        if self.long_reads:
            out["long"] = self.long_reads
        return out


def stage_dir(work_dir: Path, mode: str, stage: str) -> Path:
    d = work_dir / mode / stage
    d.mkdir(parents=True, exist_ok=True)
    return d
