"""Configuration objects for a Cerberus run.

A single CerberusConfig drives the entire pipeline. Autotune mutates a copy
of it after inspecting fastp/fastplong output; the original args from the CLI
are preserved in ``user_overrides`` so we never trample on explicit user input.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path


DEFAULT_REF_DIR = Path.home() / ".cerberus" / "refs"
DEFAULT_CACHE_DIR = Path.home() / ".cerberus" / "cache"


class Platform(str, Enum):
    AUTO = "auto"
    ILLUMINA = "illumina"
    ONT = "ont"
    PACBIO_HIFI = "pacbio-hifi"
    PACBIO_CLR = "pacbio-clr"


class ReadLengthClass(str, Enum):
    """Result of read-length autotuning. Drives parameter selection."""
    VERY_SHORT = "very_short"   # <80bp  (2x50, 2x75)
    SHORT = "short"             # 80-200bp (2x100, 2x150)
    MEDIUM = "medium"           # 200-500bp (2x250, 2x300)
    LONG = "long"               # >500bp
    VERY_LONG = "very_long"     # >5kb


@dataclass
class TunedParams:
    """Parameters resolved after auto-tuning. Applied to each stage."""
    min_length: int = 50
    min_quality: int = 20
    entropy: float = 0.7
    bbduk_k: int = 27
    bbduk_aux_enabled: bool = True
    minimap2_preset: str = "sr"
    minimap2_extra: str = ""
    bowtie2_preset: str = "--very-sensitive-local"
    bowtie2_extra: str = ""
    winnowmap_enabled: bool = False
    read_length_class: ReadLengthClass = ReadLengthClass.SHORT
    platform: Platform = Platform.ILLUMINA

    def summary(self) -> str:
        return (
            f"length-class={self.read_length_class.value} platform={self.platform.value} "
            f"min_len={self.min_length} Q={self.min_quality} entropy={self.entropy} "
            f"bbduk_k={self.bbduk_k} aux={self.bbduk_aux_enabled} "
            f"mm2={self.minimap2_preset} bt2={self.bowtie2_preset}"
        )


@dataclass
class CerberusConfig:
    # input
    r1: Path | None = None
    r2: Path | None = None
    long_input: Path | None = None
    long_mode: bool = False

    # output
    out_dir: Path = Path("cerberus_out")
    sample_id: str = "sample"

    # modes (at least one must be true)
    meta: bool = False
    profiling: bool = False
    gdpr: bool = False

    # behaviour modifiers
    platform: Platform = Platform.AUTO
    double_pass: bool = False
    fast: bool = False

    # resources
    threads: int = field(default_factory=lambda: os.cpu_count() or 4)
    memory_gb: int = 12

    # references
    ref_dir: Path = DEFAULT_REF_DIR
    cache_dir: Path = DEFAULT_CACHE_DIR
    auto_download: bool = True
    update_refs: bool = False
    kraken2_db_override: Path | None = None
    aux_refs_override: Path | None = None

    # power-user knobs (None means autotune decides)
    min_length: int | None = None
    min_quality: int | None = None
    entropy: float | None = None
    bbduk_k: int | None = None
    minimap2_args: str | None = None
    bowtie2_args: str | None = None

    # housekeeping
    keep_intermediates: bool = False
    verbose: bool = False
    quiet: bool = False
    dry_run: bool = False

    # filled in by autotune; never set by user directly
    tuned: TunedParams = field(default_factory=TunedParams)

    @property
    def modes(self) -> list[str]:
        out = []
        if self.meta:
            out.append("meta")
        if self.profiling:
            out.append("profiling")
        return out

    @property
    def work_dir(self) -> Path:
        return self.out_dir / "_work"

    @property
    def logs_dir(self) -> Path:
        return self.out_dir / "logs"

    @property
    def reports_dir(self) -> Path:
        return self.out_dir / "reports"

    def ensure_directories(self) -> None:
        for d in (self.out_dir, self.work_dir, self.logs_dir, self.reports_dir,
                  self.ref_dir, self.cache_dir):
            d.mkdir(parents=True, exist_ok=True)
