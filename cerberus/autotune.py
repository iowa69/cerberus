"""Autotune: pick the best parameters for the data Cerberus is about to process.

The contract is:

  1. ``prescan_reads()`` samples the head of the input FASTQ and reports mean
     read length and Q20/Q30 rates. This happens *before* QC, so the tuned
     values can actually reach fastp.
  2. ``autotune_from_prescan()`` (or ``autotune_from_fastp()`` when a real
     fastp report already exists) turns that into a ``TunedParams``.
  3. ``apply_user_overrides()`` lets explicit CLI flags win.
  4. Downstream stages read ``CerberusConfig.tuned``.

Earlier versions ran fastp first and tuned afterwards, which meant the tuned
``min_length``/``min_quality`` were computed too late to reach the only stage
that consumes them. The prescan removes that ordering problem.

The branching rules are deliberately tabular and testable. Any tweak to a
threshold belongs in ``_LENGTH_BUCKETS`` or ``_BASE_PARAMS`` and nowhere
else, so the behaviour stays auditable for reviewers.
"""
from __future__ import annotations

import gzip
import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import IO, Any

from cerberus.config import CerberusConfig, Platform, ReadLengthClass, TunedParams
from cerberus.utils.logger import get_logger

log = get_logger("autotune")


_LENGTH_BUCKETS: list[tuple[int, ReadLengthClass]] = [
    (80, ReadLengthClass.VERY_SHORT),
    (200, ReadLengthClass.SHORT),
    (500, ReadLengthClass.MEDIUM),
    (5_000, ReadLengthClass.LONG),
    (10**9, ReadLengthClass.VERY_LONG),
]


def classify_length(mean_length: float) -> ReadLengthClass:
    for upper, cls in _LENGTH_BUCKETS:
        if mean_length < upper:
            return cls
    return ReadLengthClass.VERY_LONG


# bbduk's ``mcf`` (minimum covered fraction) has to fall as the error rate
# rises: at 5% error only ~20% of a long read's 31-mers survive error-free,
# so the short-read default of 0.5 makes the auxiliary k-mer pass inert on
# ONT/CLR data. Long-read classes therefore get a smaller k and a much lower
# covered-fraction requirement.
_BASE_PARAMS: dict[ReadLengthClass, TunedParams] = {
    ReadLengthClass.VERY_SHORT: TunedParams(
        min_length=35, min_quality=20, entropy=0.60, bbduk_k=23,
        bbduk_aux_enabled=True, bbduk_mcf=0.5,
        minimap2_preset="sr", minimap2_extra="-k15 -w10",
        bowtie2_preset="--very-sensitive-local", bowtie2_extra="",
        entropy_window=20,
        read_length_class=ReadLengthClass.VERY_SHORT,
    ),
    ReadLengthClass.SHORT: TunedParams(
        min_length=50, min_quality=20, entropy=0.70, bbduk_k=27,
        bbduk_aux_enabled=True, bbduk_mcf=0.5,
        minimap2_preset="sr", bowtie2_preset="--very-sensitive-local",
        entropy_window=50,
        read_length_class=ReadLengthClass.SHORT,
    ),
    ReadLengthClass.MEDIUM: TunedParams(
        min_length=75, min_quality=20, entropy=0.65, bbduk_k=31,
        bbduk_aux_enabled=True, bbduk_mcf=0.5,
        minimap2_preset="sr", bowtie2_preset="--very-sensitive-local",
        entropy_window=50,
        read_length_class=ReadLengthClass.MEDIUM,
    ),
    ReadLengthClass.LONG: TunedParams(
        min_length=200, min_quality=10, entropy=0.50, bbduk_k=25,
        bbduk_aux_enabled=True, bbduk_mcf=0.10,
        minimap2_preset="map-ont", winnowmap_enabled=False,
        entropy_window=50,
        read_length_class=ReadLengthClass.LONG,
    ),
    ReadLengthClass.VERY_LONG: TunedParams(
        min_length=500, min_quality=10, entropy=0.45, bbduk_k=25,
        bbduk_aux_enabled=True, bbduk_mcf=0.05,
        minimap2_preset="map-ont", winnowmap_enabled=False,
        entropy_window=50,
        read_length_class=ReadLengthClass.VERY_LONG,
    ),
}


@dataclass
class Prescan:
    """Cheap statistics from the head of an input FASTQ."""
    reads_sampled: int = 0
    mean_length: float = 0.0
    q20_rate: float = 0.0
    q30_rate: float = 0.0
    source: str = ""

    @property
    def ok(self) -> bool:
        return self.reads_sampled > 0


def _open_text(path: Path) -> IO[str]:
    """Open a FASTQ whether or not the name advertises its compression."""
    with path.open("rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", errors="replace")
    return path.open("rt", errors="replace")


def prescan_reads(path: Path, sample_records: int = 20_000) -> Prescan:
    """Sample the first ``sample_records`` reads for length and quality.

    Compression is detected by magic bytes, not by filename, so a gzipped
    file called ``reads.fastq`` works. Any read error degrades to an empty
    Prescan rather than killing the run — the caller falls back to defaults.
    """
    total_len = 0
    q20 = q30 = qbases = 0
    n = 0
    try:
        with _open_text(path) as f:
            for i, line in enumerate(f):
                phase = i % 4
                if phase == 1:
                    total_len += len(line.rstrip("\n"))
                    n += 1
                elif phase == 3:
                    for ch in line.rstrip("\n"):
                        q = ord(ch) - 33
                        qbases += 1
                        if q >= 20:
                            q20 += 1
                            if q >= 30:
                                q30 += 1
                    if n >= sample_records:
                        break
    except (OSError, EOFError, UnicodeDecodeError) as e:
        log.warning("Could not prescan %s (%s); falling back to defaults", path, e)
        return Prescan(source=str(path))

    if n == 0:
        return Prescan(source=str(path))
    return Prescan(
        reads_sampled=n,
        mean_length=total_len / n,
        q20_rate=(q20 / qbases) if qbases else 0.0,
        q30_rate=(q30 / qbases) if qbases else 0.0,
        source=str(path),
    )


def _platform_preset(platform: Platform, length_class: ReadLengthClass) -> str:
    if platform is Platform.PACBIO_HIFI:
        return "map-hifi"
    if platform is Platform.PACBIO_CLR:
        return "map-pb"
    if platform is Platform.ONT:
        return "map-ont"
    if platform is Platform.ILLUMINA:
        return "sr"
    # AUTO: derive from length class
    return _BASE_PARAMS[length_class].minimap2_preset


def detect_platform(
    mean_len: float,
    q20_rate: float,
    q30_rate: float,
    user_choice: Platform = Platform.AUTO,
    *,
    long_mode: bool = False,
) -> Platform:
    """Infer the platform when the user said 'auto'.

    Deliberately conservative. Two rules matter:

    * ``long_mode`` is respected — a 400 bp ONT library must not be called
      Illumina merely because it sits under the length threshold.
    * HiFi is claimed only on genuinely HiFi-grade accuracy. ONT R10.4.1
      duplex reaches Q20+ but not Q30 across 85% of bases, so the Q30 test is
      what separates them; the previous ``q20 >= 0.99`` alternative filed
      modern ONT as PacBio and picked ``map-hifi`` for 1–5% error reads.
    """
    if user_choice is not Platform.AUTO:
        return user_choice

    if mean_len < 500 and not long_mode:
        return Platform.ILLUMINA
    # long-read regime
    if q30_rate >= 0.85:
        return Platform.PACBIO_HIFI
    return Platform.ONT


def detect_platform_from_fastp(report: dict[str, Any], user_choice: Platform) -> Platform:
    """Backwards-compatible wrapper reading a fastp-shaped report dict."""
    summary = _safe_dict(report, "summary")
    before = _safe_dict(summary, "before_filtering")
    return detect_platform(
        _num(before, "read1_mean_length", "read_mean_length", "read2_mean_length"),
        _num(before, "q20_rate"),
        _num(before, "q30_rate"),
        user_choice,
    )


def _safe_dict(obj: Any, key: str) -> dict:
    val = obj.get(key) if isinstance(obj, dict) else None
    return val if isinstance(val, dict) else {}


def _num(d: dict, *keys: str) -> float:
    """First numeric value among ``keys``. Missing/garbage yields 0.0."""
    for k in keys:
        v = d.get(k)
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            return float(v)
    return 0.0


def _build(mean_len: float, platform: Platform, *, origin: str) -> TunedParams:
    length_class = classify_length(float(mean_len))
    base = _BASE_PARAMS[length_class]
    tuned = replace(
        base,
        platform=platform,
        minimap2_preset=_platform_preset(platform, length_class),
        observed_mean_length=float(mean_len),
    )
    log.info("Autotune (%s): mean_len=%.0fbp class=%s platform=%s",
             origin, mean_len, length_class.value, platform.value)
    log.info("Tuned params: %s", tuned.summary())
    return tuned


def autotune_from_prescan(
    prescan: Prescan,
    *,
    user_platform: Platform = Platform.AUTO,
    long_mode: bool = False,
) -> TunedParams:
    """Produce TunedParams from a cheap pre-QC scan of the input."""
    if not prescan.ok:
        length_class = ReadLengthClass.LONG if long_mode else ReadLengthClass.SHORT
        log.warning(
            "Prescan produced no reads; defaulting to the %s profile. "
            "Pass --min-length/--min-quality/--platform to override.",
            length_class.value,
        )
        platform = user_platform if user_platform is not Platform.AUTO else (
            Platform.ONT if long_mode else Platform.ILLUMINA
        )
        return replace(
            _BASE_PARAMS[length_class],
            platform=platform,
            minimap2_preset=_platform_preset(platform, length_class),
        )

    platform = detect_platform(
        prescan.mean_length, prescan.q20_rate, prescan.q30_rate,
        user_platform, long_mode=long_mode,
    )
    return _build(prescan.mean_length, platform, origin="prescan")


def autotune_from_fastp(
    fastp_json: Path,
    *,
    user_platform: Platform = Platform.AUTO,
    long_mode: bool = False,
) -> TunedParams:
    """Read fastp's JSON output and produce a TunedParams set.

    Malformed or missing reports degrade to the default profile with a
    warning instead of aborting a multi-hour run.
    """
    try:
        with fastp_json.open() as f:
            report = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        log.warning("Could not read %s (%s); using defaults", fastp_json, e)
        report = {}

    summary = _safe_dict(report, "summary")
    before = _safe_dict(summary, "before_filtering")
    mean_len = _num(before, "read1_mean_length", "read_mean_length", "read2_mean_length")

    if mean_len <= 0:
        length_class = ReadLengthClass.LONG if long_mode else ReadLengthClass.SHORT
        log.warning(
            "fastp report %s carries no usable read length; defaulting to the %s profile",
            fastp_json.name, length_class.value,
        )
        mean_len = float(_BASE_PARAMS[length_class].min_length * 3)

    platform = detect_platform(
        mean_len, _num(before, "q20_rate"), _num(before, "q30_rate"),
        user_platform, long_mode=long_mode,
    )
    return _build(mean_len, platform, origin=fastp_json.name)


def apply_user_overrides(tuned: TunedParams, cfg: CerberusConfig) -> TunedParams:
    """User flags always win over autotune."""
    updates: dict[str, Any] = {}
    if cfg.min_length is not None:
        updates["min_length"] = cfg.min_length
    if cfg.min_quality is not None:
        updates["min_quality"] = cfg.min_quality
    if cfg.entropy is not None:
        updates["entropy"] = cfg.entropy
    if cfg.bbduk_k is not None:
        updates["bbduk_k"] = cfg.bbduk_k
    if cfg.platform is not Platform.AUTO:
        updates["platform"] = cfg.platform
        updates["minimap2_preset"] = _platform_preset(cfg.platform, tuned.read_length_class)

    if not updates:
        return tuned
    log.info("User overrides applied: %s", updates)
    return replace(tuned, **updates)


def estimate_long_read(input_path: Path, sample_records: int = 1000) -> bool:
    """True when the input really looks like long reads (mean > 1kb)."""
    return prescan_reads(input_path, sample_records=sample_records).mean_length > 1000
