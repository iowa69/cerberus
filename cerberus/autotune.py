"""Autotune: pick the best parameters for the data Cerberus is about to process.

The contract is:
  - run fastp / fastplong first and write its JSON report
  - feed that JSON to ``autotune_from_fastp()``
  - merge the result with user overrides via ``apply_user_overrides()``
  - downstream stages read ``CerberusConfig.tuned``

The branching rules are deliberately tabular and testable. Any tweak to a
threshold belongs in ``_LENGTH_BUCKETS`` or ``_PLATFORM_RULES`` and nowhere
else, so the behaviour stays auditable for reviewers.
"""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

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


_BASE_PARAMS: dict[ReadLengthClass, TunedParams] = {
    ReadLengthClass.VERY_SHORT: TunedParams(
        min_length=35, min_quality=20, entropy=0.60, bbduk_k=23,
        bbduk_aux_enabled=False,
        minimap2_preset="sr", minimap2_extra="-k15 -w10",
        bowtie2_preset="--very-sensitive-local", bowtie2_extra="--no-1mm-upfront",
        read_length_class=ReadLengthClass.VERY_SHORT,
    ),
    ReadLengthClass.SHORT: TunedParams(
        min_length=50, min_quality=20, entropy=0.70, bbduk_k=27,
        bbduk_aux_enabled=True,
        minimap2_preset="sr", bowtie2_preset="--very-sensitive-local",
        read_length_class=ReadLengthClass.SHORT,
    ),
    ReadLengthClass.MEDIUM: TunedParams(
        min_length=75, min_quality=20, entropy=0.65, bbduk_k=31,
        bbduk_aux_enabled=True,
        minimap2_preset="sr", bowtie2_preset="--very-sensitive-local",
        read_length_class=ReadLengthClass.MEDIUM,
    ),
    ReadLengthClass.LONG: TunedParams(
        min_length=200, min_quality=10, entropy=0.50, bbduk_k=31,
        bbduk_aux_enabled=True,
        minimap2_preset="map-ont", winnowmap_enabled=False,
        read_length_class=ReadLengthClass.LONG,
    ),
    ReadLengthClass.VERY_LONG: TunedParams(
        min_length=500, min_quality=10, entropy=0.45, bbduk_k=31,
        bbduk_aux_enabled=True,
        minimap2_preset="map-ont", winnowmap_enabled=True,
        read_length_class=ReadLengthClass.VERY_LONG,
    ),
}


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


def detect_platform_from_fastp(report: dict[str, Any], user_choice: Platform) -> Platform:
    """Infer platform when user said 'auto'. Heuristics are deliberately conservative."""
    if user_choice is not Platform.AUTO:
        return user_choice

    summary = report.get("summary", {})
    before = summary.get("before_filtering", {})
    mean_len = before.get("read1_mean_length") or before.get("read_mean_length") or 0
    q20_rate = before.get("q20_rate", 0)
    q30_rate = before.get("q30_rate", 0)

    if mean_len < 500:
        return Platform.ILLUMINA
    # long-read regime
    if q30_rate >= 0.85 or q20_rate >= 0.99:
        return Platform.PACBIO_HIFI
    return Platform.ONT


def autotune_from_fastp(
    fastp_json: Path,
    *,
    user_platform: Platform = Platform.AUTO,
) -> TunedParams:
    """Read fastp's JSON output and produce a TunedParams set."""
    with fastp_json.open() as f:
        report = json.load(f)

    summary = report.get("summary", {})
    before = summary.get("before_filtering", {})
    mean_len = (
        before.get("read1_mean_length")
        or before.get("read_mean_length")
        or before.get("read2_mean_length")
        or 0
    )

    length_class = classify_length(float(mean_len))
    platform = detect_platform_from_fastp(report, user_platform)

    base = _BASE_PARAMS[length_class]
    tuned = replace(
        base,
        platform=platform,
        minimap2_preset=_platform_preset(platform, length_class),
    )

    log.info(
        "Autotune: mean_len=%.0fbp class=%s platform=%s",
        mean_len, length_class.value, platform.value,
    )
    log.info("Tuned params: %s", tuned.summary())
    return tuned


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
    """Cheap pre-fastp heuristic to short-circuit long-read mode detection.

    Reads the first ``sample_records`` records, returns True if the mean
    length exceeds 1kb. Used only for the ``--long`` flag verification, not
    for parameter selection (that's autotune's job).
    """
    import gzip
    opener = gzip.open if input_path.suffix == ".gz" else open
    total = 0
    n = 0
    with opener(input_path, "rt") as f:  # type: ignore[operator]
        for i, line in enumerate(f):
            if i % 4 == 1:  # sequence line
                total += len(line.strip())
                n += 1
                if n >= sample_records:
                    break
    if n == 0:
        return False
    return (total / n) > 1000
