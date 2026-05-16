"""Unit tests for the autotune branching logic.

We do not need fastp installed; we synthesize JSON reports that exercise each
read-length bucket and each platform-detection path.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from cerberus.autotune import (
    apply_user_overrides,
    autotune_from_fastp,
    classify_length,
    detect_platform_from_fastp,
)
from cerberus.config import CerberusConfig, Platform, ReadLengthClass


def _fastp_json(mean_length: float, q20: float = 0.95, q30: float = 0.85) -> dict:
    return {
        "summary": {
            "before_filtering": {
                "read1_mean_length": mean_length,
                "q20_rate": q20,
                "q30_rate": q30,
            }
        }
    }


@pytest.mark.parametrize("length,expected", [
    (50,    ReadLengthClass.VERY_SHORT),
    (79,    ReadLengthClass.VERY_SHORT),
    (80,    ReadLengthClass.SHORT),
    (150,   ReadLengthClass.SHORT),
    (199,   ReadLengthClass.SHORT),
    (200,   ReadLengthClass.MEDIUM),
    (300,   ReadLengthClass.MEDIUM),
    (499,   ReadLengthClass.MEDIUM),
    (500,   ReadLengthClass.LONG),
    (1500,  ReadLengthClass.LONG),
    (4999,  ReadLengthClass.LONG),
    (5000,  ReadLengthClass.VERY_LONG),
    (50000, ReadLengthClass.VERY_LONG),
])
def test_classify_length_buckets(length: int, expected: ReadLengthClass) -> None:
    assert classify_length(float(length)) is expected


def test_platform_user_override_wins() -> None:
    rep = _fastp_json(mean_length=15000, q30=0.9)
    assert detect_platform_from_fastp(rep, Platform.ILLUMINA) is Platform.ILLUMINA


def test_platform_detects_illumina_from_short_reads() -> None:
    rep = _fastp_json(mean_length=150)
    assert detect_platform_from_fastp(rep, Platform.AUTO) is Platform.ILLUMINA


def test_platform_detects_hifi_from_high_q() -> None:
    rep = _fastp_json(mean_length=15000, q20=0.99, q30=0.90)
    assert detect_platform_from_fastp(rep, Platform.AUTO) is Platform.PACBIO_HIFI


def test_platform_detects_ont_from_long_low_q() -> None:
    rep = _fastp_json(mean_length=15000, q20=0.90, q30=0.50)
    assert detect_platform_from_fastp(rep, Platform.AUTO) is Platform.ONT


def test_autotune_short_reads(tmp_path: Path) -> None:
    json_path = tmp_path / "fastp.json"
    json_path.write_text(json.dumps(_fastp_json(150)))
    tuned = autotune_from_fastp(json_path, user_platform=Platform.AUTO)
    assert tuned.read_length_class is ReadLengthClass.SHORT
    assert tuned.platform is Platform.ILLUMINA
    assert tuned.minimap2_preset == "sr"
    assert tuned.bbduk_aux_enabled is True


def test_autotune_very_short_disables_aux(tmp_path: Path) -> None:
    json_path = tmp_path / "fastp.json"
    json_path.write_text(json.dumps(_fastp_json(50)))
    tuned = autotune_from_fastp(json_path, user_platform=Platform.AUTO)
    assert tuned.read_length_class is ReadLengthClass.VERY_SHORT
    assert tuned.bbduk_aux_enabled is False
    assert tuned.min_length == 35


def test_autotune_long_ont(tmp_path: Path) -> None:
    json_path = tmp_path / "fastp.json"
    json_path.write_text(json.dumps(_fastp_json(8000, q20=0.92, q30=0.55)))
    tuned = autotune_from_fastp(json_path, user_platform=Platform.AUTO)
    assert tuned.read_length_class is ReadLengthClass.VERY_LONG
    assert tuned.platform is Platform.ONT
    assert tuned.minimap2_preset == "map-ont"
    assert tuned.winnowmap_enabled is True


def test_autotune_long_hifi_picks_map_hifi(tmp_path: Path) -> None:
    json_path = tmp_path / "fastp.json"
    json_path.write_text(json.dumps(_fastp_json(15000, q20=0.99, q30=0.92)))
    tuned = autotune_from_fastp(json_path, user_platform=Platform.AUTO)
    assert tuned.platform is Platform.PACBIO_HIFI
    assert tuned.minimap2_preset == "map-hifi"


def test_user_overrides_win(tmp_path: Path) -> None:
    json_path = tmp_path / "fastp.json"
    json_path.write_text(json.dumps(_fastp_json(150)))
    tuned = autotune_from_fastp(json_path, user_platform=Platform.AUTO)

    cfg = CerberusConfig(
        meta=True,
        min_length=99,
        entropy=0.42,
        bbduk_k=21,
        platform=Platform.ONT,
    )
    overridden = apply_user_overrides(tuned, cfg)
    assert overridden.min_length == 99
    assert overridden.entropy == 0.42
    assert overridden.bbduk_k == 21
    assert overridden.platform is Platform.ONT
    assert overridden.minimap2_preset == "map-ont"
