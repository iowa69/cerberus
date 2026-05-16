"""CLI parsing tests. No external tools invoked — these only exercise argparse + config mapping."""
from __future__ import annotations

from pathlib import Path

import pytest

from cerberus.cli import _build_full_parser, _config_from_args, _parse_memory
from cerberus.config import Platform
from cerberus.orchestrator import ConfigError, validate_config


def _ns(*args: str):
    parser = _build_full_parser()
    return parser.parse_args(list(args))


def test_basic_paired_invocation(tmp_path: Path) -> None:
    r1 = tmp_path / "R1.fq.gz"
    r2 = tmp_path / "R2.fq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    ns = _ns("-r1", str(r1), "-r2", str(r2), "-o", str(tmp_path / "out"), "--meta")
    cfg = _config_from_args(ns)
    assert cfg.meta and not cfg.profiling and not cfg.gdpr
    assert cfg.r1 == r1
    assert cfg.r2 == r2


def test_all_flag_enables_three_modes(tmp_path: Path) -> None:
    r1 = tmp_path / "R1.fq.gz"
    r2 = tmp_path / "R2.fq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    ns = _ns("-r1", str(r1), "-r2", str(r2), "--all")
    cfg = _config_from_args(ns)
    assert cfg.meta and cfg.profiling and cfg.gdpr


def test_long_mode_requires_long_input(tmp_path: Path) -> None:
    ns = _ns("--long", "--meta")
    cfg = _config_from_args(ns)
    with pytest.raises(ConfigError, match="--long requires"):
        validate_config(cfg)


def test_gdpr_alone_rejected(tmp_path: Path) -> None:
    r1 = tmp_path / "R1.fq.gz"
    r2 = tmp_path / "R2.fq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    ns = _ns("-r1", str(r1), "-r2", str(r2), "--gdpr")
    cfg = _config_from_args(ns)
    with pytest.raises(ConfigError, match="--gdpr requires"):
        validate_config(cfg)


def test_fast_and_double_pass_conflict(tmp_path: Path) -> None:
    r1 = tmp_path / "R1.fq.gz"
    r2 = tmp_path / "R2.fq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    ns = _ns("-r1", str(r1), "-r2", str(r2), "--profiling", "--fast", "--double-pass")
    cfg = _config_from_args(ns)
    with pytest.raises(ConfigError, match="mutually exclusive"):
        validate_config(cfg)


def test_no_modes_rejected(tmp_path: Path) -> None:
    r1 = tmp_path / "R1.fq.gz"
    r2 = tmp_path / "R2.fq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    ns = _ns("-r1", str(r1), "-r2", str(r2))
    cfg = _config_from_args(ns)
    with pytest.raises(ConfigError, match="No mode selected"):
        validate_config(cfg)


def test_platform_choice() -> None:
    parser = _build_full_parser()
    ns = parser.parse_args(["-r1", "x", "-r2", "y", "--meta", "--platform", "ont"])
    cfg = _config_from_args(ns)
    assert cfg.platform is Platform.ONT


@pytest.mark.parametrize("inp,expected", [
    ("12G", 12),
    ("12g", 12),
    ("32GB", 32),
    ("8", 8),
])
def test_parse_memory(inp: str, expected: int) -> None:
    assert _parse_memory(inp) == expected


def test_parse_memory_default_autodetects() -> None:
    val = _parse_memory(None)
    assert val >= 4
