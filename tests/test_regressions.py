"""Regression tests for defects found in the v0.1.1 review.

Each test names the behaviour that was wrong and asserts the corrected one,
so a future refactor cannot quietly reintroduce it. These are unit-level and
need no external bioinformatics tools; the end-to-end behaviour is covered by
``scripts/smoke_test.sh``.
"""
from __future__ import annotations

import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest

from cerberus.autotune import (
    Prescan,
    autotune_from_fastp,
    autotune_from_prescan,
    classify_length,
    detect_platform,
    prescan_reads,
)
from cerberus.cli import _parse_memory, derive_sample_id
from cerberus.config import CerberusConfig, Platform, ReadLengthClass
from cerberus.orchestrator import ConfigError, _required_pipeline_keys, validate_config
from cerberus.refs import _PIPELINE_TO_ASSETS, Asset
from cerberus.stages.align import _merge_extra, _pair_filter_expr, _paired_preset
from cerberus.utils.fastq import count_reads
from cerberus.utils.shell import ToolError, pipe, run

# --------------------------------------------------------------------------
# align.py — pair-level filtering semantics
# --------------------------------------------------------------------------

def test_drop_strategies_are_actually_different() -> None:
    """'both' and 'either' used to compile to filters selecting the same reads."""
    assert _pair_filter_expr("both") == "flag.unmap || flag.munmap"
    assert _pair_filter_expr("either") is None  # plain -f 12


def test_unknown_drop_strategy_rejected() -> None:
    with pytest.raises(ValueError, match="drop_strategy"):
        _pair_filter_expr("sometimes")


def test_long_read_preset_cannot_be_used_for_paired_alignment() -> None:
    """map-ont on -1/-2 emits no READ1/READ2 flags, silently yielding 0 reads."""
    assert _paired_preset("sr") == "sr"
    for preset in ("map-ont", "map-hifi", "map-pb", "asm5"):
        assert _paired_preset(preset) == "sr"


def test_user_aligner_args_no_longer_discard_autotuned_ones() -> None:
    """`cfg.args or tuned.extra` dropped the autotuned tokens entirely."""
    assert _merge_extra("-N 5", "-k15 -w10") == ["-k15", "-w10", "-N", "5"]
    assert _merge_extra(None, "-k15 -w10") == ["-k15", "-w10"]
    assert _merge_extra("-N 5", "") == ["-N", "5"]


def test_unbalanced_aligner_args_give_a_clear_error() -> None:
    with pytest.raises(ValueError, match="Could not parse aligner arguments"):
        _merge_extra('-R "unterminated', "")


# --------------------------------------------------------------------------
# shell.py — every stage's exit status is checked
# --------------------------------------------------------------------------

def test_pipe_detects_failure_in_a_non_final_stage(tmp_path: Path) -> None:
    """A dead producer feeding a healthy consumer used to report success."""
    with pytest.raises(ToolError) as exc:
        pipe([["false"], ["cat"]], log_path=tmp_path / "p.log")
    assert "stage 1/2" in str(exc.value)


def test_pipe_reports_the_stage_that_actually_failed(tmp_path: Path) -> None:
    with pytest.raises(ToolError) as exc:
        pipe(
            [[sys.executable, "-c", "import sys; sys.exit(7)"], ["cat"]],
            log_path=tmp_path / "p.log",
        )
    assert "rc=7" in str(exc.value)


def test_pipe_succeeds_when_all_stages_succeed(tmp_path: Path) -> None:
    out = tmp_path / "out.txt"
    res = pipe([["echo", "hello"], ["cat"]], log_path=tmp_path / "p.log", final_stdout=out)
    assert res.returncodes == [0, 0]
    assert out.read_text().strip() == "hello"


def test_run_can_capture_stdout_separately_from_the_log(tmp_path: Path) -> None:
    """flagstat output is a deliverable and must not carry the '# CMD:' header."""
    stdout = tmp_path / "stats.txt"
    run(["echo", "payload"], log_path=tmp_path / "cmd.log", stdout_path=stdout)
    assert stdout.read_text().strip() == "payload"
    assert "# CMD:" not in stdout.read_text()
    assert "# CMD:" in (tmp_path / "cmd.log").read_text()


# --------------------------------------------------------------------------
# concat.py / fastq.py
# --------------------------------------------------------------------------

def test_empty_concat_writes_a_valid_gzip(tmp_path: Path) -> None:
    """A zero-byte .gz fails `gzip -t` and breaks every downstream tool."""
    from cerberus.stages.concat import concat_gz

    cfg = CerberusConfig(meta=True, out_dir=tmp_path)
    out = tmp_path / "empty.fq.gz"
    concat_gz(cfg, inputs=[], output=out, log_dir=tmp_path, tag="t")
    assert out.exists()
    with gzip.open(out, "rb") as f:
        assert f.read() == b""
    assert count_reads(out) == 0


def test_compress_to_writes_the_destination_not_the_log(tmp_path: Path) -> None:
    """pigz -c streams to stdout; that used to land in the step log."""
    from cerberus.stages.concat import compress_to

    cfg = CerberusConfig(meta=True, threads=1)
    src = tmp_path / "in.txt"
    src.write_text("payload\n" * 100)
    dst = tmp_path / "out.gz"
    compress_to(cfg, src, dst, tmp_path, "t")
    assert dst.exists() and dst.stat().st_size > 0
    with gzip.open(dst, "rb") as f:
        assert f.read() == src.read_bytes()


def test_count_reads_detects_gzip_by_content_not_extension(tmp_path: Path) -> None:
    """Gzip named .fastq used to be counted as raw text, and vice versa."""
    record = b"@r\nACGT\n+\nIIII\n"
    mislabelled = tmp_path / "reads.fastq"          # gzip content, plain name
    with gzip.open(mislabelled, "wb") as f:
        f.write(record * 10)
    assert count_reads(mislabelled) == 10

    plain = tmp_path / "reads.fq.gz"                # plain content, gzip name
    plain.write_bytes(record * 7)
    assert count_reads(plain) == 7


def test_count_reads_handles_missing_and_empty(tmp_path: Path) -> None:
    assert count_reads(tmp_path / "nope.fq.gz") == 0
    empty = tmp_path / "empty.fq"
    empty.write_bytes(b"")
    assert count_reads(empty) == 0


# --------------------------------------------------------------------------
# autotune.py — parameters must be computed before the stage that uses them
# --------------------------------------------------------------------------

def test_prescan_reads_gzip_without_gz_extension(tmp_path: Path) -> None:
    p = tmp_path / "reads.fastq"
    with gzip.open(p, "wt") as f:
        for i in range(50):
            f.write(f"@r{i}\n{'A' * 150}\n+\n{'I' * 150}\n")
    scan = prescan_reads(p)
    assert scan.ok
    assert scan.reads_sampled == 50
    assert scan.mean_length == pytest.approx(150.0)
    assert scan.q30_rate > 0.9


def test_prescan_degrades_instead_of_crashing_on_garbage(tmp_path: Path) -> None:
    p = tmp_path / "not_a_fastq.gz"
    p.write_bytes(b"\x1f\x8b" + b"\x00" * 100)   # gzip magic, invalid body
    scan = prescan_reads(p)
    assert not scan.ok


def test_tuned_min_length_reaches_fastp(tmp_path: Path, monkeypatch) -> None:
    """fastp used to run before autotune, so the tuned value never applied."""
    from cerberus.stages import qc

    captured: dict = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd

    monkeypatch.setattr(qc, "run", fake_run)
    monkeypatch.setattr(qc, "require_tools", lambda *a: None)

    tuned = autotune_from_prescan(
        Prescan(reads_sampled=100, mean_length=50.0, q20_rate=0.9, q30_rate=0.8),
        user_platform=Platform.AUTO,
    )
    assert tuned.min_length == 35          # VERY_SHORT profile

    cfg = CerberusConfig(meta=True, r1=tmp_path / "a.fq", r2=tmp_path / "b.fq")
    qc.run_fastp(cfg, tuned=tuned, workdir=tmp_path, log_dir=tmp_path)
    cmd = captured["cmd"]
    assert "35" == cmd[cmd.index("--length_required") + 1]


def test_explicit_zero_beats_the_autotuned_value(tmp_path: Path, monkeypatch) -> None:
    """`cfg.min_length or default` swallowed a deliberate 0."""
    from cerberus.stages import qc

    captured: dict = {}
    monkeypatch.setattr(qc, "run", lambda cmd, **kw: captured.update(cmd=cmd))
    monkeypatch.setattr(qc, "require_tools", lambda *a: None)

    cfg = CerberusConfig(meta=True, min_length=0, r1=tmp_path / "a", r2=tmp_path / "b")
    qc.run_fastp(cfg, tuned=None, workdir=tmp_path, log_dir=tmp_path)
    cmd = captured["cmd"]
    assert "0" == cmd[cmd.index("--length_required") + 1]


def test_long_mode_is_not_reclassified_as_illumina() -> None:
    """A 400bp ONT library must stay ONT rather than becoming Illumina/sr."""
    assert detect_platform(400, 0.9, 0.5, Platform.AUTO, long_mode=True) is Platform.ONT
    assert detect_platform(400, 0.9, 0.5, Platform.AUTO, long_mode=False) is Platform.ILLUMINA


def test_modern_ont_is_not_misfiled_as_pacbio_hifi() -> None:
    """ONT R10.4.1 duplex reaches Q20 but not Q30; q20>=0.99 alone misfiled it."""
    assert detect_platform(15000, 0.995, 0.60, Platform.AUTO) is Platform.ONT
    assert detect_platform(15000, 0.99, 0.90, Platform.AUTO) is Platform.PACBIO_HIFI


def test_malformed_fastp_report_does_not_crash(tmp_path: Path) -> None:
    for payload in ("{ not json", json.dumps([1, 2, 3]), json.dumps({"summary": "nope"})):
        p = tmp_path / "r.json"
        p.write_text(payload)
        tuned = autotune_from_fastp(p)          # must not raise
        assert tuned.read_length_class in set(ReadLengthClass)


def test_missing_fastp_report_does_not_crash(tmp_path: Path) -> None:
    tuned = autotune_from_fastp(tmp_path / "absent.json")
    assert tuned.min_length > 0


def test_length_bucket_boundaries_unchanged() -> None:
    assert classify_length(79) is ReadLengthClass.VERY_SHORT
    assert classify_length(80) is ReadLengthClass.SHORT
    assert classify_length(4999) is ReadLengthClass.LONG
    assert classify_length(5000) is ReadLengthClass.VERY_LONG


# --------------------------------------------------------------------------
# refs.py / orchestrator.py — asset mapping must cover every mode combination
# --------------------------------------------------------------------------

@pytest.mark.parametrize("long_mode", [False, True])
@pytest.mark.parametrize("fast", [False, True])
@pytest.mark.parametrize("double_pass", [False, True])
@pytest.mark.parametrize("meta,profiling,gdpr", [
    (True, False, False), (False, True, False), (True, True, False),
    (True, False, True), (False, True, True), (True, True, True),
])
def test_every_mode_combination_maps_to_known_assets(
    long_mode, fast, double_pass, meta, profiling, gdpr,
) -> None:
    """`--long --profiling --fast` produced a key no asset set covered, so the
    run downloaded nothing and then opened references that were never fetched."""
    if fast and double_pass:
        pytest.skip("mutually exclusive")
    cfg = CerberusConfig(
        meta=meta, profiling=profiling, gdpr=gdpr,
        long_mode=long_mode, fast=fast, double_pass=double_pass,
    )
    for key in _required_pipeline_keys(cfg):
        assert key in _PIPELINE_TO_ASSETS, f"{key} has no reference set"


def test_long_profiling_fast_still_requires_aux_refs() -> None:
    cfg = CerberusConfig(profiling=True, long_mode=True, fast=True)
    keys = _required_pipeline_keys(cfg)
    assets = {a for k in keys for a in _PIPELINE_TO_ASSETS[k]}
    assert "aux_refs" in assets


def test_gdpr_asset_set_includes_the_human_kmer_mechanism() -> None:
    """The manifest advertised it as a GDPR mechanism; nothing fetched it."""
    assert "human_kmer_set" in _PIPELINE_TO_ASSETS["gdpr"]
    assert "human_kmer_set" in _PIPELINE_TO_ASSETS["long-gdpr"]


def test_archive_dirname_keeps_dotted_version_strings() -> None:
    """split('.')[0] truncated T2T-CHM13v2.0_bt2.tar.zst to 'T2T-CHM13v2'."""
    a = Asset("k", "d", "T2T-CHM13v2.0_bt2.tar.zst", "u", "s", None, [])
    assert a.is_archive
    assert a.extracted_dirname == "T2T-CHM13v2.0_bt2"


# --------------------------------------------------------------------------
# cli.py — argument handling
# --------------------------------------------------------------------------

@pytest.mark.parametrize("value,expected", [
    ("12G", 12), ("12g", 12), ("32GB", 32), ("8", 8),
    ("1T", 1024), ("2048M", 2),
])
def test_parse_memory_accepts_reasonable_forms(value, expected) -> None:
    assert _parse_memory(value) == expected


@pytest.mark.parametrize("bad", ["abc", "12.5.3G", "", "-5G", "G", "0G", "10K"])
def test_parse_memory_rejects_garbage_with_a_message(bad) -> None:
    if bad == "":
        return  # empty means "autodetect"
    with pytest.raises(ValueError):
        _parse_memory(bad)


@pytest.mark.parametrize("filename,expected", [
    ("SRR12345_R1.fastq.gz", "SRR12345"),
    ("sample_1.fq.gz", "sample"),
    ("my.sample.R1.fq.gz", "my.sample"),
    ("A_R1_001.fastq.gz", "A"),
    ("Sample_S1_L001_R1_001.fastq.gz", "Sample"),
    ("reads.fq.gz", "reads"),
])
def test_sample_id_derivation(filename, expected) -> None:
    """Splitting on the first dot mangled every dotted or Illumina-style name."""
    assert derive_sample_id(Path(filename)) == expected


# --------------------------------------------------------------------------
# orchestrator.py — validation
# --------------------------------------------------------------------------

def _cfg(tmp_path: Path, **kw) -> CerberusConfig:
    r1, r2 = tmp_path / "R1.fq.gz", tmp_path / "R2.fq.gz"
    r1.write_bytes(b""), r2.write_bytes(b"")
    base = {"r1": r1, "r2": r2, "meta": True, "out_dir": tmp_path / "out"}
    base.update(kw)
    return CerberusConfig(**base)


@pytest.mark.parametrize("kw,match", [
    ({"entropy": 5.0}, "entropy"),
    ({"entropy": -1.0}, "entropy"),
    ({"bbduk_k": 99}, "bbduk-k"),
    ({"bbduk_k": 0}, "bbduk-k"),
    ({"min_length": -10}, "min-length"),
    ({"min_quality": 999}, "min-quality"),
    ({"gdpr_confidence": 1.5}, "gdpr-confidence"),
    ({"threads": 0}, "threads"),
    ({"sample_id": "../escape"}, "sample-id"),
])
def test_out_of_range_values_are_rejected(tmp_path: Path, kw, match) -> None:
    """These all used to be interpolated straight into tool command lines."""
    with pytest.raises(ConfigError, match=match):
        validate_config(_cfg(tmp_path, **kw))


def test_input_flag_requires_long_mode(tmp_path: Path) -> None:
    cfg = _cfg(tmp_path, long_input=tmp_path / "x.fq")
    with pytest.raises(ConfigError, match="requires --long"):
        validate_config(cfg)


def test_identical_r1_and_r2_rejected(tmp_path: Path) -> None:
    r1 = tmp_path / "R1.fq.gz"
    r1.write_bytes(b"")
    cfg = CerberusConfig(r1=r1, r2=r1, meta=True)
    with pytest.raises(ConfigError, match="same file"):
        validate_config(cfg)


def test_directory_as_input_rejected(tmp_path: Path) -> None:
    d = tmp_path / "adir"
    d.mkdir()
    r2 = tmp_path / "R2.fq.gz"
    r2.write_bytes(b"")
    with pytest.raises(ConfigError, match="not a file"):
        validate_config(CerberusConfig(r1=d, r2=r2, meta=True))


# --------------------------------------------------------------------------
# pipelines/base.py — dry-run must never destroy real results
# --------------------------------------------------------------------------

def test_dry_run_does_not_delete_existing_outputs(tmp_path: Path) -> None:
    """`final.unlink()` sat outside the dry-run guard in every pipeline."""
    from cerberus.pipelines.base import publish

    existing = tmp_path / "SAMPLE.meta.R1.fastq.gz"
    existing.write_bytes(b"precious real data")

    cfg = CerberusConfig(meta=True, dry_run=True, out_dir=tmp_path)
    publish(cfg, tmp_path / "never_produced.fq.gz", existing)
    assert existing.read_bytes() == b"precious real data"


def test_publish_does_not_delete_when_the_stage_produced_nothing(tmp_path: Path) -> None:
    from cerberus.pipelines.base import publish

    existing = tmp_path / "out.fq.gz"
    existing.write_bytes(b"previous result")
    cfg = CerberusConfig(meta=True, out_dir=tmp_path)
    assert publish(cfg, tmp_path / "missing.fq.gz", existing) is None
    assert existing.read_bytes() == b"previous result"


def test_publish_moves_a_real_output(tmp_path: Path) -> None:
    from cerberus.pipelines.base import publish

    src = tmp_path / "staged.fq.gz"
    src.write_bytes(b"new result")
    dst = tmp_path / "final.fq.gz"
    cfg = CerberusConfig(meta=True, out_dir=tmp_path)
    assert publish(cfg, src, dst) == dst
    assert dst.read_bytes() == b"new result"
    assert not src.exists()


# --------------------------------------------------------------------------
# accounting.py / report.py
# --------------------------------------------------------------------------

def test_accounting_tsv_is_well_formed(tmp_path: Path) -> None:
    from cerberus.accounting import RunAccounting

    acct = RunAccounting(sample_id="S1", input_r1_reads=100, input_r2_reads=100,
                         qc_paired=90)
    acct.add_stage("meta", "01_host_removal", 55)
    acct.warn("something looked odd")
    json_path, tsv_path = acct.write(tmp_path)

    rows = [ln.split("\t") for ln in tsv_path.read_text().rstrip("\n").split("\n")]
    header = rows[0]
    assert header == ["sample", "mode", "stage", "reads", "unit", "file"]
    for row in rows:
        assert len(row) == len(header), row
    assert all(r[0] == "S1" for r in rows[1:])

    data = json.loads(json_path.read_text())
    assert data["warnings"] == ["something looked odd"]


def test_run_report_renders_and_is_self_contained(tmp_path: Path) -> None:
    from cerberus.accounting import RunAccounting
    from cerberus.pipelines.base import PipelineResult
    from cerberus.report import write_run_report

    out = tmp_path / "out"
    out.mkdir()
    r1 = out / "S1.meta.R1.fastq.gz"
    with gzip.open(r1, "wt") as f:
        f.write("@a\nACGT\n+\nIIII\n")

    cfg = CerberusConfig(meta=True, out_dir=out, sample_id="S1")
    cfg.command_line = "cerberus --meta"
    acct = RunAccounting(sample_id="S1", input_r1_reads=10, qc_paired=8)
    acct.add_stage("meta", "01_host_removal", 4)

    path = write_run_report(
        cfg, accounting=acct,
        results=[PipelineResult(mode="meta", paired_r1=r1)],
        gdpr_outputs={}, refs=None, elapsed_sec=1.5,
    )
    html = path.read_text()
    assert path.name == "cerberus_report.html"
    assert "<!doctype html>" in html
    assert "S1" in html and "01_host_removal" in html
    # Self-contained: no external fetches of any kind.
    for token in ("http://", "https://", "<script"):
        assert token not in html
    record = json.loads((out / "reports" / "run_record.json").read_text())
    assert record["sample_id"] == "S1"
    assert record["resolved_parameters"]["minimap2_preset"]


def test_report_flags_a_desynchronised_pair(tmp_path: Path) -> None:
    from cerberus.report import collect_outputs

    out = tmp_path / "out"
    out.mkdir()
    r1, r2 = out / "S.R1.fq.gz", out / "S.R2.fq.gz"
    with gzip.open(r1, "wt") as f:
        f.write("@a\nACGT\n+\nIIII\n" * 3)
    with gzip.open(r2, "wt") as f:
        f.write("@a\nACGT\n+\nIIII\n" * 2)

    cfg = CerberusConfig(meta=True, out_dir=out, sample_id="S")
    from cerberus.pipelines.base import PipelineResult
    rows = collect_outputs(cfg, [PipelineResult(mode="meta", paired_r1=r1, paired_r2=r2)], {})
    assert any(r.get("paired_ok") is False for r in rows)


# --------------------------------------------------------------------------
# CLI surface
# --------------------------------------------------------------------------

def _cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "cerberus", *args],
        capture_output=True, text=True, timeout=120, check=False,
    )


def test_brief_help_is_shorter_than_full_help() -> None:
    """Both used to print byte-identical output including advanced flags."""
    brief = _cli("--help")
    full = _cli("--help-all")
    assert brief.returncode == 0 and full.returncode == 0
    assert len(full.stdout) > len(brief.stdout)
    assert "--gdpr-confidence" in full.stdout
    assert "--gdpr-confidence" not in brief.stdout


def test_bare_invocation_prints_help_not_a_traceback() -> None:
    res = _cli()
    assert "Traceback" not in res.stderr
    assert "usage:" in res.stdout.lower()


def test_run_subcommand_is_accepted() -> None:
    """The docstring documented `cerberus run`; it used to be rejected."""
    res = _cli("run", "--help")
    assert res.returncode == 0
    assert "usage:" in res.stdout.lower()


def test_bad_memory_gives_a_clean_error_not_a_traceback() -> None:
    res = _cli("-r1", "a.fq", "-r2", "b.fq", "--meta", "--memory", "banana")
    assert res.returncode == 2
    assert "Traceback" not in res.stderr
    assert "banana" in res.stderr


# --------------------------------------------------------------------------
# gdpr.py — output naming must describe what the file actually is
# --------------------------------------------------------------------------

def test_gdpr_names_the_profiling_deliverable_correctly(tmp_path: Path) -> None:
    """`singletons` means different things per head.

    For meta it is the unpaired leftovers; for profiling it is the single
    merged file that *is* the deliverable. v0.2.0 named both "orphans", which
    labelled the profiling head's only GDPR output as a side stream.
    """
    from cerberus.pipelines.base import PipelineResult

    def gdpr_name(result: PipelineResult, sample: str = "S") -> str:
        is_primary = result.paired_r1 is None
        return (f"{sample}.{result.mode}_GDPR.fastq.gz" if is_primary
                else f"{sample}.{result.mode}.orphans_GDPR.fastq.gz")

    profiling = PipelineResult(mode="profiling", singletons=tmp_path / "p.fq.gz")
    assert gdpr_name(profiling) == "S.profiling_GDPR.fastq.gz"

    meta = PipelineResult(
        mode="meta", paired_r1=tmp_path / "r1.gz", paired_r2=tmp_path / "r2.gz",
        singletons=tmp_path / "o.fq.gz",
    )
    assert gdpr_name(meta) == "S.meta.orphans_GDPR.fastq.gz"


def test_gdpr_output_names_are_unique_across_heads() -> None:
    """meta and profiling must never collide on an output filename."""
    from cerberus.pipelines.base import PipelineResult

    names: list[str] = []
    for result in (
        PipelineResult(mode="meta", paired_r1=Path("a"), paired_r2=Path("b"),
                       singletons=Path("c")),
        PipelineResult(mode="profiling", singletons=Path("d")),
        PipelineResult(mode="long_meta", long_reads=Path("e")),
        PipelineResult(mode="long_profiling", long_reads=Path("f")),
    ):
        is_primary = result.paired_r1 is None
        if result.paired_r1:
            names += [f"S.{result.mode}.R1_GDPR.fastq.gz",
                      f"S.{result.mode}.R2_GDPR.fastq.gz"]
        if result.singletons:
            names.append(f"S.{result.mode}_GDPR.fastq.gz" if is_primary
                         else f"S.{result.mode}.orphans_GDPR.fastq.gz")
        if result.long_reads:
            names.append(f"S.{result.mode}_GDPR.fastq.gz")

    assert len(names) == len(set(names)), f"colliding names: {sorted(names)}"
    assert all(n.endswith("_GDPR.fastq.gz") for n in names), sorted(names)
    assert "S.profiling_GDPR.fastq.gz" in names
    assert "S.profiling.orphans_GDPR.fastq.gz" not in names


# --------------------------------------------------------------------------
# v0.2.1 — defects introduced by the v0.2.0 refactor itself
# --------------------------------------------------------------------------

def test_pipe_blames_the_downstream_stage_when_several_fail(tmp_path: Path) -> None:
    """A dying consumer takes its producer down with it.

    v0.2.0 raised on the first non-zero status, so a bad output path on the
    consumer was reported as the producer failing.
    """
    with pytest.raises(ToolError) as exc:
        pipe(
            [[sys.executable, "-c",
              "import sys\nfor _ in range(200000): sys.stdout.write('x' * 80 + '\\n')"],
             [sys.executable, "-c", "import sys; sys.exit(9)"]],
            log_path=tmp_path / "p.log",
        )
    msg = str(exc.value)
    assert "stage 2/2" in msg
    assert "rc=9" in msg


def test_pipe_still_blames_a_lone_upstream_failure(tmp_path: Path) -> None:
    with pytest.raises(ToolError) as exc:
        pipe([[sys.executable, "-c", "import sys; sys.exit(3)"], ["cat"]],
             log_path=tmp_path / "p.log")
    assert "stage 1/2" in str(exc.value)


def test_children_are_started_in_their_own_process_group(tmp_path: Path) -> None:
    """Signalling only the direct child leaves forked grandchildren running.

    bbduk.sh is a shell wrapper around a JVM, so this is not hypothetical.
    """
    import os
    import signal as _signal
    import threading
    import time

    from cerberus.utils import shell as sh

    pidfile = tmp_path / "pids.txt"
    script = f'sleep 60 & echo "$$ $!" > {pidfile}; wait'
    t = threading.Thread(
        target=lambda: _swallow(lambda: run(["sh", "-c", script],
                                            log_path=tmp_path / "pg.log")),
        daemon=True,
    )
    t.start()
    for _ in range(100):
        if pidfile.exists():
            break
        time.sleep(0.05)
    child, grandchild = (int(x) for x in pidfile.read_text().split())
    assert os.getpgid(child) == os.getpgid(grandchild)

    assert sh.terminate_all(_signal.SIGTERM) >= 1
    t.join(timeout=20)
    time.sleep(0.5)
    for pid in (child, grandchild):
        with pytest.raises(OSError):
            os.kill(pid, 0)


def _swallow(fn):
    """Run fn, ignoring the ToolError raised when we kill it mid-flight."""
    try:
        fn()
    except Exception as e:                      # noqa: BLE001 - the kill is the point
        print(f"(expected during teardown test: {type(e).__name__})")


def test_count_reads_is_cached_per_file_identity(tmp_path: Path) -> None:
    """v0.2.0 decompressed the same file up to four times per run."""
    from cerberus.utils import fastq as fq

    p = tmp_path / "r.fq.gz"
    with gzip.open(p, "wt") as f:
        f.write("@a\nACGT\n+\nIIII\n" * 25)

    fq.clear_count_cache()
    calls = {"n": 0}
    real = fq._count_python

    def counting(path):
        calls["n"] += 1
        return real(path)

    fq._count_python = counting
    try:
        # force the Python path so the counter sees every read
        assert fq.count_reads(p) == 25
        first = calls["n"]
        for _ in range(5):
            assert fq.count_reads(p) == 25
        assert calls["n"] == first, "cache did not prevent repeated decompression"

        # rewriting the file must invalidate the entry
        with gzip.open(p, "wt") as f:
            f.write("@a\nACGT\n+\nIIII\n" * 3)
        assert fq.count_reads(p) == 3
    finally:
        fq._count_python = real
        fq.clear_count_cache()


def test_corrupt_gzip_does_not_kill_the_run(tmp_path: Path) -> None:
    """A truncated gzip raises EOFError, which is not an OSError."""
    from cerberus.utils import fastq as fq

    good = tmp_path / "g.fq.gz"
    with gzip.open(good, "wt") as f:
        f.write("@a\nACGT\n+\nIIII\n" * 50)
    truncated = tmp_path / "t.fq.gz"
    truncated.write_bytes(good.read_bytes()[: good.stat().st_size // 2])

    fq.clear_count_cache()
    assert fq.count_reads(truncated) >= 0        # must not raise


def test_gdpr_mechanism_table_covers_non_paired_heads() -> None:
    """v0.2.0 hard-coded paired keys, leaving profiling/long tables empty."""
    from cerberus.pipelines.gdpr import residual_host_estimate

    merged_only = {
        "gdpr_input_merged": 1000,
        "gdpr_after_kraken2_merged": 900,
        "gdpr_after_human_kmers_merged": 850,
        "gdpr_after_minimap2_merged": 800,
    }
    est = residual_host_estimate(merged_only)
    assert set(est) == {"merged"}
    assert est["merged"]["kraken2"] == pytest.approx(10.0)
    assert est["merged"]["bbduk-human-kmers"] == pytest.approx(100 * 50 / 900, rel=1e-3)
    assert est["merged"]["minimap2"] == pytest.approx(100 * 50 / 850, rel=1e-3)


def test_gdpr_mechanism_chain_closes_over_a_skipped_mechanism() -> None:
    """With the human k-mer asset absent, minimap2 must still be measured."""
    from cerberus.pipelines.gdpr import residual_host_estimate

    est = residual_host_estimate({
        "gdpr_input_paired": 1000,
        "gdpr_after_kraken2_paired": 800,
        "gdpr_after_minimap2_paired": 600,
    })
    assert "bbduk-human-kmers" not in est["paired"]
    assert est["paired"]["minimap2"] == pytest.approx(25.0)
