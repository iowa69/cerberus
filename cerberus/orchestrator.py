"""Top-level dispatcher. Validates config, runs QC, autotunes, runs pipelines, runs GDPR."""
from __future__ import annotations

import time
from pathlib import Path

from cerberus.accounting import RunAccounting
from cerberus.autotune import (
    apply_user_overrides,
    autotune_from_fastp,
    estimate_long_read,
)
from cerberus.config import CerberusConfig, Platform
from cerberus.pipelines.base import PipelineResult
from cerberus.pipelines.gdpr import GDPRResult, run_gdpr_for
from cerberus.pipelines.long_read import run_long_meta, run_long_profiling
from cerberus.pipelines.meta import run_meta
from cerberus.pipelines.profiling import run_profiling
from cerberus.refs import RefManager, cleanup_partial
from cerberus.stages.qc import FastpOutputs, FastplongOutputs, run_fastp, run_fastplong
from cerberus.utils.fastq import count_reads
from cerberus.utils.logger import get_logger, setup_logging

log = get_logger("orchestrator")


class ConfigError(ValueError):
    pass


def validate_config(cfg: CerberusConfig) -> None:
    if not (cfg.meta or cfg.profiling or cfg.gdpr):
        raise ConfigError(
            "No mode selected. Provide at least one of: --meta, --profiling, --all. "
            "--gdpr alone is a post-processor; it needs --meta and/or --profiling to run on."
        )
    if cfg.gdpr and not (cfg.meta or cfg.profiling):
        raise ConfigError(
            "--gdpr requires at least one of --meta or --profiling (it cleans their outputs)."
        )
    if cfg.long_mode:
        if cfg.long_input is None or not cfg.long_input.exists():
            raise ConfigError("--long requires -i FILE pointing to a long-read FASTQ.")
        if cfg.r1 or cfg.r2:
            raise ConfigError("--long is incompatible with -r1/-r2.")
    else:
        if cfg.r1 is None or cfg.r2 is None:
            raise ConfigError("Short-read mode requires both -r1 and -r2.")
        if not cfg.r1.exists() or not cfg.r2.exists():
            raise ConfigError(f"Input not found: {cfg.r1} or {cfg.r2}")
    if cfg.fast and cfg.double_pass:
        raise ConfigError("--fast and --double-pass are mutually exclusive.")
    if cfg.long_mode and cfg.meta and cfg.profiling and cfg.fast:
        log.warning("--fast on long reads only meaningfully applies to --profiling")
    if cfg.threads < 1:
        raise ConfigError("--threads must be >= 1")


def run(cfg: CerberusConfig) -> dict:
    cfg.ensure_directories()
    setup_logging(cfg.logs_dir, verbose=cfg.verbose, quiet=cfg.quiet)
    validate_config(cfg)

    cleaned_tmp = cleanup_partial(cfg.ref_dir)
    if cleaned_tmp:
        log.info("Cleaned %d leftover .tmp file(s) from %s", cleaned_tmp, cfg.ref_dir)

    log.info("Cerberus run: sample=%s threads=%s memory=%dG modes=%s",
             cfg.sample_id, cfg.threads, cfg.memory_gb,
             ",".join(cfg.modes) + (" +gdpr" if cfg.gdpr else ""))

    refs = RefManager(cfg.ref_dir, auto_download=cfg.auto_download)
    required_keys = _required_pipeline_keys(cfg)
    log.info("Required ref-asset groups: %s", required_keys)
    refs.ensure(refs.required_assets_for(required_keys))

    accounting = RunAccounting(sample_id=cfg.sample_id)
    t0 = time.time()

    if cfg.long_mode:
        result_set = _run_long(cfg, refs, accounting)
    else:
        result_set = _run_short(cfg, refs, accounting)

    if cfg.gdpr:
        for r in result_set:
            gdpr_out = run_gdpr_for(cfg, _tuned_or_default(cfg), r, refs)
            accounting.add_final(f"{r.mode}_gdpr", {
                "paired_r1": gdpr_out.paired_r1,
                "paired_r2": gdpr_out.paired_r2,
                "singletons": gdpr_out.singletons,
                "long_reads": gdpr_out.long_reads,
            })

    elapsed = time.time() - t0
    accounting.write(cfg.reports_dir)
    log.info("Cerberus finished in %.1fs", elapsed)

    return {
        "elapsed_sec": elapsed,
        "outputs": {r.mode: r.primary_output for r in result_set},
        "reports": str(cfg.reports_dir),
    }


def _required_pipeline_keys(cfg: CerberusConfig) -> list[str]:
    out: list[str] = []
    prefix = "long-" if cfg.long_mode else ""
    if cfg.meta:
        out.append(f"{prefix}meta")
    if cfg.profiling:
        out.append(f"{prefix}profiling-fast" if cfg.fast else f"{prefix}profiling")
    if cfg.gdpr:
        out.append(f"{prefix}gdpr" if cfg.long_mode else "gdpr")
    return out


def _run_short(
    cfg: CerberusConfig,
    refs: RefManager,
    accounting: RunAccounting,
) -> list[PipelineResult]:
    qc = run_fastp(
        cfg,
        tuned=None,
        workdir=cfg.work_dir / "00_qc",
        log_dir=cfg.logs_dir / "00_qc",
    )
    accounting.input_r1_reads = count_reads(cfg.r1) if not cfg.dry_run else 0  # type: ignore[arg-type]
    accounting.input_r2_reads = count_reads(cfg.r2) if not cfg.dry_run else 0  # type: ignore[arg-type]
    accounting.qc_paired = count_reads(qc.r1) if not cfg.dry_run else 0
    if qc.orphans_r1:
        accounting.qc_orphans = count_reads(qc.orphans_r1) + (
            count_reads(qc.orphans_r2) if qc.orphans_r2 else 0
        )

    tuned = autotune_from_fastp(qc.json_report, user_platform=cfg.platform)
    tuned = apply_user_overrides(tuned, cfg)
    cfg.tuned = tuned

    results: list[PipelineResult] = []
    if cfg.meta:
        meta_res = run_meta(cfg, tuned, qc, refs)
        for stage, count in meta_res.stats.items():
            accounting.add_stage("meta", stage, count)
        accounting.add_final("meta", {
            "paired_r1": meta_res.paired_r1,
            "paired_r2": meta_res.paired_r2,
            "orphans": meta_res.singletons,
        })
        results.append(meta_res)
    if cfg.profiling:
        prof_res = run_profiling(cfg, tuned, qc, refs)
        for stage, count in prof_res.stats.items():
            accounting.add_stage("profiling", stage, count)
        accounting.add_final("profiling", {"merged": prof_res.singletons})
        results.append(prof_res)

    return results


def _run_long(
    cfg: CerberusConfig,
    refs: RefManager,
    accounting: RunAccounting,
) -> list[PipelineResult]:
    if not estimate_long_read(cfg.long_input):                              # type: ignore[arg-type]
        log.warning(
            "Input mean read length appears short (<1kb). --long was requested; "
            "consider switching to short-read mode (-r1/-r2)."
        )

    qc: FastplongOutputs = run_fastplong(
        cfg, tuned=None,
        workdir=cfg.work_dir / "00_qc",
        log_dir=cfg.logs_dir / "00_qc",
    )
    accounting.input_long_reads = count_reads(cfg.long_input) if not cfg.dry_run else 0   # type: ignore[arg-type]
    accounting.qc_paired = count_reads(qc.reads) if not cfg.dry_run else 0

    tuned = autotune_from_fastp(qc.json_report, user_platform=cfg.platform)
    tuned = apply_user_overrides(tuned, cfg)
    cfg.tuned = tuned

    results: list[PipelineResult] = []
    if cfg.meta:
        r = run_long_meta(cfg, tuned, qc.reads, refs)
        for stage, count in r.stats.items():
            accounting.add_stage(r.mode, stage, count)
        accounting.add_final(r.mode, {"long_reads": r.long_reads})
        results.append(r)
    if cfg.profiling:
        r = run_long_profiling(cfg, tuned, qc.reads, refs)
        for stage, count in r.stats.items():
            accounting.add_stage(r.mode, stage, count)
        accounting.add_final(r.mode, {"long_reads": r.long_reads})
        results.append(r)
    return results


def _tuned_or_default(cfg: CerberusConfig):
    """Return cfg.tuned, or a sensible default if QC was skipped (dry-run)."""
    if cfg.tuned and cfg.tuned.read_length_class:
        return cfg.tuned
    from cerberus.config import TunedParams
    return TunedParams()
