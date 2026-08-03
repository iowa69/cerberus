"""Top-level dispatcher: validate, preflight, prescan, autotune, QC, pipelines, GDPR, report."""
from __future__ import annotations

import shutil
import time
from pathlib import Path

from cerberus.accounting import RunAccounting
from cerberus.autotune import (
    apply_user_overrides,
    autotune_from_prescan,
    estimate_long_read,
    prescan_reads,
)
from cerberus.config import CerberusConfig, TunedParams
from cerberus.pipelines.base import PipelineResult
from cerberus.pipelines.gdpr import run_gdpr_for
from cerberus.pipelines.long_read import run_long_meta, run_long_profiling
from cerberus.pipelines.meta import run_meta
from cerberus.pipelines.profiling import run_profiling
from cerberus.refs import RefManager, cleanup_partial
from cerberus.stages.qc import FastplongOutputs, run_fastp, run_fastplong
from cerberus.utils.fastq import count_reads
from cerberus.utils.logger import get_logger, setup_logging
from cerberus.utils.shell import which

log = get_logger("orchestrator")


class ConfigError(ValueError):
    pass


# Tools each selected mode will actually invoke. Checked up front so a run
# does not spend two hours in minimap2 before discovering bbduk is missing.
_TOOLS_ALWAYS = ("minimap2", "samtools")
_TOOLS_BY_MODE = {
    "qc_short": ("fastp",),
    "meta": ("bbduk.sh",),
    "profiling": ("bbduk.sh",),
    "profiling_full": ("bowtie2",),
    "gdpr": ("kraken2",),
}


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
        if not cfg.long_input.is_file():
            raise ConfigError(f"--long input is not a file: {cfg.long_input}")
        if cfg.r1 or cfg.r2:
            raise ConfigError("--long is incompatible with -r1/-r2.")
    else:
        if cfg.long_input is not None:
            raise ConfigError("-i/--input requires --long. For paired short reads use -r1/-r2.")
        if cfg.r1 is None or cfg.r2 is None:
            raise ConfigError("Short-read mode requires both -r1 and -r2.")
        for label, p in (("-r1", cfg.r1), ("-r2", cfg.r2)):
            if not p.exists():
                raise ConfigError(f"Input not found: {label} {p}")
            if not p.is_file():
                raise ConfigError(f"Input is not a file: {label} {p}")
        if cfg.r1.resolve() == cfg.r2.resolve():
            raise ConfigError("-r1 and -r2 point at the same file.")
    if cfg.fast and cfg.double_pass:
        raise ConfigError("--fast and --double-pass are mutually exclusive.")
    if cfg.threads < 1:
        raise ConfigError("--threads must be >= 1")
    if cfg.memory_gb < 1:
        raise ConfigError("--memory must be at least 1G")
    if cfg.entropy is not None and not 0.0 <= cfg.entropy <= 1.0:
        raise ConfigError(f"--entropy must be between 0.0 and 1.0 (got {cfg.entropy})")
    if cfg.bbduk_k is not None and not 1 <= cfg.bbduk_k <= 31:
        raise ConfigError(f"--bbduk-k must be between 1 and 31 (got {cfg.bbduk_k})")
    if cfg.min_length is not None and cfg.min_length < 0:
        raise ConfigError(f"--min-length must be >= 0 (got {cfg.min_length})")
    if cfg.min_quality is not None and not 0 <= cfg.min_quality <= 93:
        raise ConfigError(f"--min-quality must be between 0 and 93 (got {cfg.min_quality})")
    if not 0.0 <= cfg.gdpr_confidence <= 1.0:
        raise ConfigError(
            f"--gdpr-confidence must be between 0.0 and 1.0 (got {cfg.gdpr_confidence})"
        )
    if "/" in cfg.sample_id or cfg.sample_id in ("", ".", ".."):
        raise ConfigError(f"--sample-id must be a plain name, not a path (got {cfg.sample_id!r})")
    if cfg.out_dir.exists() and not cfg.out_dir.is_dir():
        raise ConfigError(f"--out-dir exists and is not a directory: {cfg.out_dir}")
    if cfg.fast and not cfg.profiling:
        log.warning("--fast only affects --profiling; it has no effect on the selected modes.")
    if cfg.double_pass and not cfg.profiling:
        log.warning("--double-pass only affects --profiling; it has no effect here.")


def preflight_tools(cfg: CerberusConfig) -> list[str]:
    """Return the tools this run needs that are not on PATH."""
    needed = set(_TOOLS_ALWAYS)
    if not cfg.long_mode:
        needed |= set(_TOOLS_BY_MODE["qc_short"])
    if cfg.meta:
        needed |= set(_TOOLS_BY_MODE["meta"])
    if cfg.profiling:
        needed |= set(_TOOLS_BY_MODE["profiling"])
        if not cfg.fast and not cfg.long_mode:
            needed |= set(_TOOLS_BY_MODE["profiling_full"])
    if cfg.gdpr:
        needed |= set(_TOOLS_BY_MODE["gdpr"])
    if cfg.long_mode and not (which("fastplong") or which("chopper")):
        needed.add("fastplong")
    return sorted(t for t in needed if which(t) is None)


def run(cfg: CerberusConfig) -> dict:
    validate_config(cfg)
    cfg.ensure_directories()
    setup_logging(cfg.logs_dir, verbose=cfg.verbose, quiet=cfg.quiet)

    missing = preflight_tools(cfg)
    if missing:
        raise ConfigError(
            f"Missing required tool(s) for this run: {', '.join(missing)}. "
            "Run 'cerberus doctor' for the full picture, or install the environment with "
            "'conda env create -f environment.yml'."
        )

    cleaned_tmp = cleanup_partial(cfg.ref_dir)
    if cleaned_tmp:
        log.info("Cleaned %d leftover partial download(s) from %s", cleaned_tmp, cfg.ref_dir)

    log.info("Cerberus %s: sample=%s threads=%s memory=%dG modes=%s",
             _version(), cfg.sample_id, cfg.threads, cfg.memory_gb,
             ",".join(cfg.modes) + (" +gdpr" if cfg.gdpr else ""))

    refs = RefManager(cfg.ref_dir, auto_download=cfg.auto_download,
                      kraken2_db_override=cfg.kraken2_db_override,
                      aux_refs_override=cfg.aux_refs_override)
    required_keys = _required_pipeline_keys(cfg)
    log.info("Required ref-asset groups: %s", required_keys)
    refs.ensure(refs.required_assets_for(required_keys))

    accounting = RunAccounting(sample_id=cfg.sample_id)
    t0 = time.time()

    if cfg.long_mode:
        result_set = _run_long(cfg, refs, accounting)
    else:
        result_set = _run_short(cfg, refs, accounting)

    gdpr_outputs: dict[str, dict[str, Path | None]] = {}
    if cfg.gdpr:
        for r in result_set:
            gdpr_out = run_gdpr_for(cfg, cfg.tuned, r, refs)
            paths = {
                "paired_r1": gdpr_out.paired_r1,
                "paired_r2": gdpr_out.paired_r2,
                "singletons": gdpr_out.singletons,
                "long_reads": gdpr_out.long_reads,
            }
            for stage, count in gdpr_out.stats.items():
                accounting.add_stage(f"{r.mode}_gdpr", stage, count)
            accounting.add_final(f"{r.mode}_gdpr", paths)
            gdpr_outputs[f"{r.mode}_gdpr"] = paths

    elapsed = time.time() - t0
    accounting.write(cfg.reports_dir)

    outputs: dict[str, Path | None] = {r.mode: r.primary_output for r in result_set}
    for mode, paths in gdpr_outputs.items():
        for key, p in paths.items():
            if p is not None and (cfg.dry_run or p.exists()):
                outputs[f"{mode}.{key}"] = p

    report_path = None
    try:
        from cerberus.report import write_run_report
        report_path = write_run_report(
            cfg, accounting=accounting, results=result_set,
            gdpr_outputs=gdpr_outputs, refs=refs, elapsed_sec=elapsed,
        )
    except Exception as e:  # noqa: BLE001 - a broken report must never fail a good run
        log.warning("Could not write the HTML run report: %s", e)

    if not cfg.keep_intermediates and not cfg.dry_run:
        _clean_work_dir(cfg)

    log.info("Cerberus finished in %.1fs", elapsed)

    return {
        "elapsed_sec": elapsed,
        "outputs": outputs,
        "reports": str(cfg.reports_dir),
        "report_html": str(report_path) if report_path else None,
        "warnings": accounting.warnings,
    }


def _version() -> str:
    from cerberus import __version__
    return __version__


def _clean_work_dir(cfg: CerberusConfig) -> None:
    """Remove intermediates. They can reach many times the input size."""
    work = cfg.work_dir
    if not work.exists():
        return
    try:
        freed = sum(p.stat().st_size for p in work.rglob("*") if p.is_file())
        shutil.rmtree(work)
        log.info("Removed intermediates from %s (%.1f GB freed); "
                 "use --keep-intermediates to retain them",
                 work, freed / 1e9)
    except OSError as e:
        log.warning("Could not clean %s: %s", work, e)


def _required_pipeline_keys(cfg: CerberusConfig) -> list[str]:
    out: list[str] = []
    prefix = "long-" if cfg.long_mode else ""
    if cfg.meta:
        out.append(f"{prefix}meta")
    if cfg.profiling:
        # Long-read profiling has no separate fast lane; it always needs the
        # auxiliary references. Only the short-read path has a reduced set.
        if cfg.fast and not cfg.long_mode:
            out.append("profiling-fast")
        else:
            out.append(f"{prefix}profiling")
    if cfg.gdpr:
        out.append(f"{prefix}gdpr" if cfg.long_mode else "gdpr")
    return out


def _prescan_and_tune(cfg: CerberusConfig, source: Path) -> TunedParams:
    """Sample the input before QC so tuned values can reach fastp."""
    scan = prescan_reads(source)
    if scan.ok:
        log.info("Prescan of %s: %d reads, mean %.0fbp, Q20 %.1f%%, Q30 %.1f%%",
                 source.name, scan.reads_sampled, scan.mean_length,
                 100 * scan.q20_rate, 100 * scan.q30_rate)
    tuned = autotune_from_prescan(
        scan, user_platform=cfg.platform, long_mode=cfg.long_mode,
    )
    tuned = apply_user_overrides(tuned, cfg)
    cfg.tuned = tuned
    cfg.prescan = scan
    return tuned


def _run_short(
    cfg: CerberusConfig,
    refs: RefManager,
    accounting: RunAccounting,
) -> list[PipelineResult]:
    tuned = _prescan_and_tune(cfg, cfg.r1)  # type: ignore[arg-type]

    qc = run_fastp(
        cfg,
        tuned=tuned,
        workdir=cfg.work_dir / "00_qc",
        log_dir=cfg.logs_dir / "00_qc",
    )

    if not cfg.dry_run:
        accounting.input_r1_reads = count_reads(cfg.r1)          # type: ignore[arg-type]
        accounting.input_r2_reads = count_reads(cfg.r2)          # type: ignore[arg-type]
        accounting.qc_paired = count_reads(qc.r1)
        if qc.orphans_r1:
            accounting.qc_orphans = count_reads(qc.orphans_r1) + (
                count_reads(qc.orphans_r2) if qc.orphans_r2 else 0
            )
        if accounting.input_r1_reads and not accounting.qc_paired:
            accounting.warn(
                f"fastp discarded every read ({accounting.input_r1_reads} pairs in, 0 out). "
                "Check --min-length/--min-quality and the input quality encoding."
            )

    results: list[PipelineResult] = []
    if cfg.meta:
        meta_res = run_meta(cfg, tuned, qc, refs)
        _record(accounting, meta_res)
        accounting.add_final("meta", {
            "paired_r1": meta_res.paired_r1,
            "paired_r2": meta_res.paired_r2,
            "orphans": meta_res.singletons,
        })
        results.append(meta_res)
    if cfg.profiling:
        prof_res = run_profiling(cfg, tuned, qc, refs)
        _record(accounting, prof_res)
        accounting.add_final("profiling", {"merged": prof_res.singletons})
        results.append(prof_res)

    return results


def _run_long(
    cfg: CerberusConfig,
    refs: RefManager,
    accounting: RunAccounting,
) -> list[PipelineResult]:
    if not estimate_long_read(cfg.long_input):                   # type: ignore[arg-type]
        log.warning(
            "Input mean read length appears short (<1kb). --long was requested; "
            "consider switching to short-read mode (-r1/-r2)."
        )

    tuned = _prescan_and_tune(cfg, cfg.long_input)                # type: ignore[arg-type]

    qc: FastplongOutputs = run_fastplong(
        cfg, tuned=tuned,
        workdir=cfg.work_dir / "00_qc",
        log_dir=cfg.logs_dir / "00_qc",
    )

    if not cfg.dry_run:
        accounting.input_long_reads = count_reads(cfg.long_input)  # type: ignore[arg-type]
        accounting.qc_long = count_reads(qc.reads)
        if accounting.input_long_reads and not accounting.qc_long:
            accounting.warn(
                f"Long-read QC discarded every read "
                f"({accounting.input_long_reads} in, 0 out). Check --min-length/--min-quality."
            )

    results: list[PipelineResult] = []
    if cfg.meta:
        r = run_long_meta(cfg, tuned, qc.reads, refs)
        _record(accounting, r)
        accounting.add_final(r.mode, {"long_reads": r.long_reads})
        results.append(r)
    if cfg.profiling:
        r = run_long_profiling(cfg, tuned, qc.reads, refs)
        _record(accounting, r)
        accounting.add_final(r.mode, {"long_reads": r.long_reads})
        results.append(r)
    return results


def _record(accounting: RunAccounting, result: PipelineResult) -> None:
    for stage, count in result.stats.items():
        accounting.add_stage(result.mode, stage, count)
    for w in result.warnings:
        accounting.warn(w)
