"""Quality control wrappers around fastp (short reads) and fastplong (long reads)."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import pipe, require_tools, run, which

log = get_logger("qc")


@dataclass
class FastpOutputs:
    r1: Path
    r2: Path | None
    orphans_r1: Path | None
    orphans_r2: Path | None
    json_report: Path
    html_report: Path


def _resolve(tuned: TunedParams | None, cfg_value, tuned_attr: str, fallback):
    """Precedence: explicit CLI value > autotuned value > hard fallback.

    ``cfg_value`` of 0 is a legitimate user choice, so this tests for None
    rather than using ``or``.
    """
    if cfg_value is not None:
        return cfg_value
    if tuned is not None:
        return getattr(tuned, tuned_attr)
    return fallback


def run_fastp(
    cfg: CerberusConfig,
    *,
    tuned: TunedParams | None = None,
    workdir: Path,
    log_dir: Path,
) -> FastpOutputs:
    """Run fastp on paired-end short reads. Splits orphans into separate files."""
    require_tools("fastp")
    workdir.mkdir(parents=True, exist_ok=True)

    out_r1 = workdir / "qc.R1.fq.gz"
    out_r2 = workdir / "qc.R2.fq.gz"
    orphan_r1 = workdir / "qc.unpaired_R1.fq.gz"
    orphan_r2 = workdir / "qc.unpaired_R2.fq.gz"
    json_report = workdir / "fastp.json"
    html_report = workdir / "fastp.html"

    min_len = _resolve(tuned, cfg.min_length, "min_length", 50)
    min_qual = _resolve(tuned, cfg.min_quality, "min_quality", 20)

    cmd = [
        "fastp",
        "-i", str(cfg.r1),
        "-I", str(cfg.r2),
        "-o", str(out_r1),
        "-O", str(out_r2),
        "--unpaired1", str(orphan_r1),
        "--unpaired2", str(orphan_r2),
        "--length_required", str(min_len),
        "--qualified_quality_phred", str(min_qual),
        "--trim_poly_g",
        "--trim_poly_x",
        "--detect_adapter_for_pe",
        "--correction",
        "--json", str(json_report),
        "--html", str(html_report),
        "--report_title", f"Cerberus QC — {cfg.sample_id}",
        "--thread", str(min(cfg.threads, 16)),
    ]

    run(cmd, log_path=log_dir / "fastp.log", dry_run=cfg.dry_run)

    if cfg.dry_run and not json_report.exists():
        json_report.write_text(_placeholder_json(tuned, min_len, min_qual))

    return FastpOutputs(
        r1=out_r1,
        r2=out_r2,
        orphans_r1=orphan_r1 if _has_reads(orphan_r1) else None,
        orphans_r2=orphan_r2 if _has_reads(orphan_r2) else None,
        json_report=json_report,
        html_report=html_report,
    )


@dataclass
class FastplongOutputs:
    reads: Path
    json_report: Path
    html_report: Path


def run_fastplong(
    cfg: CerberusConfig,
    *,
    tuned: TunedParams | None = None,
    workdir: Path,
    log_dir: Path,
) -> FastplongOutputs:
    """Run fastplong (or fall back to chopper) on long reads."""
    workdir.mkdir(parents=True, exist_ok=True)

    out_reads = workdir / "qc.long.fq.gz"
    json_report = workdir / "fastplong.json"
    html_report = workdir / "fastplong.html"

    min_len = _resolve(tuned, cfg.min_length, "min_length", 200)
    min_qual = _resolve(tuned, cfg.min_quality, "min_quality", 10)

    if which("fastplong"):
        cmd = [
            "fastplong",
            "-i", str(cfg.long_input),
            "-o", str(out_reads),
            "--length_required", str(min_len),
            "--mean_qual", str(min_qual),
            "--json", str(json_report),
            "--html", str(html_report),
            "--thread", str(min(cfg.threads, 8)),
        ]
        run(cmd, log_path=log_dir / "fastplong.log", dry_run=cfg.dry_run)
        if cfg.dry_run and not json_report.exists():
            json_report.write_text(_placeholder_json(tuned, min_len, min_qual))
    elif which("chopper"):
        log.warning(
            "fastplong not found; using chopper. Its report carries no per-read "
            "statistics, so autotuning falls back to the pre-QC prescan of the input."
        )
        decompress = ["cat", str(cfg.long_input)]
        if _is_gzip(cfg.long_input):
            decompress = ["pigz", "-dc", str(cfg.long_input)] if which("pigz") \
                else ["gzip", "-dc", str(cfg.long_input)]
        pipe(
            [
                decompress,
                ["chopper", "-q", str(min_qual), "-l", str(min_len),
                 "-t", str(cfg.threads)],
                ["pigz", "-c"] if which("pigz") else ["gzip", "-c"],
            ],
            log_path=log_dir / "chopper.log",
            final_stdout=out_reads,
            dry_run=cfg.dry_run,
        )
        json_report.write_text(_placeholder_json(tuned, min_len, min_qual, tool="chopper"))
    else:
        raise RuntimeError(
            "Neither fastplong nor chopper found — long-read QC cannot run. "
            "Install with: conda install -c bioconda fastplong"
        )

    return FastplongOutputs(
        reads=out_reads,
        json_report=json_report,
        html_report=html_report,
    )


def _is_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return str(path).endswith(".gz")


def _has_reads(path: Path) -> bool:
    """True when the file exists and holds at least one record.

    fastp always creates the ``--unpaired`` files, so a bare ``.exists()``
    test would send every run through a full alignment pass over nothing.
    A gzip member with no payload is ~20-30 bytes.
    """
    try:
        return path.exists() and path.stat().st_size > 40
    except OSError:
        return False


def _placeholder_json(
    tuned: TunedParams | None, min_len: int, min_qual: int, tool: str = "dry-run",
) -> str:
    """fastp-shaped stand-in used when no real report exists.

    It deliberately carries **no** read-length figure: autotune treats a
    missing length as "fall back to the prescan / defaults" rather than
    inventing a number. The previous version hard-coded 5000 bp here, which
    silently forced every chopper-based run into the VERY_LONG profile.
    """
    return json.dumps({
        "summary": {"before_filtering": {}, "after_filtering": {}},
        "filtering_result": {},
        "cerberus_placeholder": True,
        "notes": [
            f"Synthetic report written by Cerberus ({tool}).",
            f"Filters applied: min_len={min_len}, min_qual={min_qual}.",
            ("Read-length statistics are unavailable from this tool; "
             "autotuning used the pre-QC prescan instead."),
        ],
    }, indent=2)
