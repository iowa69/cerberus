"""Quality control wrappers around fastp (short reads) and fastplong (long reads)."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import require_tools, run

log = get_logger("qc")


@dataclass
class FastpOutputs:
    r1: Path
    r2: Path | None
    orphans_r1: Path | None
    orphans_r2: Path | None
    json_report: Path
    html_report: Path


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

    min_len = (tuned.min_length if tuned else None) or cfg.min_length or 50
    min_qual = (tuned.min_quality if tuned else None) or cfg.min_quality or 20

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

    return FastpOutputs(
        r1=out_r1,
        r2=out_r2,
        orphans_r1=orphan_r1 if orphan_r1.exists() else None,
        orphans_r2=orphan_r2 if orphan_r2.exists() else None,
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

    min_len = (tuned.min_length if tuned else None) or cfg.min_length or 200
    min_qual = (tuned.min_quality if tuned else None) or cfg.min_quality or 10

    from cerberus.utils.shell import which
    if which("fastplong"):
        require_tools("fastplong")
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
    elif which("chopper"):
        log.warning("fastplong not found; using chopper. JSON report will be minimal.")
        from cerberus.utils.shell import pipe
        pipe(
            [
                ["zcat", str(cfg.long_input)] if str(cfg.long_input).endswith(".gz")
                else ["cat", str(cfg.long_input)],
                ["chopper", "-q", str(min_qual), "-l", str(min_len), "-t", str(cfg.threads)],
                ["pigz", "-c"],
            ],
            log_path=log_dir / "chopper.log",
            final_stdout=out_reads,
            dry_run=cfg.dry_run,
        )
        json_report.write_text(_synthetic_long_json(out_reads, min_len, min_qual))
    else:
        raise RuntimeError("Neither fastplong nor chopper found. Install via conda env.")

    return FastplongOutputs(
        reads=out_reads,
        json_report=json_report,
        html_report=html_report,
    )


def _synthetic_long_json(reads: Path, min_len: int, min_qual: int) -> str:
    """Minimal fastp-compatible JSON for autotune when only chopper is available."""
    import json
    return json.dumps({
        "summary": {
            "before_filtering": {
                "read1_mean_length": 5000,
                "q20_rate": 0.95,
                "q30_rate": 0.7,
            },
            "after_filtering": {
                "read1_mean_length": 5000,
            },
        },
        "filtering_result": {
            "passed_filter_reads": 0,
        },
        "notes": [
            f"Synthetic JSON; chopper used with min_len={min_len}, min_qual={min_qual}",
            f"Reads: {reads.name}",
        ],
    })
