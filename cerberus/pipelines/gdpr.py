"""GDPR post-processor.

Consumes the outputs of meta and/or profiling and produces "_GDPR" variants
with zero human reads. Uses two orthogonal mechanisms in series:

  1. Kraken2 against a compact human+mammal DB (k-mer minimizer classifier)
  2. minimap2 alignment against the masked T2T-CHM13v2.0 + HLA reference,
     drop pairs where either mate maps.

The mechanisms are independent (k-mer classifier vs. alignment); a read
survives only if both passes accept it. This is what makes the output
publication-defensible.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.pipelines.base import PipelineResult, stage_dir
from cerberus.refs import RefManager
from cerberus.stages.align import minimap2_paired, minimap2_singles
from cerberus.stages.kraken import kraken2_paired, kraken2_single
from cerberus.utils.logger import get_logger

log = get_logger("pipeline.gdpr")


@dataclass
class GDPRResult:
    source_mode: str
    paired_r1: Path | None = None
    paired_r2: Path | None = None
    singletons: Path | None = None
    long_reads: Path | None = None


def run_gdpr_for(
    cfg: CerberusConfig,
    tuned: TunedParams,
    pipeline_result: PipelineResult,
    refs: RefManager,
) -> GDPRResult:
    log.info("=== Running GDPR scrub on %s output ===", pipeline_result.mode)
    mode_tag = f"gdpr_{pipeline_result.mode}"
    work = cfg.work_dir / mode_tag
    logs = cfg.logs_dir / mode_tag

    kdb_dir = refs.path_to(refs.asset("kraken2_gdpr_compact"))
    kdb = _find_kraken_db(kdb_dir)
    mm2_idx = refs.path_to(refs.asset("masked_t2t_hla_minimap2"))

    result = GDPRResult(source_mode=pipeline_result.mode)

    if pipeline_result.paired_r1 and pipeline_result.paired_r2:
        k_out = kraken2_paired(
            cfg, db=kdb,
            r1_in=pipeline_result.paired_r1,
            r2_in=pipeline_result.paired_r2,
            workdir=stage_dir(work, mode_tag, "01_kraken2"),
            log_dir=logs, tag="01_kraken2",
        )
        # Second orthogonal mechanism: alignment to masked T2T+HLA.
        # Drop pair if EITHER mate maps (most aggressive — this is the GDPR pass).
        mm2_out = minimap2_paired(
            cfg, tuned,
            index=mm2_idx,
            r1_in=k_out.cleaned_r1, r2_in=k_out.cleaned_r2,    # type: ignore[arg-type]
            workdir=stage_dir(work, mode_tag, "02_minimap2_human"),
            log_dir=logs, tag="02_minimap2_human",
            drop_strategy="either",
        )
        final_r1 = cfg.out_dir / f"{cfg.sample_id}.{pipeline_result.mode}.R1_GDPR.fastq.gz"
        final_r2 = cfg.out_dir / f"{cfg.sample_id}.{pipeline_result.mode}.R2_GDPR.fastq.gz"
        final_r1.unlink(missing_ok=True)
        final_r2.unlink(missing_ok=True)
        if not cfg.dry_run:
            mm2_out.r1.replace(final_r1)                       # type: ignore[union-attr]
            mm2_out.r2.replace(final_r2)                       # type: ignore[union-attr]
        result.paired_r1 = final_r1
        result.paired_r2 = final_r2

    if pipeline_result.singletons:
        result.singletons = _gdpr_single(
            cfg, tuned, pipeline_result.singletons, refs,
            workdir=stage_dir(work, mode_tag, "03_singletons"),
            log_dir=logs, tag_prefix="03_single",
            kdb=kdb, mm2_idx=mm2_idx,
            out_name=f"{cfg.sample_id}.{pipeline_result.mode}.GDPR.fastq.gz",
        )

    if pipeline_result.long_reads:
        result.long_reads = _gdpr_single(
            cfg, tuned, pipeline_result.long_reads, refs,
            workdir=stage_dir(work, mode_tag, "04_long"),
            log_dir=logs, tag_prefix="04_long",
            kdb=kdb, mm2_idx=mm2_idx,
            out_name=f"{cfg.sample_id}.{pipeline_result.mode}.long_GDPR.fastq.gz",
        )

    return result


def _gdpr_single(
    cfg: CerberusConfig,
    tuned: TunedParams,
    reads: Path,
    refs: RefManager,
    *,
    workdir: Path,
    log_dir: Path,
    tag_prefix: str,
    kdb: Path,
    mm2_idx: Path,
    out_name: str,
) -> Path:
    k_out = kraken2_single(
        cfg, db=kdb, reads_in=reads,
        workdir=workdir, log_dir=log_dir, tag=f"{tag_prefix}_kraken2",
    )
    mm2_out = minimap2_singles(
        cfg, tuned, index=mm2_idx, reads_in=k_out.cleaned_r1,
        workdir=workdir, log_dir=log_dir, tag=f"{tag_prefix}_minimap2_human",
    )
    final = cfg.out_dir / out_name
    final.unlink(missing_ok=True)
    if not cfg.dry_run:
        mm2_out.long_reads.replace(final)                      # type: ignore[union-attr]
    return final


def _find_kraken_db(db_dir: Path) -> Path:
    """Kraken2 expects a directory; verify hash.k2d is present."""
    if (db_dir / "hash.k2d").exists():
        return db_dir
    for sub in db_dir.iterdir():
        if sub.is_dir() and (sub / "hash.k2d").exists():
            return sub
    raise FileNotFoundError(f"No Kraken2 DB (hash.k2d) found under {db_dir}")
