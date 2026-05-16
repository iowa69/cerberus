"""Meta (assembly) pipeline — conservative, paired-end-preserving.

Stages:
  1. minimap2 vs masked-T2T+HLA: drop pairs where BOTH mates map.
  2. (orphans only) minimap2 single-end vs same reference.
  3. bbduk entropy filter on both paired and orphan streams.
  4. Final outputs: <sample>.meta.R1.fastq.gz / R2 + <sample>.meta.orphans.fastq.gz
"""
from __future__ import annotations

from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.pipelines.base import PipelineResult, stage_dir
from cerberus.refs import RefManager
from cerberus.stages.align import minimap2_paired, minimap2_singles
from cerberus.stages.entropy import entropy_paired, entropy_single
from cerberus.stages.qc import FastpOutputs
from cerberus.stages.concat import concat_gz
from cerberus.utils.fastq import count_reads
from cerberus.utils.logger import get_logger

log = get_logger("pipeline.meta")


def run_meta(
    cfg: CerberusConfig,
    tuned: TunedParams,
    qc: FastpOutputs,
    refs: RefManager,
) -> PipelineResult:
    mode = "meta"
    log.info("=== Running meta pipeline (conservative, paired) ===")
    work = cfg.work_dir / mode
    logs = cfg.logs_dir / mode

    idx = refs.path_to(refs.asset("masked_t2t_hla_minimap2"))

    # 1) paired host removal
    paired_out = minimap2_paired(
        cfg, tuned,
        index=idx,
        r1_in=qc.r1, r2_in=qc.r2,                        # type: ignore[arg-type]
        workdir=stage_dir(work, mode, "01_minimap2_paired"),
        log_dir=logs, tag="01_minimap2_paired",
        drop_strategy="both",
    )

    # 2) orphans: optional
    orphans_clean: Path | None = None
    orphan_inputs = [p for p in (qc.orphans_r1, qc.orphans_r2) if p and p.exists()]
    if orphan_inputs:
        cat_dir = stage_dir(work, mode, "00_orphans")
        merged_orphans = cat_dir / "qc.orphans.fq.gz"
        concat_gz(cfg, inputs=orphan_inputs, output=merged_orphans,
                  log_dir=logs, tag="00_orphans_cat")
        orphan_aln = minimap2_singles(
            cfg, tuned,
            index=idx,
            reads_in=merged_orphans,
            workdir=stage_dir(work, mode, "02_minimap2_orphans"),
            log_dir=logs, tag="02_minimap2_orphans",
        )
        orphans_clean = orphan_aln.long_reads

    # 3) entropy filter
    paired_entropy = entropy_paired(
        cfg, tuned,
        r1_in=paired_out.r1, r2_in=paired_out.r2,        # type: ignore[arg-type]
        workdir=stage_dir(work, mode, "03_entropy_paired"),
        log_dir=logs, tag="03_entropy_paired",
    )

    final_r1 = cfg.out_dir / f"{cfg.sample_id}.meta.R1.fastq.gz"
    final_r2 = cfg.out_dir / f"{cfg.sample_id}.meta.R2.fastq.gz"
    final_r1.unlink(missing_ok=True)
    final_r2.unlink(missing_ok=True)
    if not cfg.dry_run:
        paired_entropy.r1.replace(final_r1)
        paired_entropy.r2.replace(final_r2)                                  # type: ignore[union-attr]

    final_orphans: Path | None = None
    if orphans_clean:
        orph_entropy = entropy_single(
            cfg, tuned,
            reads_in=orphans_clean,
            workdir=stage_dir(work, mode, "04_entropy_orphans"),
            log_dir=logs, tag="04_entropy_orphans",
        )
        final_orphans = cfg.out_dir / f"{cfg.sample_id}.meta.orphans.fastq.gz"
        final_orphans.unlink(missing_ok=True)
        if not cfg.dry_run:
            orph_entropy.r1.replace(final_orphans)

    return PipelineResult(
        mode=mode,
        paired_r1=final_r1,
        paired_r2=final_r2,
        singletons=final_orphans,
        stats=_collect_counts(qc, final_r1, final_orphans, cfg),
    )


def _collect_counts(
    qc: FastpOutputs, r1: Path, orphans: Path | None, cfg: CerberusConfig
) -> dict[str, int]:
    if cfg.dry_run:
        return {}
    counts = {
        "qc_paired": count_reads(qc.r1),
        "final_paired_r1": count_reads(r1),
    }
    if orphans:
        counts["final_orphans"] = count_reads(orphans)
    return counts
