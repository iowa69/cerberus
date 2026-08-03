"""Meta (assembly) pipeline — conservative, paired-end-preserving.

Stages:
  1. minimap2 vs masked-T2T+HLA with drop_strategy="both": a pair is dropped
     only when BOTH mates map to the host. A pair with exactly one mapping
     mate is kept intact, preserving the microbial mate and the pairing an
     assembler needs, at the cost of some residual host.
  2. (orphans) minimap2 single-end vs the same reference. fastp's unpaired
     reads and the singletons from step 1 are pooled here.
  3. bbduk entropy filter on both the paired and the singleton streams.
  4. Final outputs:
       <sample>.meta.R1.fastq.gz / .R2.fastq.gz
       <sample>.meta.orphans.fastq.gz
"""
from __future__ import annotations

from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.pipelines.base import PipelineResult, StageTracker, publish, stage_dir
from cerberus.refs import RefManager
from cerberus.stages.align import minimap2_paired, minimap2_singles
from cerberus.stages.concat import concat_gz
from cerberus.stages.entropy import entropy_paired, entropy_single
from cerberus.stages.qc import FastpOutputs
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
    track = StageTracker(cfg)

    idx = refs.path_to(refs.asset("masked_t2t_hla_minimap2"))
    qc_in = track.record("00_qc_paired", qc.r1, qc.r2)

    # 1) paired host removal — conservative: drop the pair only if both mates map
    paired_out = minimap2_paired(
        cfg, tuned,
        index=idx,
        r1_in=qc.r1, r2_in=qc.r2,                        # type: ignore[arg-type]
        workdir=stage_dir(work, "01_minimap2_paired"),
        log_dir=logs, tag="01_minimap2_paired",
        drop_strategy="both",
        keep_singletons=True,
    )
    after_align = track.record("01_host_removal_paired", paired_out.r1, paired_out.r2)
    track.check_not_empty("01_host_removal_paired", qc_in, after_align)

    # 2) singleton stream: fastp orphans plus the mates rescued in step 1
    singleton_inputs = [
        p for p in (qc.orphans_r1, qc.orphans_r2, paired_out.singletons)
        if p and p.exists() and p.stat().st_size > 40
    ]
    singles_clean: Path | None = None
    if singleton_inputs:
        merged = stage_dir(work, "00_singletons") / "singletons.fq.gz"
        concat_gz(cfg, inputs=singleton_inputs, output=merged,
                  log_dir=logs, tag="00_singletons_cat")
        track.record("00_singletons_in", merged)
        aln = minimap2_singles(
            cfg, tuned,
            index=idx,
            reads_in=merged,
            workdir=stage_dir(work, "02_minimap2_singletons"),
            log_dir=logs, tag="02_minimap2_singletons",
        )
        singles_clean = aln.long_reads
        track.record("02_host_removal_singletons", singles_clean)
    else:
        log.info("No unpaired reads to process for meta.")

    # 3) entropy filter
    paired_entropy = entropy_paired(
        cfg, tuned,
        r1_in=paired_out.r1, r2_in=paired_out.r2,        # type: ignore[arg-type]
        workdir=stage_dir(work, "03_entropy_paired"),
        log_dir=logs, tag="03_entropy_paired",
    )
    track.record("03_entropy_paired", paired_entropy.r1, paired_entropy.r2)

    final_r1 = publish(cfg, paired_entropy.r1, cfg.out_dir / f"{cfg.sample_id}.meta.R1.fastq.gz")
    final_r2 = publish(cfg, paired_entropy.r2, cfg.out_dir / f"{cfg.sample_id}.meta.R2.fastq.gz")

    final_singles: Path | None = None
    if singles_clean:
        ent = entropy_single(
            cfg, tuned,
            reads_in=singles_clean,
            workdir=stage_dir(work, "04_entropy_singletons"),
            log_dir=logs, tag="04_entropy_singletons",
        )
        track.record("04_entropy_singletons", ent.r1)
        final_singles = publish(
            cfg, ent.r1, cfg.out_dir / f"{cfg.sample_id}.meta.orphans.fastq.gz"
        )

    if not cfg.dry_run:
        track.stats["final_paired_per_mate"] = count_reads(final_r1) if final_r1 else 0
        if final_singles:
            track.stats["final_orphans"] = count_reads(final_singles)

    return PipelineResult(
        mode=mode,
        paired_r1=final_r1,
        paired_r2=final_r2,
        singletons=final_singles,
        stats=track.stats,
        warnings=track.warnings,
    )
