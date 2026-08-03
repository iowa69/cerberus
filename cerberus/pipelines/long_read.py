"""Long-read variants of the pipelines (--long flag).

``--fast`` skips the auxiliary k-mer pass here, exactly as it does for short
reads. Earlier versions ignored the flag in this module while the orchestrator
still asked for a reduced reference set, so a ``--long --profiling --fast``
run downloaded nothing and then opened references that were never fetched.
"""
from __future__ import annotations

from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.pipelines.base import PipelineResult, StageTracker, publish, stage_dir
from cerberus.refs import RefManager
from cerberus.stages.align import minimap2_singles, winnowmap_singles
from cerberus.stages.entropy import entropy_single
from cerberus.stages.kmer import bbduk_kmer_single
from cerberus.utils.logger import get_logger

log = get_logger("pipeline.long")


def run_long_meta(
    cfg: CerberusConfig,
    tuned: TunedParams,
    reads: Path,
    refs: RefManager,
) -> PipelineResult:
    mode = "long_meta"
    log.info("=== Running long-read meta (conservative) ===")
    work = cfg.work_dir / mode
    logs = cfg.logs_dir / mode
    track = StageTracker(cfg)

    idx = refs.path_to(refs.asset("masked_t2t_hla_minimap2"))
    qc_in = track.record("00_qc_long", reads)

    aln = minimap2_singles(
        cfg, tuned, index=idx, reads_in=reads,
        workdir=stage_dir(work, "01_minimap2"),
        log_dir=logs, tag="01_minimap2_long_meta",
    )
    after = track.record("01_host_removal", aln.long_reads)
    track.check_not_empty("01_host_removal", qc_in, after)

    ent = entropy_single(
        cfg, tuned, reads_in=aln.long_reads,                 # type: ignore[arg-type]
        workdir=stage_dir(work, "02_entropy"),
        log_dir=logs, tag="02_entropy_long_meta",
    )
    track.record("02_entropy", ent.r1)

    final = publish(cfg, ent.r1, cfg.out_dir / f"{cfg.sample_id}.long_meta.fastq.gz")
    track.record("03_final", final)
    return PipelineResult(
        mode=mode, long_reads=final,
        stats=track.stats, warnings=track.warnings,
    )


def run_long_profiling(
    cfg: CerberusConfig,
    tuned: TunedParams,
    reads: Path,
    refs: RefManager,
) -> PipelineResult:
    mode = "long_profiling"
    log.info("=== Running long-read profiling (aggressive) ===")
    work = cfg.work_dir / mode
    logs = cfg.logs_dir / mode
    track = StageTracker(cfg)

    idx = refs.path_to(refs.asset("masked_t2t_hla_minimap2"))
    qc_in = track.record("00_qc_long", reads)

    # 1) primary aligner
    if cfg.double_pass and tuned.winnowmap_enabled:
        log.info("--double-pass + very-long reads: using winnowmap")
        primary = winnowmap_singles(
            cfg, tuned, index=idx,
            meryl_db=idx.with_suffix(".repetitive_k15.txt"),
            reads_in=reads,
            workdir=stage_dir(work, "01_winnowmap"),
            log_dir=logs, tag="01_winnowmap_long_prof",
        )
    else:
        primary = minimap2_singles(
            cfg, tuned, index=idx, reads_in=reads,
            workdir=stage_dir(work, "01_minimap2"),
            log_dir=logs, tag="01_minimap2_long_prof",
        )
    after = track.record("01_host_removal", primary.long_reads)
    track.check_not_empty("01_host_removal", qc_in, after)

    # 2) auxiliary k-mer pass (skipped under --fast, matching the short-read path)
    current = primary.long_reads
    if tuned.bbduk_aux_enabled and not cfg.fast:
        k = bbduk_kmer_single(
            cfg, tuned, ref=refs.aux_refs_path(),
            reads_in=current,                                # type: ignore[arg-type]
            workdir=stage_dir(work, "02_bbduk_kmer"),
            log_dir=logs, tag="02_bbduk_kmer_long_prof",
        )
        current = k.r1
        track.record("02_aux_kmer", current)
    elif cfg.fast:
        log.info("--fast: skipping the auxiliary k-mer pass")

    # 3) entropy
    ent = entropy_single(
        cfg, tuned, reads_in=current,                        # type: ignore[arg-type]
        workdir=stage_dir(work, "03_entropy"),
        log_dir=logs, tag="03_entropy_long_prof",
    )
    track.record("03_entropy", ent.r1)

    final = publish(cfg, ent.r1, cfg.out_dir / f"{cfg.sample_id}.long_profiling.fastq.gz")
    track.record("04_final", final)
    return PipelineResult(
        mode=mode, long_reads=final,
        stats=track.stats, warnings=track.warnings,
    )
