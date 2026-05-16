"""Long-read variants of the three pipelines (--long flag)."""
from __future__ import annotations

from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.pipelines.base import PipelineResult, stage_dir
from cerberus.refs import RefManager
from cerberus.stages.align import minimap2_singles, winnowmap_singles
from cerberus.stages.entropy import entropy_single
from cerberus.stages.kmer import bbduk_kmer_single
from cerberus.utils.fastq import count_reads
from cerberus.utils.logger import get_logger

log = get_logger("pipeline.long")


def run_long_meta(
    cfg: CerberusConfig,
    tuned: TunedParams,
    reads: Path,
    refs: RefManager,
) -> PipelineResult:
    mode = "long-meta"
    log.info("=== Running long-read meta (conservative) ===")
    work = cfg.work_dir / mode
    logs = cfg.logs_dir / mode

    idx = refs.path_to(refs.asset("masked_t2t_hla_minimap2"))
    aln = minimap2_singles(
        cfg, tuned, index=idx, reads_in=reads,
        workdir=stage_dir(work, mode, "01_minimap2"),
        log_dir=logs, tag="01_minimap2_long_meta",
    )
    ent = entropy_single(
        cfg, tuned, reads_in=aln.long_reads,                 # type: ignore[arg-type]
        workdir=stage_dir(work, mode, "02_entropy"),
        log_dir=logs, tag="02_entropy_long_meta",
    )
    final = cfg.out_dir / f"{cfg.sample_id}.long_meta.fastq.gz"
    final.unlink(missing_ok=True)
    if not cfg.dry_run:
        ent.r1.replace(final)
    return PipelineResult(
        mode=mode, long_reads=final,
        stats={"final": count_reads(final) if not cfg.dry_run else 0},
    )


def run_long_profiling(
    cfg: CerberusConfig,
    tuned: TunedParams,
    reads: Path,
    refs: RefManager,
) -> PipelineResult:
    mode = "long-profiling"
    log.info("=== Running long-read profiling (aggressive) ===")
    work = cfg.work_dir / mode
    logs = cfg.logs_dir / mode

    idx = refs.path_to(refs.asset("masked_t2t_hla_minimap2"))

    # 1) primary aligner (minimap2 by default, winnowmap on very-long if --double-pass)
    if cfg.double_pass and tuned.winnowmap_enabled:
        log.info("--double-pass + very-long reads: using winnowmap")
        meryl_db = idx.with_suffix(".meryl")
        primary = winnowmap_singles(
            cfg, tuned, index=idx, meryl_db=meryl_db, reads_in=reads,
            workdir=stage_dir(work, mode, "01_winnowmap"),
            log_dir=logs, tag="01_winnowmap_long_prof",
        )
    else:
        primary = minimap2_singles(
            cfg, tuned, index=idx, reads_in=reads,
            workdir=stage_dir(work, mode, "01_minimap2"),
            log_dir=logs, tag="01_minimap2_long_prof",
        )

    # 2) bbduk aux k-mer pass
    kmer_in = primary.long_reads
    if tuned.bbduk_aux_enabled:
        aux = refs.path_to(refs.asset("aux_refs"))
        k = bbduk_kmer_single(
            cfg, tuned, ref=aux, reads_in=primary.long_reads,        # type: ignore[arg-type]
            workdir=stage_dir(work, mode, "02_bbduk_kmer"),
            log_dir=logs, tag="02_bbduk_kmer_long_prof",
        )
        kmer_in = k.r1

    # 3) entropy
    ent = entropy_single(
        cfg, tuned, reads_in=kmer_in,                        # type: ignore[arg-type]
        workdir=stage_dir(work, mode, "03_entropy"),
        log_dir=logs, tag="03_entropy_long_prof",
    )
    final = cfg.out_dir / f"{cfg.sample_id}.long_profiling.fastq.gz"
    final.unlink(missing_ok=True)
    if not cfg.dry_run:
        ent.r1.replace(final)
    return PipelineResult(
        mode=mode, long_reads=final,
        stats={"final": count_reads(final) if not cfg.dry_run else 0},
    )
