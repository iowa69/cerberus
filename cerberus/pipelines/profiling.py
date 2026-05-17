"""Profiling pipeline — aggressive, single merged FASTQ for Kraken2/Bracken.

Standard mode:
  1. bowtie2 --very-sensitive-local vs masked-T2T+HLA (paired)
  2. bbduk k-mer pass vs aux refs (paired)
  3. bbduk entropy
  4. concat with orphans → single PROFILING.fastq.gz

Fast mode (--fast / --aligner minimap2):
  1. minimap2 sr vs masked-T2T+HLA (drop_strategy=either)
  2. bbduk entropy
  3. concat with orphans → single PROFILING.fastq.gz

Double-pass (--double-pass):
  Adds a minimap2 pre-filter before bowtie2 — disabled by default. The marginal
  yield is small; this flag exists for users with unusual host contamination.
"""
from __future__ import annotations

from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.pipelines.base import PipelineResult, stage_dir
from cerberus.refs import RefManager
from cerberus.stages.align import bowtie2_paired, minimap2_paired, minimap2_singles
from cerberus.stages.concat import concat_gz
from cerberus.stages.entropy import entropy_paired, entropy_single
from cerberus.stages.kmer import bbduk_kmer_paired, bbduk_kmer_single
from cerberus.stages.qc import FastpOutputs
from cerberus.utils.fastq import count_reads
from cerberus.utils.logger import get_logger

log = get_logger("pipeline.profiling")


def run_profiling(
    cfg: CerberusConfig,
    tuned: TunedParams,
    qc: FastpOutputs,
    refs: RefManager,
) -> PipelineResult:
    mode = "profiling"
    log.info("=== Running profiling pipeline (aggressive, single FASTQ) ===")
    work = cfg.work_dir / mode
    logs = cfg.logs_dir / mode

    mm2_idx = refs.path_to(refs.asset("masked_t2t_hla_minimap2"))

    if cfg.fast:
        log.info("--fast mode: minimap2-only path")
        primary = minimap2_paired(
            cfg, tuned,
            index=mm2_idx,
            r1_in=qc.r1, r2_in=qc.r2,                    # type: ignore[arg-type]
            workdir=stage_dir(work, mode, "01_minimap2"),
            log_dir=logs, tag="01_minimap2",
            drop_strategy="either",
        )
        paired_r1, paired_r2 = primary.r1, primary.r2
    else:
        if cfg.double_pass:
            log.info("--double-pass: minimap2 pre-filter, then bowtie2")
            pre = minimap2_paired(
                cfg, tuned,
                index=mm2_idx,
                r1_in=qc.r1, r2_in=qc.r2,                # type: ignore[arg-type]
                workdir=stage_dir(work, mode, "01_minimap2_pre"),
                log_dir=logs, tag="01_minimap2_pre",
                drop_strategy="either",
            )
            r1_after = pre.r1
            r2_after = pre.r2
        else:
            r1_after = qc.r1
            r2_after = qc.r2

        bt2_idx_dir = refs.path_to(refs.asset("masked_t2t_hla_bowtie2"))
        bt2_prefix = _bowtie2_prefix(bt2_idx_dir)
        primary = bowtie2_paired(
            cfg, tuned,
            index_prefix=bt2_prefix,
            r1_in=r1_after, r2_in=r2_after,              # type: ignore[arg-type]
            workdir=stage_dir(work, mode, "02_bowtie2"),
            log_dir=logs, tag="02_bowtie2",
        )
        paired_r1, paired_r2 = primary.r1, primary.r2

    # k-mer pass against auxiliary references (only when bbduk_aux_enabled)
    if tuned.bbduk_aux_enabled and not cfg.fast:
        aux_ref = refs.path_to(refs.asset("aux_refs"))
        kmer = bbduk_kmer_paired(
            cfg, tuned,
            ref=aux_ref,
            r1_in=paired_r1, r2_in=paired_r2,            # type: ignore[arg-type]
            workdir=stage_dir(work, mode, "03_bbduk_kmer"),
            log_dir=logs, tag="03_bbduk_kmer",
        )
        paired_r1, paired_r2 = kmer.r1, kmer.r2

    # entropy filter
    ent = entropy_paired(
        cfg, tuned,
        r1_in=paired_r1, r2_in=paired_r2,                # type: ignore[arg-type]
        workdir=stage_dir(work, mode, "04_entropy"),
        log_dir=logs, tag="04_entropy",
    )

    # orphan handling: run them through a singles-equivalent of the pipeline
    orphans_clean: Path | None = None
    orphan_inputs = [p for p in (qc.orphans_r1, qc.orphans_r2) if p and p.exists()]
    if orphan_inputs:
        merged = stage_dir(work, mode, "00_orphans") / "qc.orphans.fq.gz"
        concat_gz(cfg, inputs=orphan_inputs, output=merged, log_dir=logs, tag="00_orphan_cat")
        o_aln = minimap2_singles(
            cfg, tuned, index=mm2_idx, reads_in=merged,
            workdir=stage_dir(work, mode, "05_minimap2_orphans"),
            log_dir=logs, tag="05_minimap2_orphans",
        )
        o_kmer_in = o_aln.long_reads
        if tuned.bbduk_aux_enabled and not cfg.fast:
            o_kmer = bbduk_kmer_single(
                cfg, tuned,
                ref=refs.path_to(refs.asset("aux_refs")),
                reads_in=o_aln.long_reads,               # type: ignore[arg-type]
                workdir=stage_dir(work, mode, "06_bbduk_kmer_orphans"),
                log_dir=logs, tag="06_bbduk_kmer_orphans",
            )
            o_kmer_in = o_kmer.r1
        o_ent = entropy_single(
            cfg, tuned,
            reads_in=o_kmer_in,                          # type: ignore[arg-type]
            workdir=stage_dir(work, mode, "07_entropy_orphans"),
            log_dir=logs, tag="07_entropy_orphans",
        )
        orphans_clean = o_ent.r1

    # final concatenation: R1 + R2 + orphans into single PROFILING.fastq.gz
    final = cfg.out_dir / f"{cfg.sample_id}.profiling.fastq.gz"
    final.unlink(missing_ok=True)
    inputs = [ent.r1]
    if ent.r2:
        inputs.append(ent.r2)
    if orphans_clean:
        inputs.append(orphans_clean)
    concat_gz(cfg, inputs=inputs, output=final, log_dir=logs, tag="08_final_concat")

    return PipelineResult(
        mode=mode,
        singletons=final,
        stats={
            "qc_paired": count_reads(qc.r1) if not cfg.dry_run else 0,
            "final": count_reads(final) if not cfg.dry_run else 0,
        },
    )


def _bowtie2_prefix(bt2_dir: Path) -> Path:
    """Find the bowtie2 index prefix inside the extracted directory (recursive).

    Skips files containing ".rev." in the name — those are the auxiliary
    reverse-index files, not the primary prefix carrier.
    """
    for p in bt2_dir.rglob("*.1.bt2*"):
        if ".rev." in p.name:
            continue
        return p.with_suffix("").with_suffix("")  # strip .1.bt2 / .1.bt2l
    raise FileNotFoundError(f"No bowtie2 index found in {bt2_dir}")
