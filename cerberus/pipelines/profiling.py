"""Profiling pipeline — aggressive, single merged FASTQ for Kraken2/Bracken.

Standard mode:
  1. bowtie2 --very-sensitive-local vs masked-T2T+HLA (paired, drop the pair
     if EITHER mate maps)
  2. bbduk k-mer pass vs auxiliary references (paired)
  3. bbduk entropy filter
  4. concat with the unpaired stream -> single <sample>.profiling.fastq.gz

Fast mode (--fast):
  1. minimap2 sr vs masked-T2T+HLA (drop the pair if either mate maps)
  2. bbduk entropy
  3. concat -> single <sample>.profiling.fastq.gz

Double-pass (--double-pass):
  Adds a minimap2 pre-filter before bowtie2. The marginal yield is small;
  this flag exists for users with unusual host contamination.

R1, R2 and the unpaired reads all end up in one file, so the /1 and /2 mate
suffixes are added during the final merge. Without them the two mates of a
fragment share a read ID and any de-duplication by name silently deletes half
the data. The suffixes cannot be added earlier: bbduk's paired reader rejects
mates whose names differ.
"""
from __future__ import annotations

from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.pipelines.base import PipelineResult, StageTracker, publish, stage_dir
from cerberus.refs import RefManager
from cerberus.stages.align import bowtie2_paired, minimap2_paired, minimap2_singles
from cerberus.stages.concat import concat_gz, concat_mate_tagged
from cerberus.stages.entropy import entropy_paired, entropy_single
from cerberus.stages.kmer import bbduk_kmer_paired, bbduk_kmer_single
from cerberus.stages.qc import FastpOutputs
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
    track = StageTracker(cfg)

    mm2_idx = refs.path_to(refs.asset("masked_t2t_hla_minimap2"))
    qc_in = track.record("00_qc_paired", qc.r1, qc.r2)

    if cfg.fast:
        log.info("--fast mode: minimap2-only path")
        primary = minimap2_paired(
            cfg, tuned,
            index=mm2_idx,
            r1_in=qc.r1, r2_in=qc.r2,                    # type: ignore[arg-type]
            workdir=stage_dir(work, "01_minimap2"),
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
                workdir=stage_dir(work, "01_minimap2_pre"),
                log_dir=logs, tag="01_minimap2_pre",
                drop_strategy="either",
            )
            track.record("01_minimap2_prefilter", pre.r1, pre.r2)
            r1_after, r2_after = pre.r1, pre.r2
        else:
            r1_after, r2_after = qc.r1, qc.r2

        bt2_idx_dir = refs.path_to(refs.asset("masked_t2t_hla_bowtie2"))
        bt2_prefix = _bowtie2_prefix(bt2_idx_dir)
        primary = bowtie2_paired(
            cfg, tuned,
            index_prefix=bt2_prefix,
            r1_in=r1_after, r2_in=r2_after,              # type: ignore[arg-type]
            workdir=stage_dir(work, "02_bowtie2"),
            log_dir=logs, tag="02_bowtie2",
            drop_strategy="either",
        )
        paired_r1, paired_r2 = primary.r1, primary.r2

    after_align = track.record("01_host_removal_paired", paired_r1, paired_r2)
    track.check_not_empty("01_host_removal_paired", qc_in, after_align)

    # k-mer pass against auxiliary references
    if tuned.bbduk_aux_enabled and not cfg.fast:
        aux_ref = refs.aux_refs_path()
        kmer = bbduk_kmer_paired(
            cfg, tuned,
            ref=aux_ref,
            r1_in=paired_r1, r2_in=paired_r2,            # type: ignore[arg-type]
            workdir=stage_dir(work, "03_bbduk_kmer"),
            log_dir=logs, tag="03_bbduk_kmer",
        )
        paired_r1, paired_r2 = kmer.r1, kmer.r2
        track.record("02_aux_kmer_paired", paired_r1, paired_r2)

    ent = entropy_paired(
        cfg, tuned,
        r1_in=paired_r1, r2_in=paired_r2,                # type: ignore[arg-type]
        workdir=stage_dir(work, "04_entropy"),
        log_dir=logs, tag="04_entropy",
    )
    track.record("03_entropy_paired", ent.r1, ent.r2)

    # unpaired stream: fastp orphans through the single-end equivalent
    singles_clean: Path | None = None
    singleton_inputs = [
        p for p in (qc.orphans_r1, qc.orphans_r2, primary.singletons)
        if p and p.exists() and p.stat().st_size > 40
    ]
    if singleton_inputs:
        merged = stage_dir(work, "00_singletons") / "singletons.fq.gz"
        concat_gz(cfg, inputs=singleton_inputs, output=merged,
                  log_dir=logs, tag="00_singleton_cat")
        o_aln = minimap2_singles(
            cfg, tuned, index=mm2_idx, reads_in=merged,
            workdir=stage_dir(work, "05_minimap2_singletons"),
            log_dir=logs, tag="05_minimap2_singletons",
        )
        o_next = o_aln.long_reads
        if tuned.bbduk_aux_enabled and not cfg.fast:
            o_kmer = bbduk_kmer_single(
                cfg, tuned,
                ref=refs.aux_refs_path(),
                reads_in=o_next,                         # type: ignore[arg-type]
                workdir=stage_dir(work, "06_bbduk_kmer_singletons"),
                log_dir=logs, tag="06_bbduk_kmer_singletons",
            )
            o_next = o_kmer.r1
        o_ent = entropy_single(
            cfg, tuned,
            reads_in=o_next,                             # type: ignore[arg-type]
            workdir=stage_dir(work, "07_entropy_singletons"),
            log_dir=logs, tag="07_entropy_singletons",
        )
        singles_clean = o_ent.r1
        track.record("04_singletons", singles_clean)

    # final merge: R1 + R2 + unpaired into one FASTQ, mate-tagged so every
    # read ID in the merged file stays unique
    staged = stage_dir(work, "08_final") / "profiling.fastq.gz"
    tagged: list[tuple[Path, str]] = []
    if ent.r1 is not None:
        tagged.append((ent.r1, "/1"))
    if ent.r2 is not None:
        tagged.append((ent.r2, "/2"))
    if singles_clean is not None:
        tagged.append((singles_clean, ""))
    concat_mate_tagged(cfg, inputs=tagged, output=staged, log_dir=logs, tag="08_final_merge")
    final = publish(cfg, staged, cfg.out_dir / f"{cfg.sample_id}.profiling.fastq.gz")
    track.record("05_final_merged", final)

    return PipelineResult(
        mode=mode,
        singletons=final,
        stats=track.stats,
        warnings=track.warnings,
    )


def _bowtie2_prefix(bt2_dir: Path) -> Path:
    """Find the bowtie2 index prefix inside the extracted directory (recursive).

    Skips files containing ".rev." in the name — those are the auxiliary
    reverse-index files, not the primary prefix carrier.
    """
    for p in sorted(bt2_dir.rglob("*.1.bt2*")):
        if ".rev." in p.name:
            continue
        return p.with_suffix("").with_suffix("")  # strip .1.bt2 / .1.bt2l
    raise FileNotFoundError(
        f"No bowtie2 index (*.1.bt2) found under {bt2_dir}. "
        "The reference archive may have extracted incompletely — "
        "delete it and re-run 'cerberus fetch-refs'."
    )
