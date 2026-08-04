"""GDPR post-processor.

Consumes the outputs of meta and/or profiling and produces "_GDPR" variants
scrubbed of host reads. Three mechanisms run in series, and a read survives
only if all three accept it:

  1. **Kraken2** against a compact host (human + great apes + mouse + rat)
     database — an exact k-mer minimizer classifier.
  2. **bbduk** against a human k-mer reference — a different k-mer
     implementation over a different reference build, so a read that slips
     past Kraken2's minimizer scheme still has to survive this. Skipped when
     the asset is absent (or with ``--no-gdpr-kmer-scrub``).
  3. **minimap2** alignment against the masked T2T-CHM13v2.0 + HLA
     reference, dropping the pair if either mate aligns — an alignment-based
     mechanism rather than a k-mer one.

What this does and does not guarantee
-------------------------------------
These mechanisms reduce host content by orders of magnitude, but they are
not a proof of absence. All three ultimately derive from the same reference
assemblies, so sequence absent from those assemblies — population-specific
insertions, V(D)J recombination junctions, novel structural variants — is
invisible to every one of them. Cerberus therefore reports a *measured*
residual figure (see ``residual_host_estimate``) instead of asserting zero,
and the documentation is worded to match.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.pipelines.base import PipelineResult, StageTracker, publish, stage_dir
from cerberus.refs import RefManager
from cerberus.stages.align import minimap2_paired, minimap2_singles
from cerberus.stages.kmer import bbduk_kmer_paired, bbduk_kmer_single
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
    stats: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    mechanisms: list[str] = field(default_factory=list)


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
    track = StageTracker(cfg)

    kdb = refs.kraken2_db_path()
    mm2_idx = refs.path_to(refs.asset("masked_t2t_hla_minimap2"))
    human_kmers = refs.human_kmer_path() if cfg.gdpr_kmer_scrub else None

    mechanisms = ["kraken2", "minimap2"]
    if human_kmers is not None:
        mechanisms.insert(1, "bbduk-human-kmers")
    else:
        log.info("Human k-mer scrub disabled or asset unavailable; running 2 mechanisms.")

    result = GDPRResult(source_mode=pipeline_result.mode, mechanisms=mechanisms)

    if pipeline_result.paired_r1 and pipeline_result.paired_r2:
        before = track.record("gdpr_input_paired",
                              pipeline_result.paired_r1, pipeline_result.paired_r2)

        k_out = kraken2_paired(
            cfg, db=kdb,
            r1_in=pipeline_result.paired_r1,
            r2_in=pipeline_result.paired_r2,
            workdir=stage_dir(work, "01_kraken2"),
            log_dir=logs, tag="01_kraken2",
        )
        r1, r2 = k_out.cleaned_r1, k_out.cleaned_r2
        track.record("gdpr_after_kraken2_paired", r1, r2)

        if human_kmers is not None:
            km = bbduk_kmer_paired(
                cfg, tuned,
                ref=human_kmers,
                r1_in=r1, r2_in=r2,                        # type: ignore[arg-type]
                workdir=stage_dir(work, "02_human_kmers"),
                log_dir=logs, tag="02_human_kmers",
            )
            r1, r2 = km.r1, km.r2
            track.record("gdpr_after_human_kmers_paired", r1, r2)

        mm2_out = minimap2_paired(
            cfg, tuned,
            index=mm2_idx,
            r1_in=r1, r2_in=r2,                            # type: ignore[arg-type]
            workdir=stage_dir(work, "03_minimap2_human"),
            log_dir=logs, tag="03_minimap2_human",
            drop_strategy="either",
        )
        after = track.record("gdpr_after_minimap2_paired", mm2_out.r1, mm2_out.r2)
        _warn_on_total_loss(track, "paired", before, after)

        result.paired_r1 = publish(
            cfg, mm2_out.r1,
            cfg.out_dir / f"{cfg.sample_id}.{pipeline_result.mode}.R1_GDPR.fastq.gz")
        result.paired_r2 = publish(
            cfg, mm2_out.r2,
            cfg.out_dir / f"{cfg.sample_id}.{pipeline_result.mode}.R2_GDPR.fastq.gz")

    if pipeline_result.singletons:
        # `singletons` means different things per head: for meta it is the
        # unpaired leftovers, but for profiling it is the single merged file
        # that *is* the deliverable. Naming both "orphans" would label the
        # profiling head's only GDPR output as if it were a side stream.
        is_primary = pipeline_result.paired_r1 is None
        out_name = (
            f"{cfg.sample_id}.{pipeline_result.mode}_GDPR.fastq.gz" if is_primary
            else f"{cfg.sample_id}.{pipeline_result.mode}.orphans_GDPR.fastq.gz"
        )
        result.singletons = _gdpr_single(
            cfg, tuned, pipeline_result.singletons,
            workdir=stage_dir(work, "04_singletons"),
            log_dir=logs, tag_prefix="04_single",
            kdb=kdb, mm2_idx=mm2_idx, human_kmers=human_kmers, track=track,
            label="merged" if is_primary else "orphans",
            out_name=out_name,
        )

    if pipeline_result.long_reads:
        result.long_reads = _gdpr_single(
            cfg, tuned, pipeline_result.long_reads,
            workdir=stage_dir(work, "05_long"),
            log_dir=logs, tag_prefix="05_long",
            kdb=kdb, mm2_idx=mm2_idx, human_kmers=human_kmers, track=track,
            label="long",
            out_name=f"{cfg.sample_id}.{pipeline_result.mode}_GDPR.fastq.gz",
        )

    result.stats = track.stats
    result.warnings = track.warnings
    return result


def _gdpr_single(
    cfg: CerberusConfig,
    tuned: TunedParams,
    reads: Path,
    *,
    workdir: Path,
    log_dir: Path,
    tag_prefix: str,
    kdb: Path,
    mm2_idx: Path,
    human_kmers: Path | None,
    track: StageTracker,
    label: str,
    out_name: str,
) -> Path | None:
    before = track.record(f"gdpr_input_{label}", reads)

    k_out = kraken2_single(
        cfg, db=kdb, reads_in=reads,
        workdir=workdir, log_dir=log_dir, tag=f"{tag_prefix}_kraken2",
    )
    current = k_out.cleaned_r1
    track.record(f"gdpr_after_kraken2_{label}", current)

    if human_kmers is not None:
        km = bbduk_kmer_single(
            cfg, tuned, ref=human_kmers, reads_in=current,
            workdir=workdir, log_dir=log_dir, tag=f"{tag_prefix}_human_kmers",
        )
        current = km.r1
        track.record(f"gdpr_after_human_kmers_{label}", current)

    mm2_out = minimap2_singles(
        cfg, tuned, index=mm2_idx, reads_in=current,
        workdir=workdir, log_dir=log_dir, tag=f"{tag_prefix}_minimap2_human",
    )
    after = track.record(f"gdpr_after_minimap2_{label}", mm2_out.long_reads)
    _warn_on_total_loss(track, label, before, after)

    return publish(cfg, mm2_out.long_reads, cfg.out_dir / out_name)


def _warn_on_total_loss(track: StageTracker, label: str, before: int, after: int) -> None:
    """An empty GDPR release must never look like a clean one."""
    if before > 0 and after == 0:
        track.warnings.append(
            f"The GDPR scrub removed every {label} read ({before} in, 0 out). "
            "An empty release is not the same as a clean one — check that the Kraken2 "
            "database and the minimap2 index match the host species of this sample."
        )


# Mechanism order within one GDPR stream, as the stage keys are emitted.
_MECHANISMS = [
    ("kraken2", "gdpr_input", "gdpr_after_kraken2"),
    ("bbduk-human-kmers", "gdpr_after_kraken2", "gdpr_after_human_kmers"),
    ("minimap2", "gdpr_after_human_kmers", "gdpr_after_minimap2"),
]


def residual_host_estimate(stats: dict[str, int]) -> dict[str, dict[str, float]]:
    """Per-stream, per-mechanism percentage of reads removed, for the run report.

    Stage keys carry a stream suffix (``_paired``, ``_orphans``, ``_merged``,
    ``_long``), so this discovers the streams present rather than assuming the
    paired one — the profiling and long-read heads have no paired stream at
    all, and hard-coding the paired keys left their table empty.

    Mechanisms that did not run are skipped, and the chain closes over them:
    with the human k-mer asset absent, minimap2's "before" becomes the
    Kraken2 output rather than a missing key.
    """
    streams: dict[str, dict[str, float]] = {}
    suffixes = sorted({k[len("gdpr_input"):] for k in stats if k.startswith("gdpr_input")})

    for suffix in suffixes:
        chain: dict[str, float] = {}
        previous = stats.get(f"gdpr_input{suffix}")
        for name, _, after_key in _MECHANISMS:
            after = stats.get(f"{after_key}{suffix}")
            if after is None or previous is None:
                continue                       # mechanism did not run
            if previous > 0:
                chain[name] = round(100.0 * (previous - after) / previous, 4)
            previous = after
        if chain:
            streams[suffix.lstrip("_") or "reads"] = chain
    return streams
