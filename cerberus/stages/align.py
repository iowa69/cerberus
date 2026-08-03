"""Aligner wrappers: minimap2, bowtie2, winnowmap. Each returns unmapped read paths.

Pair-level filtering
--------------------
Host removal on paired data is a *pair-level* decision, but SAM flags are
per-record. The two strategies Cerberus offers are:

  "both"    conservative (``--meta``): drop the pair only when BOTH mates map
            to the host. A pair with exactly one mapping mate is kept intact,
            so the microbial mate survives and the pairing an assembler needs
            is preserved — at the cost of carrying its host mate along. This
            is what "retains microbial reads even at the cost of some residual
            host" actually means.

  "either"  aggressive (``--profiling``, ``--gdpr``): drop the whole pair as
            soon as EITHER mate maps. Nothing host-adjacent survives.

``samtools fastq`` cannot express "keep the pair when either mate is
unmapped" through ``-f``/``-F`` alone, because those test one record at a
time. We therefore pre-filter the BAM with a samtools filter expression
(``-e``), which can read both ``flag.unmap`` and ``flag.munmap``, and only
then convert to FASTQ. See ``_pair_filter_expr``.
"""
from __future__ import annotations

import shlex
from dataclasses import dataclass
from pathlib import Path

from cerberus.config import CerberusConfig, TunedParams
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import pipe, require_tools, run

log = get_logger("align")

# minimap2 presets that perform paired-end mapping. Every other preset treats
# the two input FASTQs as two independent single-end runs, which silently
# destroys paired output (no READ1/READ2 flags => samtools fastq routes
# everything to -0). See _paired_preset().
_PAIRED_CAPABLE_PRESETS = {"sr"}


@dataclass
class AlignOutputs:
    r1: Path | None = None
    r2: Path | None = None
    singletons: Path | None = None
    long_reads: Path | None = None
    stats: Path | None = None


def _paired_preset(preset: str) -> str:
    """Return a preset safe for paired-end alignment.

    minimap2 only emits proper paired records under ``-x sr``. Running
    ``-ax map-ont`` over ``-r1/-r2`` produces records with neither the READ1
    nor the READ2 flag set; ``samtools fastq -1/-2`` then writes nothing and
    the run "succeeds" with empty outputs. Guard against that here rather
    than letting it through.
    """
    if preset in _PAIRED_CAPABLE_PRESETS:
        return preset
    log.warning(
        "minimap2 preset %r cannot map paired-end reads; falling back to 'sr' for this "
        "stage. (Long-read presets treat -1/-2 as two single-end runs and would emit "
        "an empty paired output.)",
        preset,
    )
    return "sr"


def _pair_filter_expr(strategy: str) -> str | None:
    """samtools filter expression selecting the records to KEEP.

    Returns None when a plain flag filter suffices.
    """
    if strategy == "both":
        # Keep the pair unless both mates mapped => keep records where this
        # read or its mate is unmapped.
        return "flag.unmap || flag.munmap"
    if strategy == "either":
        # Keep only pairs where neither mate mapped.
        return None
    raise ValueError(f"Unknown drop_strategy: {strategy!r} (expected 'both' or 'either')")


def minimap2_paired(
    cfg: CerberusConfig,
    tuned: TunedParams,
    *,
    index: Path,
    r1_in: Path,
    r2_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
    drop_strategy: str = "both",
    keep_singletons: bool = False,
) -> AlignOutputs:
    """Align paired reads with minimap2 and retain the non-host fraction.

    ``keep_singletons`` captures any read whose mate was filtered away into
    its own file, instead of letting ``samtools fastq`` discard it silently.

    Read names are left untouched. Callers that merge R1 and R2 into one
    file add the /1 and /2 suffixes at merge time instead.
    """
    require_tools("minimap2", "samtools")
    workdir.mkdir(parents=True, exist_ok=True)

    bam = workdir / f"{tag}.bam"
    out_r1 = workdir / f"{tag}.unmapped.R1.fq.gz"
    out_r2 = workdir / f"{tag}.unmapped.R2.fq.gz"
    out_s = workdir / f"{tag}.unmapped.singletons.fq.gz"
    stats = workdir / f"{tag}.flagstat.txt"

    extra_tokens = _merge_extra(cfg.minimap2_args, tuned.minimap2_extra)
    preset = _paired_preset(tuned.minimap2_preset)

    minimap_cmd = [
        "minimap2",
        "-ax", preset,
        "-t", str(cfg.threads),
        "--secondary=no",
        *extra_tokens,
        str(index),
        str(r1_in),
        str(r2_in),
    ]
    view_cmd = ["samtools", "view", "-@", str(_io_threads(cfg)), "-b", "-o", str(bam), "-"]
    pipe([minimap_cmd, view_cmd],
         log_path=log_dir / f"{tag}.minimap2.log", dry_run=cfg.dry_run)

    _bam_to_paired_fastq(
        cfg, bam=bam, strategy=drop_strategy,
        out_r1=out_r1, out_r2=out_r2,
        out_singletons=out_s if keep_singletons else None,
        log_dir=log_dir, tag=tag,
    )

    run(["samtools", "flagstat", str(bam)],
        log_path=log_dir / f"{tag}.flagstat.log", stdout_path=stats, dry_run=cfg.dry_run)

    if not cfg.keep_intermediates and not cfg.dry_run:
        bam.unlink(missing_ok=True)

    return AlignOutputs(
        r1=out_r1, r2=out_r2,
        singletons=out_s if keep_singletons else None,
        stats=stats,
    )


def _bam_to_paired_fastq(
    cfg: CerberusConfig,
    *,
    bam: Path,
    strategy: str,
    out_r1: Path,
    out_r2: Path,
    out_singletons: Path | None,
    log_dir: Path,
    tag: str,
) -> None:
    """Convert a name-collated BAM to synchronised paired FASTQ."""
    expr = _pair_filter_expr(strategy)
    singleton_sink = str(out_singletons) if out_singletons else "/dev/null"

    # -n leaves read names exactly as they are. The /1 and /2 suffixes are NOT
    # added here: bbduk's paired reader rejects mates whose names differ and
    # then deadlocks rather than exiting. Callers that merge R1 and R2 into one
    # file tag the names at merge time instead (see concat.concat_mate_tagged).
    tail = [
        "-1", str(out_r1),
        "-2", str(out_r2),
        "-0", "/dev/null",
        "-s", singleton_sink,
        "-n",
    ]

    if expr is None:
        # "either": keep only pairs where neither mate mapped.
        run(
            ["samtools", "fastq", "-@", str(_io_threads(cfg)), "-f", "12", *tail, str(bam)],
            log_path=log_dir / f"{tag}.samtools_fastq.log", dry_run=cfg.dry_run,
        )
        return

    # "both": needs a pair-aware expression, so filter the BAM first.
    filt_cmd = [
        "samtools", "view", "-@", str(_io_threads(cfg)), "-b",
        "-e", expr, "-o", "-", str(bam),
    ]
    fastq_cmd = ["samtools", "fastq", "-@", str(_io_threads(cfg)), *tail, "-"]
    pipe([filt_cmd, fastq_cmd],
         log_path=log_dir / f"{tag}.samtools_fastq.log", dry_run=cfg.dry_run)


def minimap2_singles(
    cfg: CerberusConfig,
    tuned: TunedParams,
    *,
    index: Path,
    reads_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
) -> AlignOutputs:
    """Align single-end / orphan / long reads with minimap2, retain unmapped."""
    require_tools("minimap2", "samtools")
    workdir.mkdir(parents=True, exist_ok=True)

    bam = workdir / f"{tag}.bam"
    out = workdir / f"{tag}.unmapped.fq.gz"
    stats = workdir / f"{tag}.flagstat.txt"

    extra_tokens = _merge_extra(cfg.minimap2_args, tuned.minimap2_extra)

    minimap_cmd = [
        "minimap2",
        "-ax", tuned.minimap2_preset,
        "-t", str(cfg.threads),
        "--secondary=no",
        *extra_tokens,
        str(index),
        str(reads_in),
    ]
    view_cmd = ["samtools", "view", "-@", str(_io_threads(cfg)), "-b", "-o", str(bam), "-"]
    pipe([minimap_cmd, view_cmd],
         log_path=log_dir / f"{tag}.minimap2.log", dry_run=cfg.dry_run)

    # -F 0x900 drops secondary/supplementary records so a chimeric long read
    # cannot be emitted twice.
    run(
        ["samtools", "fastq", "-@", str(_io_threads(cfg)),
         "-f", "4", "-F", "0x900", "-0", str(out), "-n", str(bam)],
        log_path=log_dir / f"{tag}.samtools_fastq.log", dry_run=cfg.dry_run,
    )
    run(["samtools", "flagstat", str(bam)],
        log_path=log_dir / f"{tag}.flagstat.log", stdout_path=stats, dry_run=cfg.dry_run)
    if not cfg.keep_intermediates and not cfg.dry_run:
        bam.unlink(missing_ok=True)

    return AlignOutputs(long_reads=out, stats=stats)


def bowtie2_paired(
    cfg: CerberusConfig,
    tuned: TunedParams,
    *,
    index_prefix: Path,
    r1_in: Path,
    r2_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
    drop_strategy: str = "either",
    keep_singletons: bool = False,
) -> AlignOutputs:
    """Run bowtie2 on paired reads, extract the non-host fraction."""
    require_tools("bowtie2", "samtools")
    workdir.mkdir(parents=True, exist_ok=True)

    bam = workdir / f"{tag}.bam"
    out_r1 = workdir / f"{tag}.unmapped.R1.fq.gz"
    out_r2 = workdir / f"{tag}.unmapped.R2.fq.gz"
    out_s = workdir / f"{tag}.unmapped.singletons.fq.gz"
    stats = workdir / f"{tag}.flagstat.txt"

    preset_tokens = shlex.split(tuned.bowtie2_preset)
    extra_tokens = _merge_extra(cfg.bowtie2_args, tuned.bowtie2_extra)

    bt2_cmd = [
        "bowtie2",
        *preset_tokens,
        "-p", str(cfg.threads),
        "-x", str(index_prefix),
        "-1", str(r1_in),
        "-2", str(r2_in),
        *extra_tokens,
    ]
    view_cmd = ["samtools", "view", "-@", str(_io_threads(cfg)), "-b", "-o", str(bam), "-"]
    pipe([bt2_cmd, view_cmd],
         log_path=log_dir / f"{tag}.bowtie2.log", dry_run=cfg.dry_run)

    _bam_to_paired_fastq(
        cfg, bam=bam, strategy=drop_strategy,
        out_r1=out_r1, out_r2=out_r2,
        out_singletons=out_s if keep_singletons else None,
        log_dir=log_dir, tag=tag,
    )

    run(["samtools", "flagstat", str(bam)],
        log_path=log_dir / f"{tag}.flagstat.log", stdout_path=stats, dry_run=cfg.dry_run)
    if not cfg.keep_intermediates and not cfg.dry_run:
        bam.unlink(missing_ok=True)

    return AlignOutputs(
        r1=out_r1, r2=out_r2,
        singletons=out_s if keep_singletons else None,
        stats=stats,
    )


def winnowmap_singles(
    cfg: CerberusConfig,
    tuned: TunedParams,
    *,
    index: Path,
    meryl_db: Path,
    reads_in: Path,
    workdir: Path,
    log_dir: Path,
    tag: str,
) -> AlignOutputs:
    """winnowmap for very long reads (replaces minimap2 when --double-pass)."""
    require_tools("winnowmap", "samtools")
    if not meryl_db.exists() and not cfg.dry_run:
        raise FileNotFoundError(
            f"winnowmap needs a repetitive-k-mer file at {meryl_db}, which is not present. "
            "Build it with:\n"
            f"  meryl count k=15 output {meryl_db.with_suffix('.meryldb')} <host.fa>\n"
            f"  meryl print greater-than distinct=0.9998 {meryl_db.with_suffix('.meryldb')} "
            f"> {meryl_db}\n"
            "or drop --double-pass to use minimap2 instead."
        )
    workdir.mkdir(parents=True, exist_ok=True)

    bam = workdir / f"{tag}.bam"
    out = workdir / f"{tag}.unmapped.fq.gz"
    stats = workdir / f"{tag}.flagstat.txt"

    preset = tuned.minimap2_preset if tuned.minimap2_preset.startswith("map-") else "map-ont"
    wm_cmd = [
        "winnowmap",
        "-W", str(meryl_db),
        "-ax", preset,
        "-t", str(cfg.threads),
        str(index),
        str(reads_in),
    ]
    view_cmd = ["samtools", "view", "-@", str(_io_threads(cfg)), "-b", "-o", str(bam), "-"]
    pipe([wm_cmd, view_cmd],
         log_path=log_dir / f"{tag}.winnowmap.log", dry_run=cfg.dry_run)

    run(
        ["samtools", "fastq", "-@", str(_io_threads(cfg)),
         "-f", "4", "-F", "0x900", "-0", str(out), "-n", str(bam)],
        log_path=log_dir / f"{tag}.samtools_fastq.log", dry_run=cfg.dry_run,
    )
    run(["samtools", "flagstat", str(bam)],
        log_path=log_dir / f"{tag}.flagstat.log", stdout_path=stats, dry_run=cfg.dry_run)
    if not cfg.keep_intermediates and not cfg.dry_run:
        bam.unlink(missing_ok=True)
    return AlignOutputs(long_reads=out, stats=stats)


def _merge_extra(user_args: str | None, tuned_extra: str) -> list[str]:
    """Combine autotuned extra args with user-supplied ones.

    The user's tokens come last so they win on any repeated option, but the
    autotuned ones are no longer silently thrown away just because the user
    passed something unrelated.
    """
    tokens: list[str] = []
    for source in (tuned_extra, user_args):
        if not source:
            continue
        try:
            tokens.extend(shlex.split(source))
        except ValueError as e:
            raise ValueError(f"Could not parse aligner arguments {source!r}: {e}") from e
    return tokens


def _io_threads(cfg: CerberusConfig) -> int:
    """Threads for samtools inside a pipe.

    The aligner already claims ``cfg.threads``; giving samtools the same
    number again oversubscribes the machine, and samtools' compression
    threads saturate well before 4.
    """
    return max(1, min(4, cfg.threads))
