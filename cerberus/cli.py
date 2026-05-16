"""Cerberus command-line interface.

Two-tier help:
  cerberus --help       brief — most users only need this
  cerberus --help-all   exhaustive — every knob exposed

Subcommands:
  cerberus              (default: run the pipeline)
  cerberus fetch-refs   pre-warm references
  cerberus doctor       check installation + reference state

Design: argparse with hand-rolled subcommand dispatch so users can run
``cerberus -r1 ... -r2 ...`` without typing ``run``.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from cerberus import __version__
from cerberus.config import CerberusConfig, DEFAULT_REF_DIR, Platform


_SUBCOMMANDS = {"fetch-refs", "doctor", "run", "help-all"}


def _detect_memory_gb() -> int:
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])
                    total_gb = kb // (1024 * 1024)
                    return max(4, min(total_gb - 2, 64))
    except OSError:
        pass
    return 12


def _basic_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="cerberus",
        description="Cerberus — three-headed host removal for metagenomic data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  cerberus -r1 R1.fq.gz -r2 R2.fq.gz -o out/ --meta --profiling --gdpr\n"
            "  cerberus --long -i reads.fq.gz -o out/ --all\n"
            "  cerberus --help-all       (show every advanced flag)\n"
            "  cerberus fetch-refs       (pre-warm reference cache)\n"
            "  cerberus doctor           (validate installation)\n"
        ),
        add_help=False,
    )

    inp = p.add_argument_group("Input")
    inp.add_argument("-r1", "--reads1", type=Path, help="R1 paired-end FASTQ(.gz)")
    inp.add_argument("-r2", "--reads2", type=Path, help="R2 paired-end FASTQ(.gz)")
    inp.add_argument("--long", dest="long_mode", action="store_true",
                     help="Enable long-read mode")
    inp.add_argument("-i", "--input", type=Path, dest="long_input",
                     help="Long-read FASTQ(.gz) (with --long)")

    modes = p.add_argument_group(
        "Modes",
        "Pick one or more. --gdpr is a post-processor for the other outputs."
    )
    modes.add_argument("--meta", action="store_true",
                       help="Conservative paired R1/R2 output for assembly")
    modes.add_argument("--profiling", action="store_true",
                       help="Aggressive single-FASTQ output for taxonomic profiling")
    modes.add_argument("--gdpr", action="store_true",
                       help="Zero-host post-processed output for publication")
    modes.add_argument("--all", dest="all_modes", action="store_true",
                       help="Alias for --meta --profiling --gdpr")

    out = p.add_argument_group("Output")
    out.add_argument("-o", "--out-dir", type=Path, default=Path("cerberus_out"),
                     help="Output directory (default: %(default)s)")
    out.add_argument("-s", "--sample-id", type=str, default=None,
                     help="Sample identifier (default: derived from R1 filename)")

    res = p.add_argument_group("Resources")
    res.add_argument("-t", "--threads", type=int, default=os.cpu_count() or 4,
                     help="Threads to use (default: all CPUs)")
    res.add_argument("--memory", type=str, default=None,
                     help="Memory budget, e.g. 12G (default: autodetect)")

    tune = p.add_argument_group("Behaviour")
    tune.add_argument("--platform", type=str, default="auto",
                      choices=[p.value for p in Platform],
                      help="Sequencing platform (default: auto)")
    tune.add_argument("--fast", action="store_true",
                      help="Profiling: use minimap2-only path (faster, ~2%% less sensitive)")
    tune.add_argument("--double-pass", action="store_true",
                      help="Profiling: enable second aligner pass (slower)")

    refs = p.add_argument_group("References")
    refs.add_argument("--ref-dir", type=Path, default=DEFAULT_REF_DIR,
                      help="Reference cache (default: %(default)s)")
    refs.add_argument("--no-auto-download", dest="auto_download",
                      action="store_false", default=True,
                      help="Refuse to download missing references")
    refs.add_argument("--update-refs", action="store_true",
                      help="Force re-download of references")

    misc = p.add_argument_group("Misc")
    misc.add_argument("-v", "--verbose", action="store_true")
    misc.add_argument("-q", "--quiet", action="store_true")
    misc.add_argument("--dry-run", action="store_true",
                      help="Print commands without executing")
    misc.add_argument("-h", "--help", action="help",
                      help="Brief help")
    misc.add_argument("--help-all", action="store_true",
                      help="Show every advanced flag")
    misc.add_argument("--version", action="version", version=f"cerberus {__version__}")

    return p


def _add_advanced_args(p: argparse.ArgumentParser) -> None:
    adv = p.add_argument_group(
        "Advanced (--help-all only)",
        "Override autotuned parameters. Leave unset to let Cerberus decide."
    )
    adv.add_argument("--min-length", type=int, default=None,
                     help="Minimum read length after QC")
    adv.add_argument("--min-quality", type=int, default=None,
                     help="Minimum mean quality for fastp")
    adv.add_argument("--entropy", type=float, default=None,
                     help="bbduk entropy threshold (0.0–1.0)")
    adv.add_argument("--bbduk-k", type=int, default=None,
                     help="k-mer size for bbduk auxiliary refs")
    adv.add_argument("--minimap2-args", type=str, default=None,
                     help="Extra args appended to minimap2 invocations")
    adv.add_argument("--bowtie2-args", type=str, default=None,
                     help="Extra args appended to bowtie2 invocations")
    adv.add_argument("--kraken2-db", type=Path, default=None,
                     help="Override Kraken2 GDPR database path")
    adv.add_argument("--aux-refs", type=Path, default=None,
                     help="Override auxiliary k-mer references FASTA")
    adv.add_argument("--keep-intermediates", action="store_true",
                     help="Keep intermediate BAMs and matched-read files")


def _print_full_help() -> None:
    p = _basic_parser()
    _add_advanced_args(p)
    p.print_help()


def _build_full_parser() -> argparse.ArgumentParser:
    p = _basic_parser()
    _add_advanced_args(p)
    return p


def _parse_memory(s: str | None) -> int:
    if not s:
        return _detect_memory_gb()
    s = s.strip().upper()
    if s.endswith("G"):
        return int(s[:-1])
    if s.endswith("GB"):
        return int(s[:-2])
    return int(s)


def _config_from_args(ns: argparse.Namespace) -> CerberusConfig:
    if ns.all_modes:
        ns.meta = True
        ns.profiling = True
        ns.gdpr = True

    sample_id = ns.sample_id
    if not sample_id:
        if ns.long_input:
            sample_id = ns.long_input.name.split(".")[0]
        elif ns.reads1:
            sample_id = ns.reads1.name.split(".")[0].replace("_R1", "").replace("_1", "")
        else:
            sample_id = "sample"

    cfg = CerberusConfig(
        r1=ns.reads1,
        r2=ns.reads2,
        long_input=ns.long_input,
        long_mode=ns.long_mode,
        out_dir=ns.out_dir,
        sample_id=sample_id,
        meta=ns.meta,
        profiling=ns.profiling,
        gdpr=ns.gdpr,
        platform=Platform(ns.platform),
        double_pass=ns.double_pass,
        fast=ns.fast,
        threads=ns.threads,
        memory_gb=_parse_memory(ns.memory),
        ref_dir=ns.ref_dir,
        auto_download=ns.auto_download,
        update_refs=ns.update_refs,
        min_length=getattr(ns, "min_length", None),
        min_quality=getattr(ns, "min_quality", None),
        entropy=getattr(ns, "entropy", None),
        bbduk_k=getattr(ns, "bbduk_k", None),
        minimap2_args=getattr(ns, "minimap2_args", None),
        bowtie2_args=getattr(ns, "bowtie2_args", None),
        kraken2_db_override=getattr(ns, "kraken2_db", None),
        aux_refs_override=getattr(ns, "aux_refs", None),
        keep_intermediates=getattr(ns, "keep_intermediates", False),
        verbose=ns.verbose,
        quiet=ns.quiet,
        dry_run=ns.dry_run,
    )
    return cfg


def _run_fetch_refs(argv: list[str]) -> int:
    p = argparse.ArgumentParser(prog="cerberus fetch-refs",
                                description="Pre-warm the reference cache.")
    p.add_argument("--ref-dir", type=Path, default=DEFAULT_REF_DIR)
    p.add_argument("--update", action="store_true", help="Re-download even if present")
    ns = p.parse_args(argv)

    from cerberus.refs import RefManager
    from cerberus.utils.logger import setup_logging
    setup_logging(ns.ref_dir / "logs", verbose=True)
    rm = RefManager(ns.ref_dir, auto_download=True)
    rm.fetch_all()
    print("✓ All references present.")
    return 0


def _run_doctor(argv: list[str]) -> int:
    p = argparse.ArgumentParser(prog="cerberus doctor",
                                description="Validate installation and reference state.")
    p.add_argument("--ref-dir", type=Path, default=DEFAULT_REF_DIR)
    ns = p.parse_args(argv)

    from cerberus.refs import RefManager
    from cerberus.utils.shell import which

    print("Cerberus doctor")
    print("===============")
    print(f"Version: {__version__}")
    print(f"Python:  {sys.version.split()[0]}")
    print()

    print("External tools:")
    required = ["fastp", "minimap2", "bowtie2", "bbduk.sh", "kraken2",
                "samtools", "seqkit", "pigz"]
    optional = ["fastplong", "chopper", "winnowmap", "aria2c", "multiqc", "zstd"]
    missing_req = []
    for t in required:
        path = which(t)
        marker = "✓" if path else "✗"
        if not path:
            missing_req.append(t)
        print(f"  {marker} {t:<14s}  {path or '(missing)'}")
    print()
    print("Optional tools:")
    for t in optional:
        path = which(t)
        marker = "✓" if path else "·"
        print(f"  {marker} {t:<14s}  {path or '(missing, falls back if needed)'}")
    print()

    print(f"Reference dir: {ns.ref_dir}")
    rm = RefManager(ns.ref_dir, auto_download=False)
    problems = rm.doctor()
    if problems:
        print("Reference issues:")
        for p in problems:
            print(f"  ✗ {p}")
    else:
        print("  ✓ All references present.")
    print()

    if missing_req:
        print(f"⚠ Missing required tools: {', '.join(missing_req)}")
        print("  Install with:  conda env create -f environment.yml")
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]

    if argv and argv[0] in {"fetch-refs", "fetchrefs", "fetch_refs"}:
        return _run_fetch_refs(argv[1:])
    if argv and argv[0] in {"doctor", "check"}:
        return _run_doctor(argv[1:])

    if argv and argv[0] == "--help-all":
        _print_full_help()
        return 0

    parser = _build_full_parser()
    ns = parser.parse_args(argv)

    if ns.help_all:
        _print_full_help()
        return 0

    if ns.update_refs:
        from cerberus.refs import RefManager
        from cerberus.utils.logger import setup_logging
        setup_logging(ns.ref_dir / "logs", verbose=ns.verbose)
        RefManager(ns.ref_dir, auto_download=True).fetch_all()

    cfg = _config_from_args(ns)
    from cerberus.orchestrator import ConfigError, run as run_pipeline

    try:
        summary = run_pipeline(cfg)
    except ConfigError as e:
        parser.print_usage(sys.stderr)
        print(f"\ncerberus: error: {e}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130

    print()
    print("=" * 60)
    print(f"  Cerberus finished in {summary['elapsed_sec']:.1f}s")
    print("=" * 60)
    for mode, path in summary["outputs"].items():
        if path:
            print(f"  {mode:<18s}  {path}")
    print(f"  reports            {summary['reports']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
