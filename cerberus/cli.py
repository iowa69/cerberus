"""Cerberus command-line interface.

Two-tier help:
  cerberus --help       brief — most users only need this
  cerberus --help-all   exhaustive — every knob exposed

Subcommands:
  cerberus              (default: run the pipeline)
  cerberus run          (explicit form of the above)
  cerberus fetch-refs   pre-warm references
  cerberus doctor       check installation + reference state

Design: argparse with hand-rolled subcommand dispatch so users can run
``cerberus -r1 ... -r2 ...`` without typing ``run``.
"""
from __future__ import annotations

import argparse
import os
import re
import shlex
import signal
import sys
from pathlib import Path

from cerberus import __version__
from cerberus.config import DEFAULT_REF_DIR, CerberusConfig, Platform

_FETCH_ALIASES = {"fetch-refs", "fetchrefs", "fetch_refs"}
_DOCTOR_ALIASES = {"doctor", "check"}
_RUN_ALIASES = {"run"}


def _detect_memory_gb() -> int:
    """Memory budget in GB, respecting cgroup limits.

    A container with a 2 GB limit still reports the host's total in
    /proc/meminfo, so handing that number to ``bbduk -Xmx`` gets the JVM
    SIGKILLed instead of failing cleanly. cgroup v2 then v1 are checked first.
    """
    limits: list[int] = []
    for path in ("/sys/fs/cgroup/memory.max",
                 "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            raw = Path(path).read_text().strip()
        except OSError:
            continue
        if raw and raw != "max":
            try:
                val = int(raw)
            except ValueError:
                continue
            # cgroup v1 uses a huge sentinel for "unlimited".
            if 0 < val < (1 << 62):
                limits.append(val)

    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    limits.append(int(line.split()[1]) * 1024)
                    break
    except (OSError, ValueError, IndexError):
        pass

    if not limits:
        return 12
    total_gb = min(limits) // (1024 ** 3)
    return max(2, min(total_gb - 2, 64))


def _detect_cpus() -> int:
    """Usable CPUs, respecting CPU affinity (taskset, cgroups, schedulers)."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return os.cpu_count() or 4


def _basic_parser(*, advanced: bool = False) -> argparse.ArgumentParser:
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
                       help="Host-scrubbed output for publication")
    modes.add_argument("--all", dest="all_modes", action="store_true",
                       help="Alias for --meta --profiling --gdpr")

    out = p.add_argument_group("Output")
    out.add_argument("-o", "--out-dir", type=Path, default=Path("cerberus_out"),
                     help="Output directory (default: %(default)s)")
    out.add_argument("-s", "--sample-id", type=str, default=None,
                     help="Sample identifier (default: derived from R1 filename)")

    res = p.add_argument_group("Resources")
    res.add_argument("-t", "--threads", type=int, default=_detect_cpus(),
                     help="Threads to use (default: all usable CPUs)")
    res.add_argument("--memory", type=str, default=None,
                     help="Memory budget, e.g. 12G (default: autodetect)")

    tune = p.add_argument_group("Behaviour")
    tune.add_argument("--platform", type=str, default="auto",
                      choices=[x.value for x in Platform],
                      help="Sequencing platform (default: auto)")
    tune.add_argument("--fast", action="store_true",
                      help="Profiling: use minimap2-only path (faster, slightly less sensitive)")
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
    misc.add_argument("-h", "--help", action="store_true", dest="show_help",
                      help="Brief help")
    misc.add_argument("--help-all", action="store_true",
                      help="Show every advanced flag")
    misc.add_argument("--version", action="version", version=f"cerberus {__version__}")

    if advanced:
        _add_advanced_args(p)
    return p


def _add_advanced_args(p: argparse.ArgumentParser) -> None:
    adv = p.add_argument_group(
        "Advanced (--help-all only)",
        "Override autotuned parameters. Leave unset to let Cerberus decide."
    )
    adv.add_argument("--min-length", type=int, default=None,
                     help="Minimum read length after QC")
    adv.add_argument("--min-quality", type=int, default=None,
                     help="Minimum quality for QC (fastp: per-base; fastplong: per-read mean)")
    adv.add_argument("--entropy", type=float, default=None,
                     help="bbduk entropy threshold (0.0-1.0)")
    adv.add_argument("--bbduk-k", type=int, default=None,
                     help="k-mer size for bbduk auxiliary refs (1-31)")
    adv.add_argument("--minimap2-args", type=str, default=None,
                     help="Extra args appended to minimap2 invocations")
    adv.add_argument("--bowtie2-args", type=str, default=None,
                     help="Extra args appended to bowtie2 invocations")
    adv.add_argument("--kraken2-db", type=Path, default=None,
                     help="Override Kraken2 GDPR database directory")
    adv.add_argument("--aux-refs", type=Path, default=None,
                     help="Override auxiliary k-mer references FASTA")
    adv.add_argument("--gdpr-confidence", type=float, default=0.05,
                     help="Kraken2 --confidence for the GDPR scrub (default: %(default)s). "
                          "Low values remove more host; the database holds host taxa only.")
    adv.add_argument("--no-gdpr-kmer-scrub", dest="gdpr_kmer_scrub",
                     action="store_false", default=True,
                     help="Skip the human k-mer (bbduk) GDPR mechanism")
    adv.add_argument("--keep-intermediates", action="store_true",
                     help="Keep _work/ (intermediate FASTQ and BAMs); can be many times "
                          "the input size")


def _print_help(full: bool) -> None:
    _basic_parser(advanced=full).print_help()


def _build_full_parser() -> argparse.ArgumentParser:
    return _basic_parser(advanced=True)


_MEMORY_RE = re.compile(r"^(\d+(?:\.\d+)?)\s*(G|GB|M|MB|T|TB)?$")


def _parse_memory(s: str | None) -> int:
    """Parse a memory budget into whole GB. Raises ValueError with guidance."""
    if not s:
        return _detect_memory_gb()
    m = _MEMORY_RE.match(s.strip().upper())
    if not m:
        raise ValueError(
            f"could not understand memory budget {s!r}; use a form like 12G, 500M or 1T"
        )
    value, unit = float(m.group(1)), (m.group(2) or "G")
    gb = {"G": 1, "GB": 1, "M": 1 / 1024, "MB": 1 / 1024, "T": 1024, "TB": 1024}[unit] * value
    if gb < 1:
        raise ValueError(f"memory budget {s!r} is under 1G; Cerberus needs at least 1G")
    return int(gb)


# Illumina-style trailing fields: _S1 _L001 _R1 _001, with '.' or '-' also
# accepted as the separator (my.sample.R1 is as common as my_sample_R1).
_SAMPLE_STRIP = re.compile(
    r"([._-]S\d+)?([._-]L\d{3})?([._-]R?[12])([._-]\d{3})?$"
)


def derive_sample_id(path: Path) -> str:
    """Derive a sample name from a FASTQ filename.

    Strips known FASTQ extensions (not everything after the first dot, which
    turned ``my.sample.R1.fq.gz`` into ``my``), then removes Illumina's
    lane/read/chunk suffixes.
    """
    name = path.name
    for ext in (".fastq.gz", ".fq.gz", ".fastq.bz2", ".fq.bz2",
                ".fastq", ".fq", ".gz", ".bz2"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    stem = _SAMPLE_STRIP.sub("", name)
    stem = stem.rstrip("._-")
    return stem or name or "sample"


def _config_from_args(ns: argparse.Namespace) -> CerberusConfig:
    if ns.all_modes:
        ns.meta = True
        ns.profiling = True
        ns.gdpr = True

    sample_id = ns.sample_id
    if not sample_id:
        if ns.long_input:
            sample_id = derive_sample_id(ns.long_input)
        elif ns.reads1:
            sample_id = derive_sample_id(ns.reads1)
        else:
            sample_id = "sample"

    return CerberusConfig(
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
        gdpr_confidence=getattr(ns, "gdpr_confidence", 0.05),
        gdpr_kmer_scrub=getattr(ns, "gdpr_kmer_scrub", True),
        keep_intermediates=getattr(ns, "keep_intermediates", False),
        verbose=ns.verbose,
        quiet=ns.quiet,
        dry_run=ns.dry_run,
    )


def _run_fetch_refs(argv: list[str]) -> int:
    p = argparse.ArgumentParser(prog="cerberus fetch-refs",
                                description="Pre-warm the reference cache.")
    p.add_argument("--ref-dir", type=Path, default=DEFAULT_REF_DIR)
    p.add_argument("--update", action="store_true",
                   help="Re-download even if present and verified")
    ns = p.parse_args(argv)

    from cerberus.refs import RefManager
    from cerberus.utils.logger import setup_logging
    setup_logging(ns.ref_dir / "logs", verbose=True)
    rm = RefManager(ns.ref_dir, auto_download=True)
    newer = rm.manifest_update_available()
    if newer:
        print(f"Adopting newer packaged manifest (release {newer}).")
        rm.adopt_packaged_manifest()
    rm.fetch_all(force=ns.update)
    print("All references present.")
    return 0


def _run_doctor(argv: list[str]) -> int:
    p = argparse.ArgumentParser(prog="cerberus doctor",
                                description="Validate installation and reference state.")
    p.add_argument("--ref-dir", type=Path, default=DEFAULT_REF_DIR)
    ns = p.parse_args(argv)

    from cerberus.utils.shell import which

    print("Cerberus doctor")
    print("===============")
    print(f"Version:  {__version__}")
    print(f"Python:   {sys.version.split()[0]}")
    print(f"CPUs:     {_detect_cpus()}")
    print(f"Memory:   {_detect_memory_gb()} GB usable")
    print()

    print("External tools:")
    required = ["fastp", "minimap2", "bowtie2", "bbduk.sh", "kraken2", "samtools", "pigz"]
    optional = ["fastplong", "chopper", "winnowmap", "aria2c", "zstd"]
    missing_req = []
    for t in required:
        path = which(t)
        if not path:
            missing_req.append(t)
        print(f"  {'v' if path else 'x'} {t:<14s}  {path or '(missing)'}")
    print()
    print("Optional tools:")
    for t in optional:
        path = which(t)
        print(f"  {'v' if path else '.'} {t:<14s}  {path or '(missing, falls back if needed)'}")
    print()

    print(f"Reference dir: {ns.ref_dir}")
    problems: list[str] = []
    try:
        from cerberus.refs import RefManager
        rm = RefManager(ns.ref_dir, auto_download=False)
        problems = rm.doctor()
        newer = rm.manifest_update_available()
        if newer:
            print(f"  ! A newer reference manifest ships with this version (release {newer}). "
                  "Run: cerberus fetch-refs --update")
    except Exception as e:                               # noqa: BLE001 - diagnostics only
        print(f"  x could not read reference directory: {e}")
        problems = ["reference directory unreadable"]
    if problems:
        print("Reference issues:")
        for problem in problems:
            print(f"  x {problem}")
    else:
        print("  v All references present.")
    print()

    if missing_req:
        print(f"! Missing required tools: {', '.join(missing_req)}")
        print("  Install with:  conda env create -f environment.yml")
        return 1
    return 0


def _install_signal_handlers() -> None:
    """Tear down child processes on SIGINT/SIGTERM.

    Without this a Ctrl-C returns to the shell while a 32-thread aligner
    keeps running detached.
    """
    from cerberus.utils.shell import terminate_all

    def handler(signum, _frame):
        n = terminate_all()
        if n:
            print(f"\nStopped {n} running tool(s).", file=sys.stderr)
        raise KeyboardInterrupt

    for sig in (signal.SIGINT, signal.SIGTERM):
        try:
            signal.signal(sig, handler)
        except (ValueError, OSError):
            pass  # not on the main thread; the caller still handles KeyboardInterrupt


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]

    if argv and argv[0] in _FETCH_ALIASES:
        return _run_fetch_refs(argv[1:])
    if argv and argv[0] in _DOCTOR_ALIASES:
        return _run_doctor(argv[1:])
    if argv and argv[0] in _RUN_ALIASES:
        argv = argv[1:]
    if argv and argv[0] in {"--help-all", "help-all"}:
        _print_help(full=True)
        return 0

    # Brief help must not list advanced flags, so it is parsed by the basic
    # parser. Anything else goes through the full parser.
    brief = _basic_parser(advanced=False)
    known, _ = brief.parse_known_args(argv)
    if known.help_all:
        _print_help(full=True)
        return 0
    if known.show_help or not argv:
        _print_help(full=False)
        return 0 if argv else 1

    parser = _build_full_parser()
    ns = parser.parse_args(argv)
    if ns.show_help:
        _print_help(full=False)
        return 0

    try:
        cfg = _config_from_args(ns)
    except ValueError as e:
        print(f"cerberus: error: {e}", file=sys.stderr)
        return 2
    cfg.command_line = " ".join(shlex.quote(a) for a in ["cerberus", *argv])

    _install_signal_handlers()

    from cerberus.orchestrator import ConfigError
    from cerberus.orchestrator import run as run_pipeline
    from cerberus.refs import RefManagerError
    from cerberus.utils.shell import ToolError

    try:
        if cfg.update_refs:
            from cerberus.refs import RefManager
            from cerberus.utils.logger import setup_logging
            setup_logging(cfg.ref_dir / "logs", verbose=cfg.verbose)
            RefManager(cfg.ref_dir, auto_download=True).fetch_all(force=True)
        summary = run_pipeline(cfg)
    except ConfigError as e:
        parser.print_usage(sys.stderr)
        print(f"\ncerberus: error: {e}", file=sys.stderr)
        return 2
    except RefManagerError as e:
        print(f"\ncerberus: reference error: {e}", file=sys.stderr)
        return 3
    except ToolError as e:
        print(f"\ncerberus: {e}", file=sys.stderr)
        return 4
    except (FileNotFoundError, PermissionError, OSError) as e:
        print(f"\ncerberus: {type(e).__name__}: {e}", file=sys.stderr)
        return 5
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130

    print()
    print("=" * 66)
    print(f"  Cerberus finished in {summary['elapsed_sec']:.1f}s")
    print("=" * 66)
    for mode, path in summary["outputs"].items():
        if path:
            print(f"  {mode:<26s}  {path}")
    if summary.get("report_html"):
        print(f"  {'run report':<26s}  {summary['report_html']}")
    print(f"  {'accounting':<26s}  {summary['reports']}")
    warnings = summary.get("warnings") or []
    if warnings:
        print()
        print(f"  {len(warnings)} warning(s) — see the run report:")
        for w in warnings[:5]:
            print(f"    ! {w}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
