# Pass 4 — CLI and configuration

## Summary

The argument *grammar* is in good shape: the mode matrix (all 32 combinations of `--meta/--profiling/--gdpr/--all/--long`) is coherent, mutual exclusions are enforced, and `ConfigError` failures exit 2 with a usage block and a readable message. Everything outside that narrow `ConfigError` corridor, however, reaches the user as a raw Python traceback with exit 1 — including the two flags most likely to be typed first by a new user: `--dry-run` (crashes 100% of the time, in every mode, before it can print a single stage command) and `--memory` with any non-integer/non-`G` value. The advertised two-tier help does not exist: `cerberus --help` and `cerberus --help-all` are byte-for-byte identical, and `--help` prints a group literally headed "Advanced (`--help-all` only)". Two documented advanced flags, `--kraken2-db` and `--aux-refs`, are parsed and stored in the config but never read by any code in the package — a silent no-op with scientific consequences. Beyond that, no numeric range validation exists anywhere except `--threads`: `--entropy 5.0`, `--bbduk-k 99`, `--min-length -10`, and `--memory 0` all flow unchecked into the generated `bbduk.sh` command lines.

Environment note: the venv contains a non-editable copy of the package at `.venv/lib/python3.13/site-packages/cerberus/`, which `diff -rq` confirms is byte-identical to `/home/iowa/Desktop/cerberus/repo/cerberus/`, so tracebacks citing either path describe the same source. All line numbers below refer to the repo copy.

## Verified working

- **Mode-selection matrix (all 32 combinations)** — drove `_config_from_args` + `validate_config` over every subset of `--meta/--profiling/--gdpr/--all` × `{short, --long}`; every accept/reject is defensible and `_required_pipeline_keys` picks the right asset groups (`orchestrator.py:31-56`, `orchestrator.py:106-115`).
- **`--gdpr` alone is rejected with an actionable message** — `cerberus: error: --gdpr requires at least one of --meta or --profiling (it cleans their outputs).`, exit 2 (`orchestrator.py:37-40`).
- **`--fast` / `--double-pass` mutual exclusion** — `cerberus: error: --fast and --double-pass are mutually exclusive.`, exit 2 (`orchestrator.py:51-52`).
- **Long/short input exclusivity** — `--long -i X -r1 Y` → `--long is incompatible with -r1/-r2.`; `--long` without `-i` → `--long requires -i FILE pointing to a long-read FASTQ.`; both exit 2 (`orchestrator.py:41-45`).
- **`--threads` range check** — `--threads 0` and `--threads -1` both give `cerberus: error: --threads must be >= 1`, exit 2. It is the only numeric range check in the codebase (`orchestrator.py:55-56`).
- **argparse-level type/choice validation** — `--platform bogus` → `invalid choice: 'bogus' (choose from auto, illumina, ont, pacbio-hifi, pacbio-clr)`; `--threads abc` → `invalid int value`; `--entropy abc` → `invalid float value`; all exit 2 (`cli.py:86-94`, `cli.py:132`).
- **`ConfigError` exit path** — prints usage to stderr, then `\ncerberus: error: <msg>`, returns 2 (`cli.py:318-321`). Verified across 8 distinct rejection cases.
- **`KeyboardInterrupt` inside the pipeline** — monkeypatched `orchestrator.run` to raise it; observed `Interrupted.` on stderr and return code `130` (`cli.py:322-324`).
- **Subcommand dispatch and aliases** — `fetch-refs`/`fetchrefs`/`fetch_refs` and `doctor`/`check` all dispatch; their sub-parsers give their own `usage: cerberus fetch-refs ...` / `cerberus doctor: error: unrecognized arguments: extra_junk` (exit 2) and `--help` (exit 0) (`cli.py:291-294`).
- **`cerberus doctor` report quality and exit code** — clear ✓/✗ table for required vs optional tools, per-asset reference status, exit 1 when required tools are missing (`cli.py:236-285`); returns 0 when clean.
- **`--version`** — `cerberus 0.1.1`, exit 0, matching `pyproject.toml` and `cerberus/__init__.py`.
- **`--all` expansion** — sets `meta`, `profiling`, and `gdpr` before config construction (`cli.py:172-175`); confirmed via the matrix (`--all` → `ref-keys=['meta','profiling','gdpr']`).
- **`cfg.modes` omitting gdpr is deliberate and consistent** — `config.py:112-119` returns only meta/profiling; the sole consumer, `orchestrator.py:70`, compensates with `",".join(cfg.modes) + (" +gdpr" if cfg.gdpr else "")`, producing `modes=meta,profiling +gdpr`. A repo-wide grep for `.modes` finds no other consumer, so the omission is safe — but see F9 for a related gap.
- **`_parse_memory` happy paths** — `12G`, `12g`, `32GB`, `8`, and `None` (autodetect ⇒ `memory=27G` on this 30 GB host) all resolve correctly (`cli.py:160-168`).

## Findings

### F1. `--dry-run` crashes with an uncaught `FileNotFoundError` in every mode

- **Severity:** high
- **Location:** `cerberus/orchestrator.py:137` (short), `cerberus/orchestrator.py:181` (long) → `cerberus/autotune.py:115`
- **What:** `run()` unconditionally calls `autotune_from_fastp(qc.json_report, ...)`, which opens fastp's JSON report. Under `--dry-run`, `shell.run()` correctly short-circuits execution (`utils/shell.py:66-69`), so fastp never runs and the JSON is never written. There is no `if cfg.dry_run` guard around autotune, so the very next statement opens a file that cannot exist. `--dry-run` therefore cannot complete for *any* invocation.
- **Trigger:** any `--dry-run` run. Verified for `--meta`, `--long --meta`, and `--all`, with a fully satisfied reference dir and stub tools on `PATH`:

```
$ cerberus -r1 fx/SRR12345_R1.fastq -r2 fx/SRR12345_R2.fastq --meta --dry-run -o dr3 --ref-dir $PWD/refs
EXIT=1
[21:26:16] INFO cerberus.shell  RUN  fastp -i fx/SRR12345_R1.fastq ... --json dr3/_work/00_qc/fastp.json ...
[21:26:16] INFO cerberus.shell  DRY-RUN — skipped
Traceback (most recent call last):
  File ".../cerberus/cli.py", line 317, in main
    summary = run_pipeline(cfg)
  File ".../cerberus/orchestrator.py", line 137, in _run_short
    tuned = autotune_from_fastp(qc.json_report, user_platform=cfg.platform)
  File ".../cerberus/autotune.py", line 115, in autotune_from_fastp
    with fastp_json.open() as f:
FileNotFoundError: [Errno 2] No such file or directory: 'dr3/_work/00_qc/fastp.json'
```

  Long mode fails identically on `dr4/_work/00_qc/fastplong.json` (`orchestrator.py:181`).
- **Consequence:** the flag documented as "Print commands without executing" (`cli.py:112-113`) prints exactly one command and then dies. Users cannot preview a pipeline, CI cannot smoke-test the CLI, and reviewers cannot inspect generated command lines without a full toolchain plus ~23 GB of references. It is also the first thing most people will try after `--help`.
- **Fix:** in `_run_short`/`_run_long`, skip autotune when `cfg.dry_run` and use the existing fallback: `tuned = TunedParams() if cfg.dry_run else autotune_from_fastp(...)` — the helper `_tuned_or_default()` at `orchestrator.py:201-206` already exists for exactly this case but is only used on the GDPR path. Additionally, short-circuit reference download and `require_tools()` under `--dry-run` (see F2 triggers) so a dry run needs neither tools nor references.

### F2. No exception firewall in `main()` — every non-`ConfigError` failure is a raw traceback with exit 1

- **Severity:** high
- **Location:** `cerberus/cli.py:316-324`
- **What:** the `try` around `run_pipeline(cfg)` catches only `ConfigError` and `KeyboardInterrupt`. Every other failure — `RefManagerError`, the `RuntimeError` from `require_tools()`, `ToolError`, `FileNotFoundError`, `PermissionError`, `FileExistsError`, and the `ValueError` from `_parse_memory` (raised earlier, at `cli.py:313`, outside the `try` entirely) — propagates to the console-script wrapper and prints a stack trace. The exceptions carry good, carefully written messages; they are just buried under 8-14 lines of Python internals.
- **Trigger:** observed, all exit 1:

```
# --no-auto-download with a missing asset (message is excellent, presentation is not)
cerberus.refs.RefManagerError: Asset 'masked_t2t_hla_minimap2' missing and --no-auto-download set. Run: cerberus fetch-refs

# missing external tool
RuntimeError: Missing required tool(s): fastp. Install via the cerberus conda environment.

# manifest has no URL yet
cerberus.refs.RefManagerError: Asset 'masked_t2t_hla_minimap2' has no URL in manifest yet (awaiting first release). ...

# -o points at an existing file
FileExistsError: [Errno 17] File exists: 'afile'

# -o or --ref-dir not writable
PermissionError: [Errno 13] Permission denied: '/root/nope'
```

  Also: `KeyboardInterrupt` is only handled around `run_pipeline`. Ctrl-C during the `--update-refs` fetch (`cli.py:307-311`), `fetch-refs`, or `doctor` escapes `main()` uncaught — verified by injecting a `KeyboardInterrupt` into `RefManager.fetch_all`: *"UNCAUGHT, propagates out of main()"*. A user aborting a multi-hour 23 GB download gets a traceback.
- **Consequence:** the CLI reports "this looks like a Cerberus bug" for what are ordinary user/environment errors, so bug reports will be filed against Cerberus for missing conda packages and unwritable directories. Exit code 1 is also used indiscriminately, so wrapper scripts cannot distinguish "bad usage" (2) from "missing reference" from "tool crashed".
- **Fix:** widen the handler and assign distinct codes, e.g. `ConfigError`→2, `RefManagerError`→3, `RuntimeError`/`ToolError`→4, `OSError`→5, `KeyboardInterrupt`→130, with `--verbose` re-raising the traceback for debugging. Move `_config_from_args(ns)` (`cli.py:313`) and the `--update-refs` block inside the same guarded region.

### F3. `--kraken2-db` and `--aux-refs` are parsed, stored, documented — and never read

- **Severity:** high
- **Location:** `cerberus/cli.py:140-143` (definitions), `cerberus/cli.py:210-211` (mapped into `CerberusConfig`), `cerberus/config.py:92-93` (fields)
- **What:** both flags are plumbed into `cfg.kraken2_db_override` / `cfg.aux_refs_override`, but a repo-wide grep finds **zero** consumers outside `cli.py` and `config.py`. Every pipeline resolves these assets from the manifest instead: `refs.path_to(refs.asset("kraken2_gdpr_compact"))` (`pipelines/gdpr.py:49`) and `refs.path_to(refs.asset("aux_refs"))` (`pipelines/profiling.py:90`, `:123`, `pipelines/long_read.py:82`).
- **Trigger:** `cerberus ... --gdpr --kraken2-db /my/custom/db` — accepted without warning; the bundled `kraken2_gdpr_compact` database is used.
- **Consequence:** the most safety-critical override in the tool is a silent no-op. A user who scrubs with a stricter custom Kraken2 database believes their GDPR output was screened against it when it was screened against the default. The failure is invisible: no warning, no log line, correct-looking output files. Same for `--aux-refs`.
- **Fix:** either wire them through (`kdb_dir = cfg.kraken2_db_override or refs.path_to(refs.asset("kraken2_gdpr_compact"))`, and likewise for aux refs, with an `exists()` pre-flight check) or delete the flags. Do not leave them advertised in `--help`.

### F4. Two-tier help does not exist — `--help` is byte-identical to `--help-all`

- **Severity:** high
- **Location:** `cerberus/cli.py:1-14` (docstring claim), `cerberus/cli.py:300`, `cerberus/cli.py:114-117`
- **What:** `main()` builds the **full** parser via `_build_full_parser()` (`cli.py:300`), and `-h/--help` is an `action="help"` bound to that parser (`cli.py:114-115`), so it prints every advanced flag. The brief parser produced by `_basic_parser()` is never the one that services `--help`.
- **Trigger:**

```
$ cerberus --help > help_brief.txt; cerberus --help-all > help_all.txt; diff help_brief.txt help_all.txt
IDENTICAL (no differences)
$ wc -c help_brief.txt help_all.txt
3620 help_brief.txt
3620 help_all.txt
```

  `cerberus --help` output ends with a section headed, verbatim:

```
Advanced (--help-all only):
  Override autotuned parameters. Leave unset to let Cerberus decide.
  --min-length MIN_LENGTH ...
```

  The same happens on the error path: the bare-`cerberus` usage block lists all 24 flags including `--kraken2-db` and `--aux-refs`.
- **Consequence:** the documented UX promise ("brief — most users only need this", `cli.py:4`; README:205 `cerberus --help-all # full help with advanced flags`) is broken, and the help text contradicts itself by labelling visible options as hidden. The wall of 24 flags is exactly what the design was meant to avoid.
- **Fix:** parse in two steps — pre-scan `argv` for `--help-all`/`--help`, print from `_build_full_parser()` or `_basic_parser()` accordingly, and only then build the full parser for actual parsing (advanced flags must still *parse* even when not *displayed*). Use `argparse.SUPPRESS` as the `help=` for advanced options in the brief view, or set `parser.usage` explicitly so the error-path usage line stays short too.

### F5. `_parse_memory` raises raw `ValueError`s and silently accepts zero/negative budgets

- **Severity:** high
- **Location:** `cerberus/cli.py:160-168`, called from `cerberus/cli.py:200`
- **What:** the parser is three `int()` calls with no `try`, no unit table beyond `G`/`GB`, and no range check. It also runs at `cli.py:313`, i.e. *outside* the `try` block at `cli.py:316`, so nothing can intercept it.
- **Trigger:** observed, all exit 1 with a full traceback ending in `File ".../cli.py", line 168, in _parse_memory`:

| input | result |
|---|---|
| `--memory abc` | `ValueError: invalid literal for int() with base 10: 'ABC'` |
| `--memory 12M` | `ValueError: invalid literal for int() with base 10: '12M'` |
| `--memory 12.5G` | `ValueError: invalid literal for int() with base 10: '12.5'` |
| `--memory G` | `ValueError: invalid literal for int() with base 10: ''` |
| `--memory 1e3` | `ValueError: invalid literal for int() with base 10: '1E3'` |
| `--memory 0` | **accepted silently** → `memory_gb=0` |
| `--memory -5` | **accepted silently** → `memory_gb=-5` |

  Note `--memory abc` uppercases before reporting, so the error quotes `'ABC'` — a value the user never typed.
- **Consequence:** `12M` and `12.5G` are natural things to type (samtools, bbduk, and Java all accept them) and produce a crash instead of "memory must be an integer number of gigabytes, e.g. 12G". Worse, `0`/`-5` pass validation and are interpolated straight into `bbduk.sh -Xmx{cfg.memory_gb}g` (`stages/kmer.py:64`), producing `-Xmx0g` / `-Xmx-5g`, which the JVM rejects with an opaque error deep inside a pipeline stage, minutes or hours into a run.
- **Fix:** wrap the parse and raise `ConfigError` (or use `argparse.ArgumentTypeError` via `type=`), accept `M`/`T` and fractional values by normalising to GB, and enforce `>= 1`. Extend `tests/test_cli.py:87-94`, which currently parametrises only the four passing cases (`12G`, `12g`, `32GB`, `8`), with the malformed set above.

### F6. No range validation for any advanced numeric parameter

- **Severity:** medium
- **Location:** `cerberus/orchestrator.py:31-56` (`validate_config` checks only `threads`), `cerberus/cli.py:128-135`
- **What:** `--entropy`, `--bbduk-k`, `--min-length`, and `--min-quality` are accepted with any value argparse can coerce to `int`/`float`, then `apply_user_overrides` (`autotune.py:145-163`) copies them verbatim over the autotuned defaults.
- **Trigger:** direct probe of `validate_config` + `apply_user_overrides`:

```
entropy=5.0        -> ACCEPTED (no error)  tuned.entropy=5.0
entropy=-3.0       -> ACCEPTED (no error)  tuned.entropy=-3.0
min_length=-10     -> ACCEPTED (no error)  tuned.min_length=-10
bbduk_k=99         -> ACCEPTED (no error)  tuned.bbduk_k=99
bbduk_k=0          -> ACCEPTED (no error)  tuned.bbduk_k=0
min_quality=999    -> ACCEPTED (no error)  tuned.min_quality=999
threads=100000     -> ACCEPTED (no error)
```

  These land directly in tool invocations: `f"entropy={entropy}"` (`stages/entropy.py:43`, `:73`) and `f"k={tuned.bbduk_k}"` (`stages/kmer.py:68`).
- **Consequence:** `--entropy 5.0` is outside bbduk's documented 0.0–1.0 range (the help text at `cli.py:133` even states "0.0–1.0") and will discard every read or none; `--bbduk-k 99` exceeds bbduk's maximum k of 31 and aborts the stage; `--min-length -10` and `--min-quality 999` silently destroy the QC step. Every one of these fails late, after fastp and possibly after a multi-GB alignment.
- **Fix:** add range assertions to `validate_config` next to the existing `threads` check — `0.0 <= entropy <= 1.0`, `1 <= bbduk_k <= 31`, `min_length >= 1`, `0 <= min_quality <= 60`, `memory_gb >= 1` — each raising `ConfigError` so it exits 2 with a clean message.

### F7. `--update-refs` runs before config validation and fetches every asset

- **Severity:** medium
- **Location:** `cerberus/cli.py:307-311`
- **What:** the `--update-refs` block sits between parsing and `_config_from_args`/`run_pipeline`, so it executes before `validate_config` ever sees the invocation. It also calls `RefManager.fetch_all()` (`refs.py:225-227`), which walks the *entire* manifest rather than the assets the selected modes need.
- **Trigger:** `cerberus --update-refs` with no mode and no input files — an invocation that `validate_config` would reject instantly:

```
$ cerberus --update-refs --ref-dir $PWD/refs_up -o upx
EXIT=1
  File ".../cerberus/cli.py", line 311, in main
    RefManager(ns.ref_dir, auto_download=True).fetch_all()
  File ".../cerberus/refs.py", line 227, in fetch_all
    self.ensure(all_assets)
cerberus.refs.RefManagerError: Asset 'masked_t2t_hla_minimap2' has no URL in manifest yet ...
```

  It never reaches the `No mode selected` check. Separately, `--meta --update-refs` (which needs one asset) verified all five: `✓ masked_t2t_hla_minimap2 ... ✓ masked_t2t_hla_bowtie2 ... ✓ kraken2_gdpr_compact ... ✓ aux_refs ... ✓ human_kmer_set`.
- **Consequence:** with the shipping manifest (`data/default_manifest.json`: 7.7 + 3.8 + 11.0 + 0.9 + 0.3 GB ≈ 23.7 GB) a user who mistypes their input path or forgets `--meta` still downloads up to 23.7 GB before being told the command line was invalid — and a `--meta`-only user downloads ~16 GB they will never use.
- **Fix:** move the `--update-refs` block after `cfg = _config_from_args(ns)` and after an explicit `validate_config(cfg)` call, and scope the fetch to `refs.required_assets_for(_required_pipeline_keys(cfg))` rather than `fetch_all()`.

### F8. `--update-refs` and `fetch-refs --update` do not re-download anything

- **Severity:** medium
- **Location:** `cerberus/cli.py:106-107` (help text), `cerberus/cli.py:224` + `cli.py:230-231` (`ns.update` never read), `cerberus/refs.py:131-134`
- **What:** both "update" paths funnel into `RefManager.ensure()`, whose first action is `if self.is_satisfied(a): ... continue` (`refs.py:133-135`). There is no force/refresh parameter on `RefManager`, `fetch_all()`, or `ensure()`. In `_run_fetch_refs`, `ns.update` is parsed and then never referenced — `RefManager(ns.ref_dir, auto_download=True)` and `rm.fetch_all()` take no such argument.
- **Trigger:** against an already-satisfied ref dir:

```
$ cerberus ... --update-refs --ref-dir $PWD/refs
[21:38:44] INFO cerberus.refs  ✓ masked_t2t_hla_minimap2 present and verified      # no download
$ cerberus fetch-refs --ref-dir $PWD/refs --update
EXIT=0
✓ All references present.                                                          # no download
```
- **Consequence:** `--update-refs` promises "Force re-download of references" and `--update` promises "Re-download even if present"; neither does. The documented recovery path for a corrupted-but-hash-matching or stale asset (and the manifest-version-bump flow described in `refs.py:6-8`) silently does nothing, leaving users with no way to refresh a reference short of `rm -rf ~/.cerberus/refs`.
- **Fix:** add `force: bool = False` to `RefManager.ensure`/`fetch_all` that skips the `is_satisfied` short-circuit (and unlinks the existing target before download), then pass `ns.update` / `ns.update_refs` through. Until then, correct the two help strings.

### F9. GDPR outputs never appear in the final run summary

- **Severity:** medium
- **Location:** `cerberus/orchestrator.py:101`, printed at `cerberus/cli.py:330-333`
- **What:** `run()` returns `"outputs": {r.mode: r.primary_output for r in result_set}`, where `result_set` holds only the meta/profiling `PipelineResult`s. GDPR results (`orchestrator.py:85-93`) are recorded into `RunAccounting` and then discarded — they never enter the returned dict, so the summary loop in `cli.py:330-332` cannot print them.
- **Trigger:** `--all` run with the pipeline internals stubbed so the summary is reachable:

```
============================================================
  Cerberus finished in 0.0s
============================================================
  meta                sumdemo/S.meta.R1.fastq.gz
  profiling           sumdemo/S.profiling.fastq.gz
  reports            sumdemo/reports
```

  `S.meta.R1_GDPR.fastq.gz` and `S.profiling.GDPR.fastq.gz` (`pipelines/gdpr.py:73-74`, `:89`) are absent.
- **Consequence:** the mode the README presents as the publication-ready deliverable is the one output whose path the user is never told. They must open `reports/accounting.json` to find it — and `RunAccounting.add_final` silently drops entries whose files do not exist (`accounting.py:46-48`), so even that record can come back empty. Also, the summary label column is mis-aligned: `f"  {mode:<18s}  {path}"` vs `f"  reports            {…}"` (`cli.py:332-333`), which is a 1-space discrepancy visible in the output above.
- **Fix:** have `run()` merge GDPR outputs into the returned `outputs` dict under `f"{r.mode}_gdpr"` keys (the same keys already used for accounting at `orchestrator.py:88`), and align the `reports` label with the `:<18s` field width.

### F10. Sample-ID derivation mangles common real-world FASTQ filenames, and can produce an empty ID

- **Severity:** medium
- **Location:** `cerberus/cli.py:177-184`
- **What:** `ns.reads1.name.split(".")[0].replace("_R1", "").replace("_1", "")` truncates at the **first** dot, then strips `_R1`/`_1` as unanchored substrings anywhere in the remainder. The `.replace("_R1","")` is frequently a no-op because the dot-split already removed the token it targets, and `.replace("_1","")` is global.
- **Trigger:** exercised through the real `_config_from_args`:

```
R1 filename                      derived sample_id
SRR12345_R1.fastq.gz             'SRR12345'              OK
sample_1.fq.gz                   'sample'                OK
my.sample.R1.fq.gz               'my'                    truncated at first dot
2024.06.12_run_R1.fq.gz          '2024'                  all distinguishing info lost
A_R1_001.fastq.gz                'A_001'                 '_R1' removed mid-name
Sample_1_S1_L001_R1_001.fastq.gz 'Sample_S1_L001_001'    bcl2fastq default naming, mangled twice
sample_R2.fq.gz                  'sample_R2'             R2 token not handled
_R1.fq.gz                        ''                      EMPTY sample id
```

  Long-read derivation (`cli.py:180`) has the same dot problem: `my.sample.long.fq.gz` → `'my'`.
- **Consequence:** `A_R1_001.fastq.gz` and Illumina's standard `<Sample>_S1_L001_R1_001.fastq.gz` — the single most common filename shape in short-read metagenomics — both yield mangled IDs, and any sample name containing a dot (`Pt.01`, dated runs) collapses to its first token. Since `sample_id` names every deliverable (`pipelines/meta.py:74-75`, `pipelines/profiling.py:138`, `pipelines/gdpr.py:73-74`), batches of samples silently collide on the same output filenames and overwrite each other. The empty-ID case yields hidden dotfiles: `out/.meta.R1.fastq.gz`. Note `-s ""` is safely ignored (falsy → re-derived), but a derived-empty ID is not re-checked.
- **Fix:** strip known FASTQ suffixes rather than splitting on dots (`re.sub(r"\.(fastq|fq)(\.gz|\.bz2|\.zst)?$", "", name)`), then remove a read-token only when it is anchored (`re.sub(r"_(R?[12])(_\d+)?$", "", stem)`), and raise `ConfigError` if the result is empty. Log the derived ID at INFO so users can spot a bad derivation before the run.

### F11. Output directories and a log file are created before the config is validated

- **Severity:** medium
- **Location:** `cerberus/orchestrator.py:60-62`
- **What:** `run()` calls `cfg.ensure_directories()` and `setup_logging()` *before* `validate_config(cfg)`. `ensure_directories` (`config.py:133-136`) creates `out_dir`, `_work`, `logs`, `reports`, `ref_dir`, and `cache_dir`.
- **Trigger:** running `cerberus` with no arguments at all, in an empty directory:

```
$ cerberus
EXIT=2
cerberus: error: No mode selected. ...
$ ls
cerberus_out/
$ find cerberus_out
cerberus_out/logs      cerberus_out/reports      cerberus_out/_work
cerberus_out/logs/cerberus.log.jsonl        # 0 bytes
```
- **Consequence:** every mistyped invocation litters the working directory with a four-directory tree plus an empty JSONL log, and touches `~/.cerberus/{refs,cache}`. Users who run `cerberus` bare to see the help (a very common first action) get a directory they did not ask for, and the presence of `cerberus_out/` misleadingly suggests a run occurred. It also means a validation failure can itself fail with `PermissionError`/`FileExistsError` (see F2) before the real error is reported.
- **Fix:** call `validate_config(cfg)` first, then `ensure_directories()` and `setup_logging()`. Better still, validate in `main()` before invoking `run()` so the CLI owns argument errors end to end. Consider printing help instead of a usage error when `argv` is empty.

### F12. `-s/--sample-id` is unsanitised and can write outside `--out-dir`

- **Severity:** low
- **Location:** `cerberus/cli.py:82-83`, consumed at `cerberus/pipelines/meta.py:74-75`, `cerberus/pipelines/profiling.py:138`, `cerberus/pipelines/gdpr.py:73-74`, `cerberus/utils/fastq.py:63-65`
- **What:** the sample ID is interpolated into output paths with no check for path separators or traversal segments.
- **Trigger:**

```
-s '../escape' -o /home/iowa/Desktop/cerberus/work/p4/si2
  meta R1 -> .../work/p4/si2/../escape.meta.R1.fastq.gz
  resolved -> /home/iowa/Desktop/cerberus/work/p4/escape.meta.R1.fastq.gz   # outside -o
```

  `-s 'a/b/c'` is likewise accepted and reaches `fastp --report_title 'Cerberus QC — a/b/c'`; the resulting output path `out/a/b/c.meta.R1.fastq.gz` points into directories that were never created.
- **Consequence:** results are written outside the declared output directory (bad for containerised/HPC runs where only `-o` is mounted writable), or the run fails late with a `FileNotFoundError` when a slash-bearing ID hits a non-existent subdirectory. Not a privilege-escalation issue — the user already chose the value — but it defeats the "everything lands under `-o`" contract.
- **Fix:** validate in `validate_config`: reject IDs that are empty, contain `os.sep`/`..`, or fail `re.fullmatch(r"[A-Za-z0-9._+-]+", sample_id)`.

### F13. `_SUBCOMMANDS` is dead code, and two of its four entries are not dispatchable

- **Severity:** low
- **Location:** `cerberus/cli.py:26`, dispatch at `cerberus/cli.py:291-298`
- **What:** `_SUBCOMMANDS = {"fetch-refs", "doctor", "run", "help-all"}` is defined and never referenced — a repo-wide grep returns exactly one hit, the definition itself. The actual dispatch is three hard-coded `if argv[0] in {...}` tests that cover `fetch-refs`/`doctor` (plus aliases) and `--help-all` as the *first* argument only. Neither `run` nor bare `help-all` is handled.
- **Trigger:**

```
$ cerberus run -r1 x.fq -r2 y.fq --meta
EXIT=2   cerberus: error: unrecognized arguments: run
$ cerberus help-all
EXIT=2   cerberus: error: unrecognized arguments: help-all
```

  The module docstring at `cli.py:7-8` lists `run` under "Subcommands", and the design note at `cli.py:12-13` says the hand-rolled dispatch exists "so users can run `cerberus -r1 ... -r2 ...` without typing `run`" — which reads as though typing `run` also works. It does not.
- **Consequence:** users following the docstring, or habitually typing the explicit verb, hit an argparse error whose message ("unrecognized arguments: run") gives no hint that the verb should simply be dropped. The unused constant misleads maintainers into thinking dispatch is table-driven.
- **Fix:** drive the dispatch from `_SUBCOMMANDS` (or delete the constant), accept and strip a leading `run`, and accept bare `help-all`. Update the docstring to describe what is actually dispatched.

### F14. Several flag combinations are accepted and then silently ignored

- **Severity:** low
- **Location:** `cerberus/orchestrator.py:31-56`, `cerberus/orchestrator.py:112`
- **What:** `validate_config` checks mode presence, input presence, `--fast`/`--double-pass`, and `--threads`; nothing else. The following all pass validation and run:
  - `--meta --fast` / `--meta --double-pass` — both flags are documented as profiling-only (`cli.py:95-98`), and `_required_pipeline_keys` only consults `cfg.fast` on the profiling branch (`orchestrator.py:112`), so they are complete no-ops. Notably `--long --meta --profiling --fast` *does* get a warning (`orchestrator.py:53-54`), so the concept exists — it is just not applied to the short-read case.
  - `-i reads.fq` without `--long`, alongside valid `-r1/-r2` — the long input is silently discarded.
  - `-i reads.fq` without `--long` and without `-r1/-r2` — error is `Short-read mode requires both -r1 and -r2.`, which never mentions that the user probably meant to add `--long`.
  - `--platform illumina` with `--long`, and `--platform ont` with short reads — accepted; `apply_user_overrides` (`autotune.py:156-158`) then forces `minimap2_preset="sr"` on long reads or `map-ont` on 150 bp reads.
  - `-v` and `-q` together — accepted; `setup_logging` (`utils/logger.py:57`) resolves quiet-wins, so `-v` is silently discarded.
  - `-r1 <directory>` — passes, because only `.exists()` is checked (`orchestrator.py:49`), not `.is_file()`.
  - `Input not found: fx/SRR12345_R1.fastq or nope_R2.fq` (`orchestrator.py:50`) names both files without saying which is missing.
- **Trigger:** each verified individually; every case above reached the pipeline (exit 1 at the F1 crash point) rather than being rejected at exit 2.
- **Consequence:** users believe a flag took effect when it did not — a silent-wrong-parameters class of problem, mild here because the ignored flags are performance rather than correctness knobs, but the platform mismatch does change the aligner preset for the whole run. The `-i`-without-`--long` message actively points at the wrong fix.
- **Fix:** warn (or reject) on `--fast`/`--double-pass` without `--profiling`; reject `-i` without `--long`; warn when `--platform` conflicts with the read-length mode; make `-q` and `-v` mutually exclusive via an argparse group; use `is_file()` and report exactly which input is missing.

### F15. No pre-flight path checks, and `doctor` mutates the filesystem

- **Severity:** low
- **Location:** `cerberus/refs.py:75`, `cerberus/refs.py:79-89`, invoked from `cerberus/cli.py:271`
- **What:** `RefManager.__init__` unconditionally runs `self.ref_dir.mkdir(parents=True, exist_ok=True)` and `_load_manifest()` seeds `default_manifest.json` when absent — so `cerberus doctor`, a pure diagnostic, creates directories and writes a file. Neither `-o` nor `--ref-dir` is checked for writability or for already being a non-directory before the pipeline commits to them.
- **Trigger:**

```
$ cerberus doctor --ref-dir $PWD/docdir      # docdir did not exist
EXIT=1
created:  docdir  docdir/manifest.json       # a read-only command wrote to disk

$ cerberus doctor --ref-dir /root/forbidden/x
EXIT=1
PermissionError: [Errno 13] Permission denied: '/root/forbidden/x'   # raw traceback

$ cerberus -r1 R1 -r2 R2 --meta -o afile     # 'afile' is a regular file
EXIT=1
FileExistsError: [Errno 17] File exists: 'afile'
$ cerberus -r1 R1 -r2 R2 --meta -o /dev/null
EXIT=1
FileExistsError: [Errno 17] File exists: '/dev/null'
```
- **Consequence:** `doctor` cannot be used to inspect someone else's or a read-only reference directory — the tool that is supposed to diagnose a broken install crashes on it, and reports a misleading "missing" state for a directory it just created. A typo'd `-o` produces a `FileExistsError` traceback instead of "output path exists and is not a directory".
- **Fix:** give `RefManager` a `read_only=True` mode used by `doctor` that neither mkdirs nor seeds the manifest, and report "reference directory does not exist" as a normal finding. Add a `validate_config` check that `out_dir` is either absent or a writable directory, raising `ConfigError` (exit 2).
