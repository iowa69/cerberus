# Changelog

## v0.2.1

A second review pass — this time adversarial, against the v0.2.0 refactor
itself — found defects that the refactor had introduced. v0.2.0 was published
for less than a day and was superseded before reaching Bioconda; use v0.2.1.

### Fixed — output naming

- **The profiling head's GDPR deliverable was named `…orphans_GDPR…`.**
  `singletons` means the unpaired leftovers for `--meta` but the single merged
  file for `--profiling`, so the profiling head's only scrubbed output was
  labelled as if it were a side stream. It is now
  `<sample>.profiling_GDPR.fastq.gz`; meta's genuine leftovers keep
  `<sample>.meta.orphans_GDPR.fastq.gz`.

### Fixed — regressions introduced in v0.2.0

- **Ctrl-C left forked grandchildren running.** `terminate_all()` signalled
  only direct children, so the `bbduk.sh` wrapper died while its JVM kept
  running and writing into the work directory. Children now start in their own
  process group and the whole group is signalled.
- **Interrupting a run stalled for ten seconds.** The signal handler called
  `Popen.wait()` on a process the main thread already held the wait-lock for.
  The handler now only signals; reaping is left to the owning thread.
- **A pipe failure was blamed on the wrong tool.** v0.2.0 raised on the first
  non-zero status, but a dying consumer takes its producer down with it — so
  `minimap2 | samtools view -o /bad/path` reported minimap2 as the failure.
  The furthest-downstream failing stage is now the headline, with any upstream
  collateral listed alongside it.
- **`_line_count_pipe()` could hang forever.** It waited on the decompressor
  before draining its stderr pipe, so a tool that emitted more than a pipe
  buffer of warnings deadlocked. stderr is now drained first.
- **A truncated gzip could kill a finished run.** `count_reads()` fell through
  to the Python reader, which raises `EOFError` — not an `OSError`, so nothing
  caught it. Corrupt input is now reported and counted as zero.
- **`_work/` cleanup deleted the QC artefacts.** `fastp.html`/`.json`, every
  bbduk `stats=` file and every `flagstat.txt` went with it. They are now
  copied to `reports/qc/` before the intermediates are removed.
- **The report's GDPR mechanism table was empty for the profiling and
  long-read heads**, because it looked only for paired-stream stage keys. It
  now discovers the streams present, and closes the chain over a mechanism
  that did not run instead of dropping the rows after it.
- **The report's "removed here" column diffed across unrelated streams**,
  producing entries like "0 records → removed 8,000 → 100%". Deltas are now
  computed only within a stream.
- **The same file was decompressed up to four times per run** (58 counting
  passes over 15.6× the input size). Counts are cached on file identity
  (size and mtime), so a rewrite still invalidates the entry.

### Added

- A run now warns when the output directory holds deliverables from an earlier,
  wider run that this one did not refresh — previously a narrower re-run left
  stale files with nothing flagging them.
- Regression tests for each of the above (123 tests total).

### Packaging

- `noarch: python` restored. It was removed in v0.2.0 on the strength of a
  review finding, but the recipe already published to Bioconda uses it, as do
  600+ other pure-Python wrappers with compiled bioinformatics dependencies.
  Verified by building the package locally.
- `pyyaml`, `seqkit`, `multiqc`, `bedtools` and `chopper` dropped from the
  runtime dependencies: none is invoked by the pipeline (`fastplong` is a hard
  dependency, so chopper's fallback path is unreachable in a conda install).
  This removes multiqc's whole tree from the install closure.
- The package summary no longer claims "GDPR-compliant outputs", matching the
  wording correction made to the README in v0.2.0.

## v0.2.0 — correctness release

This release came out of a ten-pass review of v0.1.1 covering read-filtering
semantics, the subprocess layer, the reference manager, the CLI, autotuning,
pipeline composition, outputs and accounting, the GDPR claim, documentation and
packaging, and operational robustness. Several defects silently produced wrong
scientific output, so **results produced with v0.1.x should be regenerated**.

### Fixed — silent wrong results

- **`--platform ont` (or any long-read preset) on paired short reads produced
  empty output and exited 0.** minimap2 only emits paired records under
  `-x sr`; other presets treat `-1/-2` as two single-end runs, leaving no
  READ1/READ2 flags, so `samtools fastq -1/-2` wrote nothing. Paired alignment
  now forces a pair-capable preset and warns.
- **`drop_strategy` was inert: `"both"` and `"either"` selected the same
  reads.** `"both"` is documented as the conservative rule for `--meta` ("drop
  the pair only when BOTH mates map") but was implemented as `-f 12`, which
  drops the pair when *either* mate maps. `--meta` was therefore discarding
  genuine microbial mates. Pair-level filtering now uses a samtools filter
  expression (`flag.unmap || flag.munmap`) that can actually express the rule.
- **The surviving mate of a half-mapping pair was discarded** into
  `-s /dev/null`. `--meta` now keeps it as a singleton and routes it through
  the unpaired stream into `<sample>.meta.orphans.fastq.gz`.
- **`pipe()` checked only the last process's exit status.** An OOM-killed
  aligner feeding a healthy `samtools view` produced a truncated BAM and
  reported success. Every stage's status is now checked and the failing stage
  is named.
- **The merged profiling FASTQ contained duplicate read IDs** — both mates of a
  fragment shared a name, so de-duplication by name deleted half the data.
  Mate suffixes are now added during the final merge.
- **`count_reads()` detected gzip by file extension**, so a gzipped `.fastq`
  counted as 3 reads and a plain-text `.gz` counted as 0. Detection is now by
  magic bytes, and decompressor failures are reported instead of yielding a
  plausible number.
- **Autotuned `min_length`/`min_quality` never reached fastp**, because fastp
  ran first and autotune read its output. A prescan of the input now runs
  before QC, so tuned values apply. An explicit `--min-length 0` is also
  honoured rather than being swallowed by an `or` chain.

### Fixed — data loss and crashes

- **`--dry-run` crashed** (`FileNotFoundError` on a fastp report it never
  wrote) and, before crashing, **deleted the outputs of a previous real run**:
  `final.unlink()` sat outside the dry-run guard in every pipeline. Publishing
  is now atomic and never removes a file before its replacement exists.
- **`--long --profiling --fast` fetched no references at all** — it generated a
  pipeline key absent from the asset map, which `.get(key, [])` swallowed — and
  then opened references that had never been downloaded. Unknown keys now raise,
  and a test asserts every mode combination maps to a known asset set.
- **`--long --double-pass` failed with a raw traceback** on a meryl database
  that no manifest ships. winnowmap now checks for it and explains how to build
  it; `winnowmap_enabled` defaults to off.
- **`compress_to()` wrote the compressed payload into the step log** and never
  created its destination (`pigz -c` streams to stdout).
- **Interrupted archive extraction was cached as complete.** Extraction now
  goes to a `.partial` directory renamed into place only on success, with a
  marker file, and archives are removed once extracted.
- **Ctrl-C left aligners running detached.** SIGINT/SIGTERM now tear down the
  whole child process tree.
- Every failure path returns a distinct exit code (2 config, 3 reference,
  4 tool, 5 I/O, 130 interrupt) instead of a raw traceback.

### Changed — GDPR head

- **The Kraken2 confidence threshold was 0.5 against a host-only database.**
  Confidence is the fraction of a read's k-mers matching the taxon; with a
  database containing only what you want removed, a high threshold lets host
  reads escape as "unclassified". Two sequencing errors in a 150 bp read are
  enough; on ONT data at 2-5% error essentially nothing reached 0.5, making the
  mechanism a no-op. The default is now 0.05, exposed as `--gdpr-confidence`.
- **The third mechanism the manifest advertised did not exist.** The
  `human_kmer_set` asset was described as a "belt-and-braces bbduk scrub,
  orthogonal to Kraken2", was fetched by `fetch-refs`, and was read by no code.
  It is now wired into the GDPR head (disable with `--no-gdpr-kmer-scrub`).
- **"Zero detectable human reads" is no longer claimed.** There was no
  detection step, and all mechanisms derive from the same reference assemblies,
  so sequence absent from those assemblies is invisible to all of them. Cerberus
  now measures and reports per-mechanism removal, and flags a mechanism that
  removed nothing. See the README section on what the head does and does not
  guarantee.
- An empty GDPR output is no longer indistinguishable from a clean one: total
  loss raises a warning that appears in the run report.

### Added

- **`reports/cerberus_report.html` on every run** — a self-contained report
  with the resolved parameters, input prescan, per-stage read accounting with
  retention bars, verification of every output (existence, record count, gzip
  validity, pair synchronisation), the per-mechanism GDPR breakdown, warnings,
  tool versions and environment. `reports/run_record.json` carries the same
  data for methods sections.
- Per-stage read accounting for every stage in every head, including the GDPR
  passes. The TSV gained `sample` and `unit` columns so counts can be
  reconciled rather than guessed at.
- A preflight check for all tools a run will need, before any work starts.
- Reference-directory locking, so concurrent runs sharing a `--ref-dir` cannot
  delete each other's in-flight downloads.
- SHA256 verification caching keyed on size and mtime — large indices are no
  longer re-hashed on every run, but any change still forces a re-check.
- Download retries with backoff, timeouts, and aria2c resume.
- `--kraken2-db` and `--aux-refs` now actually take effect; they were parsed
  and ignored.
- cgroup-aware memory detection and affinity-aware CPU detection, so `--memory`
  and `-t` are correct inside containers and under `taskset`.
- `out/_work/` is removed at the end of a run unless `--keep-intermediates`.
- 50 regression tests, one per defect above, plus CI running them with ruff.

### Changed — behaviour worth knowing about

- Very short reads (<80 bp) no longer disable the auxiliary k-mer stage. A
  malformed QC report used to fall into this class, silently removing a
  decontamination step.
- bbduk's `mcf` and the entropy window are now tuned per read-length class.
  A fixed `mcf=0.5` made the k-mer pass statistically unreachable on
  high-error long reads; a fixed 50 bp entropy window exceeded the read length
  for 2x35 libraries.
- Platform detection no longer files ONT R10.4.1 as PacBio HiFi, and `--long`
  is respected rather than being overridden by a short mean length.
- Sample-ID derivation no longer truncates at the first dot
  (`my.sample.R1.fq.gz` gave `my`) and handles Illumina lane/chunk suffixes.
- `--help` and `--help-all` now differ; they were byte-identical.
- `cerberus run` works, as the docs always claimed.
- Out-of-range values (`--entropy 5.0`, `--bbduk-k 99`, `--memory banana`) are
  rejected with a message instead of being interpolated into a command line.
- `samtools flagstat` output is written to its own file rather than a log with
  a `# CMD:` header prepended.
- **GDPR output filenames changed.** Scripts that glob for them need updating:

  | v0.1.1 | v0.2.x |
  |---|---|
  | `<sample>.profiling.GDPR.fastq.gz` | `<sample>.profiling_GDPR.fastq.gz` |
  | `<sample>.meta.GDPR.fastq.gz` | `<sample>.meta.orphans_GDPR.fastq.gz` |
  | `<sample>.<mode>.long_GDPR.fastq.gz` | `<sample>.long_<mode>_GDPR.fastq.gz` |

  The paired names (`<sample>.<mode>.R1_GDPR.fastq.gz`) are unchanged. The
  split exists because `singletons` means different things per head: meta's
  unpaired leftovers versus profiling's single merged deliverable. Naming both
  "orphans" labelled the profiling head's only GDPR output as a side stream.
- Long-read output names use underscores consistently
  (`<sample>.long_meta.fastq.gz`).
- The conda recipe drops `noarch: python` (its runtime dependencies are
  architecture-specific) and `pyyaml` is gone from both dependency lists — it
  was never imported.

## v0.1.1

- Added `build_custom_host_ref.sh`; expanded README with reference-building
  documentation, power-user examples and the custom-host workflow.

## v0.1.0

- Initial release.
