# Changelog

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
