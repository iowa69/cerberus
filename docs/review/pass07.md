# Pass 7 — Outputs and accounting

## Summary

The accounting writer itself is mechanically sound — I constructed a `RunAccounting`, wrote it, and parsed both artifacts back: every TSV row has exactly 4 tab-separated fields, `asdict()` serialises cleanly to JSON for every field type in use, and no `Path` object ever reaches a `str`-assuming slot. The core concatenation assumption is also correct: I concatenated gzip-, pigz- and python-produced members and verified the result reads as one stream through `gzip -t`, `pigz -dc`, `zcat`, python `gzip`, `seqkit`, `bbduk` and Java's `GZIPInputStream`. Everything else in this lens is weaker than it looks. A full end-to-end run (miniature real references, real fastp/bowtie2/minimap2/bbduk/kraken2) showed that `reports/` contains only `accounting.tsv` and `accounting.json` — no fastp report, no flagstats, no tool versions, no resolved `TunedParams`, no command line — so a run is not reconstructible for a methods section; that the `reads` column silently mixes fragments and reads, making `profiling final 630` look identical to `_input input_r1 630` on a run that actually discarded half the data; that `profiling.fastq.gz` and the publication-facing `profiling.GDPR.fastq.gz` contain 630 records under only 330 unique read IDs; and that `--dry-run` into a populated output directory **deleted all 8 previous output FASTQs**. Nine of the twelve producible output filenames disagree with the README in some way, including all four long-read names and a `long_meta`/`long-meta` split within a single directory.

*Caveat on citations:* `cerberus/stages/align.py` was being concurrently rewritten in the working tree during this pass (`git status` shows ` M cerberus/stages/align.py`, +200/-66 vs HEAD, and the working copy currently raises `TypeError: run() got an unexpected keyword argument 'stdout_path'`). All `align.py` line numbers below refer to **HEAD (125a656)**; end-to-end runs used a pristine `git archive HEAD` export in `work/p7/repo_head`. The `-n` flag that drives F2 is present in both HEAD (align.py:76) and the working copy (align.py:171).

## Verified working

- **`RunAccounting.write()` TSV shape** — built a `RunAccounting` with input/QC/stage/final rows including a `Path`-valued `file` and a `None`-valued final path, wrote it, split every line on `\t`: all 16 rows returned exactly 4 fields. Rows that omit a file end in a literal trailing tab, so the field is present-but-empty rather than missing. `accounting.py:62-76`. (`work/p7/t1_accounting.py`)
- **`asdict()` → JSON round-trip** — `stages` (list of `StageCount` dataclasses) recurses correctly; `final_outputs` (`dict[str, dict]`) holds only `str`/`int` because `add_stage` (`accounting.py:40`) and `add_final` (`accounting.py:50`) both stringify `Path` at the boundary. `json.loads` of the written file reproduced every field. No `Path` is ever passed where a `str` is assumed. `accounting.py:60`.
- **Concatenated gzip members are readable as one stream** — `cat gzip.gz pigz.gz python.gz > merged.gz` (150 bytes, 3 members) decoded to the full 40 lines / 10 records under `gzip -t` (rc=0), `pigz -t` (rc=0), `zcat`, `pigz -dc`, python `gzip`, `seqkit stats` (10 seqs), `bbduk.sh` (10 reads in) and Java `GZIPInputStream`. Cerberus's own `count_reads()` returned 10. The premise in `concat.py:1` holds. (`work/p7/t3_concat.sh`)
- **`count_reads()` on a multi-member concatenation** — returns the correct record count via the `pigz -dc | wc -l` path, so the profiling final count is not corrupted by the concat. `fastq.py:11-26`.
- **`count_reads()` tolerates the zero-byte file** — short-circuits on `path.stat().st_size == 0` before invoking any decompressor, so accounting does not crash on the degenerate output. `fastq.py:15-16`.
- **R1/R2 stay synchronised through the `either` strategy** — with `-f 4` and `-0 /dev/null -s /dev/null`, `samtools fastq` routes half-surviving pairs to the discarded singleton stream, so `e.R1` and `e.R2` came out with identical name lists (4 = 4). bbduk `in1=/in2=` then consumed them as paired without complaint. The `-f 4 -F 8` variant suggested by the comment at `align.py:239` would have produced **zero** reads (`-F 8` excludes mate-unmapped, i.e. exactly the pairs you want) — good thing it was never implemented. (`work/p7/t6_sync.sh`)
- **Every real end-to-end output is a valid gzip** — `gzip -t` returned rc=0 on all 8 short-mode outputs and all 4 long-mode outputs of a real run.
- **GDPR read-delta is reconcilable** — `meta.paired_r1 300` vs `meta_gdpr.paired_r1 300` in `accounting.tsv` correctly shows the GDPR pass removed nothing on the test corpus. Of all the pipeline stages, GDPR is the only one whose removal count a reviewer can actually compute.
- **`_final` rows carry a path** — unlike stage rows, `_final` rows populate column 4, so the TSV at least ties each terminal count to a file. `accounting.py:75`.
- **Kraken2 tolerates the merged profiling file** — with a purpose-built miniature Kraken2 DB, `kraken2` classified all 12 records of a duplicate-ID merged file without error and produced a correct report (50% classified). The duplicate IDs do not break Kraken2 itself. (`work/p7/t5_kraken.sh`)
- **Long-read mode completes and writes accounting** — `--long --all` ran through chopper fallback → minimap2 → bbduk → entropy → GDPR and produced 4 outputs plus a 9-row `accounting.tsv`.

## Findings

### F1. `--dry-run` deletes the outputs of a previous real run

- **Severity:** critical
- **Location:** `cerberus/pipelines/meta.py:76-77`, `cerberus/pipelines/profiling.py:139`, `cerberus/pipelines/long_read.py:41`, `cerberus/pipelines/long_read.py:97`, `cerberus/pipelines/gdpr.py:75-76`, `cerberus/pipelines/gdpr.py:126`
- **What:** Every pipeline clears its destination with `final.unlink(missing_ok=True)` *before* the `if not cfg.dry_run:` guard that wraps the actual `replace()`/concat. In `meta.py` the unlinks are lines 76-77 and the guard is line 78; the same one-line offset appears in all six sites. `concat_gz` likewise returns early on `cfg.dry_run` (`concat.py:31-35`) — but `profiling.py:139` already deleted the target by then.
- **Trigger:** `cerberus --dry-run -o <dir>` where `<dir>` holds results from an earlier successful run — i.e. the completely ordinary "let me check what this would do before I re-run it" action.
- **Consequence:** Proved: a directory holding 8 output FASTQs from a completed `--all` run went to **0 output FASTQs** after a single `--dry-run` invocation, and `reports/accounting.tsv` was overwritten with an all-zeros stub. The CLI then printed `Cerberus finished in 0.0s` and listed `dryrun_test/S1.meta.R1.fastq.gz` and `dryrun_test/S1.profiling.fastq.gz` as the run's outputs — files it had just deleted. Irrecoverable loss of hours-to-days of compute from a flag whose entire contract is "print commands without executing" (`cli.py:112-113`). The dry run also really creates files: `kraken.py:146-147` `r1.touch()`/`r2.touch()` execute before the `cfg.dry_run` early-return in `_gzip_inplace` (`kraken.py:157`).
- **Fix:** Move each `unlink(missing_ok=True)` inside the existing `if not cfg.dry_run:` block (or guard it with `if not cfg.dry_run`). Additionally, have `orchestrator.run()` refuse to touch `cfg.out_dir` at all under `--dry-run` — write the planned command list to a scratch location instead — and suppress the "finished / here are your outputs" banner in `cli.py:326-334` when `cfg.dry_run`.

### F2. `profiling.fastq.gz` and its GDPR release file contain 100% duplicate read IDs with no way to tell R1 from R2

- **Severity:** high
- **Location:** `cerberus/pipelines/profiling.py:137-145`, `cerberus/stages/align.py:76` (HEAD; working copy `align.py:171`)
- **What:** `samtools fastq -n` is used for every extraction, which suppresses the `/1`/`/2` suffix; `samtools fastq` also drops the header comment field. `profiling.py:140-145` then concatenates `ent.r1` + `ent.r2` + orphans into one file. Both mates of a fragment therefore land in the same file under a byte-identical name.
- **Trigger:** Any `--profiling` run (the standard path, `--fast`, and `--double-pass` all reach the same concat), and any `--profiling --gdpr` run.
- **Consequence:** Measured on real pipeline output: `S1.profiling.fastq.gz` = **630 records, 330 unique IDs, 300 IDs duplicated**; `S1.profiling.GDPR.fastq.gz` = **630 records, 330 unique IDs**. Input headers `@MICRO0000 1:N:0:AA` / `@MICRO0000 2:N:0:AA` both became `@MICRO0000`, so even the Illumina `1:N:` / `2:N:` disambiguator is gone — the mates are indistinguishable. Kraken2 and Bracken themselves cope (verified against a real DB: 12 records in, 12 classified, correct report), so the headline use case survives. What does not survive: (a) `kraken2 --output` emits duplicate read IDs, so KrakenTools `extract_kraken_reads.py` and any ID-keyed join become ambiguous — verified, all 6 IDs duplicated in `out.tsv`; (b) `kraken2 --unclassified-out` inherits them — verified, 6 records / 3 unique IDs, which is exactly how the GDPR pass builds `*.GDPR.fastq.gz`; (c) any dedup step silently halves the data — `seqkit rmdup -n` reported "6 duplicated records removed" on a 12-record file; (d) `<sample>.profiling.GDPR.fastq.gz` is the file the README (line 11) tells users to release publicly, and ENA/SRA submission validators reject FASTQ with non-unique read names. Kraken2 in `--paired` mode also outperforms merged single-end classification (it joins mates' k-mers), which the merge forfeits.
- **Fix:** Use `-N` instead of `-n` for the streams destined for concatenation — verified to produce `@FRAG000/1` / `@FRAG000/2` and 12/12 unique IDs on the merged file. Keep `-n` for the `--meta` R1/R2 pair, where identical names in separate files is the convention SPAdes/MEGAHIT expect. Alternatively suffix during concat, or document loudly that the profiling output is single-end-with-duplicate-IDs and unsuitable for deposit.

### F3. No per-stage host-removal, k-mer or entropy counts exist anywhere in `reports/` — the accounting cannot reconcile reads

- **Severity:** high
- **Location:** `cerberus/pipelines/meta.py:104-115`, `cerberus/pipelines/profiling.py:150-153`, `cerberus/pipelines/long_read.py:46`, `cerberus/pipelines/long_read.py:102`, `cerberus/pipelines/gdpr.py:38-101`
- **What:** The module docstring states "A reviewer wants to see exactly how many reads survive each stage" (`accounting.py:3`). What the pipelines actually feed `add_stage` is: meta → `{qc_paired, final_paired_r1, final_orphans}`; profiling → `{qc_paired, final}`; long-meta and long-profiling → `{final}` only; GDPR → **nothing at all** (`orchestrator.py:85-93` calls only `add_final`). There is no count between QC and the terminal output for any mode.
- **Trigger:** Every run.
- **Consequence:** Verified on a real `--all` run: from `accounting.tsv` a reviewer can compute the fastp delta (630→600 pairs) and the GDPR delta (300→300), and nothing else. Meta went 600 → 300 pairs, but whether those 300 were removed by minimap2 host alignment, by the entropy filter, or split between them is unrecoverable. Profiling ran minimap2/bowtie2 + bbduk aux k-mer + entropy and reports a single number. These three stages are the tool's entire scientific claim, and the per-stage attribution a reviewer would ask for ("how many reads did host removal actually take?") does not exist. The raw material *is* produced — `samtools flagstat` files record `460 + 0 mapped (38.33%)` per alignment — but they are written to `_work/<mode>/<mode>/<stage>/*.flagstat.txt` (`align.py:48`), a directory the README describes as transient, and are never parsed into the accounting. bbduk's `stats=` files (`entropy.py:35`, `entropy.py:65`) are likewise orphaned.
- **Fix:** Have each stage return its input/output counts (bbduk already prints "Input:"/"Result:" to its log; `samtools flagstat` already gives mapped/unmapped) and call `add_stage(mode, "<NN>_<stagename>_in"/"_out", …, file=…)` for each. At minimum, copy the flagstat and bbduk stats files into `reports/` so the numbers are reachable. Cheapest correct fix: parse the counts already sitting in the logs rather than re-decompressing FASTQs (see F11).

### F4. The `reads` column silently mixes fragments and reads, making the profiling row read as "nothing was removed"

- **Severity:** high
- **Location:** `cerberus/accounting.py:62-76`, `cerberus/orchestrator.py:129-135`, `cerberus/pipelines/profiling.py:151-152`
- **What:** `input_r1_reads` counts records in R1 (= fragments). `qc_paired` counts records in `qc.R1` (= fragments). `meta final_paired_r1` counts R1 records (= fragments). But `profiling final` counts records in the **merged R1+R2+orphans** file (= reads), and `qc_orphans` sums R1 and R2 orphans (= reads). All of these are emitted under one header column named `reads`, with no unit annotation and no denominator.
- **Trigger:** Every `--profiling` run.
- **Consequence:** Verified verbatim from a real run's `accounting.tsv`:
  ```
  _input     input_r1     630
  profiling  final        630
  ```
  A reviewer reads this as "the profiling pipeline removed zero reads." The truth is that 630 input fragments (1260 reads) were reduced to 300 surviving fragments (600 reads) plus 30 orphans = 630 records. Half the data was discarded and the table shows an unchanged number. The `_input input_r1` / `profiling final` coincidence is not a fluke of my test data — it is structural: whenever roughly half the pairs survive, the merged read count equals the input fragment count.
- **Fix:** Add a `unit` column (`fragments` | `reads`) or normalise everything to reads (`input_r1 + input_r2`). Emit `profiling final_fragments` alongside `profiling final_records`. Also emit an explicit `removed` column per row so the reader never has to subtract across differing units.

### F5. Re-running into an existing `-o` directory leaves stale outputs, and overwrites the only record that would identify them

- **Severity:** medium
- **Location:** `cerberus/orchestrator.py:96`, `cerberus/config.py:133-136`
- **What:** `ensure_directories()` only `mkdir(exist_ok=True)`s. Each pipeline unlinks only *its own* target names. Nothing enumerates or clears the output directory, and `accounting.write()` unconditionally truncates `accounting.tsv` / `accounting.json`.
- **Trigger:** Any second run into the same `-o` with a different mode set — including the natural "the profiling run was wrong, let me redo just `--meta`".
- **Consequence:** Verified. Run 1: `--all` → 8 outputs. Run 2 into the same directory: `--meta` → the directory *still* contains `S1.profiling.fastq.gz`, `S1.profiling.GDPR.fastq.gz`, `S1.meta.R1_GDPR.fastq.gz`, `S1.meta.R2_GDPR.fastq.gz`, `S1.meta.GDPR.fastq.gz` from run 1, all silently stale. The new `accounting.tsv` lists only the three meta files — so the one artifact that could have flagged the other five as orphaned has been erased. Worse in long mode: long outputs use entirely disjoint filenames (`S1.long_meta.fastq.gz` etc.), so a short-mode run followed by a long-mode run into the same directory leaves a complete, plausible, wrong short-read output set alongside the real one. `logs/cerberus.log.jsonl` is opened `mode="a"` (`logger.py:61`) and accumulates both runs with no run-id or delimiter — only the `Cerberus run: … modes=…` banner and a timestamp distinguish them.
- **Fix:** Write outputs into a run-stamped subdirectory, or refuse to start when `out_dir` contains cerberus outputs not covered by the current mode set unless `--force`/`--overwrite` is given. At minimum, have `accounting` list every `<sample>.*fastq.gz` found in `out_dir` and mark ones not produced by this run as `_stale`.

### F6. No provenance: nothing records tool versions, the resolved `TunedParams`, the command line, or the cerberus version

- **Severity:** medium
- **Location:** `cerberus/orchestrator.py:96`, `cerberus/config.py:130-131`, `cerberus/autotune.py:141`
- **What:** `reports_dir` receives exactly two files. Verified on a real run: `ls reports/` → `accounting.json`, `accounting.tsv`. Nothing else. `TunedParams.summary()` exists (`config.py:52-58`) and is called in exactly one place — `log.info("Tuned params: %s", tuned.summary())` at `autotune.py:141` — as free text into `logs/cerberus.log.jsonl`, and crucially *before* `apply_user_overrides` runs (`orchestrator.py:137-138`), so even the log never shows the final resolved parameter set in one place; line 162 logs only the diff dict. `__version__` appears only in `--version` and `doctor` output (`cli.py:118`, `cli.py:247`), never in any output artifact. `sys.argv` is read at `cli.py:289` and discarded. The reference `manifest.json` with its SHA256s (`refs.py:76`) is never copied to the output.
- **Trigger:** Every run.
- **Consequence:** A run is not reproducible and cannot be written up. To state a methods section a user must reconstruct: which cerberus version, which reference bundle checksums, which fastp/minimap2/bowtie2/bbduk/kraken2/samtools versions, what `entropy`/`bbduk_k`/`min_length`/`min_quality`/`minimap2_preset` autotune actually chose, and what the invocation was. None of it is on disk. The README's Citation section (line 261) promises a preprint; a reviewer asking "what parameters did you run with?" cannot be answered from the output directory. The `--dry-run` failure mode in F1 makes this worse: the safe way to inspect the resolved parameters is destructive.
- **Fix:** Write `reports/run_manifest.json` containing `cerberus.__version__`, `" ".join(sys.argv)`, an ISO timestamp, `asdict(cfg)` (paths stringified), `asdict(cfg.tuned)` *after* overrides, the reference `manifest.json` asset keys + SHA256s, and `{tool: version}` harvested from the tools already resolved by `require_tools`/`which`. The `doctor` subcommand (`cli.py:236-285`) already enumerates the tool list — reuse it.

### F7. The README `Outputs` section does not match what the pipeline writes

- **Severity:** medium
- **Location:** `README.md:239-255`, `README.md:11`, `README.md:231`, `cerberus/stages/qc.py:39-40`, `cerberus/stages/align.py:48`, `cerberus/config.py:122`
- **What / measured mismatches**, all confirmed against a real `--all` run:
  - README lists `reports/fastp.json/html`. Actual location: `_work/00_qc/fastp.json` and `_work/00_qc/fastp.html` (`qc.py:39-40`). Never copied to `reports/`.
  - README lists `reports/*.flagstat.txt`. Actual: 7 files under `_work/<mode>/<mode>/<stage>/*.flagstat.txt` (`align.py:48`). Never copied.
  - README does not mention `_work/` at all, although it is created inside `out_dir` (`config.py:122`).
  - README:231 claims run disk is "~2× input FASTQ size … cleaned automatically unless `--keep-intermediates`". Measured: `_work` = 1.1 MB against a 68 KB input = **~16×**, with **22 intermediate `.fq.gz` files retained**. Only BAMs are removed (`align.py:84`, `:129`, `:183`, `:225`); nothing ever deletes `_work`, and `keep_intermediates` gates nothing else (`grep` shows only `align.py` BAM unlinks and two `kmer.py` matched-read paths).
  - `<sample>.meta.GDPR.fastq.gz` (the GDPR scrub of the meta orphan stream, `gdpr.py:89` with `mode="meta"`) is produced but absent from the README listing.
  - All four long-read output names are absent from the README listing.
  - README:11 tells users GDPR outputs are `<sample>.<mode>.*_GDPR.fastq.gz`. Verified: in a real output directory, `ls *_GDPR.fastq.gz` matches `S1.meta.R1_GDPR.fastq.gz` and `S1.meta.R2_GDPR.fastq.gz` but **misses** `S1.meta.GDPR.fastq.gz` and `S1.profiling.GDPR.fastq.gz` — the latter being the primary public-release artifact for the profiling head.
- **Trigger:** Reading the docs; scripting `mv *_GDPR.fastq.gz release/`.
- **Consequence:** Users script against a glob that silently omits the aggressive profiling release file; users looking for the fastp report in `reports/` don't find it; users sizing disk for a 200 GB dataset budget 2× and need considerably more; users who read "cleaned automatically" never clean `_work` and fill the volume.
- **Fix:** Regenerate the README tree from the actual code (all 12 names, listed under F8), correct the `reports/` contents, document `_work/`, correct the disk multiplier, and either implement the promised `_work` cleanup when `not keep_intermediates` or drop the claim. Change the GDPR glob in README:11 to `*GDPR.fastq.gz`, or rename outputs so a single `_GDPR` glob is correct.

### F8. Long-read outputs use two different spellings of the same mode in one directory (`long_meta` vs `long-meta`)

- **Severity:** medium
- **Location:** `cerberus/pipelines/long_read.py:24`, `cerberus/pipelines/long_read.py:40`, `cerberus/pipelines/long_read.py:57`, `cerberus/pipelines/long_read.py:96`, `cerberus/pipelines/gdpr.py:89`, `cerberus/pipelines/gdpr.py:98`
- **What:** `PipelineResult.mode` is set to `"long-meta"` / `"long-profiling"` (hyphen, `long_read.py:24`/`:57`), but the pipelines build their own filenames with a literal underscore: `f"{cfg.sample_id}.long_meta.fastq.gz"` (`long_read.py:40`) and `f"{cfg.sample_id}.long_profiling.fastq.gz"` (`long_read.py:96`). The GDPR post-processor instead interpolates `pipeline_result.mode` verbatim (`gdpr.py:98`), so it emits the hyphen form.
- **Trigger:** `cerberus --long … --gdpr` (or `--all`).
- **Consequence:** Verified — a single output directory from one `--long --all` run:
  ```
  S1.long_meta.fastq.gz
  S1.long_profiling.fastq.gz
  S1.long-meta.long_GDPR.fastq.gz
  S1.long-profiling.long_GDPR.fastq.gz
  ```
  Any glob, sample-sheet, or Snakemake/Nextflow wildcard written against one spelling misses the other. The accounting keys inherit the split too (`long-meta.long_reads` vs `long-meta_gdpr.long_reads`). The complete producible set is: `{s}.meta.R1.fastq.gz`, `{s}.meta.R2.fastq.gz`, `{s}.meta.orphans.fastq.gz`, `{s}.profiling.fastq.gz`, `{s}.meta.R1_GDPR.fastq.gz`, `{s}.meta.R2_GDPR.fastq.gz`, `{s}.meta.GDPR.fastq.gz`, `{s}.profiling.GDPR.fastq.gz`, `{s}.long_meta.fastq.gz`, `{s}.long_profiling.fastq.gz`, `{s}.long-meta.long_GDPR.fastq.gz`, `{s}.long-profiling.long_GDPR.fastq.gz`.
- **Fix:** Pick one separator and route every filename through a single helper — `fastq.output_name()` (`fastq.py:63`) already exists for precisely this and would have prevented the divergence (see F12). Suggest `long_meta`/`long_profiling` everywhere, with `PipelineResult.mode` matching.

### F9. The CLI completion block never shows the GDPR outputs, nor meta R2/orphans

- **Severity:** medium
- **Location:** `cerberus/cli.py:326-334`, `cerberus/orchestrator.py:99-103`, `cerberus/pipelines/base.py:18-20`
- **What:** `orchestrator.run()` returns `{"outputs": {r.mode: r.primary_output …}}` — one path per pipeline, and `primary_output` is `paired_r1 or singletons or long_reads` (`base.py:20`). The GDPR results are never added to the summary at all (`orchestrator.py:85-93` only feeds `accounting`), and `GDPRResult` never reaches the return value.
- **Trigger:** Any run; most visibly `--all`.
- **Consequence:** Verified on a real `--all` run that produced 8 files, the terminal showed:
  ```
    meta                out/S1.meta.R1.fastq.gz
    profiling           out/S1.profiling.fastq.gz
    reports            out/reports
  ```
  `S1.meta.R2.fastq.gz`, `S1.meta.orphans.fastq.gz` and all four GDPR files are invisible. `--gdpr` is the flagship "publication-ready, zero human reads" feature (README:11) and a user running `--all` is told nothing about where its outputs went. (Minor: the `reports` line is misaligned — `{mode:<18s}` for modes vs a hard-coded 12-space pad for `reports`, visible above.)
- **Fix:** Have `run()` return `{"outputs": {...all paths...}, "gdpr_outputs": {...}}` — the data is already in `accounting.final_outputs`, so return that structure directly — and print every non-`None` path grouped by mode. Fix the `reports` padding to `{'reports':<18s}`.

### F10. `concat_gz` writes a zero-byte `.fastq.gz` that is not a valid gzip, and can crash before creating its parent directory

- **Severity:** medium
- **Location:** `cerberus/stages/concat.py:21-26`, `cerberus/stages/concat.py:28`
- **What:** When every input is missing, `concat.py:25` does `output.write_bytes(b"")`. A zero-byte file is not a gzip stream. Separately, the early return at line 26 happens *before* `output.parent.mkdir(parents=True, exist_ok=True)` at line 28, so the degenerate path assumes the parent already exists.
- **Trigger:** Any upstream stage that produces no file (a tool that exits 0 with no output, a stage whose inputs were all filtered away and whose writer was skipped, a partially-cleaned `_work`). Reachable at `profiling.py:113` and `profiling.py:145`, and `meta.py:55`.
- **Consequence:** Verified: `gzip -t empty.fastq.gz` → `unexpected end of file`, rc=1; `pigz -t` → `skipping: empty`, rc=1; `zcat` → rc=1; `pigz -dc` → rc=1. Kraken2 pipes input through `gzip -dc` (confirmed in the kraken2 2.17.1 wrapper, lines 148-172) and does not check that subprocess's exit status, so it would silently classify zero reads rather than failing. `bbduk.sh` and `seqkit` tolerate it (both reported 0 reads). Cerberus's own `count_reads` short-circuits to 0. Net effect: a corrupt `<sample>.profiling.fastq.gz` propagates through the GDPR pass and is delivered to the user as a "successful" output, and only `gzip -t` reveals it. The codebase already knows the right answer — `kraken.py:175-179` `_write_empty_gzip()` writes a genuine 38-byte empty gzip stream with the comment "so bbduk/seqkit can read it" — `concat_gz` just doesn't use it. The mkdir ordering bug is confirmed: `concat_gz(inputs=[], output=d/"nope"/"x.fastq.gz")` raises `FileNotFoundError`, while the same call with an existing parent silently produces a 0-byte file.
- **Fix:** In the empty branch, `output.parent.mkdir(parents=True, exist_ok=True)` first, then write a valid empty gzip (reuse `kraken._write_empty_gzip`, or `gzip.open(output,'wb').close()`) if `output.suffix == ".gz"`. Better: raise instead of warning when the caller expected inputs — a missing stage output is a bug, not a benign empty result.

### F11. Accounting costs 18 full FASTQ decompressions per run, 5 of them exact duplicates, for numbers fastp already reports

- **Severity:** medium
- **Location:** `cerberus/orchestrator.py:129-135`, `cerberus/pipelines/meta.py:109-114`, `cerberus/pipelines/profiling.py:151-152`, `cerberus/accounting.py:51`
- **What:** Counts are obtained by decompressing whole files with `count_reads` rather than by reading numbers the tools already emitted. Several files are counted more than once because the pipeline's own `stats` dict and `RunAccounting.add_final` each call `count_reads` independently.
- **Trigger:** Every non-dry run.
- **Consequence:** Instrumented a real `--all` run by wrapping `count_reads`: **18 invocations over 13 distinct files**, with `qc.R1.fq.gz` decompressed **3×** (`orchestrator.py:131`, `meta.py:110`, `profiling.py:151`) and `S1.meta.R1`, `S1.meta.orphans`, `S1.profiling` decompressed **2×** each. Both raw inputs are decompressed in full (`orchestrator.py:129-130`) purely to obtain a number that `fastp.json` already contains under `summary.before_filtering.total_reads`. On a 50 GB paired library this is roughly a dozen extra full-corpus decompression passes bolted onto the pipeline for reporting alone.
- **Fix:** Read input and post-QC counts from `qc.json_report` (already parsed by `autotune_from_fastp`). Memoise `count_reads` on `(path, st_size, st_mtime_ns)`. Have `add_final` reuse the count the pipeline already put in `PipelineResult.stats` instead of recounting.

### F12. Dead code: `fastq.output_name()` is never called, and `concat.compress_to()` is never called *and* is broken

- **Severity:** low
- **Location:** `cerberus/utils/fastq.py:63-65`, `cerberus/stages/concat.py:47-59`
- **What:** `grep -rn output_name` over the whole repo (including `tests/`) returns exactly one hit: the definition. Same for `compress_to`. Beyond being unused, `compress_to`'s pigz branch is wrong: `run(["pigz","-c","-p",N,src], log_path=…)` — `run()` redirects the child's stdout into the log file (`shell.py:71-82`), so `pigz -c`'s compressed bytes go into `{tag}.pigz.log` and `dst` is never created. Verified: the function returned `/tmp/…/x.fq.gz`, `dst.exists()` was `False`, and the directory contained only `x.fq` and an 84-byte `t.pigz.log`. The fallback branch (no pigz) writes `dst` correctly, so the bug only bites where pigz *is* installed — i.e. the documented environment (`environment.yml`, `doctor` lists pigz as required).
- **Trigger:** Nothing today; a future caller of `compress_to` on a machine with pigz.
- **Consequence:** `output_name()` documents itself as "Canonical output filename layout" — had it been used it would have prevented F8's `long_meta`/`long-meta` split, so its deadness is not cosmetic. `compress_to` is a loaded gun: the next person to use it gets a silently missing output and a log file full of binary.
- **Fix:** Delete both, or fix `compress_to` (`run(..., stdout_path=dst)` once `shell.run` supports it, or use `pipe([...], final_stdout=dst)` which already does — `shell.py:117-118`) and route every output filename through `output_name()`.

### F13. Long-read accounting emits a spurious `input_r1 0` row and mislabels the long-read QC count as `qc_paired`

- **Severity:** low
- **Location:** `cerberus/accounting.py:63`, `cerberus/accounting.py:32`, `cerberus/orchestrator.py:178-179`
- **What:** `write()` emits the `_input input_r1` row unconditionally (line 63), while `input_r2` and `input_long` are conditional (lines 64-67). In long mode `input_r1_reads` is never assigned. Separately, `_run_long` stores the post-QC long-read count into the field named `qc_paired` (`orchestrator.py:179`), because that is the only survivor field available.
- **Trigger:** Any `--long` run.
- **Consequence:** Verified `accounting.tsv` from a real long run:
  ```
  _input   input_r1     0
  _input   input_long   120
  _qc      qc_paired    120
  ```
  `input_r1 0` invites the reading "R1 was empty"; `qc_paired` is meaningless for unpaired long reads. Combined with F3 (long mode records only `final`), the entire long-read accounting is three rows, one of which is noise and one of which is mislabelled.
- **Fix:** Make the `input_r1` row conditional like its siblings. Rename `qc_paired` to `qc_survived` (or emit `qc_long`), or make the row label mode-dependent.

### F14. `accounting.tsv` records whatever path string the CLI was given, and never records the sample ID

- **Severity:** low
- **Location:** `cerberus/accounting.py:52`, `cerberus/accounting.py:62-76`
- **What:** `add_final` stores `str(path)` verbatim, so a run invoked with `-o out` produces relative paths; `sample_id` is a field on `RunAccounting` (`accounting.py:28`) and appears in `accounting.json`, but `write()` never emits it into the TSV.
- **Trigger:** Relative `-o`; any attempt to collate TSVs from several samples.
- **Consequence:** Verified: `-o out` produced `_final meta.paired_r1 300 out/S1.meta.R1.fastq.gz`, interpretable only from the original working directory and broken if the results are moved or archived. A user `cat`-ing several samples' `accounting.tsv` into one table has no column identifying which sample each row belongs to — the file name is the only carrier.
- **Fix:** Store `path.resolve()` (or a path relative to `out_dir`, documented as such), and either add a `sample` column or a `# sample_id: …` comment header.

### F15. `logs/cerberus.log.jsonl` accumulates every run with no run identifier

- **Severity:** low
- **Location:** `cerberus/utils/logger.py:61`, `cerberus/utils/logger.py:13-18`
- **What:** The file handler opens `mode="a"` and the JSONL payload carries `ts`/`level`/`logger`/`msg` — no run id, no pid, no invocation marker.
- **Trigger:** Second and subsequent runs into the same `-o`.
- **Consequence:** Verified: after two runs the file contained two `Cerberus run: …` banners with nothing between them to delimit the boundary. Since this log is the *only* place the resolved tuned parameters appear (F6), a user trying to recover "what parameters produced the outputs currently on disk" has to guess which banner corresponds to the surviving files — and after a mixed-mode re-run (F5) the answer may be "several of them".
- **Fix:** Emit a `run_id` (uuid4 or timestamp) into every JSONL record and into `reports/run_manifest.json`, or rotate to `logs/<run_id>.jsonl`.
