# Pass 6 — Pipeline composition

## Summary

I enumerated all 36 supported invocations ({meta, profiling, meta+profiling} × {default, --fast, --double-pass} × {short, --long} × {±--gdpr}) by driving `orchestrator._required_pipeline_keys` and `refs.required_assets_for` directly against the real modules. Nine of the 36 cells are broken and the remainder are degraded in at least one way: the `--long --profiling --fast` combination emits a pipeline key (`long-profiling-fast`) that does not exist in `refs._PIPELINE_TO_ASSETS`, so `RefManager.ensure()` is handed an empty list and **zero reference assets are downloaded**, while `run_long_profiling` never reads `cfg.fast` at all; and the `--long --double-pass` winnowmap branch points `-W` at a `.meryl` database that exists in no manifest, no build script, and no download path — and it is the *default* branch on the shipped conda environment because `environment.yml` omits `fastplong`, forcing the chopper fallback whose synthetic JSON hard-codes a 5000 bp mean length that always classifies as `VERY_LONG` and therefore always sets `winnowmap_enabled=True`.

Beyond the long-read lanes, three cross-cutting defects affect every cell. `utils.shell.pipe()` inspects only the last process's return code, so an aligner that dies mid-stream produces a silently truncated BAM that the orchestrator accepts as success — I proved this with a live `pipe()` call. The `drop_strategy` parameter that supposedly distinguishes the "conservative" meta head from the "aggressive" profiling head cannot express a second samtools flag at all (`_filter_flags` returns a two-key dict spliced into exactly two argv slots), so `"both"` and `"either"` resolve to the same effective pair semantics and meta is not conservative. And `_work/` is never deleted by anything: for a 50 GB input under `--all` I estimate ~575 GB peak and ~430 GB left on disk, against a README claim of "~2× input FASTQ size … cleaned automatically".

On the specific questions asked: the R1/R2 files handed to bowtie2 under `--double-pass` **are** synchronised (because `-s /dev/null` discards half-mapped survivors), fastp's `--unpaired1/2` files **are** always created, so the `.exists()` gate at `qc.py:70-71` is always true and every run pays for a full orphan alignment pass over zero reads, and `PipelineResult.as_inputs_for_gdpr()` is dead code with zero call sites.

## Mode matrix

Verdict columns are for the cell *without* `--gdpr` and *with* `--gdpr`. "Assets" = what `refs.ensure()` actually downloads for that cell (proved by running `required_assets_for`). All cells inherit F3 (`pipe()` masks aligner failure), F5 (`_work/` never cleaned) and F6 (`--dry-run` deletes prior outputs); these are not repeated per row.

Asset abbreviations: **M** = `masked_t2t_hla_minimap2`, **B** = `masked_t2t_hla_bowtie2`, **A** = `aux_refs`, **K** = `kraken2_gdpr_compact`.

| # | Platform | Modes | Speed | Stage order actually executed | Assets ensured | No `--gdpr` | With `--gdpr` |
|---|---|---|---|---|---|---|---|
| 1 | short | meta | default | fastp → autotune → mm2 paired (`-f 12`) → bbduk entropy paired → move to `out/`; **always** + orphan-cat → mm2 singles → entropy single | M (+K) | DEGRADED — F4 (meta not conservative), F8 (orphan pass over zero reads) | DEGRADED — also F9 (`.meta.GDPR.fastq.gz` undocumented, normally 0 reads) |
| 2 | short | meta | `--fast` | identical to #1 — `cfg.fast` is read only at `profiling.py:49/89/120` | M (+K) | DEGRADED — F12, flag silently ignored, no warning | DEGRADED — F12 + F9 |
| 3 | short | meta | `--double-pass` | identical to #1 — `meta.py` has no `double_pass` branch | M (+K) | DEGRADED — F12, flag silently ignored | DEGRADED — F12 + F9 |
| 4 | short | profiling | default | fastp → bowtie2 `--very-sensitive-local` (`-f 12`) → bbduk k-mer (aux) → entropy → concat(R1,R2,orphans); orphans: cat → mm2 singles → bbduk → entropy | M,B,A (+K) | DEGRADED — F10 (duplicate read names in merged FASTQ), F8 | DEGRADED — output name matches README |
| 5 | short | profiling | `--fast` | fastp → mm2 paired (`-f 4`) → entropy → concat; bbduk aux skipped on both streams | M (+K) | OK (coherent; F10/F8 apply) | OK |
| 6 | short | profiling | `--double-pass` | fastp → mm2 paired pre (`-f 4`) → bowtie2 (`-f 12`) → bbduk k-mer → entropy → concat | M,B,A (+K) | DEGRADED — pre-filter and bowtie2 apply the *same* pair rule (F4), so the pass is near-redundant; pairs **are** synchronised | DEGRADED |
| 7 | short | meta+profiling | default | #1 then #4, both re-reading `qc.r1/qc.r2` from scratch | M,B,A (+K) | DEGRADED | DEGRADED — two GDPR result sets, 3 output name families (F9) |
| 8 | short | meta+profiling | `--fast` | #1 then #5 | M (+K) | DEGRADED — F12 for the meta head | DEGRADED |
| 9 | short | meta+profiling | `--double-pass` | #1 then #6 | M,B,A (+K) | DEGRADED — F12 for the meta head | DEGRADED |
| 10 | long | meta | default | chopper/fastplong → autotune → mm2 singles `-ax map-ont` **against the `-x sr` index** → entropy → `.long_meta.fastq.gz` | M (+K) | DEGRADED — F7 (wrong index preset) | DEGRADED — F7 + F9 (`.long-meta.long_GDPR.fastq.gz`) |
| 11 | long | meta | `--fast` | identical to #10 | M (+K) | DEGRADED — F7 + F12 | DEGRADED |
| 12 | long | meta | `--double-pass` | identical to #10 — `run_long_meta` has no `double_pass` branch | M (+K) | DEGRADED — F7 + F12 | DEGRADED |
| 13 | long | profiling | default | chopper → mm2 singles `map-ont` on sr index → bbduk aux → entropy → `.long_profiling.fastq.gz` | M,A (+K) | DEGRADED — F7 | DEGRADED — F7 + F9 |
| 14 | long | profiling | `--fast` | **`refs.ensure([])` — nothing downloaded.** Then the *full non-fast* path runs anyway (mm2 + bbduk aux + entropy) | **(none)** | **BROKEN** — F1: cold cache ⇒ minimap2 cannot open the index; warm cache ⇒ `--fast` silently does nothing | **BROKEN** — F1: `long-gdpr` rescues M and K, but `aux_refs` is still never ensured ⇒ bbduk fails |
| 15 | long | profiling | `--double-pass` | `winnowmap -W <ref_dir>/masked_t2t_hla.meryl` whenever `read_length_class == VERY_LONG` — which is *always* on the shipped env | M,A (+K) | **BROKEN** — F2: meryl DB exists in no manifest, no build script, and never on disk | **BROKEN** — F2 |
| 16 | long | meta+profiling | default | #10 then #13 | M,A (+K) | DEGRADED — F7 | DEGRADED — F7 + F9; this is the README quick-start `cerberus --long … --all` |
| 17 | long | meta+profiling | `--fast` | #10 (works, M came from `long-meta`) then #14 | M (+K) | **BROKEN** — F1: `aux_refs` never ensured ⇒ bbduk fails in the profiling head | **BROKEN** — F1 |
| 18 | long | meta+profiling | `--double-pass` | #10 then #15 | M,A (+K) | **BROKEN** — F2 | **BROKEN** — F2 |

`--fast` + `--double-pass` together is correctly rejected at `orchestrator.py:51-52`, so there is no fourth speed column. There are no unreachable cells other than that.

## Verified working

- **Mode-selection validation** — `--gdpr` alone, no mode at all, `--long` with `-r1/-r2`, `--long` without `-i`, and `--fast --double-pass` are all rejected with clear messages (`orchestrator.py:31-56`); exercised by `tests/test_cli.py:39-79`.
- **`_required_pipeline_keys` short-read branch** — all 18 short-read cells produce keys present in `refs._PIPELINE_TO_ASSETS`; verified by driving the real function over the full cross-product (`orchestrator.py:106-115`, `refs.py:37-45`).
- **Asset minimisation for `--fast`** — `--profiling --fast` correctly narrows to `profiling-fast` → `masked_t2t_hla_minimap2` only, and `profiling.py:89` / `profiling.py:120` correctly skip the aux k-mer pass on both the paired and orphan streams, so no unfetched asset is ever opened in that lane.
- **Bowtie2 pair synchronisation under `--double-pass`** — `minimap2_paired` passes `-1/-2/-0 /dev/null/-s /dev/null` (`align.py:68-78`); samtools `fastq` routes any read whose mate was filtered out to the `-s` file, so `pre.r1`/`pre.r2` at `profiling.py:71-72` are 1:1 and in input order. Had they desynchronised, bowtie2 aborts non-zero ("fewer reads in file specified with -1 than -2"), which `shell.run(check=True)` at `align.py:169` would surface as a `ToolError` — a loud failure, not silent corruption. (samtools is not installed on this box — `command -v samtools` → missing — so this is derived from the documented `-s`/name-collation semantics, not executed.)
- **GDPR mode dispatch shape** — `run_gdpr_for` correctly branches on which fields `PipelineResult` populated: meta fills `paired_r1/r2` + `singletons` and gets both branches; profiling fills only `singletons`; long modes fill only `long_reads` (`gdpr.py:55/83/92`, `meta.py:95-101`, `profiling.py:147-154`, `long_read.py:44-47/100-103`).
- **Kraken2 output discovery** — `_gzip_pair_outputs` globs both `<root>_1.fq` and `<root>1.fq` layouts and falls back to valid empty gzips rather than crashing (`kraken.py:129-152`, `kraken.py:175-179`); confirmed a Python-written empty gzip is 35 bytes and `count_reads` returns 0 on it.
- **bbduk empty-input guards** — `kmer.py:52` and `kmer.py:94` short-circuit on `st_size <= 40`, which does catch a 35-byte valid empty gzip, and emit valid empty gzip passthroughs so the chain does not break.
- **`concat_gz` empty handling** — filters `None`/missing inputs and writes an empty output with a warning rather than raising (`concat.py:21-26`); gzip-member concatenation is a valid gzip stream.
- **Ref `.tmp` cleanup on start** — `cleanup_partial` removes leftover aborted downloads before `ensure()` runs (`orchestrator.py:64-66`, `refs.py:238-244`).
- **`--fast`/`--double-pass` mutual exclusion** — enforced before any work is done (`orchestrator.py:51-52`), so no cell can request both lanes.
- **Autotune override precedence** — `apply_user_overrides` is applied after `autotune_from_fastp` in both `_run_short` and `_run_long` (`orchestrator.py:137-139`, `181-183`), so explicit user flags win in every cell.

## Findings

### F1. `--long --profiling --fast` requests a pipeline key that does not exist; no references are fetched and `--fast` is never honoured

- **Severity:** critical
- **Location:** `cerberus/orchestrator.py:112`, `cerberus/refs.py:37-45`, `cerberus/refs.py:109`, `cerberus/pipelines/long_read.py:50-103`
- **What:** `_required_pipeline_keys` builds `f"{prefix}profiling-fast"` with `prefix = "long-"`, producing the string `"long-profiling-fast"`. `_PIPELINE_TO_ASSETS` (`refs.py:37-45`) contains `profiling-fast` but **no** `long-profiling-fast`. `required_assets_for` looks the key up with `_PIPELINE_TO_ASSETS.get(pkey, [])` (`refs.py:109`) — an unknown key silently contributes nothing and raises no warning. Proved by running the real code: `RefManager.required_assets_for(['long-profiling-fast'])` → `[]`. Meanwhile `run_long_profiling` contains no reference to `cfg.fast` whatsoever (grep for `cfg.fast` across `cerberus/pipelines/` matches only `profiling.py:49`, `:89`, `:120`), so it runs the *full* long-profiling path — minimap2, then `bbduk` against `aux_refs` at `long_read.py:81-88`, then entropy — regardless.
- **Trigger:** `cerberus --long -i reads.fq.gz --profiling --fast` (matrix cell 14), or with `--meta` added (cell 17), or with `--gdpr` added.
- **Consequence:** On a cold reference cache `cerberus --long --profiling --fast` downloads **nothing at all** (`refs.ensure([])` is a no-op), logs `Required ref-asset groups: ['long-profiling-fast']`, then hands minimap2 a nonexistent `.mmi` path — a runtime failure inside a `pipe()` whose first-stage return code is discarded (see F3), so the run may instead continue with an empty BAM. Adding `--gdpr` pulls in `masked_t2t_hla_minimap2` and `kraken2_gdpr_compact` via the `long-gdpr` key, which masks the minimap2 failure but leaves `aux_refs` unfetched, so `bbduk.sh ref=<missing>` fails later. On a warm cache populated by an earlier non-fast run, everything "works" but `--fast` is a complete no-op: the user pays full cost and gets no warning. The pre-flight warning at `orchestrator.py:53-54` requires `cfg.meta and cfg.profiling`, so the single-mode broken case is never warned about.
- **Fix:** Add `"long-profiling-fast": ["masked_t2t_hla_minimap2"]` to `_PIPELINE_TO_ASSETS`; make `required_assets_for` raise `RefManagerError` on an unknown pipeline key instead of `.get(..., [])`; and either honour `cfg.fast` in `run_long_profiling` (skip the bbduk aux pass, matching `profiling.py:89`) or reject `--fast` with `--long --profiling` in `validate_config`.

### F2. `--long --double-pass` winnowmap branch points at a meryl database that is never downloaded or built — and it is the default branch on the shipped environment

- **Severity:** critical
- **Location:** `cerberus/pipelines/long_read.py:64-71`, `cerberus/stages/align.py:208-215`, `cerberus/data/default_manifest.json`, `environment.yml`
- **What:** `long_read.py:66` derives the winnowmap repetitive-k-mer database as `meryl_db = idx.with_suffix(".meryl")`. `idx` is `refs.path_to(refs.asset("masked_t2t_hla_minimap2"))`; running the real `RefManager` gives `<ref_dir>/masked_t2t_hla.mmi`, so `meryl_db` resolves to `<ref_dir>/masked_t2t_hla.meryl`. That file appears in **no** manifest asset (the manifest has exactly five keys: `masked_t2t_hla_minimap2`, `masked_t2t_hla_bowtie2`, `kraken2_gdpr_compact`, `aux_refs`, `human_kmer_set`), in **no** `_PIPELINE_TO_ASSETS` entry, and in **no** build script (`grep -rn meryl scripts/ README.md` → zero hits). `winnowmap_singles` passes it to `-W` unconditionally (`align.py:210`) with no existence check.
- **Trigger:** `cerberus --long -i reads.fq.gz --profiling --double-pass` whenever `tuned.winnowmap_enabled` is true. That flag is set by `_BASE_PARAMS[ReadLengthClass.VERY_LONG]` (`autotune.py:71`), reached when the autotuned mean read length is ≥ 5000. Crucially, `environment.yml` installs `chopper` and `winnowmap` but **not** `fastplong`, so every `--long` run on the documented conda environment takes the chopper fallback at `qc.py:115-129`, which writes `_synthetic_long_json` with a hard-coded `"read1_mean_length": 5000` (`qc.py:146`). I ran `classify_length(5000.0)` → `VERY_LONG` and `autotune_from_fastp` on that synthetic JSON → `winnowmap_enabled=True`. So on the shipped environment the winnowmap branch is taken **100% of the time** with `--long --profiling --double-pass`, irrespective of the actual read lengths.
- **Consequence:** `winnowmap -W <ref_dir>/masked_t2t_hla.meryl -ax map-ont …` fails immediately (cannot open the k-mer file). Because it is the first stage of a `pipe()`, its non-zero exit is discarded (F3) and `samtools view` writes a header-only BAM that exits 0 — so the run continues, `samtools fastq` yields an empty FASTQ, bbduk is skipped by the `st_size <= 40` guard, and the user receives a **zero-read `<sample>.long_profiling.fastq.gz`** with no error. Secondarily, `winnowmap_singles` hard-codes `-ax map-ont` (`align.py:211`) and ignores `tuned.minimap2_preset`, so PacBio HiFi/CLR input would be mapped with ONT parameters even if the meryl DB existed.
- **Fix:** Add a `masked_t2t_hla_meryl` asset to the manifest and to the `long-profiling` asset list, build it in `mask_t2t_hla.sh` / `build_custom_host_ref.sh` (`meryl count k=15 … && meryl print greater-than distinct=0.9998`); gate the branch on `meryl_db.exists()` with a fallback to minimap2 and a warning; drive the winnowmap preset from `tuned.minimap2_preset`; and stop the chopper fallback from fabricating a 5000 bp mean length (compute it from the actual output, or add `fastplong` to `environment.yml`).

### F3. `pipe()` only checks the last process's exit code, so an aligner crash silently yields a truncated BAM

- **Severity:** critical
- **Location:** `cerberus/utils/shell.py:132-139`; consumed at `cerberus/stages/align.py:64`, `:121`, `:169`, `:217`, and `cerberus/stages/qc.py:118`
- **What:** `pipe()` waits on every process but computes `rc = procs[-1].returncode` and raises only on that. Every alignment stage in every mode is built as `[aligner, samtools view -b -o BAM -]`, so the aligner is `procs[0]` and its return code is discarded. When the aligner dies, `samtools view` simply sees EOF and exits 0. I proved this by running the real `pipe()` with `[["bash","-c","printf 'rec1\\nrec2\\nrec3\\n'; exit 1"], ["cat"]]`: it returned `rc=0` and left the three partial records in the downstream file.
- **Trigger:** Any minimap2/bowtie2/winnowmap failure or partial failure — missing index (F1), missing meryl DB (F2), OOM kill during index load (the `.mmi` is 7.7 GB against a documented 16 GB laptop target), disk-full mid-run, or a malformed `--minimap2-args` token.
- **Consequence:** The pipeline continues with a truncated or empty BAM and produces a FASTQ that *looks* like a successful decontamination result but contains an arbitrary prefix of the data — or, for a header-only BAM, zero reads. Under `--gdpr` this silently invalidates the headline "zero detectable human reads" claim (README:11): unaligned host reads are never given the chance to be dropped, and the resulting file is published as GDPR-clean. Nothing downstream detects it, because `accounting` records the (small but non-zero) counts as if they were the real survivorship.
- **Fix:** Collect all return codes in `pipe()` and raise `ToolError` if **any** is non-zero (a SIGPIPE/`-13` on an upstream process when a downstream one exits early is the only case worth whitelisting), and add a post-alignment sanity check comparing input read count to `samtools flagstat` total.

### F4. `drop_strategy` cannot express a second samtools flag — meta's "conservative" pair rule is not implemented

- **Severity:** high
- **Location:** `cerberus/stages/align.py:230-242`, consumed at `align.py:67-78`; callers `meta.py:46`, `profiling.py:57`, `profiling.py:69`, `gdpr.py:71`
- **What:** `_filter_flags` returns a two-key dict `{"filter_flag", "filter_value"}` and `align.py:70` splices exactly those two strings into argv. The structure physically cannot emit a second flag pair, so the `-F` that the comment at `align.py:239-240` describes ("Implemented by `-f 4 -F 8` below") is never produced. `"both"` → `-f 12` (both mates unmapped). `"either"` → `-f 4` (this read unmapped), but combined with `-s /dev/null` and `-0 /dev/null` (`align.py:74-75`) any read whose mate was filtered out is classified by samtools as a singleton and discarded. The net effect of `-f 4 -s /dev/null` is therefore identical to `-f 12`: **the two strategies produce byte-identical output**. Separately, had the comment's `-f 4 -F 8` actually been emitted it would have kept *only* reads whose mate mapped — the exact inverse of the intent.
- **Trigger:** Every `--meta` run, short or long. `meta.py:46` passes `drop_strategy="both"`.
- **Consequence:** `meta.py:4` documents "drop pairs where BOTH mates map" and README:9 sells meta as "Conservative — retains microbial reads even at the cost of some residual host". The implementation drops any pair where *either* mate maps — the aggressive rule. Half-mapped pairs (a genuine microbial mate anchored to a chimeric or repeat-region host hit) are silently thrown away, and the meta and profiling heads differ only in aligner choice and the presence of the bbduk aux pass, not in pair semantics. This is the central differentiator of the product and it does not exist. It also makes the `--double-pass` minimap2 pre-filter (`profiling.py:63-70`, `"either"`) apply exactly the same rule bowtie2 applies afterwards, so the extra pass buys almost nothing beyond the aligner difference — consistent with the docstring's own "the marginal yield is small" (`profiling.py:14-17`) but for the wrong reason.
- **Fix:** Return a list of flag tokens rather than a fixed pair, and implement `"both"` as `-F 8` ∪ `-f 4` semantics correctly — i.e. keep a pair if *at least one* mate is unmapped, which requires either `samtools view -e '!flag.unmap || !flag.munmap'` piped into `samtools fastq`, or `-f 4` with the singleton mate recovered via `-s <orphans>` and re-paired. Fix the misleading comment at `align.py:236-241` in the same change.

### F5. `_work/` is never deleted; peak disk for a 50 GB input is ~575 GB, not the documented ~2×

- **Severity:** high
- **Location:** `cerberus/config.py:122-123`, `cerberus/config.py:133-136`, `cerberus/stages/align.py:84` (the only cleanup), `README.md:231`
- **What:** `--keep-intermediates` is consulted in exactly two places: BAM deletion in the four aligner wrappers (`align.py:84`, `:129`, `:183`, `:225`) and the optional bbduk `outm` matched-read files (`kmer.py:58-59`, `:105`). A repo-wide grep for `rmtree` finds only `import shutil` in `refs.py:22` and `shutil.which` in `shell.py:34` — **nothing ever removes `out/_work/`**. Every intermediate FASTQ from every stage of every mode survives the run. Only the four `.replace()` moves (`meta.py:79-80`, `meta.py:93`, `long_read.py:43`, `long_read.py:99`, `gdpr.py:78-79`, `gdpr.py:128`) actually remove a file from `_work/`, and only the terminal stage of each chain.
- **Trigger:** Any real run. The larger the input, the worse.
- **Consequence:** For a 50 GB gzipped paired input under `cerberus --all` (assuming ~90% survivorship per stage and a ~3.5× gzip→plain expansion for FASTQ):

  | `_work/` path | contents | GB |
  |---|---|---|
  | `_work/00_qc/` | `qc.R1.fq.gz` + `qc.R2.fq.gz` | ~45 |
  | `_work/meta/meta/01_minimap2_paired/` | unmapped R1+R2 | ~44 |
  | `_work/meta/meta/03_entropy_paired/` | moved out by `.replace()` | 0 |
  | `_work/profiling/profiling/02_bowtie2/` | unmapped R1+R2 | ~44 |
  | `_work/profiling/profiling/03_bbduk_kmer/` | k-mer-clean R1+R2 | ~43 |
  | `_work/profiling/profiling/04_entropy/` | entropy R1+R2 (copied, not moved, by `concat_gz`) | ~43 |
  | `_work/gdpr_meta/…/01_kraken2/` | unclassified pair, after `pigz` | ~42 (**~147 transient, uncompressed**) |
  | `_work/gdpr_profiling/…/03_singletons/` | unclassified merged, after `pigz` | ~42 (**~147 transient, uncompressed**) |
  | **`_work/` persistent total** | | **~259** |

  Plus ~170 GB of finals in `out/` (`meta.R1/R2`, `profiling`, `meta.R1/R2_GDPR`, `profiling.GDPR`) → **~429 GB left on disk after a successful run**. The transient peak adds ~147 GB because Kraken2 writes `--unclassified-out` as **plain uncompressed FASTQ** (`kraken.py:48`, `kraken.py:100`) and only gzips it afterwards in `_gzip_inplace` (`kraken.py:118`, `:152`) → **~575 GB peak, ~11.5× input**. `--double-pass` adds another ~44 GB for `01_minimap2_pre`. README:231 states "~2× input FASTQ size during processing; cleaned automatically unless `--keep-intermediates`" — wrong by roughly 6× on the persistent figure and ~12× at peak, and the "cleaned automatically" half is simply false. On the "16 GB laptop" the README targets (README:233), this fills any consumer SSD.
- **Fix:** Delete each stage directory once its successor has completed (or `shutil.rmtree(cfg.work_dir)` at the end of `orchestrator.run` when `not cfg.keep_intermediates` and the run succeeded); stream Kraken2's unclassified output through `pigz` rather than materialising it uncompressed; and correct README:231.

### F6. `--dry-run` deletes pre-existing output files without replacing them

- **Severity:** high
- **Location:** `cerberus/pipelines/meta.py:76-80`, `meta.py:91-93`, `cerberus/pipelines/profiling.py:139-145`, `cerberus/pipelines/long_read.py:41-43`, `long_read.py:97-99`, `cerberus/pipelines/gdpr.py:75-79`, `gdpr.py:126-128`
- **What:** In every pipeline the `final.unlink(missing_ok=True)` sits **outside** the `if not cfg.dry_run:` guard that protects the subsequent `.replace()`. For profiling the final write is `concat_gz`, which in dry-run returns at `concat.py:31-35` before creating anything. I demonstrated this end-to-end: wrote a file, called the exact `unlink` + `concat_gz(dry_run=True)` sequence from `profiling.py:139/145`, and the file was gone.
- **Trigger:** `cerberus … --dry-run` pointed at an `-o` directory that already holds results from a previous real run — the natural way to check "what would this invocation do differently?".
- **Consequence:** Silent destruction of finished deliverables (`.meta.R1.fastq.gz`, `.profiling.fastq.gz`, `*_GDPR.fastq.gz`) by a flag that is documented as "Print commands without executing" (`cli.py:112-113`). If `_work/` was cleaned or the machine rebooted, the data is unrecoverable without a full re-run.
- **Fix:** Move each `unlink(missing_ok=True)` inside the `if not cfg.dry_run:` block (or guard it with `if not cfg.dry_run`), and in `profiling.py` reorder so the unlink happens immediately before the real `concat_gz` write.

### F7. All long-read modes align against the short-read-preset minimap2 index

- **Severity:** high
- **Location:** `cerberus/pipelines/long_read.py:29`, `long_read.py:61`, `cerberus/refs.py:42-44`, `cerberus/data/default_manifest.json`, `README.md:174-175`, `scripts/build_custom_host_ref.sh:123-128`
- **What:** `long-meta` and `long-profiling` both resolve their index via `refs.asset("masked_t2t_hla_minimap2")`, and `_PIPELINE_TO_ASSETS` maps both long keys to that same asset. Per README:82-88 and the manifest description, `masked_t2t_hla.mmi` was built with `minimap2 -x sr -d`, which bakes `k=21, w=11` into the index. The long pipelines then invoke `minimap2 -ax map-ont` (`autotune.py:65/71` sets `minimap2_preset="map-ont"` for `LONG`/`VERY_LONG`). minimap2 does not rebuild the index from `-x`; it warns that the prebuilt index's `-k/-w` override the command line and proceeds with `k=21, w=11`. Meanwhile README:174-175 documents a separate `masked_t2t_hla.long.mmi` "for `--long` modes" and `build_custom_host_ref.sh:123-124` actually builds it — but no code path ever requests it, no manifest asset declares it, and the script only symlinks `.mmi` → `.long.mmi` in the `--platform long` case (`build_custom_host_ref.sh:126-127`), so a user who runs the default `--platform both` gets a `.long.mmi` the pipeline can never open.
- **Trigger:** Every `--long` invocation against the published reference bundle (matrix cells 10–18).
- **Consequence:** Exact-21-mer seeding on ONT reads with a 5–10% error rate is far less sensitive than the `map-ont` default of `k=15, w=10` (P(exact 21-mer) ≈ 0.11 vs ≈ 0.21 at 10% error). Host reads that `map-ont` would catch go unmapped and land in the "clean" output. For `--long --gdpr` this directly undercuts the "zero detectable human reads" claim, because the Kraken2 pass alone is explicitly acknowledged as insufficient (README:74). The degradation is silent — minimap2 emits a warning to a per-stage log file nobody reads.
- **Fix:** Add a `masked_t2t_hla_long_minimap2` asset (`masked_t2t_hla.long.mmi`) to the manifest, map `long-meta`/`long-profiling`/`long-gdpr` to it in `_PIPELINE_TO_ASSETS`, and have `long_read.py` and `gdpr.py`'s long branch select the asset key by `cfg.long_mode`.

### F8. fastp's unpaired files always exist, so every run pays for a full orphan alignment pass over zero reads

- **Severity:** medium
- **Location:** `cerberus/stages/qc.py:37-38`, `qc.py:51-52`, `qc.py:70-71`; consumers `cerberus/pipelines/meta.py:51-64`, `cerberus/pipelines/profiling.py:110-135`
- **What:** `run_fastp` always passes `--unpaired1` and `--unpaired2` with distinct paths, and fastp opens both writers at initialisation — the files are created (as valid ~20-30 byte empty gzip members) before fastp knows whether any orphan will be written. `qc.py:70-71` then decides `orphans_r1`/`orphans_r2` purely on `.exists()`, never on size or record count. Both are therefore essentially always non-`None`, so `orphan_inputs` at `meta.py:51` and `profiling.py:110` is always non-empty and the orphan branch always executes.
- **Trigger:** Every short-read run, including on perfectly paired data with zero orphans.
- **Consequence:** Wasted work, not a crash. For meta: a `concat_gz`, a full `minimap2_singles` call (which loads the 7.7 GB `.mmi` from scratch — minutes of I/O and peak RAM), a `samtools fastq`, a `samtools flagstat`, and a `bbduk` entropy pass, all over zero reads. For profiling the same plus a bbduk k-mer call (skipped by the `st_size <= 40` guard at `kmer.py:94`, but `minimap2_singles` and `entropy_single` are unguarded). Worse, the result is a zero-read `<sample>.meta.orphans.fastq.gz` presented as a deliverable, and because `PipelineResult.singletons` is non-`None`, `--gdpr` runs a **whole extra Kraken2 + minimap2 pass** on it (`gdpr.py:83-90`) and emits a zero-read `<sample>.meta.GDPR.fastq.gz`. If fastp is interrupted and leaves a truly zero-byte `.gz`, the file still passes `.exists()`; I confirmed `pigz -dc` and `zcat` both exit 1 on a zero-byte `.gz` ("empty" / "unexpected end of file"), which would surface as a bbduk/JVM error rather than a clean skip.
- **Fix:** Gate on record count, not existence: `orphans_r1=orphan_r1 if orphan_r1.exists() and count_reads(orphan_r1) > 0 else None` (or a cheap `st_size > 40` check consistent with `kmer.py:52`), and add the same guard to `minimap2_singles`/`entropy_single`.

### F9. GDPR output naming is inconsistent with itself, with the pipelines it consumes, and with the README

- **Severity:** medium
- **Location:** `cerberus/pipelines/gdpr.py:73-74`, `gdpr.py:89`, `gdpr.py:98`; compare `cerberus/pipelines/meta.py:74-75/90`, `profiling.py:138`, `long_read.py:40/96`, `README.md:245-247`
- **What:** All three GDPR name templates interpolate `pipeline_result.mode` verbatim. Generated names (verified by instantiating the f-strings):

  | source mode | GDPR output | documented? |
  |---|---|---|
  | `meta` (pair) | `<s>.meta.R1_GDPR.fastq.gz` / `.R2_GDPR.` | yes (README:245-246) |
  | `meta` (orphans) | `<s>.meta.GDPR.fastq.gz` | **no** |
  | `profiling` | `<s>.profiling.GDPR.fastq.gz` | yes (README:247) |
  | `long-meta` | `<s>.long-meta.long_GDPR.fastq.gz` | **no** |
  | `long-profiling` | `<s>.long-profiling.long_GDPR.fastq.gz` | **no** |

  Three separate inconsistencies: (a) the `.{mode}.GDPR.fastq.gz` template means "the scrubbed **orphans**" for meta but "the scrubbed **whole deliverable**" for profiling — the same suffix, two different semantics; (b) the long modes emit a hyphen (`long-meta`, from `PipelineResult.mode` at `long_read.py:24`) while the pipeline's own output for the same mode uses an underscore (`<s>.long_meta.fastq.gz`, `long_read.py:40`), so a single run produces both `SRR123.long_meta.fastq.gz` and `SRR123.long-meta.long_GDPR.fastq.gz`; (c) "long" appears twice in the long-mode names.
- **Trigger:** `--meta --gdpr` (undocumented orphan file), any `--long … --gdpr`.
- **Consequence:** Downstream globs and Snakemake/Nextflow wrappers written against README:239-255 miss three of the five GDPR outputs entirely. The hyphen/underscore split makes a single `{sample}.long_meta.*` glob fail to pick up the GDPR file it is supposed to supersede. Publication workflows keyed on "the `_GDPR` file" will pick up the near-empty meta orphan file alongside the real one.
- **Fix:** Introduce an explicit `output_slug` on `PipelineResult` (`long_meta`, `long_profiling`) separate from the internal `mode`, use `utils/fastq.py:63 output_name()` (already written for this and unused) for every final path, rename the meta orphan output to `<s>.meta.orphans_GDPR.fastq.gz`, drop the redundant `long_` infix, and document all five names in the README output tree.

### F10. The merged profiling FASTQ contains every read name twice

- **Severity:** medium
- **Location:** `cerberus/stages/align.py:76` and `align.py:178` (`-n`), consumed at `cerberus/pipelines/profiling.py:140-145`
- **What:** Both paired `samtools fastq` calls pass `-n`, which suppresses the `/1` and `/2` suffixes that samtools would otherwise append to distinguish mates. That is correct for `--meta`, whose R1 and R2 stay in separate files. But `profiling.py:140-145` then concatenates `ent.r1`, `ent.r2` and the orphans into a single `<sample>.profiling.fastq.gz`, so every fragment appears as two records carrying an identical read identifier.
- **Trigger:** Every `--profiling` run (all short-read profiling cells: 4, 5, 6, 7, 8, 9).
- **Consequence:** Kraken2 itself tolerates it, but the file is documented as a general-purpose profiling input (README:10). Anything that assumes unique identifiers — `seqkit rmdup`, BBTools dedupe, SAM/BAM round-trips, most assembly QC, and Bracken workflows that re-map — will either error or silently collapse the pair. The duplication also propagates into `<sample>.profiling.GDPR.fastq.gz` via `gdpr.py:83-90`.
- **Fix:** Drop `-n` for the profiling lane (or post-process with `seqkit replace` to append `/1`,`/2`) before the merge at `profiling.py:145`, while keeping `-n` for `--meta` where matched names are required by assemblers.

### F11. `PipelineResult.as_inputs_for_gdpr()` is dead code

- **Severity:** low
- **Location:** `cerberus/pipelines/base.py:22-32`
- **What:** A repo-wide grep for `as_inputs_for_gdpr` returns exactly one hit — the definition itself. `run_gdpr_for` (`gdpr.py:38-101`) reaches into `pipeline_result.paired_r1` / `.paired_r2` / `.singletons` / `.long_reads` directly rather than through the accessor.
- **Trigger:** N/A — never executed.
- **Consequence:** The one place that documents the contract between the pipelines and the GDPR post-processor is not the place that implements it, so the two can (and will) drift. The `dict[str, Path]` shape it returns (`"r1"`, `"r2"`, `"singletons"`, `"long"`) does not match anything `gdpr.py` consumes, so it is already stale.
- **Fix:** Either delete it, or make `run_gdpr_for` consume it so the mapping has a single definition and can be unit-tested independently of the external tools.

### F12. `--fast` and `--double-pass` are silently ignored for `--meta`, `--long --meta`, and `--long --double-pass`

- **Severity:** low
- **Location:** `cerberus/pipelines/meta.py:26-101` (no `cfg.fast` / `cfg.double_pass`), `cerberus/pipelines/long_read.py:18-47` (no `cfg.double_pass`), `cerberus/orchestrator.py:53-54`
- **What:** `cfg.fast` is read only at `profiling.py:49`, `:89`, `:120`; `cfg.double_pass` only at `profiling.py:61` and `long_read.py:64`. Neither `run_meta` nor `run_long_meta` consults either. The only guard is `orchestrator.py:53-54`, which warns solely when `cfg.long_mode and cfg.meta and cfg.profiling and cfg.fast` — a conjunction that misses `--meta --fast` (short), `--meta --double-pass` (short), `--long --meta --fast`, `--long --meta --double-pass`, and, most importantly, the genuinely broken `--long --profiling --fast` (F1).
- **Trigger:** Matrix cells 2, 3, 8, 9, 11, 12.
- **Consequence:** Users passing `--fast` expecting a faster meta run, or `--double-pass` expecting a more thorough one, get the identical default pipeline with no indication that the flag did nothing. The CLI help text for both flags does say "Profiling:" (`cli.py:95-98`), but nothing enforces it.
- **Fix:** In `validate_config`, warn (or error) whenever `--fast`/`--double-pass` is set and `--profiling` is not, and warn for `--long --meta --double-pass`; alternatively reject `--fast` outright in long mode until F1 is fixed.

### F13. Manifest, path and accounting inconsistencies in the composition layer

- **Severity:** low
- **Location:** `cerberus/data/default_manifest.json` (`human_kmer_set`), `cerberus/pipelines/base.py:36` + `meta.py:34/44`, `cerberus/orchestrator.py:99-103`, `orchestrator.py:129-131` + `meta.py:110` + `profiling.py:151`
- **What:** Three unrelated small defects in the same layer. (a) `human_kmer_set` (932 MB, `required_for: ["gdpr"]` in the manifest) appears in no `_PIPELINE_TO_ASSETS` entry and is referenced nowhere in the code — computed as `set(manifest.assets) - set(all _PIPELINE_TO_ASSETS values)` → `{'human_kmer_set'}`. It is never used by a pipeline but `fetch_all()` (`refs.py:225-227`) downloads it, so `cerberus fetch-refs` pulls ~932 MB nobody reads. (b) `stage_dir(work_dir, mode, stage)` does `work_dir / mode / stage` while callers already pass `work = cfg.work_dir / mode` (`meta.py:34/44`, `profiling.py:44`, `long_read.py:26/58`, `gdpr.py:46`), so every stage directory carries a duplicated component: `out/_work/meta/meta/01_minimap2_paired`. (c) `run()` returns `outputs` built only from `result_set` primary outputs (`orchestrator.py:101`), so the GDPR files are never printed in the CLI's final summary (`cli.py:330-332`) despite being the point of `--gdpr`. (d) `qc.r1` is fully decompressed and line-counted three times in a `--meta --profiling` run (`orchestrator.py:131`, `meta.py:110`, `profiling.py:151`) plus once per input file at `orchestrator.py:129-130` — for a 50 GB input that is several extra full passes.
- **Trigger:** (a) `cerberus fetch-refs`; (b) every run; (c) any `--gdpr` run; (d) any `--meta --profiling` run.
- **Consequence:** Wasted bandwidth and disk, confusing `_work/` layout that makes log/intermediate paths harder to reason about, GDPR deliverables invisible in the run summary, and tens of minutes of redundant decompression on large inputs.
- **Fix:** Remove `human_kmer_set` from the manifest or wire it into the GDPR asset list and a third bbduk mechanism; change `stage_dir` to `work_dir / stage` (callers already scope by mode); merge GDPR results into the returned `outputs` dict; and cache `count_reads` results on `FastpOutputs`/`PipelineResult` instead of recomputing.
