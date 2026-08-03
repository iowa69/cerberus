# Pass 5 — Autotuning

## Summary

The autotuning subsystem is well-factored, pure, and unit-tested, but roughly a third of what it computes is architecturally unreachable and several of its surviving decisions point the wrong way scientifically. The headline defect is structural: `fastp`/`fastplong` is the *only* consumer of `TunedParams.min_length`/`min_quality`, yet both call sites in `orchestrator.py` pass `tuned=None` (proven by AST over the whole tree — exactly two call sites, both `None`), so the autotuned QC thresholds are dead in every code path and fastp always hard-filters at the constants 50/20 (short) or 200/10 (long). The chicken-and-egg is never resolved by a second pass, so a 2x50 library is trimmed at `--length_required 50` before autotune ever gets to say `min_length=35`. A missing or zero read-length field in the report silently produces `VERY_SHORT`, which selects `-ax sr` and **disables the auxiliary k-mer decontamination stage entirely** — a fail-open default in a decontamination tool. Beyond that: `--dry-run` cannot complete because autotune unconditionally opens a `fastp.json` that dry-run never writes; ONT Q20+/duplex data is misclassified as PacBio HiFi; `--long` can autotune itself back to `-ax sr`; and `--no-1mm-upfront` on very-short reads is, per the Bowtie2 manual, strictly slower with no sensitivity gain.

## Verified working

- **`classify_length()` bucket boundaries** — exhaustively exercised over `[-1, 0, 79, 79.9, 80, 199, 200, 499, 500, 4999, 5000, 1e9, 1e12]`; every bucket is correctly half-open `[lower, upper)` and matches the documentation in `config.py:29-33`. `mean_length >= 1e9` falls out of the loop into the `VERY_LONG` return at `cerberus/autotune.py:39`, which is the right terminal behaviour.
- **`tuned.minimap2_preset` genuinely reaches minimap2** — read at `cerberus/stages/align.py:55` and `cerberus/stages/align.py:113`; those two functions are the entry point for all four alignment call paths (`pipelines/meta.py:40,57`, `pipelines/profiling.py:51,63,114`, `pipelines/long_read.py:30,73`, `pipelines/gdpr.py:65,121`).
- **Extra minimap2 tokens are placed *after* `-ax <preset>`** — `cerberus/stages/align.py:53-62` and `:111-119`; minimap2 honours last-wins, so `-k15 -w10` correctly overrides the `sr` preset's `-k21 -w11` rather than being silently ignored by ordering.
- **`-k15 -w10` for VERY_SHORT is directionally defensible** — smaller k and denser minimisers give more seed opportunities on a 35–50 bp read than the `sr` preset's `-k21 -w11`. The choice of *values* is arguable; the direction is right.
- **`tuned.bowtie2_preset` reaches bowtie2** — `cerberus/stages/align.py:155`, invoked from `cerberus/pipelines/profiling.py:79`.
- **`tuned.entropy` reaches bbduk with correct null-precedence** — `cerberus/stages/entropy.py:37` and `:67` use `cfg.entropy if cfg.entropy is not None else tuned.entropy`, so a legitimate `--entropy 0.0` is honoured rather than falling through a truthiness test.
- **`tuned.bbduk_k` reaches bbduk** — `cerberus/stages/kmer.py:68` and `:114`, whenever the aux stage actually runs.
- **`tuned.bbduk_aux_enabled` gates the aux stage consistently** — all three gates present and identical in intent: `cerberus/pipelines/profiling.py:89`, `:120`, `cerberus/pipelines/long_read.py:81`.
- **`tuned.winnowmap_enabled` is read and gated** — `cerberus/pipelines/long_read.py:64` correctly requires both `cfg.double_pass` and the tuned flag.
- **Explicit `--platform` overrides both platform and preset** — `cerberus/autotune.py:156-158`; verified by running `apply_user_overrides` directly: with `platform=ONT` the only changed fields are `{min_length, min_quality, platform, minimap2_preset, entropy, bbduk_k}` and `minimap2_preset` becomes `map-ont`. Matches `tests/test_autotune.py:126-127`.
- **`detect_platform_from_fastp` short-circuits on an explicit user choice** — `cerberus/autotune.py:92-93`, so autodetection can never overrule the user.
- **The `read1 → read → read2` key fallback chain is a sensible defence against fastp-vs-fastplong key naming** — `cerberus/autotune.py:120-125` and `:97`. (Its `or` semantics are a separate problem, see F3.)
- **Autotune is pure and side-effect-free**, making it genuinely auditable; `tests/test_autotune.py` covers all five length buckets and three of the platform paths.

## Findings

### F1. Autotuned `min_length` / `min_quality` are dead — fastp runs before autotune and is always passed `tuned=None`

- **Severity:** critical
- **Location:** `cerberus/orchestrator.py:123` (and `:125`), `cerberus/orchestrator.py:173` (and `:174`); consumer `cerberus/stages/qc.py:42-43`, `cerberus/stages/qc.py:98-99`
- **What:** `TunedParams.min_length` and `TunedParams.min_quality` (`cerberus/config.py:39-40`) have exactly two consumers in the entire tree, both inside `run_fastp`/`run_fastplong`:
  ```python
  min_len  = (tuned.min_length  if tuned else None) or cfg.min_length  or 50   # qc.py:42
  min_qual = (tuned.min_quality if tuned else None) or cfg.min_quality or 20   # qc.py:43
  ```
  An AST walk over every `.py` in the repo finds exactly two call sites of these functions, and **both hard-code `tuned=None`**:
  ```
  cerberus/orchestrator.py:123  run_fastp(...    tuned=None)
  cerberus/orchestrator.py:173  run_fastplong(... tuned=None)
  -> 2 call sites; all tuned=None: True
  ```
  Autotune then runs at `cerberus/orchestrator.py:137` / `:181`, i.e. strictly *after* the only thing that could have used its answer. The `tuned` parameter of both QC functions is vestigial: it is keyword-only, defaults to `None`, and is never supplied.
- **Trigger:** Every run, without exception. There is no second QC pass — grep and AST both confirm one `run_fastp` and one `run_fastplong` call site.
- **Consequence:** Ten of the values in `_BASE_PARAMS` are computed, logged as "Tuned params" (`cerberus/autotune.py:141`), written into `cfg.tuned`, and never used:
  ```
  very_short  min_length=35   min_quality=20
  short       min_length=50   min_quality=20
  medium      min_length=75   min_quality=20
  long        min_length=200  min_quality=10
  very_long   min_length=500  min_quality=10
  ```
  fastp always sees the literal `50`/`20`; fastplong always sees `200`/`10`. **For a 2x50 bp library this is a large, silent yield loss**: the reads are exactly 50 bp before trimming, and `run_fastp` unconditionally enables `--trim_poly_g`, `--trim_poly_x` and `--detect_adapter_for_pe` (`cerberus/stages/qc.py:55-57`), so *any* read that loses even one base to poly-G (NovaSeq/NextSeq two-colour chemistry — precisely the platform that produces 2x50) or adapter read-through drops below 50 and is discarded, with its mate diverted to the unpaired stream. Autotune's `min_length=35` was the correct answer and can never be applied. The log will nonetheless print `min_len=35`, so the operator believes a 35 bp floor was used. This is worse than having no autotuning: it is autotuning that *reports* a decision it did not make.
- **Fix:** Decide the read-length class *before* QC. `estimate_long_read` (`cerberus/autotune.py:166`) already implements exactly the needed sampling primitive; generalise it to `sample_mean_length(path, n=10_000) -> float`, call it on the raw `cfg.r1` / `cfg.long_input` in `_run_short`/`_run_long`, feed the result through `classify_length()` + `apply_user_overrides()`, and pass the resulting `TunedParams` into `run_fastp`/`run_fastplong`. Then re-run `autotune_from_fastp` on the resulting JSON to refine the *post*-QC decisions (preset, platform, bbduk_k, aux, entropy) which legitimately want post-trim statistics. If the two-stage design is rejected, delete `min_length`/`min_quality` from `TunedParams` and the `tuned=` parameter from both QC functions so the log stops claiming a decision that is not made.

### F2. `cerberus --dry-run` always aborts at the autotune step

- **Severity:** critical
- **Location:** `cerberus/autotune.py:115`; triggered from `cerberus/orchestrator.py:137` and `cerberus/orchestrator.py:181`
- **What:** `autotune_from_fastp` opens the report unconditionally (`with fastp_json.open() as f:`). Under `--dry-run`, `run()` short-circuits at `cerberus/utils/shell.py:66-69` and never invokes fastp, so `<workdir>/00_qc/fastp.json` is never created.
- **Trigger:** `cerberus -r1 ... -r2 ... --meta --dry-run` (or any `--dry-run` invocation that is not the chopper fallback). Reproduced end-to-end against the real `run_fastp` with `dry_run=True`:
  ```
  dry-run wrote fastp.json? False
  autotune -> FileNotFoundError: .../out/_work/00_qc/fastp.json
  ```
- **Consequence:** The advertised "print commands without executing" mode (`cerberus/cli.py:112-113`) crashes with an uncaught `FileNotFoundError` before a single pipeline command is printed — the exception is not one of the types caught in `cerberus/cli.py:316-324`, so the user gets a raw traceback and exit code 1. Every stage's `dry_run=` plumbing downstream is therefore untestable. Note the sole exception: the chopper fallback writes its synthetic JSON outside the dry-run guard (`cerberus/stages/qc.py:129`), so `--long --dry-run` *only* works when `fastplong` is absent and `chopper` is present.
- **Fix:** Guard the call: `if cfg.dry_run: tuned = TunedParams()` (or the pre-QC sampled estimate from F1) with a log line explaining the substitution, and/or make `autotune_from_fastp` return a documented default on `FileNotFoundError`. Add a dry-run smoke test to `tests/`.

### F3. A missing, null, or zero read-length silently selects VERY_SHORT and disables k-mer decontamination

- **Severity:** high
- **Location:** `cerberus/autotune.py:120-125` (and the parallel chain at `:97`), consumed at `cerberus/autotune.py:127` and `:101`
- **What:** The `or` chain terminates in a literal `0`, and `classify_length(0.0)` returns `VERY_SHORT` because the first bucket test is `0 < 80`. Four distinct malformed-report shapes all collapse to the same silent answer:
  ```
  field absent        -> class=very_short  platform=illumina  mm2=sr  aux=False  bbduk_k=23  entropy=0.6
  explicit 0          -> class=very_short  platform=illumina  mm2=sr  aux=False  bbduk_k=23  entropy=0.6
  explicit null       -> class=very_short  platform=illumina  mm2=sr  aux=False  bbduk_k=23  entropy=0.6
  whole summary gone  -> class=very_short  platform=illumina  mm2=sr  aux=False  bbduk_k=23  entropy=0.6
  ```
  Nothing is logged as anomalous; `cerberus/autotune.py:137-140` cheerfully prints `mean_len=0bp class=very_short`.
- **Trigger:** Any report where `summary.before_filtering.read1_mean_length` is absent, `null`, or `0` — a fastp/fastplong version whose key naming differs, a truncated JSON write after an OOM kill, an input FASTQ where every read was filtered, or a future third-party QC tool substituted for fastp.
- **Consequence:** This is a **fail-open default for a decontamination tool**. `VERY_SHORT` is the *only* class with `bbduk_aux_enabled=False` (`cerberus/autotune.py:45`), so the entire auxiliary k-mer host-removal stage is skipped at `cerberus/pipelines/profiling.py:89`, `:120` and `cerberus/pipelines/long_read.py:81`. A parse hiccup therefore silently removes a decontamination mechanism from a pipeline whose whole purpose is host removal, and whose GDPR mode is marketed as "publication-defensible" (`cerberus/pipelines/gdpr.py:12`). It simultaneously forces `-ax sr`, which on a long-read run is catastrophic. Note also that `refs.py:39` has already downloaded `aux_refs` for the profiling key before autotune runs, so the asset is fetched and then discarded.
- **Fix:** Distinguish "absent" from "zero". Use an explicit sentinel walk (`for key in (...): v = before.get(key); if v is not None: break`) and raise/warn loudly when no length is recoverable. Fail *safe*, not fast: on an unusable report, fall back to the `SHORT` defaults (aux enabled) and emit `log.error`, never to `VERY_SHORT`. Better still, decouple `bbduk_aux_enabled` from the length class entirely and make it an explicit `--no-aux-refs` opt-out, so no parsing accident can turn off a decontamination stage.

### F4. ONT Q20+/duplex data is misclassified as PacBio HiFi by the `q20_rate >= 0.99` clause

- **Severity:** high
- **Location:** `cerberus/autotune.py:104`
- **What:** `if q30_rate >= 0.85 or q20_rate >= 0.99: return Platform.PACBIO_HIFI`. The `or` means the `q20_rate` clause alone decides, regardless of how poor `q30_rate` is.
- **Trigger:** Modern ONT chemistry. Verified:
  ```
  ONT R10.4.1 duplex-ish (q20=0.99, q30=0.60)   -> pacbio-hifi
  high q20, poor q30      (q20=0.995, q30=0.50) -> pacbio-hifi
  typical ONT simplex     (q20=0.90, q30=0.84)  -> ont
  ```
- **Consequence:** `_platform_preset` returns `map-hifi` (`cerberus/autotune.py:79`), which reaches minimap2 at `cerberus/stages/align.py:55`/`:113`. `map-hifi` is parameterised for ~0.1–1 % error HiFi reads; applying it to ONT reads at 1–5 % error reduces host-alignment sensitivity, so human reads that `map-ont` would have caught survive into the "decontaminated" output — including the GDPR output. This is a silent contamination leak, not a crash. The threshold was clearly written against R9.4.1-era ONT (q20_rate ≈ 0.3–0.6) and has not been revisited for R10.4.1/duplex, where q20_rate ≥ 0.99 is routine.
- **Fix:** Replace the `or` with an `and` (both `q30_rate >= 0.85` *and* `q20_rate >= 0.99`) so a high q20 alone cannot outvote a poor q30 — HiFi's distinguishing property is precisely its high *q30*. Additionally, when the two clauses disagree, prefer `ONT`: `map-ont` degrades gracefully on HiFi input, whereas `map-hifi` does not degrade gracefully on ONT input. Asymmetric costs should drive an asymmetric default. Add unit tests for the duplex regime — `tests/test_autotune.py:68-70` only covers the R9-era simplex case.

### F5. `--long` can autotune itself back to Illumina and `-ax sr`

- **Severity:** high
- **Location:** `cerberus/autotune.py:101-102`; reached from `cerberus/orchestrator.py:181`
- **What:** `if mean_len < 500: return Platform.ILLUMINA` is unconditional. Neither `autotune_from_fastp` nor `detect_platform_from_fastp` receives or consults `cfg.long_mode`; only `user_platform` is threaded through.
- **Trigger:** A long-read run whose *post-QC* mean length lands under 500 bp. Verified:
  ```
  ONT run, short fragments (cfDNA/aDNA), mean 400bp -> class=medium platform=illumina mm2=-ax sr
  just under the 500bp cutoff (499)                 -> class=medium platform=illumina mm2=-ax sr
  just over (501)                                   -> class=long   platform=ont      mm2=-ax map-ont
  ```
  Realistic inputs: ONT/PacBio sequencing of cell-free DNA, ancient or formalin-degraded DNA, short-amplicon runs, or any library where fastplong's `--length_required` (which itself is the hard-coded `200`, see F1) leaves a short-fragment population.
- **Consequence:** `minimap2 -ax sr` is run on error-prone long reads. The `sr` preset assumes short, near-exact, paired fragments; on 400 bp ONT reads at 5 % error it will fail to align a large fraction of genuine host reads, so human sequence passes straight through the host-removal step into the meta/profiling/GDPR outputs. The user explicitly told the tool these are long reads via `--long`, and the tool overrode them. A one-base difference in mean length (499 vs 501) flips the aligner preset.
- **Fix:** Pass `cfg.long_mode` into `autotune_from_fastp`/`detect_platform_from_fastp` and never return `ILLUMINA` when `long_mode` is set — clamp to `ONT` as the safe long-read default. Symmetrically, warn (loudly) when short-read mode autodetects a long-read length class. Also log a warning whenever the detected platform contradicts the requested mode.

### F6. `bowtie2_extra="--no-1mm-upfront"` for VERY_SHORT reads is backwards

- **Severity:** high
- **Location:** `cerberus/autotune.py:47`; consumed at `cerberus/stages/align.py:156-157`
- **What:** VERY_SHORT is given `bowtie2_extra="--no-1mm-upfront"` on top of `--very-sensitive-local`. Per the Bowtie2 manual, by default Bowtie2 "will attempt to find either an exact or a 1-mismatch end-to-end alignment for the read *before* trying the multiseed heuristic. Such alignments can be found very quickly, and many short read alignments have exact or near-exact end-to-end alignments... This option prevents Bowtie 2 from searching for 1-mismatch end-to-end alignments before using the multiseed heuristic... **This comes at the expense of speed.**" The flag exists to make behaviour predictable when a user hand-tunes `-L`/`-N`; it is not a sensitivity knob, it is active in `--local` mode, and no preset sets or clears it (`--very-sensitive-local` expands to `-D 20 -R 3 -N 0 -L 20 -i S,1,0.50`, which is orthogonal to it).
- **Trigger:** Any 2x50 / 2x75 short-read `--profiling` run (`cerberus/pipelines/profiling.py:79`), i.e. exactly the class the flag was added for.
- **Consequence:** Strictly worse on both axes. It is *slower* (the manual says so explicitly), and it removes the fast path that finds exact/1-mismatch **end-to-end** alignments — which is the dominant alignment mode for a 35–50 bp read against the human reference. Host reads that would have been caught by that pre-pass and are not recoverable by seeded local alignment at `-L 20` now survive into the profiling output. For a host-*depletion* tool, deliberately disabling an alignment-finding pathway is the opposite of the intent, and the parameter is filed under a class whose stated goal (`config.py:29`) is to handle the hardest-to-align reads.
- **Fix:** Remove `bowtie2_extra="--no-1mm-upfront"` from the VERY_SHORT entry. If the goal was genuinely to increase short-read sensitivity, tune the seed instead — e.g. `bowtie2_extra="-L 18 -N 1"` or `"-L 15"`, since `--very-sensitive-local`'s `-L 20` seed occupies 40–57 % of a 35–50 bp read and leaves little room for a second seed. Any such change needs a benchmark on a spiked host/microbe mixture, not a table edit.

### F7. Autotune enables winnowmap for VERY_LONG, but the meryl DB it requires is never provisioned

- **Severity:** medium
- **Location:** `cerberus/autotune.py:71` (`winnowmap_enabled=True`); `cerberus/pipelines/long_read.py:64-71`; `cerberus/stages/align.py:210`; asset registry `cerberus/refs.py:37-45`
- **What:** `long_read.py:66` computes `meryl_db = idx.with_suffix(".meryl")` and passes it to `winnowmap -W <meryl_db>`. `_PIPELINE_TO_ASSETS` contains no meryl or winnowmap asset for any key — `long-profiling` maps only to `["masked_t2t_hla_minimap2", "aux_refs"]`. A repo-wide grep for `meryl` finds hits only in `align.py`, `long_read.py` and this line; no builder script, no manifest entry.
- **Trigger:** `cerberus --long -i reads.fq.gz --profiling --double-pass` on any library whose mean length ≥ 5 kb (i.e. the normal ONT/HiFi case).
- **Consequence:** The autotuned `winnowmap_enabled=True` routes the run into a stage that references a path `RefManager` never creates. `require_tools("winnowmap", ...)` passes if the binary is installed (it is in `environment.yml:22`), so the failure is a winnowmap runtime error on a missing `-W` argument, surfacing as a `ToolError` mid-run after references have been downloaded and QC has completed. This is the *only* consumer of `winnowmap_enabled`, so the field currently has no working code path at all.
- **Fix:** Either add the meryl repetitive-k-mer database as a proper asset in `refs.py` (with a `required_for` entry keyed to `long-profiling`) and gate `winnowmap_enabled` on its presence, or set `winnowmap_enabled=False` in `_BASE_PARAMS[VERY_LONG]` until the asset ships. At minimum, add an existence check in `long_read.py:66` that falls back to `minimap2_singles` with a warning rather than failing the run.

### F8. `--minimap2-args` / `--bowtie2-args` replace rather than merge the autotuned extras

- **Severity:** medium
- **Location:** `cerberus/stages/align.py:50`, `cerberus/stages/align.py:108`, `cerberus/stages/align.py:156-157`
- **What:** `extra = cfg.minimap2_args or tuned.minimap2_extra`. Merging happens here, in `align.py`, and **not** in `apply_user_overrides` — running `apply_user_overrides` with `minimap2_args="-N 5"` set leaves `tuned.minimap2_extra` untouched at `'-k15 -w10'` (the changed-field set is `{min_length, min_quality, platform, minimap2_preset, entropy, bbduk_k}`). Verified behaviour:
  ```
  cfg.minimap2_args=None   -> minimap2 extra tokens = ['-k15', '-w10']
  cfg.minimap2_args='-N 5' -> minimap2 extra tokens = ['-N', '5']
  cfg.bowtie2_args=None       -> bowtie2 extra = ['--no-1mm-upfront']
  cfg.bowtie2_args='--no-mixed' -> bowtie2 extra = ['--no-mixed']
  ```
- **Trigger:** Any user who passes `--minimap2-args` on a 2x50/2x75 library.
- **Consequence:** The autotuned `-k15 -w10` for VERY_SHORT is silently discarded the instant the user adds any unrelated minimap2 flag. A user adding `-N 5` (secondary alignment count) unknowingly reverts seeding to the `sr` preset's `-k21 -w11`, reducing sensitivity on exactly the class that needed the boost. Nothing is logged. Precedence *is* internally consistent between `align.py:50` and `align.py:156` (both `or`, replace-not-merge), but it is inconsistent with the rest of the subsystem: `apply_user_overrides` (`cerberus/autotune.py:148-157`) and `entropy.py:37`/`:67` all use `is not None`. A secondary consequence of the same `or`: `--minimap2-args ""` is falsy and silently falls back to the autotuned value rather than clearing it.
- **Fix:** Concatenate instead of replacing: `extra_tokens = shlex.split(tuned.minimap2_extra) + shlex.split(cfg.minimap2_args or "")`, with the user tokens last so last-wins still gives the user precedence on any specific flag. Log the final token list at debug level. Use `is not None` rather than truthiness so an explicit empty string can clear the autotuned extras. Apply the same to bowtie2. Document the merge semantics in `cli.py:136-139` ("appended to" currently implies merging, which is not what happens).

### F9. `bbduk_k` rises with read length while the platform error rate also rises, making the aux k-mer pass near-inert on ONT

- **Severity:** medium
- **Location:** `cerberus/autotune.py:63`, `cerberus/autotune.py:69` (`bbduk_k=31` for LONG/VERY_LONG); consumed at `cerberus/stages/kmer.py:68-69` and `:114-115` alongside the hard-coded `mcf=0.5`
- **What:** The table increases k monotonically with read length (23 → 27 → 31 → 31 → 31). But the length class also determines the platform, and long reads mean high error rates. The probability that a given k-mer is error-free is `(1-e)^k`:
  ```
  class       k    typical platform       err     P(kmer clean)
  very_short  23   Illumina               0.001         0.9773
  short       27   Illumina               0.001         0.9733
  medium      31   Illumina 2x250         0.001         0.9695
  long        31   ONT (map-ont chosen)   0.05          0.2039
  very_long   31   ONT (map-ont chosen)   0.05          0.2039
  very_long   31   PacBio HiFi            0.001         0.9695
  ```
- **Trigger:** Every `--long` `--profiling` run on ONT data (`cerberus/pipelines/long_read.py:81-88`).
- **Consequence:** With only ~20 % of 31-mers error-free on 5 %-error ONT reads, satisfying bbduk's `mcf=0.5` (at least 50 % of the read covered by matching reference k-mers) is close to impossible. The auxiliary host k-mer stage therefore runs a full pass over the data — download, I/O, `-Xmx{memory_gb}g` JVM, wall-clock — and removes essentially nothing, while its `stats` file reports a near-zero removal that reads as "the data was clean". The tuning is backwards: error-prone platforms need *smaller* k, not larger. Note this is the only knob `--bbduk-k` controls, so a user trying to fix it manually is fighting an undocumented `mcf` interaction.
- **Fix:** Key `bbduk_k` on the detected *platform*, not the length class: `k≈31` for Illumina/HiFi, `k≈19–21` for ONT/CLR. Tune `mcf` alongside it (a lower `mcf` for high-error platforms), or use `maxbadkmers`/`hdist` instead of a raw coverage fraction. Validate against a spiked host/microbe long-read mixture before committing numbers, and record the removal rate in the accounting report so an inert stage is visible.

### F10. The chopper fallback makes autotune a constant function

- **Severity:** medium
- **Location:** `cerberus/stages/qc.py:140-161` (`_synthetic_long_json`), consumed at `cerberus/stages/qc.py:129` → `cerberus/orchestrator.py:181`
- **What:** When `fastplong` is absent and `chopper` is used, the synthetic report hard-codes `read1_mean_length: 5000`, `q20_rate: 0.95`, `q30_rate: 0.7`. Feeding four different `(min_len, min_qual)` chopper settings through the real code path yields exactly one distinct `TunedParams`:
  ```
  distinct TunedParams over 4 different chopper settings: 1
    length-class=very_long platform=ont min_len=500 Q=10 entropy=0.45 bbduk_k=31 aux=True mm2=map-ont bt2=--very-sensitive-local
  ```
- **Trigger:** Any `--long` run on a machine where `fastplong` is not installed. `fastplong` is listed as *optional* in `cerberus/cli.py:254`, so this is a supported configuration.
- **Consequence:** Autotuning is a fiction on this path — it reports `mean_len=5000bp class=very_long platform=ont` for every input, whether the reads are 500 bp or 100 kb, HiFi or ONT. `5000` sits exactly on the `VERY_LONG` boundary and therefore also sets `winnowmap_enabled=True`, which combined with F7 means `chopper + --double-pass` always attempts winnowmap against a non-existent meryl DB. The `notes` field in the synthetic JSON records the chopper settings but autotune never reads it, and the warning at `cerberus/stages/qc.py:116` says only "JSON report will be minimal", not "autotuning will be disabled".
- **Fix:** Have the chopper branch compute real statistics — a sampling pass over `out_reads` (the `estimate_long_read` primitive generalised per F1) costs a fraction of a second and would populate `read1_mean_length`, `q20_rate` and `q30_rate` honestly. Failing that, make `autotune_from_fastp` detect the synthetic marker and log `log.warning("Autotune disabled: synthetic report from chopper fallback")` so the operator is not misled by the `Tuned params:` line.

### F11. `autotune_from_fastp` has no JSON validation or error handling

- **Severity:** medium
- **Location:** `cerberus/autotune.py:115-125`
- **What:** `json.load` and the subsequent `.get()` chain are entirely unguarded. Verified failure modes:
  ```
  missing file (dry-run)       -> FileNotFoundError
  truncated JSON               -> JSONDecodeError: Expecting value: line 1 column 12
  empty file (0 bytes)         -> JSONDecodeError: Expecting value: line 1 column 1
  not JSON at all              -> JSONDecodeError
  JSON is a top-level list     -> AttributeError: 'list' object has no attribute 'get'
  summary is a string          -> AttributeError: 'str' object has no attribute 'get'
  ```
- **Trigger:** Truncated write after an OOM kill of fastp, a disk-full condition, a `--dry-run` (see F2), or a fastp fork whose schema differs.
- **Consequence:** An unhandled `JSONDecodeError`/`AttributeError` propagates all the way to `cerberus/cli.py:317`, which catches only `ConfigError` and `KeyboardInterrupt` — the user gets a raw traceback with no indication that the problem is a corrupt QC report, after the expensive reference download and QC have already completed. The `AttributeError` variants are particularly opaque.
- **Fix:** Wrap the load in `try/except (OSError, json.JSONDecodeError) as e` and raise a domain exception naming the file and the likely cause. Type-check the result (`if not isinstance(report, dict): ...`) and each nested `.get()` before descending. Combined with F3's fix, resolve to a safe documented default (`SHORT`, aux enabled) plus `log.error` rather than propagating.

### F12. `estimate_long_read` sniffs compression by filename suffix only, and crashes on any mismatch

- **Severity:** medium
- **Location:** `cerberus/autotune.py:174`; called from `cerberus/orchestrator.py:167`
- **What:** `opener = gzip.open if input_path.suffix == ".gz" else open`. `Path.suffix` returns only the final component, so the check is a pure filename convention. Verified:
  ```
  plain .fastq                 -> True
  proper .fastq.gz             -> True
  gzip content, .fastq name    -> UnicodeDecodeError: 'utf-8' codec can't decode byte 0x8b in position 1
  gzip content, .bgz name      -> UnicodeDecodeError
  gzip content, .gzip name     -> UnicodeDecodeError
  plain text, .gz name         -> BadGzipFile: Not a gzipped file (b'@r')
  ```
  bgzip output is gzip-compatible and works *if* named `.gz`; it fails on the common `.bgz` naming. `.bz2`, `.xz` and `.zst` are all silently treated as plain text.
- **Trigger:** `cerberus --long -i reads.fq.bgz ...` or any gzip-compressed file not named `.gz`. This runs at `cerberus/orchestrator.py:167`, *before* QC and before any output is produced, so it aborts the run immediately.
- **Consequence:** An uncaught `UnicodeDecodeError` from a function whose sole purpose is to emit an advisory warning kills the run. That is a hard failure caused by a purely cosmetic check. Two secondary issues with the same function: (a) its warning at `cerberus/orchestrator.py:168-171` is advisory only — nothing changes, so a user who genuinely fed short reads to `--long` proceeds anyway with fastplong's hard-coded 200 bp filter; (b) there is no symmetric check in `_run_short`, so long reads fed to `-r1/-r2` produce no warning at all.
- **Fix:** Sniff the magic bytes (`\x1f\x8b` for gzip, `\x28\xb5\x2f\xfd` for zstd, `BZh` for bzip2) rather than the extension, and wrap the whole function in `try/except Exception` returning `False` with a `log.debug` — it is a heuristic, and a heuristic must never abort a run. Add the mirrored check to `_run_short`. Reuse the same sniffer for the F1 pre-QC sampling and for the `zcat`/`cat` decision at `cerberus/stages/qc.py:120-121`, which has the identical suffix-only bug.

### F13. `--min-length 0` and `--min-quality 0` are silently ignored

- **Severity:** medium
- **Location:** `cerberus/stages/qc.py:42-43`, `cerberus/stages/qc.py:98-99`
- **What:** The `or` chain treats a legitimate `0` as absent. Verified:
  ```
  --min-length None  --min-quality None  -> fastp gets (50, 20)
  --min-length 0     --min-quality 0     -> fastp gets (50, 20)
  --min-length 1     --min-quality 1     -> fastp gets (1, 1)
  ```
- **Trigger:** A user who wants to disable length or quality filtering entirely — a normal request when Cerberus is used purely as a decontamination step in front of a pipeline that does its own QC.
- **Consequence:** `--min-quality 0` (intent: no quality filtering) silently becomes 20; `--min-length 0` silently becomes 50. The CLI accepts the value, `apply_user_overrides` at `cerberus/autotune.py:148-151` correctly uses `is not None` and dutifully records it in `tuned`, and the log prints the user's `0` — but fastp receives the default. The user's explicit instruction is contradicted with a log line that says it was honoured. Note that `--entropy 0.0` *does* work (`entropy.py:37` uses `is not None`), so the behaviour is inconsistent across flags.
- **Fix:** Replace both `or` chains with explicit `is not None` resolution. Optionally validate at the CLI (`cerberus/cli.py:128-131`) that the values are `>= 0` and document that `0` disables the filter.

### F14. `Platform.PACBIO_CLR` is unreachable from autodetection

- **Severity:** medium
- **Location:** `cerberus/autotune.py:80-81`; detection logic `cerberus/autotune.py:101-106`
- **What:** Exhaustive sweep over `mean_len ∈ {0, 100, 499, 500, 1000, 15000, 100000} × q20 ∈ {0, 0.5, 0.9, 0.98, 0.99, 1.0} × q30 ∈ {0, 0.5, 0.84, 0.85, 0.99, 1.0}`:
  ```
  platforms reachable over exhaustive grid: ['illumina', 'ont', 'pacbio-hifi']
  PACBIO_CLR reachable from --platform auto: False
  ```
  `detect_platform_from_fastp` has exactly three return paths and none of them is `PACBIO_CLR`. The `map-pb` branch at `cerberus/autotune.py:80-81` executes only via `apply_user_overrides` when the user passes `--platform pacbio-clr` explicitly.
- **Trigger:** Any PacBio CLR run without an explicit `--platform` flag.
- **Consequence:** CLR (≈85–88 % accuracy, 10–30 kb) falls through to `ONT` — `q30 < 0.85` and `q20 < 0.99` — and gets `-ax map-ont`. The practical damage is modest since `map-ont` and `map-pb` have similar error tolerances, but `README.md:13` states the tool "Autotunes its parameters from the data — you do not need to know the read length, platform, or sensible thresholds", and `README.md:112` documents `--platform pacbio-clr` as an *override*. For one of the four advertised platforms, autodetection is structurally impossible, and there is no way for the user to discover this short of reading the source. The `_platform_preset` CLR branch is untested — `tests/test_autotune.py` has no CLR case.
- **Fix:** Either add a CLR detection rule (CLR is distinguishable from ONT primarily by its very low `q20_rate` at long lengths, though the separation is genuinely weak) or, more honestly, document at `cerberus/autotune.py:91` and in `README.md:13` that CLR requires an explicit `--platform` flag and log a hint when a long low-quality library is classified as ONT. Add a `--platform pacbio-clr` test asserting `map-pb`.

### F15. `min_quality` conflates two different statistics between fastp and fastplong

- **Severity:** low
- **Location:** `cerberus/stages/qc.py:54` vs `cerberus/stages/qc.py:111`; table values at `cerberus/autotune.py:44-70`; help text at `cerberus/cli.py:130-131`
- **What:** For short reads the value is passed as `--qualified_quality_phred`, which is fastp's *per-base* threshold defining a "qualified" base — it filters a read only in combination with fastp's `--unqualified_percent_limit` (default 40 %). For long reads the same value is passed as fastplong's `--mean_qual`, which is a genuine per-read mean quality. `cli.py:131` describes the flag as "Minimum mean quality for fastp", which is wrong for the fastp branch.
- **Trigger:** Reading the tuning table or the CLI help; any attempt to reason about why `min_quality` is 20 for Illumina and 10 for long reads.
- **Consequence:** The `20` vs `10` contrast in `_BASE_PARAMS` is not a like-for-like comparison — a `--qualified_quality_phred 20` with a 40 % unqualified allowance is far more permissive than a `--mean_qual 20` would be, and `--mean_qual 10` on ONT is roughly "keep everything". The table reads as a calibrated sensitivity gradient but is comparing two different quantities. (Given F1, none of these values currently reach either tool, so this is latent rather than active.)
- **Fix:** Split into two fields (`fastp_qualified_phred`, `fastplong_mean_qual`) or at minimum document the distinction inline in `_BASE_PARAMS` and correct the `--min-quality` help string. If `--unqualified_percent_limit` matters for the short-read decision, tune it explicitly rather than relying on fastp's default.

### F16. `bbduk_k` for VERY_SHORT — and the user's `--bbduk-k` — are unreachable

- **Severity:** low
- **Location:** `cerberus/autotune.py:44-45`; gates at `cerberus/pipelines/profiling.py:89`, `:120`, `cerberus/pipelines/long_read.py:81`
- **What:** `bbduk_k` is read only at `cerberus/stages/kmer.py:68` and `:114`, inside `bbduk_kmer_paired`/`bbduk_kmer_single`, which are called only when `tuned.bbduk_aux_enabled` is true. `VERY_SHORT` sets `bbduk_aux_enabled=False` *and* `bbduk_k=23`, so the 23 can never reach bbduk.
- **Trigger:** Any 2x50/2x75 library, and any run that falls into `VERY_SHORT` via F3.
- **Consequence:** Cosmetic on its own — but a user on a 2x50 library who passes `--bbduk-k 19` to try to improve host removal gets no error, no warning, and no effect: `apply_user_overrides` records the value, the log line prints `bbduk_k=19`, and the stage that would consume it never runs. The same is true for `--aux-refs`, whose `cfg.aux_refs_override` is stored at `cerberus/cli.py:211` but grep shows it is never read anywhere in `cerberus/`.
- **Fix:** Warn in `apply_user_overrides` when `cfg.bbduk_k is not None` but the resolved class has `bbduk_aux_enabled=False`. Remove the misleading `bbduk_k=23` from the VERY_SHORT entry or give VERY_SHORT a working aux configuration (see F3's fix). Separately, wire up or remove `aux_refs_override`.

### F17. The entropy table is non-monotonic and `entropywindow`/`entropyk` are not tuned alongside it

- **Severity:** low
- **Location:** `cerberus/autotune.py:44-70`; consumers `cerberus/stages/entropy.py:37,43-45` and `:67,73-75`
- **What:** The thresholds are `very_short=0.6 → short=0.7 → medium=0.65 → long=0.5 → very_long=0.45` — neither increasing nor decreasing in read length; `SHORT` is the strictest, and `MEDIUM` is looser than `SHORT` despite having more sequence context to estimate entropy from. There is no comment explaining the peak at `SHORT`. Meanwhile `entropywindow=50` and `entropyk=5` are hard-coded in both entropy functions and are not part of `TunedParams` at all.
- **Trigger:** Every run; most visibly on VERY_SHORT, where the 50 bp sliding window is ≥ the read length, so bbduk effectively computes a single whole-read entropy rather than a windowed maximum — a different statistic in kind from what a 150 bp or 10 kb read gets.
- **Consequence:** The threshold and the window are not comparable across classes, so the apparent "gradient" in the table does not correspond to a consistent stringency. A 0.60 whole-read entropy cut on a 35 bp read is not looser than a 0.70 windowed cut on a 150 bp read — it is a different filter. Related documentation drift: `README.md:86` describes "bbmask with entropy=0.7, window=80, plus k-mer repeat masking (kr=5 minlen=40 mincount=4)", but the code runs `bbduk.sh` with `entropywindow=50 entropyk=5` and no masking at all.
- **Fix:** Either make the table monotonic with a documented rationale, or add `entropy_window`/`entropy_k` to `TunedParams` and scale the window with read length (e.g. `min(50, read_length)`) so the statistic is comparable across classes. Reconcile `README.md:86` with `cerberus/stages/entropy.py:43-45`. Any change to these thresholds needs an empirical false-positive/false-negative measurement, not a table edit.

### F18. Dead autotune plumbing: override branches, `_tuned_or_default`, and `tuned.platform`

- **Severity:** low
- **Location:** `cerberus/autotune.py:148-151`; `cerberus/orchestrator.py:201-206`; `cerberus/config.py:50`
- **What:** Three pieces of vestigial machinery, each verified:
  1. `apply_user_overrides` copies `cfg.min_length`/`cfg.min_quality` into `tuned` (`cerberus/autotune.py:148-151`), but per F1 nothing reads `tuned.min_length`. The user's values reach fastp through the *separate* `cfg.min_length` term at `cerberus/stages/qc.py:42`, so these two branches are pure no-ops. `tests/test_autotune.py:123` asserts `overridden.min_length == 99`, testing a value that has no consumer.
  2. `_tuned_or_default` (`cerberus/orchestrator.py:201-206`) guards with `if cfg.tuned and cfg.tuned.read_length_class:` — a dataclass instance is always truthy and `read_length_class` is always a non-empty enum member. Verified: on a fresh `CerberusConfig` that never went through autotune, `bool(cfg.tuned)` is `True`, `bool(cfg.tuned.read_length_class)` is `True`, and the function returns `cfg.tuned is r` → `True`. The "QC was skipped (dry-run)" fallback at `:205-206` can never execute, and its docstring documents behaviour that does not exist.
  3. `tuned.platform` (`cerberus/config.py:50`) is written by autotune and printed by `summary()`, but a grep for `tuned.platform` across `cerberus/` returns nothing outside `autotune.py` itself. It influences the run only indirectly, via `minimap2_preset`, which is computed inside autotune. No stage ever branches on it.
- **Trigger:** Always.
- **Consequence:** No runtime misbehaviour, but the code reads as though these paths matter, which is how F1 and F2 stayed invisible. `_tuned_or_default`'s docstring in particular asserts dry-run protection that does not exist, and its presence at `cerberus/orchestrator.py:87` suggests the dry-run gap (F2) was noticed and mis-fixed. `tests/test_autotune.py` passing while `min_length` is unreachable shows the tests validate the table, not the wiring.
- **Fix:** Delete the dead branches, or make them live. `_tuned_or_default` should check a real sentinel (e.g. `cfg.tuned is None` with `tuned: TunedParams | None = None` in `CerberusConfig`) so the dry-run fallback actually engages — which would also fix F2. Add an integration test that asserts the *command line* fastp receives, not just the contents of `TunedParams`; that single test would have caught F1, F2 and F13.
