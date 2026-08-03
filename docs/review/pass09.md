# Pass 9 — Documentation and packaging

## Summary

The README is well written and reads as if the tool were more finished than it is: roughly a quarter of its concrete claims do not match the code. The headline defect is that the advertised two-tier help does not exist — `cerberus --help` and `cerberus --help-all` produce byte-identical output (verified by `diff`, empty), because `main()` builds the *full* parser before argparse's `--help` action fires. Two documented power-user flags (`--kraken2-db`, `--aux-refs`) are parsed into the config and then never read by any stage, so they silently do nothing; and `fastplong` — the tool the README's architecture diagram names for long reads — is declared in neither `environment.yml` nor the conda recipe, so every `--long` run falls back to `chopper` and autotunes off a **hardcoded synthetic JSON**, making the "autotunes from the data" and "PacBio HiFi" claims false for long reads. Packaging is otherwise sound: the wheel builds cleanly and correctly ships `data/default_manifest.json`, the entry point works, the recipe's `test.commands` pass, and the recipe's `sha256` provably matches the real `v0.1.1` GitHub tarball — but `pyyaml` is a phantom dependency (zero `yaml` imports anywhere) and `noarch: python` is wrong for a package whose runtime deps are platform-specific bioconda binaries. Repository hygiene is thin: no CI, no CONTRIBUTING/CHANGELOG/CITATION, the 34 passing tests are never run automatically, and none of the upstream reference-data sources (T2T-CHM13, IPD-IMGT/HLA, Ensembl, RefSeq, NCBI taxonomy) are acknowledged or licence-checked anywhere despite the project redistributing derivatives of all of them.

## README claim audit

| Claim | Location in README | Status | Note |
|---|---|---|---|
| `<sample>.meta.R1/R2/orphans.fastq.gz` | :9, :241-243 | OK | `meta.py:74,75,90` |
| `<sample>.profiling.fastq.gz` | :10, :244 | OK | `profiling.py:138` |
| `<sample>.<mode>.*_GDPR.fastq.gz` | :11, :245-247 | OK | `gdpr.py:73,74,89` |
| Works on Illumina PE | :13 | OK | `autotune.py:83-84`, `_run_short` |
| Works on ONT | :13 | OK | `autotune.py:82`, auto-detected at `autotune.py:106` |
| Works on PacBio HiFi | :13 | WRONG | Preset exists (`autotune.py:78-79`) but auto-detection depends on a real fastp JSON that long-read runs never get — see F3 |
| Works on PacBio CLR | :13 | WRONG | Reachable *only* via explicit `--platform pacbio-clr`; `detect_platform_from_fastp` (`autotune.py:90-106`) can never return `PACBIO_CLR`. Contradicts "you do not need to know the platform" |
| "Autotunes its parameters from the data" | :13 | WRONG | True for short reads; for `--long` the JSON is fabricated (`qc.py:140-161`) — F3 |
| First run downloads ~22 GB | :38 | OK | manifest sizes sum to 23.5 GB (`default_manifest.json:12,24,35,46,57`) |
| `cerberus fetch-refs` | :38, :202 | OK | `cli.py:220-233` |
| `cerberus doctor` | :38, :203 | OK | `cli.py:236-285`, but see F16 |
| Architecture: `fastp / fastplong` | :59 | WRONG | `fastplong` in neither `environment.yml` nor `meta.yaml` — F3 |
| profiling → "bowtie2-vsl (or minimap2 fast)" | :62-63 | OK | `profiling.py:49-86` |
| GDPR = Kraken2 + minimap2 | :68-74 | OK | `gdpr.py:56-72` |
| bbmask `entropy=0.7 window=80 kr=5 minlen=40 mincount=4` | :86 | OK | `mask_t2t_hla.sh:51-54` |
| minimap2 asm5 viral, ≥85% id, `bedtools merge -d 50` | :87 | OK | `mask_t2t_hla.sh:60-66` |
| Kraken2 DB: human/chimp/gorilla/mouse/rat via `kraken:taxid` headers | :92-94 | OK | `build_kraken2_gdpr.sh:44-59` |
| aux refs = Ensembl 113 ncRNA + NC_012920.1 | :96-98 | OK | `curate_aux_refs.sh:28-65` |
| `--platform {auto,illumina,ont,pacbio-hifi,pacbio-clr}` | :109-112, :217 | OK | `config.py:19-24`, `cli.py:92-94` |
| `--fast` / `--double-pass` | :119, :122, :218-219 | OK | `profiling.py:49,61`; mutual exclusion at `orchestrator.py:51` |
| "Per-stage knobs (everything under `--help-all`)" | :125 | WRONG | Advanced group is printed by plain `--help` too — F1 |
| `--min-length` / `--min-quality` / `--entropy` / `--bbduk-k` | :129-135 | OK | `cli.py:128-135`, consumed at `autotune.py:145-163` |
| `--minimap2-args` / `--bowtie2-args` | :138-140 | OK | `align.py:50-51,156-157` |
| `--memory 24G` = "bbduk/Kraken2 memory cap" | :143, :216 | WRONG | Reaches bbduk only (`kmer.py:64,110`, `entropy.py:40,70`); `kraken2` is invoked with no memory flag (`kraken.py:52-63`) |
| `--keep-intermediates` keeps BAMs/FASTQs | :146, :231 | WRONG | Only BAMs are conditionally kept (`align.py:84`); intermediate FASTQs are *always* kept — F6 |
| `cerberus --help-all` shows every flag | :149, :205 | WRONG | Identical to `--help` — F1 |
| `build_custom_host_ref.sh -i … --mask -t 16` | :159-163 | OK | `build_custom_host_ref.sh:64-118` |
| Script is invocable as `bash scripts/…` | :159, :183-190 | WRONG | `scripts/` is not in the wheel or the conda package (verified wheel listing); only a git clone has it |
| `masked_t2t_hla.mmi` built by `minimap2 -x sr -d` | :174 | WRONG | Custom builder uses `-x sr` (`build_custom_host_ref.sh:114`) but the *published* bundle uses the default preset (`mask_t2t_hla.sh:81`) |
| `masked_t2t_hla.long.mmi` used for `--long` modes | :175 | WRONG | Never read; `long_read.py:29,61` use the `masked_t2t_hla_minimap2` asset — F12 |
| `masked_t2t_hla_bt2/…` from `bowtie2-build` | :176 | OK | `build_custom_host_ref.sh:118`; found by `profiling.py:163` |
| `manifest.json` makes RefManager skip download/verify | :177 | OK | `build_custom_host_ref.sh:166-180`; `refs.py:127-129` |
| No `--kraken-db` ⇒ `--gdpr` unavailable | :193 | OK | `refs.py:141-145` raises, though the message says "awaiting first release" |
| `--all` = `--meta --profiling --gdpr` | :212 | OK | `cli.py:171-175` |
| `--ref-dir`, `--no-auto-download` | :220-221 | OK | `cli.py:101-105` |
| `--kraken2-db` / `--aux-refs` (implied by ":125 everything") | help text, :125 | WRONG | Parsed and stored, never read — F2 |
| RAM peak ~9 GB on `--gdpr` | :229 | WRONG | Kraken2 loads a "~14 GB extracted" DB (:94) with no `--memory-mapping` — F10 |
| Disk (refs) ~13 GB extracted | :230 | WRONG | Contradicts both ":38 ~22 GB" and ":94 ~14 GB" for the Kraken2 DB alone — F10 |
| Run disk ~2× input, "cleaned automatically" | :231 | WRONG | Nothing ever removes `out/_work/` — F6 |
| "Designed for 16 GB laptops" | :233 | WRONG | Not supported by the DB size the README itself states — F10 |
| `reports/accounting.tsv` + `.json` | :249-250 | OK | `accounting.py:56-58` |
| `reports/fastp.json/html` | :251 | WRONG | Written to `out/_work/00_qc/` (`qc.py:39-40`) — F5 |
| `reports/*.flagstat.txt` | :252 | WRONG | Written to per-stage `_work` dirs (`align.py:48,106,153,206`) — F5 |
| `logs/*.log` one per stage + JSONL run log | :253-254 | OK | `shell.py:64,71`; `logger.py:50` |
| Long-read output filenames | — | MISSING | `<sample>.long_meta.fastq.gz` / `.long_profiling.fastq.gz` (`long_read.py:40,96`) documented nowhere |
| `<sample>.meta.GDPR.fastq.gz` (orphan GDPR) | — | MISSING | Produced by `gdpr.py:83-90`, absent from the output tree |
| `human_k27.fa.gz` asset (932 MB) | — | UNDOCUMENTED | Downloaded by `fetch-refs`, checked by `doctor`, used by nothing — F17 |
| `--update-refs`, `-s/--sample-id`, `-v/-q/--dry-run` | — | UNDOCUMENTED | Real flags (`cli.py:106-113`) absent from the Commands section |
| Upstream data licences / how-to-cite for T2T, IMGT/HLA, Ensembl, RefSeq | :259-261 | MISSING | Only Cerberus itself is cited — F15 |

## Verified working

- **Conda recipe `sha256` is genuine** — downloaded `https://github.com/iowa69/cerberus/archive/refs/tags/v0.1.1.tar.gz`; `sha256sum` = `c17944783acded3876b5d49a572ac77371b8e1737646b05f8b4c87f5d8d149bf`, exactly `conda-recipe/meta.yaml:10`. URL, tag `v{{version}}`, and version `0.1.1` are self-consistent with `pyproject.toml:7` and `cerberus/__init__.py:3`.
- **Wheel ships the packaged manifest** — built `cerberus_mg-0.1.1-py3-none-any.whl` into `work/p9/dist/`; `zipfile` listing contains `cerberus/data/default_manifest.json`, so `refs.py:87` `resources.files("cerberus.data")` resolves. `pyproject.toml:43-44` works (via both `package-data` and the `data/__init__.py` at `cerberus/data/__init__.py`).
- **Console entry point** — `pyproject.toml:34` → wheel `entry_points.txt` `cerberus = cerberus.cli:main`; installed and ran successfully.
- **Recipe `test:` block passes** — `cerberus --version` → `cerberus 0.1.1` (rc 0), `cerberus --help` (rc 0), `import cerberus` (rc 0), matching `conda-recipe/meta.yaml:44-48`.
- **All eight shell scripts are syntactically valid and strict** — `bash -n` passes on every file under `scripts/` and `scripts/build_refs/`; `grep -L "set -euo pipefail"` returns nothing, i.e. all eight carry it (`build_custom_host_ref.sh:16`, `link_local_refs.sh:6`, `smoke_test.sh:8`, `mask_t2t_hla.sh:17`, `build_kraken2_gdpr.sh:14`, `curate_aux_refs.sh:19`, `run_all_builds.sh:5`, `run_remaining.sh:4`).
- **Zenodo token handling is clean** — read from `CERBERUS_ZENODO_TOKEN` only (`zenodo_upload.py:75-79`, `zenodo_publish.py:29-33`), sent as a Bearer header, never printed, never written to the cache file (`zenodo_upload.py:105-107` persists only `deposition_id`/`bucket_url`/`created_at`).
- **`environment.yml` and the recipe's `run:` list are pin-consistent** — every shared package carries the identical lower bound (`environment.yml:9-24` vs `conda-recipe/meta.yaml:25-41`); the recipe correctly omits `curl`/`tar`/`pip`.
- **Test suite passes** — `pytest -q` → 34 passed, covering `_parse_memory`, mode validation, and every autotune branch (`tests/test_cli.py`, `tests/test_autotune.py`).
- **`--platform` choices exactly match the README** — `cli.py:92-94` renders `{auto,illumina,ont,pacbio-hifi,pacbio-clr}` in `--help`, identical to `README.md:217`.
- **Reference-build parameters match their prose** — `README.md:86` bbmask flags are verbatim `mask_t2t_hla.sh:51-54`; `README.md:87` viral masking matches `mask_t2t_hla.sh:60-71`; `README.md:98` aux-ref biotypes match `curate_aux_refs.sh:34-51`.
- **MIT `LICENSE` present and shipped** — `LICENSE:1-21`, included in the wheel as `cerberus_mg-0.1.1.dist-info/licenses/LICENSE` via `pyproject.toml:11`.
- **`.gitignore` correctly excludes build noise** — `git ls-files` shows 47 tracked files with no `cerberus_mg.egg-info/`, `__pycache__/`, `.ruff_cache/` or `.pytest_cache/` leakage; `git status --porcelain` is clean.

## Findings

### F1. Two-tier help does not exist: `--help` and `--help-all` are byte-identical
- **Severity:** high
- **Location:** `cerberus/cli.py:300` (with `cerberus/cli.py:114`, `cerberus/cli.py:123-126`)
- **What:** `main()` calls `_build_full_parser()`, which is `_basic_parser()` **plus** `_add_advanced_args()`. The `-h/--help` action is registered on that same object at `cli.py:114`, so argparse prints the full parser — advanced group included. `_print_full_help()` (`cli.py:148-151`) builds the identical parser, so `--help-all` prints exactly the same text.
- **Trigger:** `cerberus --help`.
- **Consequence:** The central UX promise of `cli.py:2-4` and `README.md:125, 149, 204-205` is not delivered. Verified: `cerberus --help > a; cerberus --help-all > b; diff a b` produces no output. Worse, the group is literally captioned `Advanced (--help-all only):` while appearing in brief help, and all nine advanced flags plus their metavars are dumped into the brief usage line — the opposite of the "most users only need this" intent.
- **Fix:** Have `main()` parse with the basic parser for `-h` purposes: register `-h` with `action="store_true"` and dispatch manually, or parse with `_basic_parser()` first and only build the full parser once no help was requested. Simplest: add the advanced group to a parser used solely for parsing, and give `-h` a custom action that prints `_basic_parser().format_help()`.

### F2. `--kraken2-db` and `--aux-refs` are parsed, stored, and silently ignored
- **Severity:** high
- **Location:** `cerberus/cli.py:140-143` → `cerberus/cli.py:210-211` → `cerberus/config.py:92-93`
- **What:** Both flags populate `cfg.kraken2_db_override` / `cfg.aux_refs_override`. A repo-wide grep for those attribute names returns only the definition and the assignment — no consumer. `gdpr.py:49` unconditionally does `refs.path_to(refs.asset("kraken2_gdpr_compact"))`; `profiling.py:90` and `long_read.py:82` unconditionally use `refs.asset("aux_refs")`.
- **Trigger:** `cerberus … --gdpr --kraken2-db /my/custom/db`.
- **Consequence:** Silent wrong behaviour, not an error. A user pointing at a custom Kraken2 database gets the bundled human/mammal DB and has no way to tell from the logs. For a tool whose selling point is publication-defensible host removal, silently substituting the reference is the worst possible failure mode. `--aux-refs` has the same shape. Both are additionally implied by `README.md:125` ("everything under `--help-all`").
- **Fix:** In `gdpr.py:49` use `cfg.kraken2_db_override or refs.path_to(...)`; in `profiling.py:90` and `long_read.py:82` use `cfg.aux_refs_override or refs.path_to(...)`. Until then, make `orchestrator.validate_config` raise `ConfigError` if either override is set.

### F3. `fastplong` is in no dependency list, so every long-read run autotunes off a fabricated JSON
- **Severity:** high
- **Location:** `environment.yml:9-24` and `conda-recipe/meta.yaml:25-41` (absent), consumed at `cerberus/stages/qc.py:102`, fallback at `cerberus/stages/qc.py:115-129`, fabrication at `cerberus/stages/qc.py:140-161`
- **What:** `grep -n fastplong environment.yml conda-recipe/meta.yaml` returns nothing. `run_fastplong` therefore always takes the `elif which("chopper")` branch, and writes a **synthetic** fastp report with `read1_mean_length: 5000`, `q20_rate: 0.95`, `q30_rate: 0.7` hardcoded (`qc.py:145-148`).
- **Trigger:** `cerberus --long -i reads.fq.gz -o out/ --all` on a stock `conda env create -f environment.yml` install.
- **Consequence:** `autotune_from_fastp` reads constants, not data. `classify_length(5000)` always yields `VERY_LONG` (`autotune.py:26-32`: `5000 < 5000` is false), and `detect_platform_from_fastp` always yields `ONT` (`autotune.py:104-106`: 0.7 < 0.85 and 0.95 < 0.99). A PacBio HiFi sample in `auto` mode silently gets `map-ont` instead of `map-hifi`, and `winnowmap_enabled` is forced on. This falsifies `README.md:13` ("Autotunes its parameters from the data"), `README.md:59` (the diagram naming fastplong), and the HiFi platform-support claim. `cerberus doctor` lists `fastplong` as merely optional (`cli.py:254`), reinforcing the misconception.
- **Fix:** Add `fastplong` to `environment.yml` and to the recipe's `run:` list; move it to the `required` list in `cli.py:252-253` for long-read use. Separately, when the chopper fallback is taken, either compute real length/quality statistics from the FASTQ (`autotune.estimate_long_read` already streams it) or refuse to autotune and require `--platform`.

### F4. PacBio CLR is unreachable by autodetection despite the "no need to know the platform" claim
- **Severity:** medium
- **Location:** `cerberus/autotune.py:90-106` (vs `README.md:13`, `README.md:112`)
- **What:** `detect_platform_from_fastp` has exactly three outcomes — `ILLUMINA` (mean < 500), `PACBIO_HIFI` (high Q), `ONT` (otherwise). `Platform.PACBIO_CLR` (`config.py:24`) is only ever produced by the user typing `--platform pacbio-clr`, which then routes to `map-pb` via `autotune.py:80-81` and `apply_user_overrides` (`autotune.py:156-158`).
- **Trigger:** Any PacBio CLR input run in the default `--platform auto`.
- **Consequence:** CLR reads are aligned with `map-ont` (or, per F3, always `map-ont`). CLR has ~10-15% error versus ONT's profile; the preset difference materially changes host-removal sensitivity — precisely the thing the GDPR head must not get wrong. The README lists CLR as a first-class supported platform on line 13 with no caveat.
- **Fix:** Either extend the heuristic (CLR is distinguishable from HiFi by q30 and from ONT by length/accuracy distribution) or amend `README.md:13` and `:217` to state that CLR requires an explicit `--platform pacbio-clr`.

### F5. The documented output tree is wrong: `reports/` holds only the accounting files
- **Severity:** medium
- **Location:** `README.md:248-252` vs `cerberus/stages/qc.py:39-40` and `cerberus/stages/align.py:48,106,153,206`
- **What:** `fastp.json` / `fastp.html` are written into `workdir`, which the orchestrator sets to `cfg.work_dir / "00_qc"` = `out/_work/00_qc/` (`orchestrator.py:126-127`). `*.flagstat.txt` are written into each stage's `_work` subdirectory. Only `accounting.json` and `accounting.tsv` reach `reports/` (`accounting.py:56-58`, called from `orchestrator.py:96`).
- **Trigger:** Any run; confirmed by executing the pipeline until it aborted on a missing binary — `out/reports/` was created and empty.
- **Consequence:** Users and reviewers following the README look in `reports/` for the QC report and alignment statistics and find nothing. Because the flagstat files live under `_work/`, they are exactly the artefacts a reviewer would want and exactly the ones buried in the scratch tree.
- **Fix:** Copy or symlink `fastp.json`/`fastp.html` and the `*.flagstat.txt` files into `cfg.reports_dir` at the end of `orchestrator.run()`, or correct the README tree.

### F6. "Cleaned automatically unless `--keep-intermediates`" — nothing is ever cleaned
- **Severity:** medium
- **Location:** `README.md:231`; `cerberus/stages/align.py:84-85, 129-130, 183-184, 225-226`
- **What:** A repo-wide grep for `rmtree` returns zero hits. The only cleanup is `bam.unlink(missing_ok=True)` in the four aligner wrappers. Every intermediate FASTQ — `01_minimap2.unmapped.R{1,2}.fq.gz`, `03_bbduk_kmer.*`, `04_entropy.*`, the Kraken2 unclassified outputs — persists under `out/_work/` forever. `--keep-intermediates` only additionally preserves BAMs (`align.py:84`) and bbduk matched-read files (`kmer.py:58-59, 105`).
- **Trigger:** Any successful run.
- **Consequence:** On `--all`, the meta head alone writes 3 intermediate FASTQ pairs and profiling writes 4 more, plus GDPR's Kraken2 and minimap2 stages. Real disk use is closer to 8-10× the input, not the documented "~2× … cleaned automatically". On the "16 GB laptop" the README targets, a modest 50 GB FASTQ pair silently fills the disk.
- **Fix:** Add `shutil.rmtree(cfg.work_dir, ignore_errors=True)` guarded by `not cfg.keep_intermediates` at the end of `orchestrator.run()` (after `accounting.write`), or correct the README.

### F7. `noarch: python` is wrong for a package whose runtime deps are architecture-specific binaries
- **Severity:** medium
- **Location:** `conda-recipe/meta.yaml:14` with `conda-recipe/meta.yaml:24-41`
- **What:** The build is declared `noarch: python` — one artefact installed on every platform — while `run:` requires `fastp`, `minimap2`, `bowtie2`, `bbmap`, `kraken2`, `samtools`, `seqkit`, `chopper`, `winnowmap`, `bedtools`, `aria2`, `multiqc`. These are compiled bioconda packages with no `win-64` builds and, for several, no `osx-arm64` build. The project's own metadata contradicts the recipe: `pyproject.toml:18` declares `Operating System :: POSIX :: Linux`.
- **Trigger:** `conda install -c bioconda cerberus-mg` on Windows or Apple Silicon — the package resolves as installable (noarch) and then fails to solve its dependencies.
- **Consequence:** bioconda's linter rejects noarch recipes with platform-restricted run dependencies; this will block the merge the README's install line at `:25` is waiting on. Even if merged, the package advertises availability on platforms where it cannot work.
- **Fix:** Drop `noarch: python`, add `skip: true  # [not linux]` (or `[win]`), and let the recipe build per-platform. Since `build.script` is a plain `pip install`, no other change is needed.

### F8. Both Zenodo scripts hardcode absolute paths that do not exist
- **Severity:** medium
- **Location:** `scripts/zenodo_upload.py:24` and `scripts/zenodo_publish.py:18`
- **What:** `BUILD_DIR = Path("/home/iowa/Desktop/cerberus/scripts/build_refs/build")` and `MANIFEST_PATH = Path("/home/iowa/Desktop/cerberus/cerberus/data/default_manifest.json")`. The repository actually lives at `/home/iowa/Desktop/cerberus/repo/`, so both paths are missing the `repo/` component. Verified: both report `MISSING` on a filesystem check.
- **Trigger:** `CERBERUS_ZENODO_TOKEN=… python scripts/zenodo_upload.py`, exactly as the module docstring at `zenodo_upload.py:5` documents.
- **Consequence:** `zenodo_upload.py` exits at `:202` with `Missing asset: …` for the first file; `zenodo_publish.py` raises `FileNotFoundError` at `:67` *after* having already POSTed the irreversible `actions/publish` call at `:55` — the deposition is published and the DOI minted, but the manifest is never updated and the run aborts with a traceback. These scripts cannot work on any machine, including the author's.
- **Fix:** Derive both from `Path(__file__).resolve().parents[1]`, and make them overridable by env var/CLI argument.

### F9. `zenodo_publish.py` never writes `sha256` — the only checksum RefManager verifies
- **Severity:** medium
- **Location:** `scripts/zenodo_publish.py:81-89` (vs `cerberus/refs.py:127-129, 158-165`)
- **What:** The publisher sets `entry["url"]`, `entry["size_bytes"]`, and `entry["md5"]`, taking the md5 from Zenodo's `checksum` field. It never touches `entry["sha256"]`. But `RefManager.is_satisfied` verifies `asset.sha256` (`refs.py:127-128`) and `_download` verifies `asset.sha256` (`refs.py:158`); nothing anywhere reads `md5`.
- **Trigger:** Following the documented release flow (`zenodo_upload.py` → `zenodo_publish.py`).
- **Consequence:** After a release, `sha256` is left at whatever it was. If that is `""` or `"PENDING"`, `refs.py:129` returns `True` unconditionally and `refs.py:166-168` merely logs a warning — integrity verification of a 22 GB download is silently disabled, on the exact code path that protects the "zero human reads" guarantee. The current `default_manifest.json` has real sha256 values, so they were filled in by hand outside this script; the next release will not be. The script also hardcodes `manifest["release"] = "0.1.0"` at `:68` while the package is `0.1.1`.
- **Fix:** Compute the sha256 locally during `zenodo_upload.py` (it already streams every file) and have `zenodo_publish.py` write it; derive `release` from `cerberus.__version__`.

### F10. Memory and disk figures are internally contradictory and understate Kraken2 RAM
- **Severity:** medium
- **Location:** `README.md:229-233` (vs `README.md:38`, `README.md:94`, `cerberus/stages/kraken.py:52-63`)
- **What:** Three mutually inconsistent numbers: `:38` "~22 GB" bundle, `:94` "Final DB is ~14 GB extracted", `:230` "Disk (refs) ~13 GB extracted" for *all* references. Summing `default_manifest.json` (`:12,24,35,46,57`) gives 23.5 GB **compressed**; extracted is materially larger. Separately, `:229` claims "~9 GB peak on `--gdpr` (Kraken2)", but `kraken.py:52-63` passes no `--memory-mapping`, so Kraken2 loads the whole `hash.k2d` into RAM — i.e. the ~14 GB the README itself states.
- **Trigger:** A user sizing a machine from the README, then running `--gdpr` on a 16 GB laptop.
- **Consequence:** OOM on the machine class the README explicitly targets ("Designed for 16 GB laptops. Tested on 4-core/16 GB", `:233`). The disk figure understates the requirement by more than 2×, on top of F6.
- **Fix:** Recompute all three figures from the actual extracted sizes; either pass `--memory-mapping` to Kraken2 when `cfg.memory_gb` is below the DB size, or state the true RAM requirement.

### F11. `pyyaml` is a phantom dependency and `package-data` globs a non-existent file type
- **Severity:** medium
- **Location:** `pyproject.toml:27` and `pyproject.toml:44` (also `environment.yml:10`, `conda-recipe/meta.yaml:27`)
- **What:** `grep -rn "import yaml\|from yaml\|yaml\." --include=*.py` over the whole tree returns nothing. The only occurrence of "yml" in Python is a help string at `cli.py:283`. `[tool.setuptools.package-data]` also lists `data/*.yaml`; `cerberus/data/` contains only `default_manifest.json` and `__init__.py`.
- **Trigger:** Any install.
- **Consequence:** Every consumer — pip, conda, the environment file — pulls a dependency the code never uses, expanding the install footprint and the security surface for no benefit, and (for the bioconda recipe) adding a needless constraint to the solve. Separately, `tqdm` (`pyproject.toml:26`) is a *hard* dependency but is only ever imported inside `try/except ImportError` with a `None` fallback (`refs.py:188-190`, `hashing.py:8-10`) — it belongs in an optional extra.
- **Fix:** Delete `pyyaml` from `pyproject.toml:27`, `environment.yml:10` and `conda-recipe/meta.yaml:27`; drop `data/*.yaml` from `pyproject.toml:44`; move `tqdm` to an extra such as `progress = ["tqdm>=4.66"]`.

### F12. `masked_t2t_hla.long.mmi` is documented, built, and never used; the published `.mmi` preset contradicts the docs
- **Severity:** medium
- **Location:** `README.md:174-175` vs `cerberus/refs.py:37-45`, `cerberus/pipelines/long_read.py:29,61`, `scripts/build_custom_host_ref.sh:124`, `scripts/build_refs/mask_t2t_hla.sh:81`
- **What:** The README's custom-reference table states `masked_t2t_hla.long.mmi` is "Used for `--long` modes". No code path references `.long.mmi`: `_PIPELINE_TO_ASSETS["long-meta"]` and `["long-profiling"]` both map to `masked_t2t_hla_minimap2`, whose filename is `masked_t2t_hla.mmi` (`default_manifest.json:9`). The table also claims `masked_t2t_hla.mmi` is built by `minimap2 -x sr -d` — true for `build_custom_host_ref.sh:114`, but the *published* bundle is built by `minimap2 -d` with the default preset (`mask_t2t_hla.sh:81`).
- **Trigger:** `bash scripts/build_custom_host_ref.sh -i mouse.fa -o ./refs` then `cerberus --long --ref-dir ./refs`.
- **Consequence:** The user waits for a second full index build that is then ignored; long reads are aligned against the short-read `-x sr` index. Conversely, for the published bundle, short reads are aligned against a default-preset index — minimap2's indexing parameters (k, w) are baked into the `.mmi` and override the `-ax sr` preset, so `-ax sr` does not actually give short-read indexing behaviour. Either way the documented mapping is not what happens.
- **Fix:** Add a distinct `masked_t2t_hla_minimap2_long` asset key wired into the `long-*` entries of `_PIPELINE_TO_ASSETS`, or delete the `.long.mmi` row from `README.md:175` and stop building it. Reconcile the preset used by `mask_t2t_hla.sh:81` with `README.md:174`.

### F13. `smoke_test.sh` cannot run anywhere but the author's machine, and its inputs are deleted by the build it depends on
- **Severity:** medium
- **Location:** `scripts/smoke_test.sh:20` (also `scripts/build_refs/run_all_builds.sh:11`, `scripts/build_refs/run_remaining.sh:8`)
- **What:** `source /home/iowa/miniconda3/etc/profile.d/conda.sh` under `set -euo pipefail` (`:8`) aborts immediately on any machine without that exact path. `run_all_builds.sh:22,26,30` and `run_remaining.sh:19,23` additionally hardcode `df -h /home/iowa`. Compounding this, `smoke_test.sh:28` reads `$BUILD/masked_t2t_hla/_work/t2t.fa.gz` and `:38-39` reads `$BUILD/masked_t2t_hla/_work/virus.fa.gz`, but `mask_t2t_hla.sh:110` does `rm -rf "$WORK"` — that directory never survives a completed build.
- **Trigger:** `bash scripts/smoke_test.sh` on a fresh clone.
- **Consequence:** The only end-to-end test in the repository is unrunnable, which is presumably why nothing runs it (F14). The `||` fallback at `:32` substitutes `$REF_DIR/human_k27.fa.gz` — an asset no pipeline uses (F17) and which only exists after `fetch-refs`, so on a normal install even the fallback fails. The 6 output names it asserts (`:75-82`) are correct, so the check itself is sound; only the plumbing is broken.
- **Fix:** Replace the hardcoded `source` with `command -v conda >/dev/null || { echo "activate the cerberus env first" >&2; exit 1; }`, or honour `$CONDA_PREFIX`. Point read simulation at a small committed fixture or download it on demand rather than at a deleted `_work/` directory. Replace `df -h /home/iowa` with `df -h .`.

### F14. No CI, no CONTRIBUTING/CHANGELOG/CITATION, and the test suite is never run automatically
- **Severity:** medium
- **Location:** repository root (`/home/iowa/Desktop/cerberus/repo/`); `pyproject.toml:30-31`
- **What:** `.github/` does not exist (verified). `git ls-files` lists 47 files: no `CONTRIBUTING.md`, `CHANGELOG.md`, `CITATION.cff`, `CODE_OF_CONDUCT.md`, or issue templates. A `dev` extra declaring `pytest`, `pytest-cov` and `ruff` exists at `pyproject.toml:31` but nothing invokes it; `[tool.ruff]` (`:46-48`) and `[tool.pytest.ini_options]` (`:50-52`) are configured yet unenforced.
- **Trigger:** Any push or pull request.
- **Consequence:** The 34 tests pass today by luck of the working tree; nothing prevents a regression, and nothing would have caught F1 (a two-line `diff` in CI). `pyproject.toml:37-38` advertises a public Issues URL with no template. For a tool asking reviewers to trust a "zero human reads" guarantee, the absence of any automated verification is itself a credibility problem.
- **Fix:** Add `.github/workflows/ci.yml` running `ruff check .` and `pytest` on 3.10-3.12, plus a `python -m build` + `twine check` job. Add `CHANGELOG.md` (the git log already reads like one) and `CITATION.cff`.

### F15. Upstream reference-data licences and citations are acknowledged nowhere
- **Severity:** medium
- **Location:** `README.md:78-98` and `README.md:259-261`; `scripts/zenodo_upload.py:57-64`
- **What:** Cerberus builds and redistributes derivatives of T2T-CHM13v2.0 (`mask_t2t_hla.sh:24`), IPD-IMGT/HLA (`mask_t2t_hla.sh:25`), RefSeq viral (`mask_t2t_hla.sh:26`), Ensembl 113 ncRNA (`curate_aux_refs.sh:28`), NCBI taxonomy plus chimp/gorilla/mouse/rat assemblies (`build_kraken2_gdpr.sh:26,55-59`), and the README names UHGG at `:90`. The build scripts document the URLs but no file states the sources' licences or citations. The Zenodo deposition is stamped `"license": "cc-by-4.0"` with a single creator (`zenodo_upload.py:58,64`) and no upstream attribution in its description (`:39-56`). The README's "Citation" section (`:259-261`) cites only Cerberus.
- **Trigger:** Publishing the reference bundle; any user citing the pipeline in a paper.
- **Consequence:** Two distinct problems. (a) Attribution: IPD-IMGT/HLA, Ensembl and the T2T consortium all request citation, and reviewers of a host-removal paper will ask which reference build was used. (b) Licence compatibility: relicensing a derived work as CC BY 4.0 requires that every input permit it — IPD-IMGT/HLA is distributed under CC BY-ND terms, under which redistributing a masked and re-indexed derivative needs explicit checking before the deposition is made public. This has not visibly been done.
- **Fix:** Add a `REFERENCES.md` (or a README section) listing each source with its URL, version, licence and canonical citation; mirror it into the Zenodo `description` and `related_identifiers`; verify the IPD-IMGT/HLA terms before asserting CC BY 4.0 on the bundle; add a `CITATION.cff` and a "How to cite" section covering both Cerberus and its data sources.

### F16. `--dry-run` requires the full external toolchain and fails with a raw traceback
- **Severity:** low
- **Location:** `cerberus/stages/qc.py:32` (with `cerberus/utils/shell.py:66`, `cerberus/cli.py:316-324`)
- **What:** `run_fastp` calls `require_tools("fastp")` before any command is built, whereas the `dry_run` short-circuit lives inside `shell.run` at `:66`. So `--dry-run` cannot get past QC without the binaries installed. `require_tools` raises `RuntimeError` (`shell.py:41`), which `main()` does not catch — it only handles `ConfigError` and `KeyboardInterrupt`.
- **Trigger:** `cerberus -r1 a.fq.gz -r2 b.fq.gz -o out --all --dry-run` on a machine without the conda env. Verified: full Python traceback ending in `RuntimeError: Missing required tool(s): fastp`.
- **Consequence:** The flag documented as "Print commands without executing" (`cli.py:113`) cannot serve its main purpose — previewing a run before provisioning. The unhandled traceback also affects the normal path: any mid-run missing tool or `ToolError` surfaces as a stack trace rather than a message. Similarly `cerberus doctor` is the documented way to find missing tools (`README.md:38`) but the pipeline itself does not fail gracefully.
- **Fix:** Skip `require_tools` when `cfg.dry_run` is set (in all stage wrappers), and add `except (RuntimeError, ToolError)` to `cli.py:316-324` to print a clean message and return a non-zero code.

### F17. `human_kmer_set` (932 MB) is downloaded and health-checked but used by nothing, and its description is false
- **Severity:** low
- **Location:** `cerberus/data/default_manifest.json:52-62`; `cerberus/refs.py:37-45, 225-235`; `scripts/build_refs/curate_aux_refs.sh:87`
- **What:** The manifest describes the asset as a "Human-specific 27-mer set for GDPR belt-and-braces bbduk scrub". `curate_aux_refs.sh:87` is `cp "$T2T" "$HUMAN27"` — it is a byte copy of the whole T2T-CHM13v2.0 genome, not a k-mer set. No entry in `_PIPELINE_TO_ASSETS` references `human_kmer_set`, and the GDPR pass uses minimap2, not bbduk (`gdpr.py:65-72`). `zenodo_upload.py:51-52` repeats the false description publicly.
- **Trigger:** `cerberus fetch-refs`, which iterates every manifest key (`refs.py:226`); or `cerberus doctor`, which does the same (`refs.py:231`).
- **Consequence:** `fetch-refs` downloads 932 MB nobody needs. Worse, a user who runs the pipeline normally (which fetches only required assets, `orchestrator.py:73-75`) and then runs `cerberus doctor` is told `missing or corrupt: human_kmer_set` on a perfectly healthy install — undermining the command the README recommends for validating an installation (`:38`).
- **Fix:** Remove the asset from the manifest, or have `fetch_all`/`doctor` operate on the union of `_PIPELINE_TO_ASSETS` rather than all manifest keys. Correct the description in both the manifest and the Zenodo metadata.

### F18. `run_all_builds.sh` and `run_remaining.sh` exit non-zero after a successful build
- **Severity:** low
- **Location:** `scripts/build_refs/run_all_builds.sh:36-38`, `scripts/build_refs/run_remaining.sh:29-31`
- **What:** The final statement is `for f in <globs>; do [ -f "$f" ] && echo … && cat "$f"; done`. If the last glob matches nothing it stays literal, `[ -f ]` fails, and the loop's exit status — hence the script's — is 1. Verified: `bash -c 'set -euo pipefail; for f in build/*/manifest_fragment.json; do [ -f "$f" ] && …; done'` returns rc=1.
- **Trigger:** Any partial build, e.g. re-running `run_remaining.sh` when `build/aux_refs/` has not been produced.
- **Consequence:** A multi-hour reference build that fully succeeded reports failure, and any CI or wrapper checking `$?` treats it as broken. `run_all_builds.sh:36` also globs `build/*/manifest_fragment.json`, which never matches the aux-refs fragment (`curate_aux_refs.sh:79` writes `aux_refs_manifest_fragment.json`), so that fragment is silently never displayed.
- **Fix:** Append `|| true` to the loop body, or restructure as `for f in …; do [ -f "$f" ] || continue; …; done` and add an explicit `exit 0`. Fix the glob to `build/*/*manifest_fragment.json`.

### F19. `.gitignore` patterns that can silently swallow real source files
- **Severity:** low
- **Location:** `.gitignore:9`, `.gitignore:36-47`, `.gitignore:54-55`
- **What:** Verified with `git check-ignore -v`: `tests/data/tiny.fastq.gz` → ignored by `:38 *.fastq.gz`; `cerberus/build/x.py` and `docs/build/index.html` → ignored by `:9 build/`; `notes.log` → `:54 *.log`; `logs/keep.md` → `:55 logs/`.
- **Trigger:** Adding an integration-test fixture, a nested `build/` package, or documentation under `logs/`.
- **Consequence:** Committing a small FASTQ/BAM fixture — precisely what this project needs to make `smoke_test.sh` runnable in CI (F13, F14) — fails silently: `git add` reports nothing and the file never lands in the repo. The `*.mmi`/`*.bt2`/`*.k2d` rules (`:44-47`) have the same effect on any tiny toy index.
- **Fix:** Add negations for test fixtures, e.g. `!tests/data/**`, and scope the artefact patterns to the directories where they actually appear (`/build/`, `cerberus_out/`, `scripts/build_refs/build/`).

### F20. Assorted dead code and documentation drift
- **Severity:** low
- **Location:** `cerberus/pipelines/profiling.py:9`; `cerberus/cli.py:12-13,26`; `cerberus/stages/kraken.py:21`; `cerberus/pipelines/long_read.py:24,40`
- **What:** (a) `profiling.py:9` documents a mode selector `--aligner minimap2`; no such flag exists anywhere — a grep for "aligner" finds only this docstring, an unrelated `--double-pass` help string, and two prose mentions in the README. (b) `cli.py:26` defines `_SUBCOMMANDS = {"fetch-refs", "doctor", "run", "help-all"}` which is never referenced; `cli.py:12-13` says the dispatch exists "so users can run … without typing `run`", implying `cerberus run -r1 …` works — it does not, argparse rejects `run` as an unrecognised argument. (c) `kraken.py:21` `GDPR_DROP_TAXA` is unused. (d) `long_read.py:24` sets `mode = "long-meta"` (hyphen) while `:40` writes `{sample}.long_meta.fastq.gz` (underscore); since `gdpr.py:98` interpolates `pipeline_result.mode`, the long-read GDPR output is named `{sample}.long-meta.long_GDPR.fastq.gz`, mixing both conventions in one filename.
- **Trigger:** Reading the module docs; running `cerberus run …`; running `--long --gdpr`.
- **Consequence:** Individually cosmetic, collectively they signal drift between the docs and the code and make the module docstrings unreliable as a specification — the same class of defect as F1 and F2 but without functional impact. The filename inconsistency in (d) will break any downstream glob written against either convention.
- **Fix:** Delete the `--aligner` mention from `profiling.py:9`; remove `_SUBCOMMANDS` and `GDPR_DROP_TAXA` or wire them up; pick one separator for long-read mode strings and filenames; document the long-read outputs in the README's output tree.
