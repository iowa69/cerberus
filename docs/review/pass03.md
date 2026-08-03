# Pass 3 — Reference manager

## Summary

The reference manager's happy path (short-read `--meta` / `--profiling` / `--gdpr`) maps assets correctly, and SHA256 verification genuinely catches corrupt non-archive files. Everything else is weaker than it looks. `_PIPELINE_TO_ASSETS` (`cerberus/refs.py:37`) is missing the `long-profiling-fast` key that `_required_pipeline_keys` (`cerberus/orchestrator.py:112`) can generate, so `cerberus --long --profiling --fast` ensures **zero** assets and then immediately opens a reference that was never downloaded. Worse, archives are exempted from hash verification entirely by an early return in `is_satisfied` (`cerberus/refs.py:123-124`): extraction is non-atomic, so an interrupted or silently-truncated Kraken2 DB is permanently reported as "present and verified" by both `ensure()` and `doctor()` — verified empirically with a 9-byte `hash.k2d`. That directly undermines the GDPR mode's "zero detectable human reads" claim, which is the tool's headline guarantee. The re-hash cost the brief asked me to quantify turned out to be **bounded at 7.7 GB (~3.5–7.0 s measured here, ~64 s on an HDD-backed laptop)** rather than the ~23 GB implied — precisely *because* the archive short-circuit skips the other 15 GB, so the thing that saves the time is the same bug that makes the cache unsafe.

## Verified working

- **Short-read asset mapping is complete and correct** — I enumerated all 36 valid `--meta`/`--profiling`/`--gdpr` × `--fast`/`--double-pass` × `--long` combinations against `_PIPELINE_TO_ASSETS` (`cerberus/refs.py:37-45`) and against the assets each pipeline function actually opens. All 18 short-read combinations ensure a superset of what is used.
- **`profiling-fast` is a genuinely correct subset** — `run_profiling` under `cfg.fast` uses only the minimap2 index (`cerberus/pipelines/profiling.py:47`); bowtie2 (`:77`) and the aux k-mer pass (`:89-90`) are both correctly gated behind `not cfg.fast`, matching the single-asset mapping at `cerberus/refs.py:40`.
- **`--double-pass` (short read) needs no extra assets** — the minimap2 pre-filter at `cerberus/pipelines/profiling.py:61-70` reuses `mm2_idx` already ensured by the `profiling` key (`cerberus/refs.py:39`).
- **GDPR asset mapping matches the code** — `run_gdpr_for` opens exactly `kraken2_gdpr_compact` and `masked_t2t_hla_minimap2` (`cerberus/pipelines/gdpr.py:49,51`), matching both `gdpr` and `long-gdpr` entries (`cerberus/refs.py:41,44`).
- **SHA256 verification does catch corrupt non-archive files** — wrote garbage over a `masked_t2t_hla.mmi` in a scratch ref dir; `is_satisfied` returned `False` via `verify_sha256` (`cerberus/refs.py:127-128`) and re-queued the download.
- **An *empty* extracted dir is correctly rejected** — `any(target.iterdir())` (`cerberus/refs.py:124`) returns `False`, so a zero-progress extraction does re-download. Only *partial* extractions slip through.
- **`.tar.zst` vs `.tar.gz` branch selection is correct** — probed every manifest filename through the `_extract` predicate (`cerberus/refs.py:210-223`): both shipped `.tar.zst` assets route to the zstd pipe branch, and `aux_refs.fa.gz` / `human_k27.fa.gz` are correctly classified `is_archive=False` (`cerberus/refs.py:59-60`) so they are never extracted.
- **`.tmp` naming and atomic rename are sound** — `target.with_suffix(target.suffix + ".tmp")` (`cerberus/refs.py:150`) produces correct names for all five manifest filenames (`.mmi.tmp`, `.tar.zst.tmp`, `.fa.gz.tmp`), and `tmp.rename(target)` (`cerberus/refs.py:170`) is a same-filesystem POSIX rename, so a truncated download is never visible under the final name.
- **Hash mismatch deletes the bad file** — `cerberus/refs.py:160-165` unlinks the `.tmp` before raising, so a corrupted download cannot be promoted.
- **`extracted_dirname` agrees with `_extract`'s output dir** — both use `split(".")[0]` (`cerberus/refs.py:64` and `:208`), so `path_to()` and the extractor never disagree for the shipped manifest.
- **`tarfile.extractall(filter="data")` works on the installed interpreter** — Python 3.13.9, `tarfile.data_filter` present, `filter` in the `extractall` signature.

## Findings

### F1. `--long --profiling --fast` generates a pipeline key that maps to no assets — the run downloads nothing

- **Severity:** critical
- **Location:** `cerberus/orchestrator.py:112`, `cerberus/refs.py:37-45`
- **What:** `_required_pipeline_keys` builds the profiling key as `f"{prefix}profiling-fast"` where `prefix = "long-"` when `cfg.long_mode`. That yields `"long-profiling-fast"`, which is **not** a key in `_PIPELINE_TO_ASSETS` (keys are `gdpr, long-gdpr, long-meta, long-profiling, meta, profiling, profiling-fast`). `required_assets_for` uses `.get(pkey, [])` (`cerberus/refs.py:109`), so the unknown key silently contributes an empty list — no warning, no error.
- **Trigger:** Any of the 4 combinations containing `--long --profiling --fast`. I enumerated them:
  - `--long --profiling --fast` → ensured assets: **NONE**
  - `--long --profiling --fast --gdpr` → ensures only `kraken2_gdpr_compact`, `masked_t2t_hla_minimap2`
  - `--long --meta --profiling --fast` → ensures only `masked_t2t_hla_minimap2`
  - `--long --meta --profiling --fast --gdpr` → ensures only `kraken2_gdpr_compact`, `masked_t2t_hla_minimap2`
- **Consequence:** In the first case `refs.ensure([])` (`cerberus/orchestrator.py:75`) is a no-op — on a fresh machine Cerberus prints "Required ref-asset groups: ['long-profiling-fast']", downloads nothing at all, runs fastplong, then `run_long_profiling` calls `refs.path_to(refs.asset("masked_t2t_hla_minimap2"))` (`cerberus/pipelines/long_read.py:61`) and hands minimap2 a path that does not exist. Compounding this: `run_long_profiling` **never reads `cfg.fast`** — the function body (`cerberus/pipelines/long_read.py:50-103`) has no `fast` branch — so `--fast` does not actually simplify the long-read path. It still runs the aux k-mer pass whenever `tuned.bbduk_aux_enabled`, which is `True` for both long classes (`cerberus/autotune.py:64,70`). So in **all four** combinations `aux_refs` is used but never ensured, even the ones that look partially fine. With `--no-auto-download` the guard at `cerberus/refs.py:136-140` never fires either, because there is nothing in the list to check.
- **Fix:** Two independent defects. (a) Make an unmapped pipeline key a hard error, not a silent empty list: in `required_assets_for` raise `RefManagerError(f"No asset mapping for pipeline key {pkey!r}")` instead of `.get(pkey, [])` at `cerberus/refs.py:109` — this converts every future mapping gap from a mid-run crash into a startup error. (b) Add `"long-profiling-fast": ["masked_t2t_hla_minimap2", "aux_refs"]` to `_PIPELINE_TO_ASSETS`, and either honour `cfg.fast` inside `run_long_profiling` or make `_required_pipeline_keys` emit plain `long-profiling` when `cfg.long_mode` (since `--fast` is a no-op there). Add a unit test that asserts every key `_required_pipeline_keys` can emit, over the full flag product, exists in `_PIPELINE_TO_ASSETS`.

### F2. A half-extracted or truncated archive is permanently reported as "present and verified"

- **Severity:** critical
- **Location:** `cerberus/refs.py:121-129` (short-circuit at `:123-124`), `cerberus/refs.py:206-221` (non-atomic extraction)
- **What:** `is_satisfied` returns `target.is_dir() and any(target.iterdir())` for archives and returns **before** reaching the SHA256 branch at `:127-128`. So archive assets are hash-verified exactly once — inside `_download` (`:158-165`) — and never again. Extraction is not atomic: `_extract` does `out_dir.mkdir(parents=True, exist_ok=True)` (`:209`) and then writes members directly into the final location, with no temp dir and no completion marker.
- **Trigger:** Any interruption during extraction of the 11 GB Kraken2 DB or the 3.8 GB bowtie2 index — Ctrl-C, OOM kill, laptop suspend, or a full disk (this test machine has 22 G free on `/` against a ~42 GB steady-state requirement, see F11). Also any silent truncation per F3.
- **Consequence:** Verified empirically. I created `kraken2_gdpr_compact/hash.k2d` containing 9 bytes, with no `taxo.k2d` or `opts.k2d`:
  - `is_satisfied()` → `True`
  - `ensure()` → returns OK, no download, no hash check
  - `doctor()` → does **not** list it among problems
  - `_find_kraken_db()` (`cerberus/pipelines/gdpr.py:132-139`) → accepts the directory, because it only tests `(db_dir / "hash.k2d").exists()` with no size or sibling-file check
  There is no recovery path: the state is self-perpetuating across every subsequent run. Because this is the Kraken2 DB used for the GDPR scrub, the consequence is a `_GDPR.fastq.gz` that silently under-filters human reads while the README (`README.md:11`) advertises "**Zero detectable human reads** via dual orthogonal mechanisms". The second mechanism (minimap2) is hash-verified, so the failure is partial and therefore invisible — output is produced, read counts look plausible, and nothing errors. `cerberus doctor`, which `README.md:38` presents as "validate the installation", cannot detect it.
- **Fix:** Extract to a sibling temp dir (`<name>.partial-<pid>`) and `os.rename` it into place only after the extractor exits 0 — rename of a directory is atomic on POSIX. Write a `.cerberus-complete` stamp file inside the extracted dir recording the source archive's SHA256, and have `is_satisfied` require that stamp (and, cheaply, that it matches the manifest hash) rather than `any(target.iterdir())`. Separately, harden `_find_kraken_db` to require `hash.k2d`, `taxo.k2d` and `opts.k2d` and a non-trivial `hash.k2d` size. Offer `cerberus doctor --deep` that re-hashes the retained tarballs.

### F3. `pipe()` only checks the last command's exit code, so a corrupt `.tar.zst` "extracts successfully"

- **Severity:** high
- **Location:** `cerberus/utils/shell.py:132-139`, used by `cerberus/refs.py:217-221`
- **What:** `pipe()` waits on every process (`for p in procs: p.wait()`) but then evaluates only `rc = procs[-1].returncode`. Intermediate exit codes are discarded entirely.
- **Trigger:** `zstd -dc <archive> | tar -x -C <out>` where `zstd` fails — a truncated download that nonetheless passed SHA256 is impossible, but a disk read error, an out-of-space condition mid-decompress, or a zstd frame error on a filesystem-corrupted tarball all cause `zstd` to exit non-zero after emitting a partial stream. `tar` happily extracts the partial stream and exits 0.
- **Consequence:** Confirmed empirically — `pipe([["sh","-c","echo partial; exit 1"], ["cat"]])` returned `returncode=0` and raised no `ToolError`. `_extract` therefore returns normally, `_download` returns normally, and the half-populated directory is then locked in as "satisfied" forever by F2. This is the most likely real-world route into the F2 state, and it is completely silent.
- **Fix:** In `pipe()`, collect all return codes and raise `ToolError` for the first non-zero one, e.g. `for cmd, p in zip(cmds, procs): if p.returncode != 0: raise ToolError(cmd, p.returncode, log_path)`. Note SIGPIPE (`-13`) on upstream stages should be tolerated when a downstream stage exits early, but here the downstream is `tar`, which consumes the whole stream. Additionally pass `--long=31 -c` / use `tar --zstd` directly (GNU tar ≥1.31) to remove the pipe entirely.

### F4. Long-read modes align against the short-read (`-x sr`) minimap2 index; the documented long-read index does not exist

- **Severity:** high
- **Location:** `cerberus/refs.py:42-44`, `cerberus/pipelines/long_read.py:29,61`, `README.md:174-176`
- **What:** `README.md:175` documents `masked_t2t_hla.long.mmi`, "built by `minimap2 -d` (default preset)", "Used for `--long` modes". No such asset exists in `cerberus/data/default_manifest.json`, and the string `long.mmi` appears nowhere in the codebase (grepped `*.py`, `*.json`). `_PIPELINE_TO_ASSETS` maps both `long-meta` and `long-profiling` to `masked_t2t_hla_minimap2` (`cerberus/refs.py:42-43`), which `README.md:174` states is built with `minimap2 -x sr -d` and whose manifest description confirms is the short-read index.
- **Trigger:** Any `--long` run against the published Zenodo bundle.
- **Consequence:** Autotune selects `minimap2_preset="map-ont"` for both long classes (`cerberus/autotune.py:65,71`) and `align.py:113` passes it as `-ax map-ont`. minimap2 ignores a preset's `-k`/`-w` when given a prebuilt index (it emits `[WARNING] Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index`), so ONT/PacBio reads are seeded with the `sr` index's short-read k/w. Host depletion sensitivity on long noisy reads drops substantially — which again means residual human reads in a `--long --gdpr` output that claims zero. The user sees only a minimap2 stderr warning buried in a per-stage log. The README/manifest divergence means users who follow the docs will assume the correct index is in use.
- **Fix:** Either publish `masked_t2t_hla.long.mmi` as a distinct manifest asset and add `masked_t2t_hla_minimap2_long` to the `long-*` mappings, or — cheaper and arguably better — drop the prebuilt `.mmi` for long modes and align against the masked FASTA so the preset's indexing parameters actually apply. Until then, correct `README.md:174-176`. Add a startup assertion that warns loudly when a `map-*` preset is paired with an `sr`-built index.

### F5. `--long --profiling --double-pass` on very-long reads requires a meryl DB that is in no manifest

- **Severity:** high
- **Location:** `cerberus/pipelines/long_read.py:66`, `cerberus/refs.py:43`
- **What:** The winnowmap branch derives its repetitive-k-mer DB as `meryl_db = idx.with_suffix(".meryl")`, i.e. `<ref_dir>/masked_t2t_hla.meryl`. That path is not a manifest asset, is not in `_PIPELINE_TO_ASSETS["long-profiling"]`, is never downloaded, and is never built locally. It is passed to winnowmap as `-W` (`cerberus/stages/align.py:210`).
- **Trigger:** `cerberus --long --profiling --double-pass` when autotune classifies the input as `VERY_LONG`, which sets `winnowmap_enabled=True` (`cerberus/autotune.py:71`). Both conditions are checked at `cerberus/pipelines/long_read.py:64`. This is the documented use of `--double-pass` for long reads.
- **Consequence:** winnowmap is invoked with `-W /…/refs/masked_t2t_hla.meryl` pointing at a non-existent path and fails. The resulting `ToolError` is not a `ConfigError`, and `cli.main` catches only `ConfigError` and `KeyboardInterrupt` (`cerberus/cli.py:316-324`), so the user gets a raw traceback after having already paid for QC and a 7.7 GB download. Note this fails *late*, not at the ensure step, because the meryl DB is not modelled as an asset at all.
- **Fix:** Add the meryl DB as a first-class manifest asset (`masked_t2t_hla_meryl`, a `.tar.zst`) and include it in `_PIPELINE_TO_ASSETS["long-profiling"]` / `["long-profiling-fast"]`. Alternatively, gate the winnowmap branch on `meryl_db.exists()` and fall back to minimap2 with a warning. Either way `--double-pass` must not be able to reach winnowmap without the DB being ensured.

### F6. Every invocation re-hashes the 7.7 GB minimap2 index before doing any work

- **Severity:** medium
- **Location:** `cerberus/refs.py:121-129`, `cerberus/utils/hashing.py:16-30`, called from `cerberus/orchestrator.py:75`
- **What:** `ensure()` calls `is_satisfied()` for every required asset on every run, and for non-archive assets that means a full `sha256_file()` read. There is no mtime/size/inode cache and no stamp file. Note the brief's premise that ~23 GB is re-hashed is **not** what happens: the two `.tar.zst` assets (3.79 GB + 11.05 GB) short-circuit at `cerberus/refs.py:123-124` and are never re-hashed. The recurring cost is `masked_t2t_hla.mmi` (7,697,109,690 B) plus `aux_refs.fa.gz` (331,131 B, negligible).
- **Trigger:** Every single `cerberus` invocation, including `--dry-run`. Also `cerberus doctor` and `cerberus fetch-refs` on a complete install, which additionally hash `human_k27.fa.gz` (932,696,125 B) for a total of 8.63 GB.
- **Consequence:** Benchmarked on this machine (NVMe Gen4, `sha_ni` + `avx2`, 29 GB RAM) by hashing a 1 GiB file through `cerberus.utils.hashing.sha256_file`:

  | condition | measured throughput | cost for the 7.697 GB index |
  |---|---|---|
  | pure sha256, in-RAM buffer | 2309–2501 MB/s | 3.1 s |
  | warm page cache, via `sha256_file()` | 2072–2264 MB/s | **3.5 s** |
  | cold page cache (`POSIX_FADV_DONTNEED`) | 1100 MB/s | **7.0 s** |
  | raw disk read, `dd iflag=direct` | 2.7 GB/s | (I/O is not the limit here) |

  Cold-cache is the realistic case: `README.md:233` targets 16 GB laptops and `README.md:229` notes ~9 GB peak RAM on the GDPR step, so a 7.7 GB index will not stay resident between runs. Extrapolating the cold-cache figure to the stated target hardware: a SATA SSD at ~550 MB/s gives **~14 s**, a 5400/7200 rpm or USB-3 HDD at ~120 MB/s gives **~64 s**, on every run. On a CPU without SHA-NI the hash itself becomes the limit — `sha512` on this box (not hardware-accelerated) runs at 1417 MB/s, and older mobile CPUs are considerably slower, so ~15–25 s is a fair estimate even on fast storage. This is pure overhead paid before the first read is processed, and it is most painful in exactly the iterate-on-a-small-sample workflow. It is a genuine annoyance rather than a blocker on modern hardware, and the honest framing is that it is *bounded only because of F2* — verifying the archives too, as correctness requires, would take the per-run cost to ~23 GB (~21 s cold here, minutes on an HDD).
- **Fix:** Cache verification results. Write `<ref_dir>/.verified.json` mapping asset key → `(size, mtime_ns, sha256, verified_at)`; skip re-hashing when size and mtime match. Re-verify on an explicit `cerberus doctor --deep` or when the manifest hash changes. This makes full verification affordable for archives too (F2), turning the current trade-off into a non-issue. A cheap interim step: skip the hash when `st_size != manifest.size_bytes` is already decisive, and always check size first.

### F7. Manifest versioning and the "user is prompted" behaviour are entirely unimplemented; a stale manifest is never refreshed

- **Severity:** medium
- **Location:** `cerberus/refs.py:8` (docstring claim), `cerberus/refs.py:79-89` (implementation), `cerberus/cli.py:307-311`
- **What:** The module docstring states: "If a newer version is shipped in the package, the user is prompted (unless `--update-refs`)." Nothing implements this. `_load_manifest` seeds the packaged default only when `manifest.json` is absent (`:80-82`) and otherwise loads whatever is on disk. The manifest's own `schema_version` and `release` fields are never read anywhere in `cerberus/` — grepping shows the only reader is `scripts/zenodo_publish.py:68`, which *writes* it. `--update-refs` does not update the manifest at all; `cerberus/cli.py:307-311` just calls `fetch_all()` against the **existing on-disk manifest**.
- **Trigger:** Upgrading the Cerberus package after any prior run has created `~/.cerberus/refs/manifest.json`.
- **Consequence:** The on-disk manifest is pinned forever at whatever version first ran. A new release that fixes a bad checksum, moves to a new Zenodo record, or adds an asset will be silently ignored — users keep validating against stale hashes and downloading from stale URLs, and a newly-added asset key raises `RefManagerError: Unknown asset key` (`cerberus/refs.py:94`) from deep inside a pipeline. `--update-refs` is the natural thing a user reaches for and it does not help. There is no `cerberus refs --refresh-manifest` either.
- **Fix:** On load, compare the on-disk `release`/`schema_version` against the packaged default; if the package is newer, either prompt (as documented) or auto-merge new/changed asset entries while preserving user edits, and make `--update-refs` re-seed the manifest before calling `fetch_all()`. If prompting is not wanted, delete the claim at `cerberus/refs.py:8` rather than leaving a docstring that describes software that does not exist.

### F8. `_urllib()` has no timeout, no retry, no resume, and no HTTP error handling

- **Severity:** medium
- **Location:** `cerberus/refs.py:186-204`
- **What:** `urllib.request.urlopen(url)` at `:194` is called with no `timeout=` argument, so it inherits the global default of `None` — an unbounded block. There is no retry loop, no `Range`-header resume, no `User-Agent`, and no handling of `HTTPError`/`URLError`. The downloaded length is never compared against the `Content-Length` read at `:195`.
- **Trigger:** The fallback path taken on every machine without `aria2c` (`cerberus/refs.py:153-156`) — which is the default, since `aria2c` is not in `environment.yml`'s critical path for most users. Any flaky network, Zenodo 502/503 (common under load), or captive portal.
- **Consequence:** A stalled TCP connection hangs Cerberus indefinitely with a frozen tqdm bar and no timeout — the user must Ctrl-C. A transient 503 raises a raw `urllib.error.HTTPError` that escapes `RefManagerError`, escapes `orchestrator.run`, and escapes `cli.main` (which catches only `ConfigError`/`KeyboardInterrupt`, `cerberus/cli.py:316-324`), so the user sees an unhandled traceback rather than a "could not reach Zenodo, retry with…" message. Because there is no resume, a failure at 10.9 GB of the 11.0 GB Kraken2 DB discards all of it: `cleanup_partial` deletes the `.tmp` at the start of the next run (`cerberus/orchestrator.py:64`), so the full download restarts from zero. On a slow link this can make the tool effectively uninstallable.
- **Fix:** Pass `timeout=60` to `urlopen`; wrap the whole download in a retry loop (3–5 attempts, exponential backoff) that catches `URLError`/`HTTPError`/`socket.timeout` and re-raises as `RefManagerError` with the URL and status. Implement resume by stat-ing the existing `.tmp` and sending `Range: bytes=N-`, accepting 206. Verify the written byte count against `Content-Length` before handing off to the hash check. Catch `RefManagerError` in `cli.main` and print it without a traceback.

### F9. `cleanup_partial()` misses every aria2c artifact, and `_aria2c` cannot resume or retry

- **Severity:** medium
- **Location:** `cerberus/refs.py:238-244` (glob at `:241`), `cerberus/refs.py:174-184`
- **What:** `cleanup_partial` globs only `*.tmp`. aria2c writes a control file named `<target>.aria2` — i.e. `masked_t2t_hla.mmi.tmp.aria2` — which does not match `*.tmp`. The log path is `dst.with_suffix(".aria2.log")` (`:176`), which for `dst = masked_t2t_hla.mmi.tmp` resolves to `masked_t2t_hla.mmi.aria2.log` (`with_suffix` replaces `.tmp`), also unmatched. `_extract` similarly leaves `<name>.tar.extract.log` (`cerberus/refs.py:220`).
- **Trigger:** Any interrupted or failed aria2c download.
- **Consequence:** Verified empirically — seeding six realistic leftovers, `cleanup_partial` removed 2 and left 4: `masked_t2t_hla.mmi.tmp.aria2`, `kraken2_gdpr_compact.tar.zst.tmp.aria2`, `masked_t2t_hla.mmi.aria2.log`, `kraken2_gdpr_compact.tar.extract.log`. The stale `.aria2` control files are the harmful ones: `cleanup_partial` deletes the data `.tmp` but leaves the control file describing it, so the next `aria2c` invocation finds a control file referencing a file that no longer exists and either errors or silently restarts, with the "cleaned N leftover .tmp file(s)" log (`cerberus/orchestrator.py:65-66`) implying the directory is clean when it is not. Separately, `_aria2c` passes neither `--continue=true` nor `--max-tries`/`--retry-wait`, so aria2c's main advantages (resume, retry) are unused — `--allow-overwrite=true` actively forces a restart-from-scratch. A failed aria2c raises `ToolError` (not `RefManagerError`) from `run()`, which escapes to a traceback exactly as in F8.
- **Fix:** Broaden the glob to `*.tmp`, `*.tmp.aria2`, `*.aria2`, and `*.aria2.log`, or track downloads in a dedicated `<ref_dir>/.partial/` directory that can be cleared wholesale. Add `--continue=true --max-tries=5 --retry-wait=5 --timeout=60` to the aria2c argv and drop `--allow-overwrite=true` in favour of resume. Consider `--checksum=sha-256=<hash>` so aria2c verifies inline. Wrap the `run()` call so `ToolError` becomes `RefManagerError`.

### F10. `human_kmer_set` is dead weight: 0.93 GB downloaded by `fetch-refs`, used by nothing, and reported as a fault by `doctor`

- **Severity:** medium
- **Location:** `cerberus/data/default_manifest.json:52-62`, `cerberus/refs.py:37-45`, `cerberus/refs.py:225-227`, `cerberus/refs.py:229-235`
- **What:** The manifest declares `human_kmer_set` (`human_k27.fa.gz`, 932,696,125 B) with `"required_for": ["gdpr"]` and a description promising a "belt-and-braces bbduk scrub (orthogonal to Kraken2)". But `_PIPELINE_TO_ASSETS["gdpr"]` is `["kraken2_gdpr_compact", "masked_t2t_hla_minimap2"]` — no k-mer set — and `cerberus/pipelines/gdpr.py` never calls `refs.asset("human_kmer_set")` or any bbduk stage; grepping `human_kmer_set|human_k27` across `cerberus/` finds only the manifest entry. Meanwhile `fetch_all()` iterates `self.manifest["assets"]` unconditionally (`:226`), so it **does** download it, and `doctor()` iterates the same (`:231`), so it **does** flag it.
- **Trigger:** `cerberus fetch-refs`, `cerberus --update-refs`, or `cerberus doctor`.
- **Consequence:** Every user pays 0.93 GB of download and disk for a file no code path opens. Users who let assets download lazily (the default) and then run `cerberus doctor` are told `missing or corrupt: human_kmer_set (human_k27.fa.gz)` — confirmed empirically — on a completely healthy installation, which is exactly the kind of false alarm that trains people to ignore `doctor`. More seriously, the manifest description and `README.md:74`'s "dual orthogonal mechanisms" framing imply a third defence that is not wired up; a reviewer reading the manifest would reasonably believe the bbduk k-mer scrub runs during `--gdpr`.
- **Fix:** Decide which is true. If the bbduk scrub is intended, wire it into `run_gdpr_for` and add `human_kmer_set` to the `gdpr`/`long-gdpr` mappings. If it is deferred (as `README.md:90` hints for the related UHGG masking), remove it from `default_manifest.json` or mark it `"optional": true` and have `fetch_all()`/`doctor()` skip optional assets by default.

### F11. Downloaded tarballs are never deleted after extraction; the README's disk figure is ~3× low

- **Severity:** medium
- **Location:** `cerberus/refs.py:170-172`, `README.md:230`
- **What:** `_download` renames the `.tmp` into place and calls `_extract` (`:171-172`), but nothing ever unlinks the archive. Both `.tar.zst` assets therefore persist alongside their extracted directories forever, since `is_satisfied` only looks at the directory.
- **Trigger:** Normal first-run behaviour.
- **Consequence:** Steady-state disk after `cerberus fetch-refs`: retained downloads 23.47 GB (7.697 + 3.792 + 11.049 + 0.933 + 0.0003) **plus** extracted content — the Kraken2 DB is ~14 GB extracted per `README.md:94`, and a bowtie2 index of the human genome is roughly 4–5 GB — for a total around **42 GB**. `README.md:230` states "~13 GB extracted in `~/.cerberus/refs/` (one-time)", and `README.md:233` says "Designed for 16 GB laptops". A user who provisions ~15–20 GB on the strength of the README will run out of space partway through extracting the 11 GB Kraken2 DB — which, per F2/F3, then leaves a half-extracted DB that is permanently considered valid. This machine (22 G free on `/`) could not complete `fetch-refs`. The two defects compound: the under-documented disk requirement is a likely *cause* of the silent-corruption state.
- **Fix:** Unlink the archive after a verified extraction (gated behind a `--keep-archives` flag if re-extraction without re-download is wanted — see also the related waste that deleting an extracted dir forces a full re-download rather than a re-extract, because `_download` is the only path to `_extract`). Check `shutil.disk_usage(ref_dir).free` against the sum of `size_bytes` plus an extraction multiplier before starting, and fail early with a clear message. Correct `README.md:230` to the real steady-state figure.

### F12. The manifest's `required_for` field is parsed but never used, and contradicts `_PIPELINE_TO_ASSETS`

- **Severity:** low
- **Location:** `cerberus/refs.py:56,102`, `cerberus/data/default_manifest.json`
- **What:** `Asset.required_for` is populated from the manifest (`:102`) and then read by nothing — grep for `required_for` across `cerberus/` returns only the dataclass field and its assignment. The authoritative mapping is the hardcoded `_PIPELINE_TO_ASSETS`, and the two disagree: the manifest says `masked_t2t_hla_minimap2` is required for `["meta", "profiling-fast"]`, but the code also uses it for `profiling`, `gdpr`, `long-meta`, `long-profiling` and `long-gdpr`; the manifest says `human_kmer_set` is required for `["gdpr"]`, but the code never uses it (F10). `README.md:174-177` reproduces the manifest's (wrong) view.
- **Trigger:** Reading the manifest or README to understand which assets a mode needs.
- **Consequence:** Two sources of truth that have already drifted. Anyone maintaining the asset list will update the manifest and see no effect, or will trust `required_for` when deciding what to publish. This is the same class of drift that produced F1 and F4.
- **Fix:** Delete `_PIPELINE_TO_ASSETS` and derive the mapping by inverting `required_for` at load time, so the manifest is the single source of truth — then fix the manifest's entries to match reality. Add a test asserting the inverted map covers every key `_required_pipeline_keys` can emit (see F1).

### F13. `tarfile.extractall(filter="data")` raises `TypeError` on Python 3.10.0–3.10.11 and 3.11.0–3.11.3, which `requires-python` permits

- **Severity:** low
- **Location:** `cerberus/refs.py:212`, `pyproject.toml:10`
- **What:** `pyproject.toml:10` declares `requires-python = ">=3.10"`. The `filter=` keyword on `TarFile.extractall` was added in Python 3.12 and back-ported only to 3.11.4, 3.10.12 and 3.9.17 (the CVE-2007-4559 / PEP 706 fix). On any permitted interpreter older than those patch releases, `extractall(out_dir, filter="data")` raises `TypeError: extractall() got an unexpected keyword argument 'filter'`. Verified present on the installed interpreter (3.13.9, `tarfile.data_filter` exists, `filter` in the signature), so this is latent rather than active here.
- **Trigger:** Extracting a `.tar`, `.tar.gz` or `.tar.xz` asset on Python 3.10.0–3.10.11 / 3.11.0–3.11.3. Currently unreachable with the shipped manifest, since both archives are `.tar.zst` and route to the zstd branch — so this only bites custom manifests (`README.md:165-177` documents `build_custom_host_ref.sh` producing local ref dirs) or a future `.tar.gz` asset.
- **Consequence:** A confusing `TypeError` from inside the reference manager on a supported interpreter, after a multi-GB download has already completed.
- **Fix:** Raise `requires-python` to `>=3.11.4` (or `>=3.12`), or guard: `kw = {"filter": "data"} if hasattr(tarfile, "data_filter") else {}` and pass `**kw`, logging a warning about unsanitised extraction on the fallback.

### F14. `extracted_dirname` truncates at the first dot, so any dotted stem breaks the cache

- **Severity:** low
- **Location:** `cerberus/refs.py:62-64`, mirrored at `cerberus/refs.py:208`
- **What:** `self.filename.split(".")[0]` takes everything before the *first* dot, not the archive stem. Correct for all five shipped filenames, but a natural name like `T2T-CHM13v2.0_bt2.tar.zst` yields `T2T-CHM13v2` — silently dropping part of the name.
- **Trigger:** Adding a manifest asset (or a custom-built ref per `README.md:165`) whose filename contains a dot before the extension. Version numbers in filenames are the obvious case, and the project's own descriptions use `T2T-CHM13v2.0` throughout.
- **Consequence:** Two assets whose names differ only after the first dot would collide into the same extraction directory and silently overwrite each other; a single asset merely gets an unexpected directory name. `path_to()` and `_extract` use the same expression so they stay consistent, which is why this is low rather than higher.
- **Fix:** Strip known archive suffixes explicitly, e.g. `re.sub(r"\.tar(\.(gz|zst|xz))?$", "", self.filename)`, and factor the single implementation so `path_to` and `_extract` cannot diverge.

### F15. `--kraken-db` / `--aux-refs` overrides are parsed into config and then ignored

- **Severity:** low
- **Location:** `cerberus/config.py:92-93`, `cerberus/cli.py:210-211`
- **What:** `kraken2_db_override` and `aux_refs_override` are defined on `CerberusConfig` and populated from the CLI, but no pipeline reads them — grep finds only the definition and the assignment. `run_gdpr_for` unconditionally resolves the Kraken2 DB through `refs.asset("kraken2_gdpr_compact")` (`cerberus/pipelines/gdpr.py:49`), and `run_profiling` / `run_long_profiling` unconditionally resolve `refs.asset("aux_refs")` (`cerberus/pipelines/profiling.py:90,123`, `cerberus/pipelines/long_read.py:82`).
- **Trigger:** Passing `--kraken-db` or `--aux-refs`.
- **Consequence:** The flags are accepted without complaint and have no effect; Cerberus downloads the 11 GB bundled DB anyway. This matters for the custom-reference workflow the README promotes at `README.md:189-193`, which explicitly shows `--kraken-db` being used to reuse a published DB alongside a custom host reference — that documented invocation does not work.
- **Fix:** Honour the overrides at the point of use, e.g. `kdb_dir = cfg.kraken2_db_override or refs.path_to(refs.asset("kraken2_gdpr_compact"))`, and drop the corresponding asset from the required list when an override is supplied so the download is skipped. Validate the override path exists at config time.
