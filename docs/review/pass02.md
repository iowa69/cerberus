# Pass 2 — Subprocess and I/O layer

## Summary

The subprocess layer is safe against shell injection (nothing anywhere uses `shell=True`) and the single-process `run()` path is solid, but the multi-process `pipe()` path has a defect that undermines the scientific validity of the whole pipeline: it checks only the **last** process's return code, so a failure in any earlier stage is silently discarded. Every host-depletion alignment in Cerberus is `minimap2 | samtools view`, and `samtools view` happily exits 0 on a truncated stream — I demonstrated an OOM-killed aligner (rc 137) producing a 70-byte truncated BAM with `pipe()` returning 0 and the pipeline continuing. A second proven defect is `concat.compress_to()`, where `pigz -c` writes its gzip stream to **the log file** instead of the destination (I recovered the full 140 KB payload out of `test.pigz.log`, and the destination was never created); it is currently unreachable dead code, which is the only reason it has not caused a data loss incident. `count_reads()` compounds the problem downstream by detecting gzip from the filename suffix alone and ignoring the decompressor's exit status, so truncated and mislabelled files silently produce undercounts or plain zero in the user-facing accounting report.

## Verified working

- **No shell interpolation anywhere** — `grep -rn "shell=True" --include=*.py` over the repo returns zero hits. `run()` accepts a `str` only to `shlex.split` it into an argv list (`cerberus/utils/shell.py:57-60`) and always execs a list; `pipe()` only ever accepts `list[list[str]]` (`cerberus/utils/shell.py:90`). User-supplied `--minimap2-args`/`--bowtie2-args` cannot inject shell metacharacters.
- **`run()` log-header ordering is correct** — the `# CMD:` line is written *and flushed* before the child is spawned (`cerberus/utils/shell.py:72-73`), so the header cannot land after the child's output. Verified: `run(["sh","-c","echo CHILD-STDOUT; echo CHILD-STDERR >&2"])` produced exactly `"# CMD: ...\nCHILD-STDOUT\nCHILD-STDERR\n"`.
- **`run()` return-code checking is correct** — `check=True` by default and the non-zero path raises `ToolError` with a `shlex.quote`d command plus the log path (`cerberus/utils/shell.py:84-85`, `cerberus/utils/shell.py:20-23`).
- **`pipe()` does not leak file descriptors on the success path** — each intermediate read-end is closed immediately after being handed to the next child (`cerberus/utils/shell.py:128-129`), and the last stage's `p.stdout` is always `None` because its stdout is `DEVNULL` or a file object (`cerberus/utils/shell.py:117-120`). Measured: 200 consecutive `pipe()` calls with `final_stdout` left the process fd count unchanged at 4.
- **`final_stdout` data is fully on disk when `pipe()` returns** — the parent's file object is only used to hand a fd to the child; the parent never buffers payload bytes, and `p.wait()` completes before the return. Verified by reading `'partial-data\n'` out of `final_stdout` immediately after `pipe()` returned.
- **SIGPIPE behaviour is normal** — `Popen`'s default `restore_signals=True` resets SIGPIPE to `SIG_DFL` in children, so producers die rather than hang when a consumer exits early. Verified: `pipe([["yes"],["head","-c","10"]])` terminated promptly and wrote exactly 10 bytes.
- **Last-stage failures *are* detected** — `pipe([["true"],["false"]])` raised `ToolError rc=1`, and a fully corrupt `.tar.zst` raised `ToolError rc=2` because `tar` itself rejected the input (`cerberus/utils/shell.py:137-139`).
- **No deadlock risk in the sequential `p.wait()` loop** — stderr for every stage goes to a real file, never a pipe (`cerberus/utils/shell.py:125`), so no child can block writing to an unread stderr pipe while the parent waits on an earlier process (`cerberus/utils/shell.py:132-133`).
- **`require_tools()` / `which()` are correct and actionable** — missing binaries are collected and reported together with an install hint (`cerberus/utils/shell.py:33-44`).
- **`hashing.py` is clean** — `sha256_file` streams in 1 MiB chunks and never loads a reference into memory, with the tqdm bar closed in a `finally` (`cerberus/utils/hashing.py:21-29`); `verify_sha256` compares case-insensitively (`cerberus/utils/hashing.py:35`). Critically, `refs.py` hashes the `.tmp` file and only then renames (`cerberus/refs.py:160-170`), so a verification failure can never leave a corrupt file at the final path.
- **`concat_gz` streaming is sound** — 4 MiB chunked copy with both handles context-managed (`cerberus/stages/concat.py:37-41`), and `cerberus/pipelines/profiling.py:139` unlinks the destination first so a re-run cannot append to a stale output.
- **`dry_run` never execs** — both `run()` and `pipe()` short-circuit before any `Popen` and write a marker log (`cerberus/utils/shell.py:66-69`, `cerberus/utils/shell.py:103-105`).

## Findings

### F1. `pipe()` checks only the final process's return code; every upstream failure is silently swallowed
- **Severity:** critical
- **Location:** `cerberus/utils/shell.py:132-139`
- **What:** The wait loop `for p in procs: p.wait()` discards each process's status, and only `procs[-1].returncode` is tested. Any non-zero exit — or fatal signal — from stages `0 .. n-2` is invisible to the caller. There is no `pipefail` equivalent.
- **Trigger:** Reproduced directly against the real function:
  - `pipe([["false"],["true"]])` → **no exception, `returncode=0`**
  - `pipe([["false"],["cat"]])` → **no exception, `returncode=0`**
  - `pipe([["sh","-c","echo x"],["sh","-c","exit 7"],["cat"]])` (middle stage fails) → **no exception, `returncode=0`**
  - `pipe([["sh","-c","echo partial-data; exit 3"],["cat"]], final_stdout=out)` → rc 0 and `out` contains the truncated `'partial-data\n'`.

  Realistic call-site reproduction using shims on `PATH` for the `cerberus/stages/align.py:64` pattern: a `minimap2` that emits a valid `@HD`/`@SQ` header plus one record, prints `[ERROR] failed to load index / out of memory` to stderr, then exits **137** (OOM-kill), piped into a `samtools view -b -o aln.bam -` that consumes stdin and exits 0. Result: `pipe() returned rc=0`, a **70-byte truncated BAM was written**, and the real error sat unread in `<tag>.minimap2.log`.

  Same class at `cerberus/stages/qc.py:117-128`: `pipe([["zcat", missing_file], ["chopper", ...], ["pigz","-c"]], final_stdout=qc.long.fq.gz)`. `zcat` failed with `No such file or directory`; `pigz -c` still exited 0 and emitted a **valid 20-byte empty gzip member**; `pipe()` returned rc 0; `count_reads()` on the result returned **0**.

  Also at `cerberus/refs.py:217-221`: a `.tar.zst` with valid content followed by trailing garbage makes `zstd -dc` exit 1 (`unsupported format`) *after* `tar` has already seen the end-of-archive blocks and exited 0 → `pipe()` rc 0, reference directory accepted as good. (`cerberus/refs.py:166-168` only *warns* when the manifest carries no SHA256, so this is the unguarded path.)
- **Consequence:** Silent scientific corruption. A host-depletion alignment that dies from OOM, a missing index, a full disk, or a bad `--minimap2-args` token produces a truncated or empty BAM; `samtools fastq` then emits a tiny/empty unmapped FASTQ; the run completes "successfully" and `accounting.json` reports a plausible-looking but wrong depletion rate. The long-read chopper path can silently deliver a zero-read output file that is a structurally valid gzip. Nothing in the run surfaces the failure — the evidence exists only in a log file no code reads.
- **Fix:** Collect all statuses and fail on the first non-zero, while tolerating SIGPIPE on non-final stages (which is legitimate when a consumer stops early, e.g. `tar` finishing before `zstd`):
  ```python
  rcs = [p.wait() for p in procs]
  for i, (rc, cmd) in enumerate(zip(rcs, cmds)):
      is_last = i == len(cmds) - 1
      if rc == 0:
          continue
      if not is_last and rc == -signal.SIGPIPE:   # consumer closed early: benign
          continue
      raise ToolError(cmd, rc, log_path)
  ```
  Additionally, `cerberus/stages/qc.py:118` should not let `pigz -c` mask an upstream failure — verify `out_reads` is non-empty after the call.

### F2. `compress_to()` writes pigz's compressed output into the log file; the destination is never created
- **Severity:** high
- **Location:** `cerberus/stages/concat.py:52-53`
- **What:** `pigz -c` writes the gzip stream to **stdout**, and `run()` redirects the child's stdout into `log_path` (`cerberus/utils/shell.py:75-77`). `dst` appears nowhere in the argv, so it is never written. The function then `return dst` unconditionally (`cerberus/stages/concat.py:59`), reporting success for a file that does not exist. The no-pigz fallback branch (`cerberus/stages/concat.py:55-58`) *does* write `dst` correctly, so behaviour silently diverges based on whether `pigz` happens to be installed. That branch also ignores `cfg.dry_run` and writes for real during a dry run.
- **Trigger:** Ran it for real with a 140,000-byte `reads.fq` source and `CerberusConfig(threads=2)`:
  ```
  returned:    /.../d2/reads.fq.gz
  dst exists?  False
  log file test.pigz.log: 492 bytes
    first 80 bytes: b'# CMD: pigz -c -p 2 /.../reads.fq\n\x1f\x8b\x08\x08...'
  Decompressed payload from the LOG FILE: 140000 bytes; matches src? True
  ```
  The entire compressed payload landed in `<tag>.pigz.log`, prefixed by the `# CMD:` text header — so the "log" is not readable text *and* not a valid gzip file either.
- **Consequence:** Any caller silently loses the data it thought it compressed and gets back a path to a non-existent file, while the log directory fills with binary garbage. **Reachability: `compress_to` has no callers** — `grep -rn "compress_to"` across the repo matches only its definition and its own `ValueError` string. It is dead code today; wiring it into any pipeline would immediately cause silent data loss.
- **Fix:** Either delete the function, or write to the destination properly and honour `dry_run` in both branches:
  ```python
  run(["pigz", "-c", "-p", str(cfg.threads), str(src)],
      log_path=dst, dry_run=cfg.dry_run)   # still wrong: log header pollutes output
  ```
  is not sufficient because `run()` always writes the `# CMD:` header. Use `pigz -p N -c src > dst` via an explicit `subprocess.run(stdout=dst_fh, stderr=log_fh)`, or simply `["pigz", "-f", "-p", str(cfg.threads), "-k", str(src)]` and rename — and add a `stdout_path` parameter to `run()` so binary payloads never share a handle with the log.

### F3. `_line_count_pipe()` ignores the decompressor's exit status — truncated and corrupt files yield silent wrong counts
- **Severity:** high
- **Location:** `cerberus/utils/fastq.py:29-36`
- **What:** `p1.wait()` is called but its return code is discarded, and `p2`'s status from `communicate()` is discarded too. Whatever `wc -l` saw is returned as fact.
- **Trigger:** Ran against real `pigz`:
  - Truncated gzip (first 50% of a valid 10,000-read `.fq.gz`, i.e. an interrupted write or partial download). `pigz -dc` alone: `rc=1`, `corrupted -- incomplete deflate data`, 228,023 bytes recovered. **`count_reads()` returned 4981** for a 10,000-read file, with no error and no warning.
  - A `.fq.gz` that is not gzip at all (e.g. an HTML 404 page saved under a `.gz` name). `pigz -dc` alone: `rc=1`, `unrecognized format`. **`count_reads()` returned 0**, no exception.
  - Direct: `_line_count_pipe(["sh","-c","echo a; echo b; exit 9"])` returned `2` with no error raised.
- **Consequence:** `accounting.json` / `accounting.tsv` — the user-facing read-accounting deliverable produced from these counts (`cerberus/accounting.py:44-52`, `cerberus/orchestrator.py:129-135`) — reports fabricated numbers. A half-written intermediate reads as a ~50% depletion rate; a corrupt file reads as 100% depletion. Both are indistinguishable from a real biological result.
- **Fix:** Check both statuses and raise:
  ```python
  out, _ = p2.communicate()
  rc1, rc2 = p1.wait(), p2.returncode
  if rc1 not in (0, -signal.SIGPIPE) or rc2 != 0:
      raise ToolError(cmd, rc1 or rc2, None)
  ```

### F4. `count_reads()` detects gzip by filename suffix only, ignoring the magic-byte sniffer that already exists
- **Severity:** high
- **Location:** `cerberus/utils/fastq.py:18`
- **What:** `is_gz = path.suffix == ".gz"` drives the entire dispatch (`cerberus/utils/fastq.py:19-26`). The module already contains `is_gzipped()` (`cerberus/utils/fastq.py:53-60`), which correctly falls back to reading the `\x1f\x8b` magic bytes — and it is never called from `count_reads`. Suffix matching also misses `.bgz`, `.fq.bz2`, and `.fastq.zst`.
- **Trigger:** Measured on 1000-read files:
  - Gzip *content* under a plain name (`b.fastq`, produced by many upstream tools and by `cp x.fq.gz x.fastq`): `is_gzipped()` correctly returns `True`, but `count_reads()` takes the `wc -l` branch (`cerberus/utils/fastq.py:23-25`) and counts `\n` **bytes inside the compressed stream** → **returned 3** instead of 1000.
  - Plain FASTQ text under a `.gz` name (`c.fq.gz` — common when a tool is told to write `.gz` but does not compress): `pigz -dc` fails with `unrecognized format`, and combined with F3 → **returned 0** instead of 1000.
- **Consequence:** Read counts off by three orders of magnitude, or zero, with no diagnostic. Because `count_reads` is the sole source of every number in the accounting report, a single mislabelled input poisons the entire summary.
- **Fix:** Replace line 18 with `is_gz = is_gzipped(path)` and add an explicit `raise ValueError` when the sniffed format is one the function cannot handle, rather than falling through to `wc -l`.

### F5. `pipe()` leaks a file descriptor and orphans running child processes when a later `Popen` raises
- **Severity:** high
- **Location:** `cerberus/utils/shell.py:112-135`
- **What:** The `finally` block closes only `stderr_log`. If constructing stage *i* fails (missing binary → `FileNotFoundError`, `PermissionError`, `OSError`), stages `0 .. i-1` are already running and are never terminated or reaped, and `prev_stdout` (the parent's read end of stage `i-1`'s pipe) is never closed — so the still-running producer never receives SIGPIPE and keeps running until it blocks on a full pipe buffer, forever. The `final_stdout` file object opened at `cerberus/utils/shell.py:118` is likewise never explicitly closed (CPython refcounting reclaims it, but that is incidental, not guaranteed).
- **Trigger:** `pipe([["sh","-c","sleep 30; echo done"], ["/nonexistent/binary"]])`:
  ```
  raised: FileNotFoundError [Errno 2] ... '/nonexistent/binary'
  fds before=4 after=5 (leak=1)
  orphan stage-0 still running?: 2927965 sh -c sleep 30; echo done
                                 2927967 sleep 30
  ```
  The orphan survived the exception and had to be killed manually.
- **Consequence:** A typo'd tool name or a partially installed conda environment leaves a `minimap2` process pinned to N threads and tens of GB of index in RAM after Cerberus has already exited with a traceback — on a shared HPC node this holds the allocation hostage. Repeated failures leak one fd each; in a long-lived embedding process this eventually hits `EMFILE`.
- **Fix:** Wrap the launch loop so failures clean up:
  ```python
  except BaseException:
      if prev_stdout is not None:
          prev_stdout.close()
      for p in procs:
          p.kill()
      for p in procs:
          p.wait()
      raise
  finally:
      stderr_log.close()
      if isinstance(stdout, io.IOBase):
          stdout.close()
  ```

### F6. `pipe()` has no timeout, so a stalled stage hangs the run indefinitely
- **Severity:** medium
- **Location:** `cerberus/utils/shell.py:90-97`, `cerberus/utils/shell.py:132-133`
- **What:** `run()` exposes and honours a `timeout` parameter (`cerberus/utils/shell.py:54`, `cerberus/utils/shell.py:80`), but `pipe()` has no such parameter at all — verified: `signature(pipe) == (cmds, *, log_path, final_stdout=None, cwd=None, dry_run=False)`. `p.wait()` is unbounded.
- **Trigger:** Any stage that blocks — a reference index on a stalled NFS mount, an aligner stuck on a malformed FASTQ, a tool waiting on a tty. All four aligner invocations (`cerberus/stages/align.py:64`, `:121`, `:169`, `:217`) and both decompression pipes (`cerberus/stages/qc.py:118`, `cerberus/refs.py:218`) go through `pipe()`.
- **Consequence:** The whole run hangs with no diagnostic and no way to bound wall-clock time; on a scheduler the job is eventually killed at the walltime limit with no partial results and no error message.
- **Fix:** Add `timeout: int | None = None` to `pipe()`, use `deadline = time.monotonic() + timeout` and `p.wait(timeout=max(0, deadline - time.monotonic()))`, and on `TimeoutExpired` kill/reap every process before raising.

### F7. `setup_logging()` destroys the host application's logging configuration and forces the root logger to DEBUG
- **Severity:** medium
- **Location:** `cerberus/utils/logger.py:52-54`
- **What:** `root.handlers.clear()` removes every handler on the root logger, including ones Cerberus did not install, without calling `.close()` on them. `root.setLevel(logging.DEBUG)` then globally lowers the threshold for *all* libraries in the process. Cerberus ships a public package (`cerberus.orchestrator.run` is importable), so this is not private CLI behaviour — and it is invoked from three places (`cerberus/cli.py:229`, `cerberus/cli.py:310`, `cerberus/orchestrator.py:61`).
- **Trigger:** Ran a host-application simulation:
  ```
  host captured:                            'HOST: before cerberus\n'
  host captured after setup_logging():      'HOST: before cerberus\n'   <- handler destroyed
  root level now: DEBUG
  jsonl contains third-party debug line?:   True
  ```
  A `logging.getLogger("urllib3.connectionpool").debug("GET /secret?token=...")` record was written into `cerberus.log.jsonl` purely because Cerberus forced root to DEBUG.
- **Consequence:** Any workflow manager, notebook, or service that embeds Cerberus loses its own logging the moment a run starts. The forced DEBUG level also means third-party libraries (`urllib3`, `botocore`, `asyncio`) dump internal traffic into the JSONL — noise at best, and at worst request URLs containing download tokens are persisted to disk in the reference-cache log directory.
- **Fix:** Attach handlers to the `cerberus` logger, not root, and stop clearing other people's handlers:
  ```python
  root = logging.getLogger("cerberus")
  for h in list(root.handlers):
      h.close(); root.removeHandler(h)
  root.setLevel(logging.DEBUG)
  root.propagate = False
  ```
  Also add `logging.getLogger("cerberus").addHandler(logging.NullHandler())` at package import so library consumers get silence by default.

### F8. `ToolError` from `pipe()` always blames the last command, even when an earlier one failed
- **Severity:** medium
- **Location:** `cerberus/utils/shell.py:139`
- **What:** `raise ToolError(cmds[-1], rc, log_path)` hard-codes the final command, so the exception message names the wrong tool whenever the failure originates upstream (once F1 is fixed and upstream failures actually raise, this becomes the primary source of misdirection).
- **Trigger:** `pipe([["false"],["false"]])` → `ToolError.cmd == ['false']` — indistinguishable from the case where only the last stage failed. With the `zstd`/`tar` pipe, a `zstd` decompression failure is reported to the user as *"Command failed (rc=2): tar -x -C ..."*.
- **Consequence:** Operators debug the wrong tool. In the reference-extraction path this sends people looking at `tar`/disk permissions when the actual problem is a corrupt download.
- **Fix:** Track the failing index and pass `cmds[i]` (see the F1 fix snippet), and include the full pipeline string in the message.

### F9. `--minimap2-args` / `--bowtie2-args` are unvalidated and spliced into argv positions where they can silently redirect or reinterpret I/O
- **Severity:** medium
- **Location:** `cerberus/stages/align.py:51`, `cerberus/stages/align.py:58`, `cerberus/stages/align.py:156`, `cerberus/stages/align.py:166`
- **What:** There is **no shell-injection risk** (confirmed: zero `shell=True` hits repo-wide, everything is an argv list). The risk is argument-level. Three concrete problems: (a) `shlex.split()` raises an uncaught `ValueError` on unbalanced quotes; (b) for minimap2 the user tokens are inserted *before* the positional index/reads arguments (`cerberus/stages/align.py:58`), so an option that takes a value can consume the index path, and `-o FILE` diverts the SAM away from the pipe entirely; (c) for bowtie2 the tokens are appended *after* `-1`/`-2` (`cerberus/stages/align.py:166`), so a stray bare token is parsed as an unpaired-reads file. Users can also override the pipeline's own invariants (`--secondary=no` at `cerberus/stages/align.py:57` is defeated by `--secondary=yes`).
- **Trigger:** Verified `shlex.split` behaviour:
  ```
  '-k 21 -w 11'      -> ['-k', '21', '-w', '11']            ok
  "-x map-ont'"      -> ValueError: No closing quotation     UNCAUGHT -> traceback
  '-o /tmp/pwned.sam'-> ['-o', '/tmp/pwned.sam']
  ```
  With `--minimap2-args "-o /tmp/out.sam"` the constructed argv is
  `minimap2 -ax sr -t 8 --secondary=no -o /tmp/out.sam INDEX R1.fq R2.fq` — minimap2 writes to the file, the pipe carries nothing, `samtools view` writes an empty BAM, and (per F1) `pipe()` returns 0.
- **Consequence:** An unbalanced quote crashes with a raw traceback instead of a usage error. A well-meaning `-o` or a stray token produces an empty/garbage alignment that F1 then reports as success — the two defects compose into completely silent wrong output.
- **Fix:** Wrap the `shlex.split` calls in `try/except ValueError` and re-raise as a CLI usage error; reject a denylist of tokens that conflict with the pipeline's contract (`-o`, `-a`, `--secondary`, `-t`/`-p`, `-1`, `-2`, `-x` for bowtie2, `-S`); and for minimap2 place `extra_tokens` immediately after the fixed flags but assert none of them are value-taking options left dangling before the positionals.

### F10. Read accounting re-decompresses the same files up to three times each — 14 full passes per run
- **Severity:** medium
- **Location:** `cerberus/orchestrator.py:129-135`, `cerberus/pipelines/meta.py:110-114`, `cerberus/pipelines/profiling.py:151-152`, `cerberus/accounting.py:44-52`
- **What:** `count_reads()` is a full streaming decompression with no caching or memoisation. For a short-read run with `--meta --profiling` the same paths are counted repeatedly: `qc.r1` is decompressed at `cerberus/orchestrator.py:131`, again at `cerberus/pipelines/meta.py:110`, and again at `cerberus/pipelines/profiling.py:151` — **three full passes over the same file**. Each pipeline's final output is counted once in `_collect_counts` and again in `accounting.add_final` (`cerberus/accounting.py:51`) — **two passes each**. Total: 14 `count_reads()` invocations, 14 full decompressions.
- **Trigger:** Measured on a synthetic 600,000-read, 150 bp FASTQ (194 MB uncompressed / 94 MB gzipped): one `pigz -dc | wc -l` pass took 0.55 s ≈ **353 MB/s of uncompressed throughput**.
- **Consequence:** Scaled to a realistic 50 GB uncompressed WGS FASTQ set, 14 passes cost **~33 minutes of pure re-decompression** on top of the actual analysis, and saturate the I/O path for the whole duration. On a shared filesystem this is the difference between a 2-hour and a 2.5-hour job for numbers that are already known.
- **Fix:** Memoise on `(resolved_path, st_size, st_mtime_ns)`:
  ```python
  @functools.lru_cache(maxsize=256)
  def _count_cached(key): ...
  def count_reads(path):
      st = path.stat()
      return _count_cached((str(path.resolve()), st.st_size, st.st_mtime_ns))
  ```
  Better still, take counts from the tools that already computed them (fastp's JSON `filtering_result.passed_filter_reads`, `samtools flagstat`) instead of recounting.

### F11. `count_reads()` floor-divides by 4, silently absorbing truncated records and missing trailing newlines
- **Severity:** low
- **Location:** `cerberus/utils/fastq.py:20`, `cerberus/utils/fastq.py:22`, `cerberus/utils/fastq.py:25`, `cerberus/utils/fastq.py:45`
- **What:** Every branch returns `lines // 4` with no check that the line count is a multiple of 4, and `wc -l` counts newline *characters*, not lines.
- **Trigger:** Measured:
  - 100 complete records plus a 3-line truncated record → `count_reads()` returned **100**, silently discarding evidence that the file is malformed.
  - 100 complete records with the final newline stripped (very common for files produced by `head`/manual truncation) → `wc -l` reports 399 → **`count_reads()` returned 99**, losing a whole record.
- **Consequence:** Truncated FASTQ files pass through the pipeline undetected; downstream paired-end sync checks (`detect_paired_orphans`, `cerberus/utils/fastq.py:48-50`) compare two silently-rounded numbers and can report "no orphans" for a genuinely desynchronised pair.
- **Fix:** Compute `lines` first and `log.warning` (or raise) when `lines % 4 != 0`; for the `wc -l` branch, add 1 when the file's last byte is not `\n`.

### F12. `run()` prepends a `# CMD:` line to `log_path`, which for `samtools flagstat` is a user-facing report file
- **Severity:** low
- **Location:** `cerberus/utils/shell.py:72`, used as a report sink at `cerberus/stages/align.py:81-82`, `cerberus/stages/align.py:128`, `cerberus/stages/align.py:182`, `cerberus/stages/align.py:224`
- **What:** `run()` unconditionally writes `# CMD: <argv>\n` into the file it is given. Four call sites pass `log_path=stats` where `stats` is `<tag>.flagstat.txt` — an output documented as a user deliverable at `README.md:252`.
- **Trigger:** Verified: `run(["sh","-c","printf '100 + 0 in total\\n42 + 0 mapped\\n'"], log_path=stats)` produced
  `"# CMD: sh -c 'printf ...'\n100 + 0 in total\n42 + 0 mapped\n"`.
- **Consequence:** Cosmetic today (nothing in the repo parses `.flagstat.txt` — `AlignOutputs.stats` is only stored, never read), but any standard flagstat parser (MultiQC, `samtools`-aware scripts) will choke on or misread the first line. In `--dry-run` the same files contain only `# DRY-RUN`.
- **Fix:** Add a `header: bool = True` parameter to `run()` and pass `header=False` for report sinks, or give `run()` a separate `stdout_path` distinct from `log_path` and route flagstat there.

### F13. `pipe(dry_run=True)` does not create `final_stdout`, and the chopper branch writes a report file during dry runs
- **Severity:** low
- **Location:** `cerberus/utils/shell.py:103-105`, `cerberus/stages/qc.py:126-129`
- **What:** The dry-run short-circuit writes the marker log and returns, but never touches `final_stdout`. At `cerberus/stages/qc.py:129`, `json_report.write_text(...)` is executed unconditionally — it is outside any `cfg.dry_run` guard — so a dry run creates a real `fastplong.json` on disk describing an output that does not exist.
- **Trigger:** Verified: `pipe([["echo","x"],["cat"]], final_stdout=out, dry_run=True)` → `rc=0`, `out.exists() == False`. The returned `FastplongOutputs.reads` then points at a missing path.
- **Consequence:** Dry runs are not side-effect-free and do not exercise the same code shape as real runs; a dry run that "passes" gives false confidence about the chopper fallback path.
- **Fix:** In dry-run mode, `final_stdout.parent.mkdir(parents=True, exist_ok=True)` and `final_stdout.touch()` so downstream `exists()` checks behave consistently; move `cerberus/stages/qc.py:129` inside a `if not cfg.dry_run:` guard.

### F14. Log timestamps are taken at format time, not event time, and the two handlers disagree about timezone
- **Severity:** low
- **Location:** `cerberus/utils/logger.py:14`, `cerberus/utils/logger.py:39`
- **What:** Both formatters call `datetime.now(...)` inside `format()` and ignore `record.created`. The JSONL formatter uses UTC while the console formatter uses naive local time with no offset marker, so the same event is stamped with two different times and neither is authoritative.
- **Trigger:** Created a `LogRecord`, slept 1.2 s, then formatted it:
  ```
  record.created : 19:24:26 (UTC)
  JSONL 'ts'     : 2026-08-03T19:24:27.582742+00:00   <- 1.2s late
  console line   : [21:24:27] INFO    x  hello        <- local time, no tz
  ```
- **Consequence:** JSONL timestamps drift under any formatting delay (large messages, slow disk, handler contention), so they cannot be used to time stages. Correlating the human-readable console transcript against `cerberus.log.jsonl` requires knowing the machine's UTC offset, which is recorded nowhere — a real problem for post-hoc debugging of runs on remote/HPC machines in a different timezone.
- **Fix:** Use the record's own timestamp in both formatters: `datetime.fromtimestamp(record.created, timezone.utc)` for the JSONL, and `.astimezone()` for the console so it prints local time *with* an offset.
