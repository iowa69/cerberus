# Pass 1 — Read-filtering semantics

## Summary

The core filtering logic is a thin wrapper around `samtools fastq` flag filters, and the flag
arithmetic does not do what the code says it does. `_filter_flags()` returns `-f 12` for
`drop_strategy="both"` and `-f 4` for `"either"`, but because every paired `samtools fastq`
invocation also passes `-s /dev/null`, the two strategies emit **exactly the same reads**: a pair
survives only when *neither* mate aligns. The `drop_strategy` parameter is therefore inert, and the
`meta` pipeline — advertised in the README and in its own docstring as the conservative,
"drop only if BOTH mates map" head — is silently as aggressive as the GDPR head. Pair
synchronisation itself is sound on every path I traced (samtools' `-1/-2` + `-s` routing, BBDuk's
`removeifeitherbad=t`, and kraken2's joint `--paired` classification all keep R1/R2 in lockstep),
and secondary/supplementary alignments cannot duplicate reads because they are mapped records that
the `-f 4`/`-f 12` filter discards before samtools' name-grouping runs. The most dangerous defects
are silent-empty-output modes: a non-`sr` minimap2 preset on paired input produces zero surviving
reads with exit code 0, and `pipe()` only checks the last process in the chain, so an OOM-killed
aligner yields a truncated dataset that the pipeline reports as success. Nothing in `tests/` covers
any of these four stage modules.

## Verified working

- **Pair synchronisation on every paired alignment path** — `samtools fastq` is always given both
  `-1/-2` *and* `-s /dev/null` (`cerberus/stages/align.py:72-76` and `align.py:176-178`), so a read
  whose mate was filtered out is routed to the singleton stream and discarded rather than being
  appended to R1 or R2. I enumerated all 10 relevant flag combinations
  (`/home/iowa/Desktop/cerberus/work/p1/flag_math.py`): for both `-f 12` and `-f 4`, records are
  written to `-1`/`-2` only when both mates of a name group pass, so R1 and R2 always have equal
  record counts in the same order.
- **`-f 12` on the bowtie2 path is the canonical idiom and matches its comment** —
  `align.py:175` (`"-f", "12",  # both unmapped`) retains only pairs where flags 0x4 and 0x8 are
  both set, i.e. neither mate aligned to the masked T2T+HLA reference. This is the standard
  host-removal recipe and is correct for the aggressive profiling head.
- **Secondary/supplementary alignments cannot duplicate reads** — 0x100 and 0x800 records are by
  definition *mapped* (0x4 clear), so the `-f 4` / `-f 12` filter removes them before
  `samtools fastq` groups records by name; my enumeration confirms `keep(0x1|0x800|0x40, f=4)` and
  `keep(0x1|0x100|0x40, f=4)` are both False. `--secondary=no` (`align.py:56`, `align.py:116`)
  additionally suppresses minimap2 secondaries. No `-F 0x900` is needed for the *unmapped*
  selection (see F7 for why its absence is still a latent hazard).
- **No spurious reverse-complementing** — `samtools fastq` reverse-complements records carrying
  0x10, but every surviving record here has 0x4 set and unmapped records never carry 0x10, so
  sequences in `*.unmapped.R*.fq.gz` are in original read orientation.
- **BBDuk paired invocations are pair-preserving** — both `kmer.py:65-66` (`in1/in2` → `out1/out2`)
  and `entropy.py:41-42` use BBDuk's paired mode, where the default `removeifeitherbad=t` discards
  both mates together. R1/R2 cannot desynchronise here (BBDuk's default `ordered=f` may reorder
  *pairs* relative to the input, but each mate pair is written as a unit).
- **kraken2 `--paired` keeps mates together** — with `--paired` (`kraken.py:56`) kraken2 classifies
  the concatenated mate pair as a single unit and writes both mates or neither to the
  `--unclassified-out` template (`kraken.py:48,58`), so the R1/R2 handed to the second GDPR
  mechanism are in sync.
- **The `#` template and the discovery glob agree** — kraken2 substitutes `#` → `_1`/`_2`, and
  `_gzip_pair_outputs` globs `f"{name_root}*1.fq"` / `*2.fq` (`kraken.py:136-139`), which matches
  `<tag>.kraken_unclass_1.fq` / `_2.fq`. The `kraken_class` root cannot cross-match the
  `kraken_unclass` root, so `keep_classified=True` does not confuse the two.
- **Empty-gzip emission is valid** — `_write_empty_gzip` (`kraken.py:175-179`) and
  `_empty_paired_passthrough` (`kmer.py:30-33`) write real gzip members (25–53 bytes depending on
  the embedded filename), which decompress cleanly; `count_reads` returns 0 for them
  (verified in `/home/iowa/Desktop/cerberus/work/p1/empty_probe.py`).
- **Single-end paths use the right output slot** — `align.py:125` and `align.py:221` pass
  `-f 4 -0 <out>`; records with neither 0x40 nor 0x80 are designated READ_OTHER, which is exactly
  what `-0` receives, so long-read/orphan streams are written in full.
- **Unknown strategies fail loudly** — `_filter_flags` raises `ValueError` on an unrecognised
  strategy (`align.py:242`) rather than silently defaulting.
- **Concatenation is stream-legal** — `concat_gz` (`concat.py:37-41`) byte-concatenates gzip
  members, which is a valid multi-member gzip stream readable by minimap2 (zlib `gzread`), BBDuk
  and kraken2. The merged profiling FASTQ and the orphan merges are therefore well-formed.
- **Entropy thresholds are scaled by read class and the invocation is a legal ref-free BBDuk call** —
  0.60/0.70/0.65/0.50/0.45 across the length buckets (`autotune.py:43-73`) fed to
  `entropy={...}` (`entropy.py:43`); BBDuk accepts entropy filtering with no `ref=`, so this stage
  runs as a pure filter.

## Findings

### F1. `drop_strategy` is inert — `"both"` and `"either"` select identical read sets, and `meta` is not conservative

- **Severity:** critical
- **Location:** `cerberus/stages/align.py:230-242` (with `align.py:34,38-41`, `cerberus/pipelines/meta.py:46`, `README.md:9,61`)
- **What:** `_filter_flags("both")` returns `-f 12`; `_filter_flags("either")` returns `-f 4`.
  Because the caller always adds `-s /dev/null` (`align.py:75`), a read that passes `-f 4` but whose
  mate does not becomes a *singleton* and is written to `/dev/null`. The pair-level outcome table
  (`work/p1/flag_math.py`) is identical for the two strategies:
  both-unmapped → `R1+R2`; one-mate-mapped → dropped; both-mapped → dropped. So both strategies
  implement "drop the pair if **either** mate maps". The documented semantics of `"both"`
  ("only drop reads where both mates map", `align.py:39`; "drop pairs where BOTH mates map",
  `meta.py:5`; "Conservative — retains microbial reads even at the cost of some residual host",
  `README.md:9`) is implemented nowhere in the codebase.
- **Trigger:** any `--meta` run. Concretely: a pair where R1 is a genuine microbial fragment and R2
  overlaps a segmental-duplication or HLA region of the masked T2T reference well enough to align.
- **Consequence:** the whole pair is deleted from `<sample>.meta.R1.fastq.gz`/`R2`. For a
  host-heavy sample (biopsy, saliva) where 1–5 % of pairs are half-mapping, `meta` loses 1–5 % of
  the pairs it promises to keep, and it loses precisely the reads that sit at
  host/microbe homology boundaries — exactly the reads assembly needs to close those contigs. The
  three "heads" that justify the tool's existence collapse to two distinct filters, not three.
- **Fix:** make the two strategies actually differ. For a true "drop only if both mates map",
  select on the *pair* rather than the record: e.g.
  `samtools view -@ N -u -e '(flag & 4) || (flag & 8)' in.bam | samtools fastq -1 ... -2 ... -s <orphans> -n -`
  (keep any pair with at least one unmapped mate), or equivalently
  `samtools view -b -F 2 -f 12` for "both" versus the current `-f 12` for "either" after renaming
  them honestly. At minimum, change `_filter_flags` to return a *list* of tokens (see F6) and add a
  unit test asserting the two strategies produce different token lists.

### F2. A non-`sr` minimap2 preset on paired input silently produces empty R1/R2

- **Severity:** critical
- **Location:** `cerberus/stages/align.py:53-79` (preset from `cerberus/autotune.py:77-87`, `autotune.py:156-158`)
- **What:** minimap2 only performs paired-end ("fragment") alignment when the preset enables
  `--frag=yes`, which only `-x sr` does. With `map-ont` / `map-hifi` / `map-pb`, the two FASTQ
  arguments (`align.py:60-61`) are processed as two independent single-end read sets: no record
  carries 0x1/0x40/0x80/0x8. `samtools fastq -f 12` then matches nothing, and `-f 4` matches
  records that have neither READ1 nor READ2 designated — which are routed to `-0 /dev/null`
  (`align.py:74`).
- **Trigger:** two independent routes, neither of which is validated in `orchestrator.validate_config`.
  (a) `cerberus -r1 R1.fq.gz -r2 R2.fq.gz --meta --platform ont` — `apply_user_overrides`
  (`autotune.py:156-158`) sets `minimap2_preset="map-ont"` for short paired input.
  (b) Paired input with mean read length ≥ 500 bp: `classify_length` → `LONG`, and
  `detect_platform_from_fastp` (`autotune.py:101-106`) returns `PACBIO_HIFI`/`ONT` for
  `mean_len >= 500`, so the preset becomes `map-hifi`/`map-ont` for genuinely paired data.
- **Consequence:** `<sample>.meta.R1.fastq.gz` and `R2` are valid but contain **zero reads**, and
  the run exits 0. `accounting.tsv` records `final_paired_r1 = 0` while the flagstat file next to it
  shows millions of unmapped reads — nothing cross-checks the two. In profiling `--fast` the same
  path yields an empty `<sample>.profiling.fastq.gz`.
- **Fix:** in `minimap2_paired`, assert the preset supports fragment mode
  (`if tuned.minimap2_preset != "sr": raise ValueError(...)` or force-append `--frag=yes`), and in
  `validate_config` reject `--platform ont|pacbio-*` combined with `-r1/-r2` unless `--long`.
  Additionally assert `count_reads(out_r1) == count_reads(out_r2)` and warn when the output is empty
  but the BAM contained unmapped records.

### F3. `pipe()` only checks the last process, so an aligner crash yields a silently truncated dataset

- **Severity:** high
- **Location:** `cerberus/utils/shell.py:132-139` (used by `align.py:64,121,169,217`)
- **What:** `for p in procs: p.wait()` then `rc = procs[-1].returncode`. Only `samtools view`'s exit
  status is checked; minimap2/bowtie2/winnowmap return codes are discarded.
- **Trigger:** minimap2 is OOM-killed part-way through loading/streaming (the masked T2T `.mmi`
  needs ~8–10 GB and nothing enforces `cfg.memory_gb` on minimap2), or bowtie2 aborts on a
  corrupt index after emitting the SAM header and some alignments.
- **Consequence:** `samtools view` sees a clean EOF after N records, writes a valid BAM containing
  the first N alignments, and exits 0. The pipeline continues, `samtools fastq` extracts the
  unmapped subset of that truncated BAM, and the final FASTQ contains a silently truncated
  fraction of the library — with exit code 0 and no warning. Reads never emitted by the aligner are
  simply absent, so a user can lose 90 % of a sample and only notice from the accounting numbers.
- **Fix:** check every element of `procs`: `bad = [(c, p.returncode) for c, p in zip(cmds, procs) if p.returncode != 0]`
  and raise `ToolError` on the first failure. Also set `pipefail`-equivalent semantics and consider
  comparing `samtools flagstat` "total" against the FASTQ input count.

### F4. Half-mapped pairs are thrown away instead of being kept as orphans

- **Severity:** high
- **Location:** `cerberus/stages/align.py:74-75` and `align.py:177` (`-0 /dev/null -s /dev/null`); unused field at `align.py:19`
- **What:** every paired path sends both the singleton stream and the READ_OTHER stream to
  `/dev/null`. `AlignOutputs.singletons` (`align.py:19`) exists but is never populated by any
  function in the module, and no caller ever reads it.
- **Trigger:** any pair where exactly one mate aligns to host. In `--meta` this is the *defining*
  case for a conservative filter.
- **Consequence:** the non-host mate — a read the pipeline has explicitly judged to be non-host —
  is deleted rather than routed to `<sample>.meta.orphans.fastq.gz`, which already exists as an
  output slot and is already fed by the fastp orphan stream. Assemblers lose usable single-end
  coverage; the loss is invisible because nothing counts what lands in `/dev/null`.
- **Fix:** write `-s <workdir>/<tag>.unmapped.singletons.fq.gz` (and `-0` to a real file), populate
  `AlignOutputs.singletons`, and in `meta.py` concatenate it into the orphan stream before the
  entropy pass. At the very least, count what would be discarded and log it.

### F5. `--confidence 0.5` against a host-only Kraken2 DB inverts the intended safety margin

- **Severity:** high
- **Location:** `cerberus/stages/kraken.py:62` and `kraken.py:113`
- **What:** kraken2's confidence score is (k-mers assigned to the taxon's clade) / (total k-mers in
  the query), where k-mers with no DB hit count in the denominator. The `--confidence` knob exists
  to suppress *false positives* in a multi-taxon DB. Here the DB contains only human/chimp/gorilla/
  mouse/rat (`README.md:94`), so any hit is host by construction and the only possible error is a
  false *negative*. Requiring 50 % of a read's k-mers to hit host makes false negatives the default
  failure mode.
- **Trigger:** (computed in `work/p1/flag_math.py`, k=35) a 150 bp human read with 2 well-spaced
  sequencing errors retains 46/116 host k-mers → confidence 0.40 → **unclassified**. With
  `--paired`, mates are scored jointly: a pair with one perfectly human mate and one microbial mate
  scores exactly 0.50 (borderline), and a *single* mismatch anywhere in the human mate drops it to
  0.35 → the whole pair, including the intact human mate, is written to `--unclassified-out`.
  100 bp reads escape with a single mismatch (0.47).
- **Consequence:** the first of the two "orthogonal mechanisms" that make the GDPR output
  "publication-defensible" (`README.md:74`) leaks a large class of human reads — every human read
  carrying ≥2 errors, every human read straddling a masked/absent region, and every human/microbial
  chimeric pair. The claim of "**Zero detectable human reads**" (`README.md:11`) then rests entirely
  on the second (minimap2) pass, i.e. on one mechanism, not two.
- **Fix:** use `--confidence 0.0` (or ≤0.05) for a host-only exclusion DB — any host hit should
  condemn the read — and document the choice. If a non-zero value is wanted for a chimera-tolerant
  mode, expose it as a CLI flag and default it to 0 for `--gdpr`. Also consider
  `--minimum-hit-groups 2` rather than a fractional threshold.

### F6. The `-f 4 -F 8` comment is both unimplemented and semantically inverted

- **Severity:** medium
- **Location:** `cerberus/stages/align.py:236-241`
- **What:** the comment states "Implemented by `-f 4 -F 8` below", but the function returns
  `{"filter_flag": "-f", "filter_value": "4"}` — a single flag/value pair. The `dict[str, str]`
  return type structurally cannot express a second filter, and the caller
  (`align.py:71`) splices exactly two tokens.
- **Trigger:** a future maintainer "fixing" the code to match the comment.
- **Consequence:** `-f 4 -F 8` means "read unmapped **AND** mate NOT unmapped", i.e. it keeps only
  the unmapped mates of half-mapped pairs and **excludes every both-unmapped pair** — the exact
  inverse of the intent. My enumeration confirms: with `-f 4 -F 8`, `both unmapped → dropped` and
  `R1 unmapped / R2 mapped → singleton`. Applied literally, `<sample>.meta.R*.fastq.gz` would become
  empty (everything routed to `-s /dev/null`) while the pipeline reported success.
- **Fix:** delete the misleading comment, change the return type to `list[str]` (e.g.
  `["-f", "12"]` / `["-F", "0x900", "-f", "4"]`) so multi-flag filters are expressible, and add a
  unit test pinning the exact token list for each strategy.

### F7. Name-collation is required by `samtools fastq` but never guaranteed, asserted, or documented

- **Severity:** medium
- **Location:** `cerberus/stages/align.py:63-79`, `align.py:168-181`
- **What:** the samtools manual states the input to `samtools fastq` must be name-collated;
  the code pipes the aligner straight into `samtools view -b` (no `samtools collate`, no comment
  noting the requirement). It happens to hold today because minimap2 emits all records of a
  fragment consecutively and bowtie2 writes both mates of a pair as one unit — but nothing in the
  code encodes or checks that invariant, and `cfg.minimap2_args` / `cfg.bowtie2_args`
  (`align.py:50-51`, `align.py:156-157`) inject arbitrary user tokens into the aligner command.
- **Trigger:** any change that breaks adjacency — inserting a coordinate `samtools sort` for
  flagstat/depth purposes, a future aligner or a `--reorder`-less multi-writer, or user args that
  alter output grouping.
- **Consequence:** because `-s /dev/null` is set, a non-collated BAM does not fail loudly: every
  record becomes its own name group, every group looks like a singleton, and **100 % of reads are
  written to `/dev/null`**. R1 and R2 are both empty (not desynchronised), the pipeline exits 0, and
  the only symptom is a zero in `accounting.tsv`.
- **Fix:** insert `samtools collate -u -O` (or `samtools sort -n -u`) between the aligner and
  `samtools fastq`, or document the adjacency invariant at the call site and assert
  `count_reads(out_r1) == count_reads(out_r2) > 0` (when the BAM is non-empty) after the call.
  Adding `-F 0x900` would also make the filter robust to any future path that emits secondaries.

### F8. Kraken2's missing-output fallback fabricates a zero-read result (and writes files during `--dry-run`)

- **Severity:** medium
- **Location:** `cerberus/stages/kraken.py:141-152`
- **What:** if the `_1.fq`/`_2.fq` glob finds nothing, the code logs a warning and `touch()`es two
  empty files, which are then gzipped and returned as the cleaned reads. The `touch()` calls are not
  guarded by `cfg.dry_run` (contrast `_gzip_inplace`, which returns early at `kraken.py:157`).
- **Trigger:** any kraken2 build whose `#` substitution differs from `_1`/`_2` (the docstring at
  `kraken.py:3-6` explicitly anticipates version variation), a full disk during kraken2's output
  write, or a stale/renamed template. Also fires on every `--dry-run` invocation, since no
  kraken2 output exists.
- **Consequence:** `<sample>.<mode>.R1_GDPR.fastq.gz` / `R2_GDPR` are empty; the subsequent
  minimap2 pass runs on nothing and the run reports success. A "zero human reads" deliverable that
  is actually "zero reads" is indistinguishable from a successful scrub except by reading
  `accounting.json`. In dry-run mode the function additionally creates real files inside the work
  directory.
- **Fix:** raise `FileNotFoundError` listing `sorted(workdir.iterdir())` instead of fabricating
  placeholders; only tolerate missing files when kraken2's report shows zero input reads. Guard the
  `touch()` path with `if cfg.dry_run: return ...`.

### F9. `mcf=0.5` makes the auxiliary k-mer pass a no-op for long and error-prone reads

- **Severity:** medium
- **Location:** `cerberus/stages/kmer.py:69-70` and `kmer.py:115-116` (invoked for long reads at `cerberus/pipelines/long_read.py:81-88`)
- **What:** `mcf` (mincoveredfraction) requires ≥50 % of a read's **bases** to be covered by
  reference k-mers before the read is called a match. The value is hard-coded and independent of
  read length and platform error rate, while `bbduk_k` rises to 31 for the LONG/VERY_LONG classes
  (`autotune.py:62-73`).
- **Trigger:** ONT long reads. Simulation (`work/p1/mcf_sim.py`, 10 kb read that is 100 % host):
  covered fraction is 0.503 at 5 % error with k=31 (a coin flip against the 0.5 threshold) and
  0.122 at 10 % error — i.e. **no match**. Separately, any read that is only partially host (a 10 kb
  read with a 2 kb human mitochondrial segment) tops out at 0.20 coverage and can never match.
- **Consequence:** in `run_long_profiling` and in the long-read GDPR path the bbduk aux stage
  consumes runtime and produces an output file identical to its input, while the logs and the
  architecture diagram (`README.md:64`) present it as an active host-removal mechanism. Human
  mitochondrial and rRNA-derived long reads pass through untouched.
- **Fix:** scale the criterion with read length — use `mkf`/`minkmerhits` or an absolute
  `mincovered`-style bases threshold (e.g. "≥300 covered bases") for long reads, lower `k` to
  19–21 for ONT, and skip/disable the aux pass explicitly (with a log line) when
  `read_length_class` is LONG/VERY_LONG rather than pretending it ran.

### F10. Entropy window is fixed at 50 bp while `min_length` can be 35 bp, and the entropy stats file records nothing

- **Severity:** medium
- **Location:** `cerberus/stages/entropy.py:44-45` and `entropy.py:74-75` (vs `cerberus/autotune.py:44`)
- **What:** `entropywindow=50` and `entropyk=5` are hard-coded for every read-length class. For the
  VERY_SHORT class, fastp's `--length_required` is 35 (`autotune.py:44`, `qc.py:53`), so reads
  shorter than the entropy window reach BBDuk. BBDuk's behaviour for reads shorter than
  `entropywindow` is version-dependent (whole-read fallback vs. skip) and is not pinned by any test
  here, so the effective cutoff for 2×50 bp libraries is unknown.
- **Trigger:** `cerberus -r1 R1.fq.gz -r2 R2.fq.gz --all` on a 2×50 bp library (autotune →
  VERY_SHORT, `entropy=0.60`, min_length=35).
- **Consequence:** either short reads bypass the low-complexity filter entirely (host satellite/
  poly-A reads survive into a "GDPR" output) or they are all removed (silent loss of the shortest
  reads). Neither outcome is observable, because `stats={stats}` on a ref-free BBDuk call records
  only reference-match counts — the returned `EntropyOutputs.stats` file carries no information
  about how many reads the entropy filter removed, and no stage count is added to `accounting`.
- **Fix:** derive the window from the tuned class (`entropywindow=min(50, max(20, min_length))`) or
  pass `entropywindow` through `TunedParams`; and record pre/post read counts around the entropy
  stage via `RunAccounting.add_stage` so the filter's effect is auditable.

### F11. The "empty input" guard is a magic byte-size threshold that mis-fires in both directions

- **Severity:** low
- **Location:** `cerberus/stages/kmer.py:52` and `kmer.py:94`
- **What:** `if not r1_in.exists() or r1_in.stat().st_size <= 40:` treats the input as empty. The
  size of an empty gzip member depends on the embedded FNAME field, i.e. on the producer's
  *filename length*: measured (`work/p1/gzname.py`) 25 bytes for `a.fq.gz` but **50 bytes** for
  `03_bbduk_kmer.kmerclean.R1.fq.gz` and 53 bytes for `02_minimap2_human.unmapped.R1.fq.gz` — both
  above the threshold. Conversely a real single-record file can be as small as 45 bytes. Only `r1_in`
  is inspected; `r2_in` is never checked. The equivalent guard is absent from `entropy.py`, so the
  empty files this function fabricates are handed straight to another BBDuk call anyway.
- **Trigger:** an empty stage output written by Python's `gzip` (i.e. `_write_empty_gzip` in
  `kraken.py:175` or `_empty_paired_passthrough` itself) with a long stage-tag filename.
- **Consequence:** the guard fails to fire, BBDuk is invoked on a zero-read input, and whatever
  BBDuk does with that (error or silent pass-through) determines the run's fate; conversely a
  genuinely tiny but non-empty library can be discarded wholesale and replaced with empty outputs
  plus a `# bbduk skipped` stats file.
- **Fix:** replace the size heuristic with `count_reads(r1_in) == 0` (already available in
  `cerberus/utils/fastq.py:11`), check both mates, and apply the same guard consistently in
  `entropy.py` — or drop the guard entirely and let BBDuk handle empty input.

### F12. `-n` yields identical R1/R2 read names, and the profiling concat emits duplicate IDs

- **Severity:** low
- **Location:** `cerberus/stages/align.py:76`, `align.py:178` (`-n`) and `cerberus/pipelines/profiling.py:140-145`
- **What:** `-n` suppresses the `/1` and `/2` suffixes, so mates share an identical name; the
  profiling pipeline then concatenates R1, R2 and orphans into a single FASTQ, producing a file in
  which every ID appears at least twice.
- **Trigger:** every `--profiling` run.
- **Consequence:** harmless for Kraken2/Bracken, but breaks any downstream step that assumes unique
  IDs (seqkit rmdup, read-level joins back to the kraken2 `--output` TSV, most dedup tools), and
  makes the merged file unusable as evidence for per-read provenance in a GDPR audit.
- **Fix:** drop `-n` (letting samtools append `/1`, `/2`) for the profiling path, or add
  `-N`/rename during the concat so the merged file has unique identifiers.
