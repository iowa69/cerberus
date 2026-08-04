# Cerberus v0.1.1 — ten-pass review

Ten passes examined the codebase from ten different angles, each with its own
lens and no knowledge of the others' conclusions. **146 findings: 15 critical,
37 high, 60 medium, 34 low.**

GitHub renders the Markdown below directly. The `.html` files are the same
reports with severity styling and cross-links — clone the repo and open
`index.html` to read them that way.

| Pass | Lens | Report |
|---|---|---|
| 1 | Read-filtering semantics — SAM flag arithmetic, pair synchronisation | [md](pass01.md) · [html](pass01.html) |
| 2 | Subprocess and I/O layer — error propagation, file descriptors | [md](pass02.md) · [html](pass02.html) |
| 3 | Reference manager — download, verification, extraction, caching | [md](pass03.md) · [html](pass03.html) |
| 4 | CLI and configuration — argument handling, validation, exit codes | [md](pass04.md) · [html](pass04.html) |
| 5 | Autotuning — do the tuned parameters reach the tools that use them? | [md](pass05.md) · [html](pass05.html) |
| 6 | Pipeline composition — the full mode × flag matrix | [md](pass06.md) · [html](pass06.html) |
| 7 | Outputs and accounting — what the user actually receives | [md](pass07.md) · [html](pass07.html) |
| 8 | Adversarial audit of the "zero detectable human reads" claim | [md](pass08.md) · [html](pass08.html) |
| 9 | Documentation and packaging — claim-by-claim against the code | [md](pass09.md) · [html](pass09.html) |
| 10 | Operational robustness — memory, threads, disk, crash safety | [md](pass10.md) · [html](pass10.html) |

## How the findings were verified

Findings were checked against a running pipeline, not just read off the source.
Miniature but real references — a minimap2 index, a bowtie2 index, a Kraken2
database and auxiliary FASTA, all built from synthetic genomes — let the whole
tool execute end-to-end on reads whose origin is encoded in their names, so
host removal could be measured rather than argued about.

That step earned its keep. It disproved an early hypothesis about pair
desynchronisation (samtools' singleton routing was quietly protecting the
pairing) and confirmed the more serious bug underneath: the "conservative" and
"aggressive" filters were selecting identical reads, so `--meta` was discarding
the microbial reads it promised to keep.

## The three that mattered most

1. **A long-read minimap2 preset on paired input wrote zero reads and exited 0.**
   minimap2 only emits paired records under `-x sr`; other presets leave no
   READ1/READ2 flags, so `samtools fastq -1/-2` produced nothing.
2. **The conservative and aggressive filters were the same filter.** On a
   fixture of read pairs with one host mate and one microbial mate, meta now
   keeps 500 pairs and profiling drops all 500 — a difference that did not
   previously exist.
3. **A dead process in a pipe reported success.** Only the last command's exit
   status was checked, so an OOM-killed aligner produced a truncated dataset
   that looked like a biological result.

## Outcome

Everything here is fixed in v0.2.0, alongside 50 regression tests — one per
defect, each documenting the wrong behaviour it prevents — and an end-to-end
smoke test that asserts host removal by read provenance rather than by exit
code. The same fixture now runs in CI as `scripts/smoke_test.sh`.

See [`../../CHANGELOG.md`](../../CHANGELOG.md) for the complete fix list.
