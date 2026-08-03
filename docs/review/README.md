# Cerberus v0.1.1 — ten-pass review

Open [`index.html`](index.html) in a browser for the rendered report; the
`pass*.md` files are the raw findings.

Ten passes examined the codebase from ten different angles, each with its own
lens and no knowledge of the others' conclusions:

| Pass | Lens |
|---|---|
| 1 | Read-filtering semantics — SAM flag arithmetic, pair synchronisation |
| 2 | Subprocess and I/O layer — error propagation, file descriptors |
| 3 | Reference manager — download, verification, extraction, caching |
| 4 | CLI and configuration — argument handling, validation, exit codes |
| 5 | Autotuning — do the tuned parameters reach the tools that use them? |
| 6 | Pipeline composition — the full mode × flag matrix |
| 7 | Outputs and accounting — what the user actually receives |
| 8 | Adversarial audit of the "zero detectable human reads" claim |
| 9 | Documentation and packaging — claim-by-claim against the code |
| 10 | Operational robustness — memory, threads, disk, crash safety |

Findings were then verified against a running pipeline. Miniature but real
references — a minimap2 index, a bowtie2 index, a Kraken2 database and
auxiliary FASTA, all built from synthetic genomes — let the whole tool execute
end-to-end on reads whose origin is encoded in their names, so host removal
could be measured rather than argued about.

That step earned its keep: it corrected an early hypothesis about pair
desynchronisation (samtools' singleton routing was quietly protecting the
pairing) and confirmed the more serious finding underneath it — that the
conservative and aggressive filters were selecting identical reads.

Everything found here is fixed in v0.2.0. See [`../../CHANGELOG.md`](../../CHANGELOG.md).
The same fixture now runs in CI as `scripts/smoke_test.sh`.
