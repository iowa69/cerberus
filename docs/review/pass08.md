# Pass 8 — GDPR claim audit

## Summary

**No — the claim is not supported by the code as written.** Two mechanisms do run in series and the samtools filtering is (by accident) strict, so the pipeline is a genuinely aggressive host scrubber; but "**Zero detectable human reads** via dual orthogonal mechanisms" (`README.md:11`) fails on all three of its load-bearing words. *Zero*: nothing in the codebase ever measures residual host content in a `*_GDPR.fastq.gz` — there is no positive control, no spike-in, no re-scan, and no assertion, so the claim is unfalsifiable by construction and an empty output file is reported identically to a clean one. *Detectable*: `--confidence 0.5` (`cerberus/stages/kraken.py:62,113`) against a database that contains **only** the thing being removed inverts Kraken2's cost matrix — measured below, it raises the per-read miss rate of mechanism 1 from ~0.00–0.01% to 0.36–5.6% for Illumina, and to **~100% for ONT ≥2% error and all PacBio CLR**, where mechanism 1 becomes a complete no-op and only minimap2 votes. *Orthogonal*: a third mechanism (`human_kmer_set`, `cerberus/data/default_manifest.json:52`, explicitly labelled "orthogonal to Kraken2", `required_for: ["gdpr"]`, 932 MB) is declared in the shipped manifest and **never referenced by a single line of pipeline code**; and the two mechanisms that do exist share a common blind spot — human sequence absent from CHM13v2.0 — that neither can see.

The defensible statement the code supports is "aggressive dual-pass host depletion with no residual measurement", not "zero human reads". The README as written invites a depositor to treat the output as anonymised data under GDPR Art. 4(1)/Recital 26 on the strength of a guarantee the software does not test.

## Verified working

- **Mechanisms are correctly composed in series** — `gdpr.py:65-72` feeds `k_out.cleaned_r1/cleaned_r2` (Kraken2's unclassified stream) into minimap2, and `gdpr.py:117-124` does the same for singles/long. A read must therefore pass *both* filters to reach the output. The "and" logic the README claims is genuinely implemented.
- **A read whose MATE aligned to human does NOT survive** — empirically verified against samtools 1.23.1 with a synthetic BAM containing all four mapping states. With `-f 4` plus `-1/-2/-s /dev/null` (`align.py:67-79`), pairs where exactly one mate maps are reclassified as singletons by `bam2fq` and routed to `-s`, i.e. `/dev/null`. Only the both-mates-unmapped pair reached `R1`/`R2`. `-f 4` and `-f 12` produced **byte-identical** paired output. This is the strict behaviour the GDPR pass needs. (It is arrived at by accident, not design — see F7.)
- **Supplementary/secondary alignments cannot leak** — `minimap2 --secondary=no` (`align.py:56,116`) plus `samtools fastq`'s default exclusion of 0x100/0x800; a supplementary record implies a mapped primary, which fails `-f 4` regardless.
- **The Kraken2 DB is built from UNMASKED T2T-CHM13v2.0** — `build_kraken2_gdpr.sh:55` pulls the raw NCBI FASTA and `:48` passes `--no-masking` to `kraken2-build --add-to-library`, suppressing dustmasker. This is the single best design decision in the whole GDPR path: the satellite / low-entropy / viral-homologous sequence that `mask_t2t_hla.sh:51-71` N-masks out of the minimap2 index is *retained* in the Kraken2 DB, so on that axis the two mechanisms genuinely are decorrelated. (F1 then throws most of that benefit away.)
- **The Kraken2 DB is multi-mammal on purpose, and the pipeline keeps only `--unclassified-out`** — `build_kraken2_gdpr.sh:55-59` adds chimp/gorilla/mouse/rat alongside human. Because the retained stream is the *unclassified* one (`kraken.py:58,109`), a human read whose LCA resolves to Homininae or Euarchontoglires is still classified and therefore still dropped. Widening the DB strictly increases removal here; it cannot cause a leak.
- **IPD-IMGT/HLA is concatenated into the reference** — `mask_t2t_hla.sh:44-46`. The MHC is the most polymorphic human locus and the one most likely to evade a single-haplotype reference; including the allele database is a real, correctly-motivated mitigation, and the sequences are normal-complexity so `bbmask` largely leaves them intact.
- **The GDPR pass consumes already entropy-filtered input** — `meta.py:67-80` and `profiling.py:101-106` run bbduk `entropy=0.7 entropywindow=50 entropyk=5` before `run_gdpr_for` sees the data (`orchestrator.py:85-87`). Low-complexity human reads — precisely the class that evades alignment — are therefore largely removed *before* the GDPR pass. This is a real (if incidental) third line of defence, and it is why the missing post-GDPR entropy re-filter is not itself a finding.
- **The Kraken2 DB directory is validated, not assumed** — `_find_kraken_db` (`gdpr.py:132-139`) checks for `hash.k2d` at the top level and one level down, and raises `FileNotFoundError` rather than proceeding with a bad path.
- **Tool failures abort the run** — `shell.run()` uses `check=True` and raises `ToolError` on non-zero exit (`utils/shell.py:84-85`). A kraken2 or samtools failure in the GDPR pass terminates the run rather than emitting a partial release. (Exception: `pipe()` — see F9.)
- **Reference integrity is hash-verified for the published bundle** — `refs.py:121-129,158-165` verifies SHA256 on both download and every subsequent `is_satisfied` check, so a corrupted or substituted `.mmi`/Kraken2 DB is caught. (Not true for custom ref dirs, where `sha256` is `""` — see F5.)
- **`minimap2_singles` uses plain `-f 4`, which is correct for unpaired data** — `align.py:125`. For orphans and long reads there is no mate, so retaining "this read unmapped" is exactly right.

## Leak routes

Ordered by likelihood of contributing at least one human read to a released `*_GDPR.fastq.gz`.

1. **Divergent short read from an N-masked reference region.** Kraken2 declines to classify it because <50% of its 35-mers are in the DB (`kraken.py:62`); minimap2 has no seeds because the source locus was N-masked by `mask_t2t_hla.sh:51-53`. Requires both, but the two are *positively correlated*: masked regions are satellite/repeat regions, which are exactly where a read diverges most from CHM13. Order-of-magnitude: P(≥2 mismatches per 150 bp) = 0.36–5.6% (measured, F1) × masked fraction of the assembly (bbmask at `entropy=0.7 window=80` typically N-masks ~8–12% of the human genome) ⇒ roughly **3×10⁻⁴–7×10⁻³ of human reads**, i.e. 300–7000 residual human reads per million human reads entering the GDPR pass, before the correlation penalty.
2. **Human sequence absent from CHM13v2.0 entirely.** In neither the `.mmi` nor the Kraken2 DB, so *both* mechanisms fail with probability ≈1. CHM13 is one haploid hydatidiform mole of predominantly European ancestry; the HPRC pangenome adds ~119 Mb of non-syntenic sequence across 47 assemblies. Concretely: population-specific insertions and non-reference HERV/LINE insertions; reads spanning novel SV breakpoints; **somatically recombined IG/TCR V(D)J junctions** (present in any blood, gut-biopsy or lymphoid metagenome, absent from every germline reference, and highly identifying); and, if the tool is ever pointed at metatranscriptomic data, human reads spanning exon-exon junctions (`-ax sr` does not do spliced alignment and the junction k-mers are not in the DB).
3. **Long-read mode: mechanism 1 is a no-op.** For `--platform ont` at ≥2% error and for all PacBio CLR, expected Kraken2 confidence is 0.49 / 0.011 respectively (measured, F1) — below the hardcoded 0.5, so essentially *every* human read is left unclassified and written to `--unclassified-out`. `--gdpr` on long reads is therefore a single-mechanism minimap2 filter, and any human read minimap2 misses (masked region, satellite, chimeric, subtelomeric) goes straight to the release. This route alone falsifies "dual orthogonal mechanisms" for two of the four advertised platforms.
4. **Custom non-human host ref-dir, a configuration the README explicitly recommends.** `README.md:188-191` tells the user to combine `--ref-dir <mouse>` with `--kraken-db <human GDPR DB>`. `gdpr.py:51` then unconditionally resolves `masked_t2t_hla_minimap2`, which in that ref dir is the **mouse** index (`build_custom_host_ref.sh:114,173`). Mechanism 2 aligns human reads against mouse; human removal rests entirely on Kraken2 at confidence 0.5. No warning is emitted anywhere.
5. **Deliberately deleted human sequence: viral-homology masking.** `mask_t2t_hla.sh:60-71` maps RefSeq viral genomes onto the reference with `-x asm5` and `bedtools maskfasta`-N's every interval at ≥85% identity. This removes *human* sequence — endogenous retroviral LTRs and internal regions, integrated herpesvirus — from the aligner's view by construction. Reads from those loci are invisible to minimap2 with probability 1; they leak whenever route 1's Kraken2 condition also holds.
6. **Chimeric / mixed pair in the paired GDPR path.** Kraken2 `--paired` (`kraken.py:56`) scores the pair as one unit — hit and total k-mer counts are summed across both mates — so a pristine human R1 paired with a microbial R2 scores ≈116/232 = **0.500**, sitting exactly on the threshold; any single error in R1 tips it below and *both* mates are written to `--unclassified-out`. Mechanism 1 fails for exactly the pairs where it matters most. Mechanism 2 catches it only if the human mate aligns — i.e. unless it is also a route-1/route-5 read, in which case the pair is written to `.R1_GDPR`/`.R2_GDPR`.
7. **Reads at or near the k-mer floor.** Kraken2 is built with `--kmer-len 35` (`build_kraken2_gdpr.sh:64`). A 35 bp read yields exactly one 35-mer and cannot satisfy Kraken2's default `--minimum-hit-groups 2` (not overridden anywhere), so mechanism 1 is a no-op for it; anything shorter yields zero k-mers. `autotune.py:44` sets `min_length=35` for the VERY_SHORT class, so fastp will pass such reads through. `-ax sr` on a 35 bp read is also at its sensitivity limit.
8. **User-injected aligner flags reach the GDPR pass.** `cfg.minimap2_args` (`README.md:137-140`) overrides the tuned value for *every* minimap2 invocation including the GDPR one (`align.py:50-51`), `shlex.split` with no allow-list. A flag intended to speed up `--profiling` silently desensitises the publication-release filter.
9. **(Not a leak, but a corrupted release presented as a clean one.)** `_gzip_pair_outputs`'s fallback branch (`kraken.py:141-147`) uses `Path.touch()`, which does **not** truncate. Verified: with R1 present and the R2 glob missing, R1 keeps its real reads while R2 is emitted empty — a desynchronised pair handed to `minimap2_paired`. Combined with F9, the run can complete "successfully" with an empty or truncated `_GDPR` file, which the accounting records as a normal zero count.

## Findings

### F1. `--confidence 0.5` against a host-only database inverts Kraken2's error tradeoff

- **Severity:** critical
- **Location:** `cerberus/stages/kraken.py:62` (paired), `cerberus/stages/kraken.py:113` (single)
- **What:** Kraken2's confidence score is the fraction of a query's k-mers whose LCA falls within the clade rooted at the candidate taxon, divided by the total number of non-ambiguous k-mers in the query (`Q = L − k + 1`, here `k = 35` per `build_kraken2_gdpr.sh:64`). Kraken2 walks the highest-weighted root-to-leaf path from leaf toward root and returns the most specific node scoring ≥ threshold; **if no node on the path — including root — reaches the threshold, the query is reported unclassified** and is written to `--unclassified-out`, i.e. it is *kept* by Cerberus.

  Because this DB contains only human, chimp, gorilla, mouse and rat, root's clade is the entire database, so root's score is simply `hit_kmers / total_kmers`. `--confidence 0.5` therefore reduces to the requirement **"at least half of this read's 35-mers must exist somewhere in the mammalian k-mer set."** It performs no taxonomic discrimination whatsoever — there is nothing in the DB to discriminate against — and is a pure sensitivity tax.

  The threshold exists to suppress false-positive *assignments* in a mixed-taxon DB. In host exclusion the cost matrix is inverted: a false positive costs one microbial read; a false negative leaks personal genetic data into a public archive. The correct value for a host-only exclusion DB is **`--confidence 0`** (accept any hit), optionally with `--minimum-hit-groups 1` for short reads. It is not user-tunable — there is no CLI flag, only the two hardcoded literals.

- **Trigger:** Every `--gdpr` invocation. Measured impact (simulated k-mer survival, `k=35`, mismatches = sequencing error + divergence from CHM13, 60k reads/cell):

  | scenario | P(unclassified) @ conf 0.0 | @ conf 0.05 | @ conf 0.5 |
  |---|---|---|---|
  | 150 bp, Q30, 1 SNV/1000 | 0.0000% | 0.0000% | **0.3550%** |
  | 150 bp, Q25, 1 SNV/700 (non-European donor) | 0.0050% | 0.0117% | **5.6033%** (480×) |
  | 100 bp trimmed, Q30 | 0.0050% | 0.0167% | **3.7700%** (226×) |
  | 75 bp trimmed, Q30 | 0.0400% | 0.0917% | **3.9333%** (43×) |

  Deterministically, for a 150 bp read (`Q = 116`): 0 mismatches → conf 1.000 (classified); 1 → 0.698 (classified); **2 → 0.397 (unclassified, leaks)**; 3 → 0.259; 4 → 0.069. **Two well-separated mismatches relative to CHM13 are sufficient to make a 150 bp human read invisible to mechanism 1.** Note the ancestry gradient: a donor further from CHM13's predominantly-European haplotype leaks ~16× more than a matched one.

  For long reads, expected confidence ≈ `(1 − e)^35`:

  | platform | error | E[confidence] | outcome |
  |---|---|---|---|
  | PacBio HiFi | 0.1–0.5% | 0.97–0.84 | classified |
  | ONT R10.4 | 1% | 0.703 | classified (marginal) |
  | ONT R10.4 | 2% | **0.493** | **every human read unclassified** |
  | ONT legacy | 5% | **0.166** | **every human read unclassified** |
  | PacBio CLR | 12–15% | **0.011–0.003** | **every human read unclassified** |

  `autotune.py:65,71` selects `map-ont` for the LONG/VERY_LONG classes and `_platform_preset` honours `--platform ont`/`pacbio-clr`, so these are supported, advertised configurations (`README.md:13,110-112`).

- **Consequence:** For Illumina, mechanism 1's miss rate rises 40–480× — from statistically zero to 0.36–5.6% of human reads — leaving minimap2 as the effective sole filter for that fraction and making the two mechanisms' failures correlated rather than independent. For ONT ≥2% and all PacBio CLR, mechanism 1 removes **nothing**, and the "a read survives only if both a k-mer classifier and an aligner agree it isn't host" argument at `README.md:74` is simply false. The model above is optimistic: it assumes every error-free 35-mer is present in the DB, ignoring indels, non-reference alleles, and quality-trimmed fragments.
- **Fix:** Set `--confidence 0` for both calls, and add `--minimum-hit-groups 1`. Expose the value as `--gdpr-kraken-confidence` (default 0) with a loud warning above 0, and record the value used in `reports/`. Do not treat the resulting microbial-read loss as a regression: for the GDPR head the whole point is that false positives are cheap and false negatives are not. If microbial loss is unacceptable, that is an argument for a separate `--gdpr-strict`/`--gdpr-balanced` split, not for a silently insensitive default.

### F2. "Zero detectable human reads" with no detection step anywhere in the codebase

- **Severity:** critical
- **Location:** `README.md:11`; `cerberus/pipelines/gdpr.py:38-101`; `cerberus/orchestrator.py:85-93`; `cerberus/accounting.py:43-54`
- **What:** Confirmed by exhaustive reading: **there is no verification step.** `run_gdpr_for` writes the final files and returns; `orchestrator.py:88-93` calls only `accounting.add_final(...)`, which records `{path, reads, size_bytes}` per output. There are **zero `add_stage` calls for the GDPR pass**, so `accounting.tsv` contains no row for how many reads Kraken2 removed, how many minimap2 removed, or what fraction of input survived. There is no residual-host scan, no spike-in positive control, no assertion that the output is non-empty, and no check that the Kraken2 DB even contains taxid 9606 (`_find_kraken_db` at `gdpr.py:132-139` verifies only that `hash.k2d` exists — a DB whose human library was silently dropped at build time would pass and would classify nothing). `GDPR_DROP_TAXA` at `kraken.py:21` is declared and never referenced by any code path.
- **Trigger:** Every run. Also: if `_gzip_pair_outputs` takes its fallback branch (`kraken.py:141-147`) or `pipe()` swallows a minimap2 failure (F9), the GDPR output is empty and `add_final` records `reads: 0` — **byte-identical in the report to a genuine 100%-host sample.** A reviewer reading `accounting.tsv` cannot distinguish "we removed all human reads" from "the pipeline silently produced nothing".
- **Consequence:** The headline claim is unfalsifiable by the tool that makes it. "Detectable" is undefined — detectable by what, at what sensitivity, at what depth? A depositor has no per-run evidence to hand a reviewer, IRB or DPO beyond "the software says so", and the one failure mode that would be most catastrophic to miss (a silently empty or silently unfiltered release) is indistinguishable from success in the accounting.
- **Fix:** Four additions, in priority order.
  1. **Spike-in positive control (the only thing that measures per-run sensitivity).** Inject N labelled reads (e.g. 10,000) simulated from a human assembly *not* in the reference — an HPRC/1KG sample of different ancestry to CHM13 — into the QC'd stream with a recognisable name prefix; track them by name; assert 0 survive to `_GDPR`; report the residual count and the implied per-million rate in `accounting.tsv`. Strip them from the final output.
  2. **Residual scan of the final file.** Re-run Kraken2 at `--confidence 0`, or bbduk `k=27` against `human_k27.fa.gz` — which is already downloaded and sitting unused on disk (F3) — over each `*_GDPR.fastq.gz`; fail the run, or at minimum emit a prominent non-zero warning, if any hit is found.
  3. **Per-stage GDPR accounting.** `gdpr_input`, `gdpr_kraken2_removed`, `gdpr_minimap2_removed`, `gdpr_final` for each source mode, plus an explicit `WARNING: GDPR output has 0 reads` guard so empty ≠ clean.
  4. **Provenance block in `reports/`.** kraken2 and minimap2 versions, the confidence value used, and the SHA256 of the DB and `.mmi` actually consumed — without this the release is not reproducible or auditable.

### F3. `human_kmer_set` — a third "orthogonal" mechanism that is declared, downloaded, and never used

- **Severity:** high
- **Location:** `cerberus/data/default_manifest.json:52-62`; `cerberus/refs.py:41`; `scripts/build_refs/curate_aux_refs.sh:83-88`
- **What:** The shipped manifest declares an asset described as *"Human-specific 27-mer set for GDPR belt-and-braces bbduk scrub (orthogonal to Kraken2)"* with `required_for: ["gdpr"]` and `size_bytes: 932696125`. Grep across the entire package: `human_kmer_set` and `human_k27` appear **only** in `default_manifest.json`, in `scripts/zenodo_*.py`, and in the two ref-dir generator shell scripts. No pipeline, stage or ref-resolution code ever requests it. `_PIPELINE_TO_ASSETS["gdpr"]` (`refs.py:41`) lists only `kraken2_gdpr_compact` and `masked_t2t_hla_minimap2`; `required_assets_for` (`refs.py:105-114`) iterates that hardcoded map and **ignores the manifest's `required_for` field entirely**, which is therefore purely decorative.

  Two further problems. First, the description is false twice over: `curate_aux_refs.sh:83-88` shows `human_k27.fa.gz` is a **verbatim `cp` of the full T2T-CHM13v2.0 genome FASTA** — not a k-mer set, not "human-specific", not filtered in any way. Second, `build_custom_host_ref.sh:160-162` symlinks `human_k27.fa.gz` → `aux_refs.fa.gz` (which may itself be an empty gzip placeholder, `:157`) with the comment *"not strictly needed since Cerberus uses minimap2 for the GDPR second-mechanism, but the manifest still lists it"* — so the authors know it is dead and kept the manifest entry anyway.
- **Trigger:** `cerberus fetch-refs` (`refs.py:225-227`) and `cerberus doctor` (`refs.py:229-235`) both iterate *all* manifest assets, so every user downloads 932 MB that nothing reads, and `doctor` reports a hard "missing or corrupt" problem for an asset no pipeline needs.
- **Consequence:** The manifest — a shipped, machine-readable artifact — asserts that the GDPR pass has a bbduk k-mer mechanism orthogonal to Kraken2. It does not. Worse, this is the mechanism that would have closed the largest hole in the design: bbduk `k=27` against the **unmasked** full genome would catch exactly the reads that route 1 and route 5 leak (masked satellite, viral-homology-masked HERV loci), and at `k=27` it tolerates ~2 mismatches per 150 bp read far better than Kraken2 at `k=35`/`conf=0.5`. The asset is already on the user's disk. Not wiring it up is the single largest missed mitigation in the audit.
- **Fix:** Either wire it up — add `human_kmer_set` to `_PIPELINE_TO_ASSETS["gdpr"]` and insert a `bbduk_kmer_single`/`bbduk_kmer_paired` stage after minimap2 in `run_gdpr_for` — or delete the asset from the manifest and the README's mechanism story. Do not ship a manifest that claims a mechanism the code lacks. Independently: make `required_assets_for` derive from the manifest's `required_for` field so the two can never diverge again, and fix the description to say what the file actually is.

### F4. The two mechanisms are not orthogonal for sequence absent from CHM13

- **Severity:** high
- **Location:** `cerberus/pipelines/gdpr.py:4-12` (docstring), `README.md:74`; `scripts/build_refs/mask_t2t_hla.sh:51-71`; `scripts/build_refs/build_kraken2_gdpr.sh:55`
- **What:** Both mechanisms derive from the same source assembly. Credit where due: the Kraken2 DB uses the *unmasked* NCBI FASTA (`build_kraken2_gdpr.sh:55` + `--no-masking` at `:48`) while the minimap2 index uses the *masked* one, so the masking axis genuinely is decorrelated. But masking is not the only axis, and the residual shared blind spot is **any human sequence not present in CHM13v2.0 at all**, which is invisible to both by construction:
  - **Population-specific sequence.** CHM13 is a single haploid hydatidiform mole of predominantly European ancestry. The HPRC pangenome adds ~119 Mb of non-syntenic sequence across 47 assemblies; non-reference insertions and novel SV alleles in a given donor are in neither the `.mmi` nor the k-mer DB. This compounds with F1's measured ancestry gradient (5.6% vs 0.36% miss rate).
  - **Somatic recombination.** IG/TCR V(D)J junctions exist in no germline reference. Any blood, gut-biopsy or lymphoid-tissue metagenome contains them; they are human, identifying, and unreachable by both mechanisms.
  - **Deliberately deleted human sequence.** `mask_t2t_hla.sh:60-71` N-masks every interval where a RefSeq viral genome aligns at ≥85% identity — that is human ERV/LTR sequence removed from the aligner's view on purpose. Mechanism 2 cannot see it at all; mechanism 1 sees it only when F1's confidence gate is cleared.
  - **Non-genomic human sequence.** Spliced exon-exon junctions (`-ax sr` is not a spliced preset, and the junction k-mers are not in a genomic DB) if the tool is ever pointed at metatranscriptomes.
- **Trigger:** Every run, at a rate set by donor ancestry, tissue type and library composition.
- **Consequence:** "Orthogonal" is claimed as the basis for publication-defensibility (`README.md:74`, `gdpr.py:10-12`). Two mechanisms drawn from one assembly are *complementary within that assembly's coverage* and *jointly blind outside it*. The residual is not bounded by the product of two independent miss rates, as the claim implies; for these classes it is bounded by neither.
- **Fix:** State the shared limitation explicitly in the README. Substantively: add a pangenome component (HPRC minigraph-cactus non-reference sequence, or a decoy set) to at least one of the two references so the mechanisms stop sharing a coverage boundary; and wire up F3's unmasked k-mer pass, which at minimum closes the masked-region half of the problem.

### F5. Custom `--ref-dir` silently reduces `--gdpr` to a single mechanism

- **Severity:** high
- **Location:** `cerberus/pipelines/gdpr.py:51`; `scripts/build_custom_host_ref.sh:114,160-177`; `README.md:188-191`
- **What:** `run_gdpr_for` unconditionally resolves `refs.asset("masked_t2t_hla_minimap2")` for mechanism 2. In a ref dir produced by `build_custom_host_ref.sh`, that key maps to the **custom host's** index (`:173`, `filename: "$NAME.mmi"`). The README explicitly documents combining a non-human host with the published human Kraken2 DB (`README.md:188-191`: *"Reuse the published human Kraken2 GDPR DB alongside a custom host"*). In that configuration mechanism 2 aligns reads against mouse (or plant, or whatever), so human removal depends entirely on Kraken2 at confidence 0.5. Nothing checks that the minimap2 index and the Kraken2 DB describe the same organism, and nothing warns.

  Secondary defects in the same path: the generated manifest sets every `sha256` to `""`, so `is_satisfied` (`refs.py:127-129`) short-circuits to "exists ⇒ fine" with no integrity checking at all; `masked_t2t_hla_minimap2` is declared `required_for: ["meta","profiling-fast"]` with no `"gdpr"`, further confirming `required_for` is ignored; and `README.md:174-176` documents a `masked_t2t_hla.long.mmi` "used for `--long` modes" that **no code path references** — `refs.py:42-44` maps every long pipeline back to `masked_t2t_hla_minimap2`, so a `--platform both` custom ref dir runs `-ax map-ont` against a `-x sr` (k=21,w=11) index.
- **Trigger:** Any `--gdpr` run with `--ref-dir` pointing at a `build_custom_host_ref.sh` output — the exact workflow the README recommends.
- **Consequence:** A user follows documented instructions and receives a file named `*_GDPR.fastq.gz`, produced by a code path the README says uses dual orthogonal mechanisms, when in fact only one mechanism ever looked for human sequence — and that one is F1-degraded. For ONT/CLR input, *no* mechanism looks for human sequence.
- **Fix:** Record the organism/reference identity of both GDPR assets in the manifest and refuse (or emit a hard warning and set a flag in `accounting.json`) when the Kraken2 DB and the minimap2 index do not match. Better: give the GDPR pass its own asset key (`gdpr_human_minimap2`) that always resolves to the human index regardless of `--ref-dir` host, so a custom-host run genuinely has two human-facing mechanisms.

### F6. `--kraken2-db` and `--aux-refs` are parsed, stored, and never read

- **Severity:** high
- **Location:** `cerberus/cli.py:140-142`; `cerberus/config.py:92-93`; `cerberus/pipelines/gdpr.py:49-50`
- **What:** `--kraken2-db` and `--aux-refs` are declared as advanced CLI options and threaded into `CerberusConfig.kraken2_db_override` / `aux_refs_override` (`cli.py:210-211`). Grep across the package: those two attributes appear **only** in `config.py` and `cli.py`. No pipeline reads them. `gdpr.py:49-50` unconditionally resolves the bundled `kraken2_gdpr_compact`.
- **Trigger:** Any user passing `--kraken2-db /path/to/my/db` — for example to substitute a more comprehensive human DB, a differently-built DB, or a validated in-house one for a regulated release.
- **Consequence:** The user believes they have redirected the GDPR pass's primary classifier and receives output from the bundled DB instead, with no error, warning, or log line. This is a silently-ignored override on the single most GDPR-critical control surface the CLI exposes. If the substitution was made to satisfy an institutional requirement, the resulting release is non-compliant and the operator has no way to know.
- **Fix:** Honour both overrides in `run_gdpr_for` and `run_profiling`/`run_long_profiling`, or delete the flags. If honoured, log the resolved path and its checksum into `reports/`. A flag that silently does nothing on a compliance-critical path is worse than no flag.

### F7. `_filter_flags("either")` was never implemented; `drop_strategy` is a dead parameter

- **Severity:** medium
- **Location:** `cerberus/stages/align.py:230-242` (esp. `:236-241`), `cerberus/stages/align.py:67-79`, `cerberus/pipelines/gdpr.py:71`, `cerberus/pipelines/meta.py:46`
- **What:** The comment block at `align.py:236-240` reasons its way to `-f 4 -F 8` and then ends *"Implemented by `-f 4 -F 8` below"* — but `:241` returns only `{"filter_flag": "-f", "filter_value": "4"}`, and the call site at `:70-71` splices exactly one flag/value pair into the command. `-F 8` was never emitted. There is no test for any of this (F11).

  Empirically verified against samtools 1.23.1 with a synthetic BAM covering all four mapping states: `-f 4` and `-f 12` produce **identical** `-1`/`-2` output (pairs where one mate maps are demoted to singletons by `bam2fq` and sent to `-s /dev/null`, `align.py:76`). So `drop_strategy` has **no effect**: `"both"` (`meta.py:46`) and `"either"` (`gdpr.py:71`, `profiling.py:57,69`) are the same filter.
- **Trigger:** Every paired alignment in every pipeline.
- **Consequence:** Two separate problems. (a) For the GDPR pass the outcome is *correct by accident* — the strict "drop the pair if either mate maps" behaviour is what is wanted, and F1's leak routes are not made worse. But it is undefended: change `-s /dev/null` to a real file (a plausible "let's keep singletons" edit) and pairs with a human mate immediately start leaking into the release. (b) For `--meta` the outcome is *wrong*: `meta.py:4` and `align.py:39` both document `--meta` as conservative — "only drop reads where both mates map" — and `README.md:9` sells that conservatism as the head's entire reason to exist ("retains microbial reads even at the cost of some residual host"). It is not implemented; `--meta` is exactly as aggressive as `--profiling`, and the three-head differentiation at `README.md:44-49` is partly fiction.

  Additionally, the whole construction depends on an unguarded invariant: `bam2fq` requires **name-collated** input. Verified — with mates non-adjacent, all four records were routed to the singleton file and both `-1`/`-2` outputs were **empty**. minimap2 currently preserves pair adjacency so this does not fire, but nothing asserts it, and the failure mode is a silently empty release (see F2).
- **Fix:** Make the two strategies actually differ and match their docstrings — `"either"` ⇒ `-f 12` (explicit, self-documenting, no reliance on singleton routing); `"both"` ⇒ retain pairs unless both mates map, which requires a real implementation (e.g. `-f 4` writing singletons to a kept file, or a name-keyed pass) rather than the current no-op. Add unit tests over synthetic BAMs asserting exactly which records reach which file for both strategies. Delete the misleading comment at `:236-240`.

### F8. An empty output is indistinguishable from a clean one, and a rerun can silently blank a good result

- **Severity:** medium
- **Location:** `cerberus/stages/kraken.py:129-152` (fallback at `:141-147`), `cerberus/stages/kraken.py:155-172` (esp. `:159-161`); `cerberus/accounting.py:43-54`
- **What:** When the `--unclassified-out` globs miss, the fallback logs a warning and `touch()`es placeholders. Three verified behaviours:
  - `Path.touch()` **does not truncate**. With R1 present and only the R2 glob missing, R1 retains its 3 real reads while R2 is emitted as a valid empty gzip — a **desynchronised pair** handed straight to `minimap2_paired`.
  - With neither present, both outputs are empty gzips and the pipeline proceeds normally to produce an empty release.
  - `_gzip_inplace` at `:159-161`: if `src` does not exist it **overwrites** `dst` with an empty gzip. Verified — a directory containing a previously-good `*_1.fq.gz` (5 reads) came back with 0 reads. Any re-entry into this branch destroys prior good output rather than reusing it.

  The failure is `log.warning` only; it does not raise, does not set a flag, and does not reach `accounting`. Combined with F2, the resulting `accounting.tsv` row is `_final  meta_gdpr.paired_r1  0  <path>` — exactly what a genuine 100%-host sample produces.
- **Trigger:** Any kraken2 build/version whose `--unclassified-out` naming escapes `{root}*1.fq` / `{root}*[12].fastq`; a partial write or full disk; or re-entry into the branch on a rerun.
- **Consequence:** Directly answers the audit question: **no, the failure is not distinguishable from success in the accounting.** A silently empty release would be reported as a successful zero-human-read release, and a desynchronised pair is passed to an aligner that has no pairing guard.
- **Fix:** Raise instead of warning — a missing Kraken2 output is a pipeline failure, not a recoverable condition. If a placeholder path is genuinely needed, `write_bytes(b"")` rather than `touch()`, never overwrite an existing non-empty `.gz`, and record a `gdpr_placeholder_emitted` flag in `accounting.json` so the condition is visible downstream. Add a read-count equality assertion on R1/R2 before every `minimap2_paired` call.

### F9. `pipe()` checks only the last process's exit code

- **Severity:** medium
- **Location:** `cerberus/utils/shell.py:132-139`
- **What:** `pipe()` waits on every process but inspects `procs[-1].returncode` only. In `minimap2 | samtools view -b -o bam -`, a minimap2 failure (index load OOM against the 7.7 GB `.mmi`, corrupt index, bad user-injected `--minimap2-args`) closes the pipe cleanly; `samtools view` reads EOF, writes a valid empty/short BAM, and exits 0. `pipe()` returns success.
- **Trigger:** Any minimap2 (or bowtie2, or winnowmap) failure in any pipeline, including the GDPR pass. Most plausible on the 16 GB laptops the README targets (`README.md:233`), where the actual reference footprint is ~23.5 GB downloaded (measured from the manifest) versus the ~13 GB extracted claimed at `README.md:230`.
- **Consequence:** In the GDPR pass, mechanism 2 silently does nothing and produces an empty or truncated BAM; the run reports success. Feeds directly into F2 and F8 — an empty release is emitted and recorded as a zero count.
- **Fix:** Collect all return codes and raise on the first non-zero, naming the failing stage. Use `set -o pipefail` semantics explicitly. Additionally assert the output BAM is non-empty and has a header before running `samtools fastq` on it.

### F10. The shipped manifest and README misstate how the GDPR reference was built

- **Severity:** medium
- **Location:** `cerberus/data/default_manifest.json:8`, `:53`; `README.md:90`, `:174-176`; `scripts/build_refs/mask_t2t_hla.sh:7-10,81`
- **What:** Three provenance statements in shipped artifacts contradict the build scripts:
  - `default_manifest.json:8` describes the minimap2 index as *"masked against UHGG bacterial and RefSeq viral pan-genomes"*. `mask_t2t_hla.sh:7-10` and `README.md:90` both state plainly that UHGG masking is **deferred to v0.2** and that v0.1 performs only viral masking. The manifest is the machine-readable record of what the GDPR-critical reference *is*, and it is wrong.
  - `README.md:174` states `masked_t2t_hla.mmi` is built by `minimap2 -x sr -d`. `mask_t2t_hla.sh:81` builds it with plain `minimap2 -d` (default k=15, w=10). At alignment time `-ax sr` is passed (`align.py:55`), and minimap2 will emit `[WARNING] Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index` — the index's k/w win. Here that is benign-to-favourable for sensitivity (k=15/w=10 seeds more densely than sr's k=21/w=11), but it is undocumented, untested, and the reverse mismatch does bite in the custom-ref long-read path (F5).
  - `README.md:175` documents `masked_t2t_hla.long.mmi` as "Used for `--long` modes". No code path references it; `refs.py:42-44` maps every long pipeline to `masked_t2t_hla_minimap2`.

  Also: `default_manifest.json:53`'s description of `human_kmer_set` is false on both counts (F3), and `README.md:230`'s "~13 GB extracted" understates the measured 23.47 GB of downloads (plus ~14 GB extracted Kraken2 DB, plus retained archives — `refs.py:171-172` never deletes them post-extraction).
- **Trigger:** Any reviewer or DPO auditing the provenance of the reference used to produce a public release.
- **Consequence:** The provenance chain for the artifact that governs the GDPR claim cannot be trusted from the documentation; it has to be reverse-engineered from the build scripts, and the two disagree. For a tool whose selling point is *publication-defensibility*, an incorrect reference description in the shipped manifest is a direct hit on the value proposition.
- **Fix:** Generate the manifest descriptions from the build scripts rather than by hand, or at minimum reconcile all three statements. Remove the `masked_t2t_hla.long.mmi` row from `README.md:175` or add the code path that uses it. Correct the disk-footprint table.

### F11. Zero test coverage for the entire GDPR code path

- **Severity:** medium
- **Location:** `tests/` (only `test_autotune.py`, `test_cli.py`)
- **What:** The test suite covers autotune parameter tables and CLI argument parsing. There is **no test** for `pipelines/gdpr.py`, `stages/kraken.py`, or `stages/align.py`. The only GDPR references in `tests/` are three argument-validation assertions (`test_cli.py:25,37,47-54`). `_filter_flags`, `_gzip_pair_outputs`, `_gzip_inplace` and `_find_kraken_db` are all pure, trivially testable functions with no test.
- **Trigger:** N/A — a standing gap.
- **Consequence:** F7 (a comment describing flags the code does not emit) and F8 (`touch()` not truncating; `_gzip_inplace` blanking a good file) are precisely the class of defect a handful of unit tests would have caught immediately — I found both in minutes with ~30 lines of test code. The highest-stakes module in the project is the only one with no coverage at all.
- **Fix:** Unit-test `_filter_flags` against its docstring; test `_gzip_pair_outputs` for all four presence/absence combinations and assert it never silently desynchronises or blanks existing data; test `_find_kraken_db` for both layouts and the raising case; add an end-to-end `--dry-run` test asserting the exact kraken2 and samtools command lines (which would pin `--confidence` and make F1 a visible, reviewable decision rather than a buried literal).

### F12. The profiling GDPR pass has no pair-level information

- **Severity:** low
- **Location:** `cerberus/pipelines/profiling.py:137-145`; `cerberus/pipelines/gdpr.py:83-90`
- **What:** `run_profiling` concatenates R1 + R2 + orphans into one FASTQ before returning (`PipelineResult(mode="profiling", singletons=final)`). `run_gdpr_for` therefore takes the `singletons` branch and runs `kraken2_single` / `minimap2_singles`. Every read is judged alone.
- **Trigger:** Every `--profiling --gdpr` run.
- **Consequence:** The paired GDPR path gets "guilt by association" for free — if either mate looks human, both are dropped. The profiling path cannot, so a human read whose mate is obviously human but which is itself from a masked region survives on its own merits. Mechanism strength is strictly lower than for `--meta --gdpr`, and the README makes no such distinction. Inherent to the merged-output design rather than a bug, but it should not be silently equated with the paired path.
- **Fix:** Document the difference, or preserve `/1`-`/2` name suffixes through the concatenation so the GDPR pass can group mates by base name and apply the same either-mate rule.

## GDPR / regulatory framing

Short answer: **"zero human reads" is not a defensible claim to make in software documentation, and no tool can make it.**

- Under GDPR Art. 4(1) and Recital 26, the anonymisation determination is made by the **controller** — the depositor — in the context of the specific dataset, the recipients, and the means reasonably likely to be used for re-identification. A pipeline cannot discharge that assessment on the user's behalf, and a README that says "zero detectable human reads" invites exactly that substitution.
- The residual that matters is small. Human identifiability from sequence is not a function of coverage: on the order of 30–80 common SNPs uniquely identify an individual, and Homer et al. (2008) showed presence-inference from aggregate allele frequencies alone. A few hundred residual human reads in a public metagenome — well within the range F1's measured miss rates imply — are enough to support re-identification against a reference genotype panel. Residual host content is a re-identification-risk quantity, not a QC nuisance.
- "**Detectable**" is doing all the work in the sentence and is undefined. Detectable by whom, with what method, at what sensitivity, at what depth? Because the tool contains no detector (F2), the claim cannot be tested by the person relying on it. A claim that cannot fail is not a guarantee.
- Comparable tools do not make absolute claims. Hostile and the ENA/EBI human-read-removal guidance both report measured removal rates and residual estimates rather than asserting zero. `README.md:47` even critiques Kraken2-only filters for not satisfying *"'0 human reads' reviewers"* — then makes the absolute claim itself, on a Kraken2 configuration tuned in the wrong direction.

**Caveats the README must carry**, replacing the current line 11:

1. A **measured** removal figure from a spiked benchmark, with the residual expressed per million human reads, the donor ancestry used, and the platform — not an absolute.
2. An explicit statement that the two mechanisms share the CHM13v2.0 coverage boundary and that human sequence absent from that assembly (population-specific insertions, V(D)J junctions, novel SVs) is removed by neither (F4).
3. An explicit statement that on ONT and PacBio CLR the Kraken2 mechanism contributes essentially nothing at the current confidence setting, so long-read GDPR output is single-mechanism (F1).
4. An explicit statement that `--ref-dir` with a non-human host disables the human-facing minimap2 mechanism (F5).
5. A clear allocation of responsibility: the user remains the controller; this tool performs depletion, not anonymisation; independent verification before public release is required; and controlled-access deposition (EGA/dbGaP) remains the appropriate route where residual host risk is not acceptable.
6. The `--confidence` value used, surfaced in `reports/` and in the README, so the sensitivity/specificity tradeoff is a visible, reviewable choice rather than a buried literal at `kraken.py:62`.
