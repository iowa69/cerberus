# Cerberus

**Three-headed host removal for metagenomic data.**

Cerberus is an opinionated, all-in-one host-decontamination pipeline that produces three publication-ready outputs from a single run:

| Head | Output | Use case |
|---|---|---|
| **meta** | `<sample>.meta.R1.fastq.gz` + `R2` + orphans | Assembly (SPAdes, MEGAHIT). Conservative — retains microbial reads even at the cost of some residual host. |
| **profiling** | `<sample>.profiling.fastq.gz` | Taxonomic profiling (Kraken2, Bracken). Aggressive — single merged FASTQ, host removed hard. |
| **gdpr** | `<sample>.<mode>*_GDPR.fastq.gz` | Public release. Host removed by **three mechanisms in series** (Kraken2 + human k-mers + minimap2), with the residual measured and reported. See [What the GDPR head does and does not guarantee](#what-the-gdpr-head-does-and-does-not-guarantee). |

Works on **Illumina paired-end short reads**, **ONT long reads**, **PacBio HiFi**, and **PacBio CLR**. Autotunes its parameters from the data — you do not need to know the read length, platform, or sensible thresholds.

- GitHub: <https://github.com/iowa69/cerberus>
- Reference data: [Zenodo DOI 10.5281/zenodo.20258069](https://doi.org/10.5281/zenodo.20258069)
- License: MIT

---

## Quick start

```bash
# Install
conda install -c bioconda cerberus-mg            # once the recipe is merged

# Or from source
git clone https://github.com/iowa69/cerberus
cd cerberus
conda env create -f environment.yml
conda activate cerberus

# Run it
cerberus -r1 R1.fq.gz -r2 R2.fq.gz -o out/ -t 8 --all
cerberus --long -i reads.fq.gz -o out/ -t 8 --all
```

First run downloads the reference bundle (~23 GB compressed, ~24 GB on disk after extraction; archives are deleted once extracted) from Zenodo into `~/.cerberus/refs/`. Pre-warm with `cerberus fetch-refs`; validate the installation with `cerberus doctor`.

---

## Why three outputs?

Existing tools force you to pick one tradeoff:
- **Hostile** is precise — keeps microbes, but leaves residual host. Bad for GDPR-restricted public release.
- **Kneaddata** is aggressive — strips host, but kills 10× more microbes than Hostile. Bad for assembly.
- **Kraken2 host filter** is k-mer based — catches what alignment misses, but doesn't satisfy "0 human reads" reviewers asking for mechanism diversity.

Cerberus runs the right pipeline for each downstream task — **from raw reads each time, not iteratively** — and emits all three outputs in a single invocation. No more rerunning your decontamination tool because you now want to publish what you assembled.

---

## Architecture

```
INPUT (paired-end OR long reads)
   │
   ▼
fastp / fastplong  ──→  autotune (read-length & platform)
   │
   ├──── meta ─────────► minimap2 (conservative) → entropy → R1/R2/orphans
   │
   ├──── profiling ────► bowtie2-vsl (or minimap2 fast)
   │                       → bbduk k-mer (aux refs)
   │                       → entropy
   │                       → single merged FASTQ
   │
   └──── gdpr ─────────► consumes meta/profiling outputs
                          → Kraken2 (compact human+mammal DB)
                          → minimap2 (masked T2T+HLA, different mechanism)
                          → *_GDPR.fastq.gz
```

**Why three mechanisms for GDPR?** Kraken2 alone has false negatives at the read level (it is a minimizer-based classifier optimised for taxonomy, not exclusion). Cerberus follows it with a bbduk pass against a human k-mer reference and then an independent alignment pass, so a read survives only if all three agree it is not host.

Note that the GDPR Kraken2 database contains **host taxa only**. Every classification is therefore a hit on the thing being removed, and a high `--confidence` simply lets host reads escape as "unclassified". The default is deliberately low (`--gdpr-confidence 0.05`); raise it only if you are knowingly trading host removal for retained microbial signal.

---

## How the human reference was built

The default reference bundle (Zenodo: [10.5281/zenodo.20258069](https://doi.org/10.5281/zenodo.20258069)) was produced by three scripts in [`scripts/build_refs/`](scripts/build_refs/):

### `mask_t2t_hla.sh` — masked T2T-CHM13v2.0 + IPD-IMGT/HLA

1. Download **T2T-CHM13v2.0** (NCBI `GCA_009914755.4`) and **IPD-IMGT/HLA** (EBI `hla_gen.fasta`).
2. Concatenate into a single FASTA.
3. **bbmask** with entropy=0.7, window=80, plus k-mer repeat masking (`kr=5 minlen=40 mincount=4`) to N-out low-complexity and tandem repeat regions.
4. **minimap2 (asm5)** RefSeq viral genomes against the masked reference; convert hits (≥85% identity) to BED with `bedtools merge -d 50`; **bedtools maskfasta** to N-out viral-homologous regions.
5. Build a **minimap2** index (`.mmi`) and a **bowtie2** index for the final masked FASTA. Pack the bowtie2 set into `tar.zst`.

Microbial-homology masking against UHGG (bacterial pan-genome) is deferred to v0.2; the bbduk auxiliary k-mer pass in the profiling stage compensates for residual bacterial-like host regions.

### `build_kraken2_gdpr.sh` — compact Kraken2 DB for `--gdpr`

Coverage: human (T2T-CHM13v2.0), chimp, gorilla, mouse, rat. Built using the `>kraken:taxid|N|...` header-injection technique so we don't need NCBI's 7 GB accession-to-taxid maps. Final DB is ~14 GB extracted (~11 GB compressed).

### `curate_aux_refs.sh` — auxiliary k-mer references for `--profiling`

Human host-decoy ncRNA from Ensembl 113: rRNA, snRNA, snoRNA, scaRNA, miRNA, Y_RNA, vaultRNA, ribozyme, Mt_rRNA, Mt_tRNA. Plus the human mitochondrion (NC_012920.1). Concatenated into a single gzipped multi-FASTA used by the bbduk k-mer pass in the profiling pipeline.

---

## Power-user usage

### Override platform autodetection

Cerberus runs fastp first, reads the JSON, and picks parameters. You can override:

```bash
cerberus -r1 R1.fq.gz -r2 R2.fq.gz --all --platform illumina
cerberus --long -i reads.fq.gz --all --platform ont
cerberus --long -i reads.fq.gz --all --platform pacbio-hifi
cerberus --long -i reads.fq.gz --all --platform pacbio-clr
```

### Force the profiling pipeline lane

```bash
# minimap2-only (3× faster, ~2% less sensitive than the standard path)
cerberus -r1 R1.fq.gz -r2 R2.fq.gz --profiling --fast

# bowtie2 + minimap2 stacked (slightly more aggressive — for edge cases)
cerberus -r1 R1.fq.gz -r2 R2.fq.gz --profiling --double-pass
```

### Per-stage knobs (everything under `--help-all`)

```bash
# Tighter QC for low-quality short reads
cerberus -r1 R1.fq.gz -r2 R2.fq.gz --all --min-length 75 --min-quality 25

# Looser entropy filter
cerberus -r1 R1.fq.gz -r2 R2.fq.gz --all --entropy 0.55

# Custom k for the aux k-mer pass
cerberus -r1 R1.fq.gz -r2 R2.fq.gz --profiling --bbduk-k 31

# Inject aligner flags directly
cerberus -r1 R1.fq.gz -r2 R2.fq.gz --profiling \
  --minimap2-args "-N 5 --secondary=no" \
  --bowtie2-args "--score-min L,0,-0.4"

# Resource caps
cerberus -r1 R1.fq.gz -r2 R2.fq.gz --all -t 16 --memory 24G

# Keep all intermediate BAMs/FASTQs for debugging
cerberus -r1 R1.fq.gz -r2 R2.fq.gz --all --keep-intermediates -v
```

Run `cerberus --help-all` to see every flag.

---

## Custom host reference (non-human, plant, mouse, etc.)

If you're working with a non-human host, build a Cerberus-compatible reference directory from any FASTA with one command:

```bash
# Mouse, fully indexed for both short and long reads, with masking
bash scripts/build_custom_host_ref.sh \
  -i GRCm39.fa.gz \
  -o ~/cerberus_refs/mouse \
  --mask \
  -t 16

# Use it
cerberus -r1 R1.fq.gz -r2 R2.fq.gz -o out/ --all \
  --ref-dir ~/cerberus_refs/mouse
```

`build_custom_host_ref.sh` produces:

| File | Built by | Used for |
|---|---|---|
| `masked_t2t_hla.mmi` | `minimap2 -x sr -d` | `--meta`, `--profiling --fast` (short reads) |
| `masked_t2t_hla_bt2/...` | `bowtie2-build` | `--profiling` standard path |
| `manifest.json` | (generated) | makes `RefManager` skip the download/verify step |

Options:

```bash
# Short-read indexes only (skip the long-read .mmi)
bash scripts/build_custom_host_ref.sh -i ref.fa -o ./refs --platform short

# Long-read setup only (skip bowtie2 + sr-preset .mmi)
bash scripts/build_custom_host_ref.sh -i ref.fa -o ./refs --platform long

# Reuse the published human Kraken2 GDPR DB alongside a custom host
bash scripts/build_custom_host_ref.sh -i mouse.fa -o ./refs \
  --kraken-db ~/.cerberus/refs/kraken2_gdpr_compact
```

If you don't supply `--kraken-db`, the `--gdpr` mode won't be available with that custom ref dir — Kraken2 needs a built DB and that step requires taxonomy setup beyond this script's scope.

---

## Commands

```
cerberus -r1 FILE -r2 FILE -o DIR [MODES] [OPTIONS]
cerberus --long -i FILE -o DIR [MODES] [OPTIONS]
cerberus fetch-refs            # download references
cerberus doctor                # validate installation
cerberus --help                # brief help
cerberus --help-all            # full help with advanced flags
```

### Required: at least one mode
- `--meta` — paired output for assembly
- `--profiling` — single FASTQ for Kraken2
- `--gdpr` — post-process selected modes for publication
- `--all` — alias for `--meta --profiling --gdpr`

### Common options
- `-t, --threads` — default: all CPUs
- `--memory NG` — bbduk/Kraken2 memory cap (default: autodetect)
- `--platform {auto,illumina,ont,pacbio-hifi,pacbio-clr}` — default: `auto`
- `--fast` — profiling: minimap2-only path
- `--double-pass` — profiling: pre-filter with minimap2 before bowtie2
- `--ref-dir PATH` — override `~/.cerberus/refs/`
- `--no-auto-download` — refuse to download missing references

---

## Memory & disk

| Resource | Need |
|---|---|
| RAM (peak) | ~14 GB on the `--gdpr` step (Kraken2 loads the whole database). 4-6 GB during alignment. `--memory` is honoured by bbduk and is detected from cgroup limits, so it is correct inside containers. |
| Disk (refs) | ~24 GB in `~/.cerberus/refs/` (one-time). Downloaded archives are removed after extraction; set `CERBERUS_KEEP_ARCHIVES=1` to retain them. |
| Disk (run) | Intermediates in `out/_work/` can reach **10-15× the input FASTQ size** while running (Kraken2 writes uncompressed FASTQ). They are deleted when the run finishes unless `--keep-intermediates` is set. Size `out/` accordingly. |

Designed for 16 GB laptops. Tested on 4-core/16 GB.

---

## Review

v0.2.0 is the result of a ten-pass review of v0.1.1 that found and fixed
several defects producing silently wrong scientific output — including two
filter strategies that compiled to the same filter, and a platform flag that
emptied the output while exiting 0. The findings, and how each was verified by
running the pipeline, are in [`docs/review/`](docs/review/) (open
`index.html`). The fixes are listed in [`CHANGELOG.md`](CHANGELOG.md).

**Results produced with v0.1.x should be regenerated.**

---

## What the GDPR head does and does not guarantee

The `--gdpr` head runs three host-removal mechanisms in series, and a read is
published only if all three accept it:

| # | Mechanism | Kind | Reference |
|---|---|---|---|
| 1 | Kraken2 (`--gdpr-confidence`, default 0.05) | exact k-mer minimizer classifier | compact host DB: human, great apes, mouse, rat |
| 2 | bbduk against a human k-mer reference | a second, independent k-mer implementation | `human_k27.fa.gz` |
| 3 | minimap2 alignment, dropping the pair if **either** mate aligns | alignment, not k-mers | masked T2T-CHM13v2.0 + IPD-IMGT/HLA |

**What this does not prove.** All three mechanisms ultimately derive from the
same reference assemblies. Human sequence that is absent from those assemblies
is invisible to every one of them — population-specific insertions, V(D)J
recombination junctions, novel structural variants, and regions removed by the
bacterial/viral masking applied when the index was built. No combination of
reference-based filters can establish that a file contains zero human reads.

Cerberus therefore **measures and reports** what each mechanism removed rather
than asserting a zero. Every run writes `reports/cerberus_report.html` with a
per-mechanism breakdown; a mechanism that removed nothing is flagged, because
that usually means it could not act rather than that the sample was clean.

**If you are releasing data publicly**, treat this as one control among
several: confirm the removal rates in the run report are what you expect for
your sample type, keep `reports/run_record.json` with the release, and follow
your institution's own review process. Cerberus reduces host content by orders
of magnitude; it does not discharge your legal obligations under GDPR or any
other framework.

---

## Outputs

```
out/
├── <sample>.meta.R1.fastq.gz        # if --meta
├── <sample>.meta.R2.fastq.gz
├── <sample>.meta.orphans.fastq.gz
├── <sample>.profiling.fastq.gz      # if --profiling
├── <sample>.meta.R1_GDPR.fastq.gz      # if --gdpr (per source head)
├── <sample>.meta.R2_GDPR.fastq.gz
├── <sample>.meta.orphans_GDPR.fastq.gz  # meta's unpaired leftovers
├── <sample>.profiling_GDPR.fastq.gz    # profiling's scrubbed deliverable
├── reports/
│   ├── cerberus_report.html  # full run report: parameters, per-stage
│   │                         # accounting, output verification, warnings
│   ├── run_record.json       # the same, machine-readable, for methods sections
│   ├── accounting.tsv        # per-stage read counts (reviewer-friendly)
│   └── accounting.json       # same, machine-readable
└── logs/
    └── *.log                 # one per stage, plus a JSONL run log

# --long produces <sample>.long_meta.fastq.gz / <sample>.long_profiling.fastq.gz
# (and <sample>.long_meta_GDPR.fastq.gz with --gdpr).
```

---

## Citation

Pre-print pending. For now: `Cerberus (v0.2.0). https://github.com/iowa69/cerberus  DOI: 10.5281/zenodo.20258069`

### Reference data provenance

The reference bundle is derived from the following sources. If you publish
results produced with Cerberus, please cite them as well as Cerberus itself:

| Source | Used for | Reference |
|---|---|---|
| T2T-CHM13v2.0 | host genome for the minimap2/bowtie2 indices | Nurk *et al.* (2022) *Science* 376:44-53 |
| IPD-IMGT/HLA | HLA allele sequences added to the host reference | Barker *et al.* (2023) *Nucleic Acids Res* 51:D1053 |
| UHGG | bacterial pan-genome used to mask the host reference | Almeida *et al.* (2021) *Nat Biotechnol* 39:105-114 |
| NCBI RefSeq viral | viral pan-genome used to mask the host reference | O'Leary *et al.* (2016) *Nucleic Acids Res* 44:D733 |
| Ensembl | human rRNA/ncRNA in the auxiliary k-mer references | Ensembl release 113 |
| NCBI Taxonomy | taxonomy backbone of the Kraken2 GDPR database | Schoch *et al.* (2020) *Database* baaa062 |

Cerberus itself is MIT-licensed. The reference data retains the licence of its
upstream sources; check those terms before redistributing the bundle.
