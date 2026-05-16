# Cerberus

**Three-headed host removal for metagenomic data.**

Cerberus is an opinionated, all-in-one host-decontamination pipeline that produces three publication-ready outputs from a single run:

| Head | Output | Use case |
|---|---|---|
| **meta** | `<sample>.meta.R1.fastq.gz` + `R2` + orphans | Assembly (SPAdes, MEGAHIT). Conservative — retains microbial reads even at the cost of some residual host. |
| **profiling** | `<sample>.profiling.fastq.gz` | Taxonomic profiling (Kraken2, Bracken). Aggressive — single merged FASTQ, host removed hard. |
| **gdpr** | `<sample>.<mode>.*_GDPR.fastq.gz` | Public release. **Zero detectable human reads** via dual orthogonal mechanisms (Kraken2 + bbduk k-mer). |

Works on **Illumina paired-end short reads**, **ONT long reads**, **PacBio HiFi**, and **PacBio CLR**. Autotunes its parameters from the data — you do not need to know the read length, platform, or sensible bbduk thresholds.

---

## Quick start

```bash
# Conda install (after bioconda submission)
conda install -c bioconda cerberus-mg

# Or from source
git clone https://github.com/iowa69/cerberus
cd cerberus
conda env create -f environment.yml
conda activate cerberus

# Standard short-read invocation
cerberus -r1 R1.fq.gz -r2 R2.fq.gz -o out/ -t 8 --meta --profiling --gdpr

# Long reads — single input file, --long flag
cerberus --long -i reads.fq.gz -o out/ -t 8 --all

# Just want everything at maximum quality
cerberus -r1 R1.fq.gz -r2 R2.fq.gz --all -o out/
```

On the first run, Cerberus downloads its references (~13 GB) from Zenodo into `~/.cerberus/refs/`. Use `cerberus fetch-refs` to pre-warm the cache; `cerberus doctor` to validate your installation.

---

## Why three outputs?

Existing tools force you to pick one tradeoff:
- **Hostile** is precise — keeps microbes, but leaves residual host. Bad for GDPR-restricted public release.
- **Kneaddata** is aggressive — strips host, but kills 10× more microbes than Hostile. Bad for assembly.
- **Kraken2 host filter** is k-mer based — catches what alignment misses, but doesn't satisfy mechanism-diversity reviewers asking for "0 human reads."

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
                          → bbduk (human-only 27-mer set, different k)
                          → *_GDPR.fastq.gz
```

**Why two mechanisms for GDPR?** Kraken2 alone has false negatives at the read level (it's a minimizer-based classifier optimized for taxonomy, not exclusion). bbduk's exact k-mer match catches what Kraken2 misses; using a different k value makes the two truly orthogonal.

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
- `--fast` — profiling: minimap2-only path (~3× faster, ~2% less sensitive)
- `--double-pass` — profiling: pre-filter with minimap2 before bowtie2 (rare cases)
- `--ref-dir PATH` — override `~/.cerberus/refs/`
- `--no-auto-download` — refuse to download missing references

### Advanced (--help-all)
Power-user overrides: `--min-length`, `--min-quality`, `--entropy`, `--bbduk-k`, `--minimap2-args`, `--bowtie2-args`, `--kraken2-db`, `--aux-refs`, `--keep-intermediates`.

---

## Memory & disk

| Resource | Need |
|---|---|
| RAM (peak) | ~9 GB on `--gdpr` step (Kraken2). 4-6 GB during alignment. |
| Disk (refs) | ~13 GB in `~/.cerberus/refs/` (one-time). |
| Disk (run) | ~2× input FASTQ size during processing; cleaned automatically unless `--keep-intermediates`. |

Designed for 16 GB laptops. Tested on 4-core/16 GB.

---

## Outputs

A single Cerberus run produces:

```
out/
├── <sample>.meta.R1.fastq.gz        # if --meta
├── <sample>.meta.R2.fastq.gz
├── <sample>.meta.orphans.fastq.gz
├── <sample>.profiling.fastq.gz      # if --profiling
├── <sample>.meta.R1_GDPR.fastq.gz   # if --gdpr (per source mode)
├── <sample>.meta.R2_GDPR.fastq.gz
├── <sample>.profiling.GDPR.fastq.gz
├── reports/
│   ├── accounting.tsv     # per-stage read counts (reviewer-friendly)
│   ├── accounting.json    # same, machine-readable
│   ├── fastp.json/html
│   └── *.flagstat.txt
└── logs/
    └── *.log              # one per stage, plus a JSONL run log
```

---

## License

MIT — see [LICENSE](LICENSE).

## Citation

Pre-print pending. For now: `Cerberus (v0.1.0). https://github.com/iowa69/cerberus`
