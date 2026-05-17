#!/usr/bin/env bash
# Curate the auxiliary k-mer reference used by Cerberus's profiling pipeline.
#
# v0.1 composition (sources that we have reliable public URLs for):
#   - Human ncRNA from Ensembl 113: rRNA, Mt_rRNA, Mt_tRNA, snRNA, snoRNA,
#     scaRNA, miRNA, Y_RNA, vaultRNA, ribozyme  (everything that could
#     introduce host-derived k-mers into a metagenomic sample)
#   - Human mitochondrion (NC_012920.1 from NCBI eutils)
#
# Deferred to v0.2 (URLs not currently reliable):
#   - RepeatMasker hg38 low-copy subset
#   - RNAcentral human active fragment
# Both contribute marginal value; the rRNA + ncRNA + mtDNA already captures
# the dominant false-positive sources for alignment-based host removal.
#
# Also builds: human_k27.fa.gz = full T2T-CHM13v2.0, used by the GDPR pass
# bbduk run as a second orthogonal mechanism after Kraken2.

set -euo pipefail

OUT="${OUT:-build/aux_refs}"
mkdir -p "$OUT"

log() { printf '\033[36m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*"; }

# ----- 1. Ensembl ncRNA → extract host-decoy biotypes -----
log "Downloading Ensembl 113 human ncRNA"
ENS_URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"
ENS="$OUT/_ensembl_ncrna.fa.gz"
[ -f "$ENS" ] || curl -L --fail -o "$ENS" "$ENS_URL"

log "Extracting host-decoy biotypes (rRNA/snRNA/snoRNA/miRNA/Y_RNA/etc.)"
NCRNA_FA="$OUT/_ncrna.fa"
zcat "$ENS" | awk '
  BEGIN { keep = 0 }
  /^>/ {
    keep = 0
    if ($0 ~ /gene_biotype:rRNA/) keep = 1
    else if ($0 ~ /gene_biotype:Mt_rRNA/) keep = 1
    else if ($0 ~ /gene_biotype:Mt_tRNA/) keep = 1
    else if ($0 ~ /gene_biotype:snRNA/) keep = 1
    else if ($0 ~ /gene_biotype:snoRNA/) keep = 1
    else if ($0 ~ /gene_biotype:scaRNA/) keep = 1
    else if ($0 ~ /gene_biotype:miRNA/) keep = 1
    else if ($0 ~ /gene_biotype:Y_RNA/) keep = 1
    else if ($0 ~ /gene_biotype:vault_RNA/) keep = 1
    else if ($0 ~ /gene_biotype:ribozyme/) keep = 1
    else if ($0 ~ /transcript_biotype:rRNA/) keep = 1
  }
  { if (keep) print }
' > "$NCRNA_FA"
NCRNA_COUNT=$(grep -c '^>' "$NCRNA_FA" || echo 0)
log "Captured $NCRNA_COUNT host-decoy ncRNA sequences"

# ----- 2. Mitochondrion -----
log "Downloading human mitochondrion (NC_012920.1)"
MT="$OUT/_mt.fa"
curl -L --fail \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text" \
  -o "$MT"

# ----- 3. Concatenate -----
log "Concatenating into single aux_refs.fa.gz"
FINAL="$OUT/aux_refs.fa.gz"
cat "$NCRNA_FA" "$MT" | gzip > "$FINAL"
FINAL_RECORDS=$(zcat "$FINAL" | grep -c '^>' || echo 0)
log "aux_refs.fa.gz contains $FINAL_RECORDS records ($(du -h "$FINAL" | cut -f1))"

# ----- 4. Aux refs manifest fragment -----
{
  echo "{"
  echo "  \"aux_refs\": {"
  echo "    \"filename\": \"aux_refs.fa.gz\","
  printf "    \"sha256\": \"%s\",\n" "$(sha256sum "$FINAL" | awk '{print $1}')"
  printf "    \"size_bytes\": %d,\n" "$(stat -c%s "$FINAL")"
  printf "    \"record_count\": %d\n" "$FINAL_RECORDS"
  echo "  }"
  echo "}"
} > "$OUT/aux_refs_manifest_fragment.json"

# ----- 5. Build human k27 fasta for the GDPR bbduk pass -----
log "Building human-only k27 reference fasta (full T2T-CHM13v2.0)"
HUMAN27="$OUT/human_k27.fa.gz"
T2T_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz"
T2T="$OUT/_t2t.fa.gz"
[ -f "$T2T" ] || curl -L --fail -o "$T2T" "$T2T_URL"
cp "$T2T" "$HUMAN27"

{
  echo "{"
  echo "  \"human_kmer_set\": {"
  echo "    \"filename\": \"human_k27.fa.gz\","
  printf "    \"sha256\": \"%s\",\n" "$(sha256sum "$HUMAN27" | awk '{print $1}')"
  printf "    \"size_bytes\": %d\n" "$(stat -c%s "$HUMAN27")"
  echo "  }"
  echo "}"
} > "$OUT/human_kmer_set_manifest_fragment.json"

# ----- 6. Cleanup -----
rm -f "$NCRNA_FA" "$MT" "$ENS" "$T2T"

log "Done."
ls -lh "$OUT"
