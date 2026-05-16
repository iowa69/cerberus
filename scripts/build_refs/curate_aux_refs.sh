#!/usr/bin/env bash
# Curate the auxiliary k-mer reference used by Cerberus's profiling pipeline.
#
# Composition:
#   - Human rRNA from Ensembl 113 ncRNA fasta (rRNA biotype only)
#   - Human mitochondrion (RefSeq NC_012920.1)
#   - Common low-copy human repeats (from RepeatMasker hg38 lib, filtered)
#   - Human ribosomal RNA / sno / tRNA decoys from RNAcentral
#
# Output: a single gzipped multi-fasta used by bbduk.sh ref=...

set -euo pipefail

OUT="${OUT:-build/aux_refs}"
mkdir -p "$OUT"

log() { printf '\033[36m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*"; }

# ----- 1. Ensembl ncRNA → extract rRNA -----
log "Downloading Ensembl 113 human ncRNA"
ENS_URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"
ENS="$OUT/_ensembl_ncrna.fa.gz"
[ -f "$ENS" ] || curl -L --fail -o "$ENS" "$ENS_URL"

log "Extracting rRNA biotype only"
RRNA="$OUT/_rrna.fa"
zcat "$ENS" | awk '
  /^>/ { keep = ($0 ~ /gene_biotype:rRNA/ || $0 ~ /transcript_biotype:rRNA/) }
  { if (keep) print }
' > "$RRNA"

# ----- 2. Mitochondrion -----
log "Downloading human mitochondrion (NC_012920.1)"
MT="$OUT/_mt.fa"
curl -L --fail \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text" \
  -o "$MT"

# ----- 3. Low-copy repeats (subset of RepeatMasker hg38 lib) -----
# We exclude high-copy repeats (Alu, L1) that would also match microbial sequences.
log "Building low-copy repeat subset"
REP_URL="https://www.repeatmasker.org/species/hg38_dfam_lib.fa.gz"
REP="$OUT/_repeats_raw.fa.gz"
[ -f "$REP" ] || curl -L --fail -o "$REP" "$REP_URL" || \
  log "WARNING: RepeatMasker download failed; producing aux refs without repeat masker subset"
LOWCOPY="$OUT/_lowcopy_repeats.fa"
: > "$LOWCOPY"
if [ -s "$REP" ]; then
  zcat "$REP" | awk '
    BEGIN { keep = 0 }
    /^>/ {
      keep = !(/Alu|L1|L2|LINE|SINE/);
    }
    { if (keep) print }
  ' > "$LOWCOPY"
fi

# ----- 4. RNAcentral human ncRNA decoys -----
log "Downloading RNAcentral human active fragment"
RNAC="$OUT/_rnacentral.fa.gz"
RNAC_URL="https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/by-database/refseq_active.fa.gz"
[ -f "$RNAC" ] || curl -L --fail -o "$RNAC" "$RNAC_URL" || \
  log "WARNING: RNAcentral download failed; aux refs will lack tRNA/snoRNA decoys"

# ----- 5. Concatenate -----
log "Concatenating into single aux_refs.fa.gz"
FINAL="$OUT/aux_refs.fa.gz"
{
  cat "$RRNA"
  cat "$MT"
  cat "$LOWCOPY"
  [ -s "$RNAC" ] && zcat "$RNAC" | awk '/^>HSAP/,/^>[^H]/' || true
} | gzip > "$FINAL"

# ----- 6. Manifest -----
{
  echo "{"
  echo "  \"aux_refs\": {"
  echo "    \"filename\": \"aux_refs.fa.gz\","
  printf "    \"sha256\": \"%s\",\n" "$(sha256sum "$FINAL" | awk '{print $1}')"
  printf "    \"size_bytes\": %d\n" "$(stat -c%s "$FINAL")"
  echo "  }"
  echo "}"
} > "$OUT/aux_refs_manifest_fragment.json"

# ----- 7. Also build a human-only 27-mer set for GDPR bbduk pass -----
log "Building human-only 27-mer fasta for GDPR pass"
HUMAN27="$OUT/human_k27.fa.gz"
T2T_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz"
T2T="$OUT/_t2t.fa.gz"
[ -f "$T2T" ] || curl -L --fail -o "$T2T" "$T2T_URL"
zcat "$T2T" | gzip > "$HUMAN27"

{
  echo "{"
  echo "  \"human_kmer_set\": {"
  echo "    \"filename\": \"human_k27.fa.gz\","
  printf "    \"sha256\": \"%s\",\n" "$(sha256sum "$HUMAN27" | awk '{print $1}')"
  printf "    \"size_bytes\": %d\n" "$(stat -c%s "$HUMAN27")"
  echo "  }"
  echo "}"
} > "$OUT/human_kmer_set_manifest_fragment.json"

# ----- 8. Cleanup -----
rm -f "$RRNA" "$MT" "$LOWCOPY" "$REP" "$RNAC" "$T2T" "$ENS"

log "Done."
ls -lh "$OUT"
