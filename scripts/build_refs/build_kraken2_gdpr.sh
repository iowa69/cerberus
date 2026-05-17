#!/usr/bin/env bash
# Build a compact Kraken2 DB for Cerberus's GDPR pass.
#
# Coverage: human (T2T-CHM13v2.0), chimp, gorilla, mouse, rat. The full
# Kraken2 "Standard" DB is ~50GB; this slim variant fits in <12GB RAM and
# only needs to distinguish "any mammal" from "everything else" — that's
# all the GDPR pass requires.

set -euo pipefail

THREADS="${THREADS:-$(nproc 2>/dev/null || echo 8)}"
OUT="${OUT:-build/kraken2_gdpr_compact}"
mkdir -p "$OUT"

log() { printf '\033[36m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*"; }

# ----- 1. Build taxonomy -----
log "Downloading Kraken2 taxonomy"
kraken2-build --download-taxonomy --use-ftp --threads "$THREADS" --db "$OUT"

# ----- 2. Add reference genomes -----
log "Adding human (T2T-CHM13v2.0)"
kraken2-build --add-to-library /dev/stdin --no-masking --db "$OUT" --threads "$THREADS" \
  < <(curl -sL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz | gunzip -c)

declare -A GENOMES=(
  [chimp]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/858/775/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz"
  [gorilla]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/281/585/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna.gz"
  [mouse]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz"
  [rat]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.fna.gz"
)

for name in "${!GENOMES[@]}"; do
  log "Adding $name"
  kraken2-build --add-to-library /dev/stdin --no-masking --db "$OUT" --threads "$THREADS" \
    < <(curl -sL "${GENOMES[$name]}" | gunzip -c)
done

# ----- 3. Build DB -----
log "Building Kraken2 DB (this is the long step)"
kraken2-build --build --threads "$THREADS" --db "$OUT" --kmer-len 35 --minimizer-len 31

# ----- 4. Compact / clean -----
log "Cleaning intermediate files"
kraken2-build --clean --db "$OUT"

# ----- 5. Package -----
log "Packaging"
( cd "$(dirname "$OUT")" && tar -cf - "$(basename "$OUT")" \
    | zstd -19 -T0 -o "$(basename "$OUT").tar.zst" )

# ----- 6. Manifest -----
ARCHIVE="$(dirname "$OUT")/$(basename "$OUT").tar.zst"
{
  echo "{"
  echo "  \"kraken2_gdpr_compact\": {"
  echo "    \"filename\": \"$(basename "$ARCHIVE")\","
  printf "    \"sha256\": \"%s\",\n" "$(sha256sum "$ARCHIVE" | awk '{print $1}')"
  printf "    \"size_bytes\": %d\n" "$(stat -c%s "$ARCHIVE")"
  echo "  }"
  echo "}"
} > "$(dirname "$OUT")/kraken2_gdpr_manifest_fragment.json"

log "Done. Archive: $ARCHIVE"
