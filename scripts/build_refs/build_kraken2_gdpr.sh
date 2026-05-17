#!/usr/bin/env bash
# Build a compact Kraken2 DB for Cerberus's GDPR pass.
#
# Coverage: human (T2T-CHM13v2.0), chimp, gorilla, mouse, rat. The full
# Kraken2 Standard DB is ~50 GB; this slim variant fits in <12 GB RAM and
# only needs to distinguish "any mammal" from "everything else" — which is
# all the GDPR pass requires.
#
# Strategy: download the NCBI taxonomy dump (~80 MB) directly, skip the
# 7 GB accession-to-taxid maps by injecting taxonomy IDs into the FASTA
# headers using the `>kraken:taxid|TAXID|...` convention that kraken2-build
# understands natively.

set -euo pipefail

THREADS="${THREADS:-$(nproc 2>/dev/null || echo 8)}"
DB="${DB:-build/kraken2_gdpr_compact}"
mkdir -p "$DB/taxonomy"

log() { printf '\033[36m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*"; }

# ----- 1. Taxonomy dump -----
if [ ! -f "$DB/taxonomy/nodes.dmp" ]; then
  log "Downloading NCBI taxdump (~80 MB)"
  curl -L --fail -o "$DB/taxonomy/taxdump.tar.gz" \
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
  tar -xzf "$DB/taxonomy/taxdump.tar.gz" -C "$DB/taxonomy"
  rm "$DB/taxonomy/taxdump.tar.gz"
fi

# kraken2-build expects these files to exist even when using kraken:taxid headers
touch "$DB/taxonomy/nucl_gb.accession2taxid" "$DB/taxonomy/nucl_wgs.accession2taxid"

# ----- 2. Add reference genomes with embedded taxonomy IDs -----
add_genome() {
  local name="$1" taxid="$2" url="$3"
  if [ -f "$DB/library/added/.${name}.done" ]; then
    log "Skipping $name (already added)"
    return
  fi
  log "Adding $name (taxid=$taxid)"
  local tmp_fa="$DB/_${name}.tagged.fna"
  curl -sL --fail "$url" | gunzip -c \
    | awk -v tid="$taxid" '
        /^>/ { print ">kraken:taxid|" tid "|" substr($0, 2); next }
        { print }
      ' > "$tmp_fa"
  kraken2-build --add-to-library "$tmp_fa" --no-masking \
    --db "$DB" --threads "$THREADS"
  rm "$tmp_fa"
  mkdir -p "$DB/library/added"
  touch "$DB/library/added/.${name}.done"
}

add_genome human   9606  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz"
add_genome chimp   9598  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/858/775/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz"
add_genome gorilla 9593  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/281/585/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna.gz"
add_genome mouse   10090 "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz"
add_genome rat     10116 "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.fna.gz"

# ----- 3. Build the DB -----
log "Building Kraken2 DB (k=35, minimizer=31; ~30-45 min CPU)"
kraken2-build --build --threads "$THREADS" --db "$DB" \
  --kmer-len 35 --minimizer-len 31

# ----- 4. Cleanup intermediate libraries -----
log "Pruning intermediates with kraken2-build --clean"
kraken2-build --clean --db "$DB"

# ----- 5. Package -----
log "Packaging compact DB as .tar.zst"
ARCHIVE="$(dirname "$DB")/$(basename "$DB").tar.zst"
( cd "$(dirname "$DB")" && tar -cf - "$(basename "$DB")" \
    | zstd -19 -T0 -o "$(basename "$ARCHIVE")" )

# ----- 6. Manifest -----
{
  echo "{"
  echo "  \"kraken2_gdpr_compact\": {"
  echo "    \"filename\": \"$(basename "$ARCHIVE")\","
  printf "    \"sha256\": \"%s\",\n" "$(sha256sum "$ARCHIVE" | awk '{print $1}')"
  printf "    \"size_bytes\": %d\n" "$(stat -c%s "$ARCHIVE")"
  echo "  }"
  echo "}"
} > "$(dirname "$DB")/kraken2_gdpr_manifest_fragment.json"

log "Done. Archive: $ARCHIVE ($(du -h "$ARCHIVE" | cut -f1))"
