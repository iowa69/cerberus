#!/usr/bin/env bash
# Build the masked T2T-CHM13v2.0 + IPD-IMGT/HLA reference for Cerberus.
#
# Strategy:
#   1. Download T2T-CHM13v2.0 from NCBI and IPD-IMGT/HLA gen.fasta from EBI.
#   2. Concatenate.
#   3. Mask repeats with bbmask (default windowed entropy).
#   4. Mask regions that match bacterial/viral pan-genomes by aligning
#      representative UHGG bacterial reps + RefSeq virus reps to T2T with
#      minimap2 (asm5) and converting hits to BED → bedtools maskfasta.
#   5. Build minimap2 (.mmi) and bowtie2 (.bt2) indexes.
#   6. Emit SHA256 sums and a manifest snippet.
#
# Run from a directory with ~80 GB free. Output: ./build/masked_t2t_hla/
#
# Tools required: curl, gunzip, samtools, minimap2, bowtie2, bbmask.sh,
# bedtools, sha256sum.

set -euo pipefail

THREADS="${THREADS:-$(nproc 2>/dev/null || echo 8)}"
OUT="${OUT:-build/masked_t2t_hla}"
WORK="${WORK:-$OUT/_work}"
mkdir -p "$OUT" "$WORK"

T2T_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz"
HLA_URL="https://github.com/ANHIG/IMGTHLA/raw/Latest/hla_gen.fasta"

# Representative bacterial set: UHGG species reps (gut-derived but broadly
# applicable). Final fasta is large (~6 GB compressed); the script below uses
# a reduced "first 5000 reps" subset for masking speed; full set is optional.
UHGG_REPS_URL="https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/uhgg_genomes_reps.fa.gz"
VIRUS_URL="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"

log() { printf '\033[36m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*"; }

# ----- 1. Download -----
log "Downloading T2T-CHM13v2.0"
[ -f "$WORK/t2t.fa.gz" ] || curl -L --fail -o "$WORK/t2t.fa.gz" "$T2T_URL"
log "Downloading IPD-IMGT/HLA"
[ -f "$WORK/hla.fa" ] || curl -L --fail -o "$WORK/hla.fa" "$HLA_URL"
log "Downloading UHGG bacterial reps"
[ -f "$WORK/uhgg.fa.gz" ] || curl -L --fail -o "$WORK/uhgg.fa.gz" "$UHGG_REPS_URL"
log "Downloading RefSeq viral genomes"
[ -f "$WORK/virus.fa.gz" ] || curl -L --fail -o "$WORK/virus.fa.gz" "$VIRUS_URL"

# ----- 2. Concatenate T2T + HLA -----
log "Concatenating T2T + HLA"
COMBINED="$WORK/t2t_hla.fa"
if [ ! -f "$COMBINED" ]; then
  gunzip -c "$WORK/t2t.fa.gz" > "$COMBINED"
  cat "$WORK/hla.fa" >> "$COMBINED"
fi
samtools faidx "$COMBINED"

# ----- 3. Mask repeats with bbmask -----
log "Masking repeats with bbmask (entropy + tandem)"
ENTROPY_MASKED="$WORK/t2t_hla.entropy_masked.fa"
if [ ! -f "$ENTROPY_MASKED" ]; then
  bbmask.sh in="$COMBINED" out="$ENTROPY_MASKED" \
    entropy=0.7 window=80 maskrepeats=t mintandem=10 \
    threads="$THREADS" -Xmx32g
fi

# ----- 4. Mask bacterial/viral hits -----
log "Mapping bacterial reps onto reference (asm5)"
BAC_PAF="$WORK/bacterial_hits.paf"
[ -f "$BAC_PAF" ] || minimap2 -x asm5 -t "$THREADS" "$ENTROPY_MASKED" "$WORK/uhgg.fa.gz" > "$BAC_PAF"

log "Mapping viral reps onto reference (asm5)"
VIR_PAF="$WORK/viral_hits.paf"
[ -f "$VIR_PAF" ] || minimap2 -x asm5 -t "$THREADS" "$ENTROPY_MASKED" "$WORK/virus.fa.gz" > "$VIR_PAF"

log "Converting hits to BED and merging"
HITS_BED="$WORK/microbial_hits.bed"
awk -v OFS='\t' '$10/$11 >= 0.85 {print $6, $8, $9}' "$BAC_PAF" "$VIR_PAF" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -d 50 -i - > "$HITS_BED"

log "Masking microbial-homologous regions in reference"
FINAL_FA="$OUT/masked_t2t_hla.fa"
bedtools maskfasta -fi "$ENTROPY_MASKED" -bed "$HITS_BED" -fo "$FINAL_FA"
samtools faidx "$FINAL_FA"

# ----- 5. Build indexes -----
log "Building minimap2 index (.mmi)"
minimap2 -d "$OUT/masked_t2t_hla.mmi" "$FINAL_FA"

log "Building bowtie2 index"
mkdir -p "$OUT/masked_t2t_hla_bt2"
bowtie2-build --threads "$THREADS" "$FINAL_FA" "$OUT/masked_t2t_hla_bt2/masked_t2t_hla"

log "Packaging bowtie2 index"
( cd "$OUT" && tar -cf - masked_t2t_hla_bt2 | zstd -19 -T0 -o masked_t2t_hla_bt2.tar.zst )
rm -rf "$OUT/masked_t2t_hla_bt2"

# ----- 6. Hashes + manifest snippet -----
log "Computing SHA256 and manifest fragment"
{
  echo "{"
  echo "  \"masked_t2t_hla_minimap2\": {"
  echo "    \"filename\": \"masked_t2t_hla.mmi\","
  printf "    \"sha256\": \"%s\",\n" "$(sha256sum "$OUT/masked_t2t_hla.mmi" | awk '{print $1}')"
  printf "    \"size_bytes\": %d\n" "$(stat -c%s "$OUT/masked_t2t_hla.mmi")"
  echo "  },"
  echo "  \"masked_t2t_hla_bowtie2\": {"
  echo "    \"filename\": \"masked_t2t_hla_bt2.tar.zst\","
  printf "    \"sha256\": \"%s\",\n" "$(sha256sum "$OUT/masked_t2t_hla_bt2.tar.zst" | awk '{print $1}')"
  printf "    \"size_bytes\": %d\n" "$(stat -c%s "$OUT/masked_t2t_hla_bt2.tar.zst")"
  echo "  }"
  echo "}"
} > "$OUT/manifest_fragment.json"

log "Done. Outputs:"
ls -lh "$OUT"
log "Manifest fragment: $OUT/manifest_fragment.json"
