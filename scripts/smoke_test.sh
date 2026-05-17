#!/usr/bin/env bash
# End-to-end smoke test for Cerberus.
#
# 1. Simulate tiny FASTQs: ~5000 reads from T2T human, ~5000 from viral RefSeq.
# 2. Run cerberus --all against the local reference cache.
# 3. Verify the three output files exist and are non-empty.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$REPO_DIR/scripts/build_refs/build"
TEST_DIR="$REPO_DIR/scripts/build_refs/build/_smoke"
REF_DIR="${CERBERUS_REF_DIR:-$HOME/.cerberus/refs}"

mkdir -p "$TEST_DIR"

log() { printf '\033[36m[smoke] %s\033[0m %s\n' "$(date +%H:%M:%S)" "$*"; }

# ----- 0. Activate env -----
source /home/iowa/miniconda3/etc/profile.d/conda.sh
conda activate cerberus

# ----- 1. Simulate reads -----
if [ ! -f "$TEST_DIR/mix.R1.fq.gz" ]; then
  log "Simulating reads with bbmap randomreads.sh"

  # 3000 paired human reads
  randomreads.sh ref="$BUILD/masked_t2t_hla/_work/t2t.fa.gz" \
    out1="$TEST_DIR/human.R1.fq.gz" out2="$TEST_DIR/human.R2.fq.gz" \
    reads=3000 length=150 paired=t illuminanames=t \
    -Xmx2g 2>/dev/null || \
  randomreads.sh ref="$REF_DIR/human_k27.fa.gz" \
    out1="$TEST_DIR/human.R1.fq.gz" out2="$TEST_DIR/human.R2.fq.gz" \
    reads=3000 length=150 paired=t illuminanames=t \
    -Xmx2g

  # 3000 paired viral reads — use RefSeq viral as a stand-in for "microbial"
  if [ -f "$BUILD/masked_t2t_hla/_work/virus.fa.gz" ]; then
    VIRUS_REF="$BUILD/masked_t2t_hla/_work/virus.fa.gz"
  else
    # Re-download tiny viral subset if work dir is gone
    VIRUS_REF="$TEST_DIR/virus.fa.gz"
    [ -f "$VIRUS_REF" ] || curl -L --fail -o "$VIRUS_REF" \
      "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"
  fi

  randomreads.sh ref="$VIRUS_REF" \
    out1="$TEST_DIR/virus.R1.fq.gz" out2="$TEST_DIR/virus.R2.fq.gz" \
    reads=3000 length=150 paired=t illuminanames=t \
    -Xmx2g

  log "Combining into mix (50:50 human:viral)"
  cat "$TEST_DIR/human.R1.fq.gz" "$TEST_DIR/virus.R1.fq.gz" > "$TEST_DIR/mix.R1.fq.gz"
  cat "$TEST_DIR/human.R2.fq.gz" "$TEST_DIR/virus.R2.fq.gz" > "$TEST_DIR/mix.R2.fq.gz"
fi

ls -lh "$TEST_DIR/mix.R1.fq.gz" "$TEST_DIR/mix.R2.fq.gz"

# ----- 2. Run cerberus --all -----
log "Running cerberus --all"
OUT="$TEST_DIR/cerberus_out"
rm -rf "$OUT"
cerberus \
  -r1 "$TEST_DIR/mix.R1.fq.gz" \
  -r2 "$TEST_DIR/mix.R2.fq.gz" \
  -o "$OUT" \
  -s smoke \
  --all \
  --threads 4 \
  --memory 8G \
  --ref-dir "$REF_DIR"

# ----- 3. Verify outputs -----
log "Verifying outputs"
EXPECTED=(
  "smoke.meta.R1.fastq.gz"
  "smoke.meta.R2.fastq.gz"
  "smoke.profiling.fastq.gz"
  "smoke.meta.R1_GDPR.fastq.gz"
  "smoke.meta.R2_GDPR.fastq.gz"
  "smoke.profiling.GDPR.fastq.gz"
)
FAIL=0
for f in "${EXPECTED[@]}"; do
  p="$OUT/$f"
  if [ -s "$p" ]; then
    n=$(zcat "$p" | wc -l)
    reads=$((n / 4))
    printf "  \033[32m✓\033[0m %-32s  %s reads\n" "$f" "$reads"
  else
    printf "  \033[31m✗\033[0m %-32s  MISSING or EMPTY\n" "$f"
    FAIL=$((FAIL + 1))
  fi
done

log "Accounting report:"
cat "$OUT/reports/accounting.tsv" 2>/dev/null || echo "  (no accounting.tsv produced)"

if [ "$FAIL" -gt 0 ]; then
  log "FAIL: $FAIL output(s) missing"
  exit 1
fi
log "PASS: all 6 expected outputs present"
