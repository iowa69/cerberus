#!/usr/bin/env bash
# Run only kraken2_gdpr + mask_t2t_hla (aux_refs already done).

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

source /home/iowa/miniconda3/etc/profile.d/conda.sh
conda activate cerberus

LOG_DIR="$SCRIPT_DIR/build/_run_logs"
mkdir -p "$LOG_DIR"

ts() { date +'%Y-%m-%d %H:%M:%S'; }
banner() { printf '\n\033[1;36m=== [%s] %s ===\033[0m\n\n' "$(ts)" "$*"; }

banner "1/2 Build Kraken2 GDPR compact DB (~5 GB downloads, ~45 min)"
bash build_kraken2_gdpr.sh 2>&1 | tee "$LOG_DIR/02_kraken2_gdpr.log"
df -h /home/iowa | tail -1

banner "2/2 Build masked T2T+HLA indexes (~5 GB downloads, ~1-2 h)"
bash mask_t2t_hla.sh 2>&1 | tee "$LOG_DIR/03_mask_t2t_hla.log"
df -h /home/iowa | tail -1

banner "DONE."
ls -lh build/aux_refs/ build/kraken2_gdpr_compact.tar.zst build/masked_t2t_hla/ 2>/dev/null

banner "Manifest fragments:"
for f in build/*/manifest_fragment.json build/*manifest_fragment.json build/aux_refs/*_fragment.json; do
    [ -f "$f" ] && echo "=== $f ===" && cat "$f"
done
