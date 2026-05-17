#!/usr/bin/env bash
# Run all three Cerberus reference builds in sequence with unified logging.
# Each script cleans its own _work/ on success; failure short-circuits the chain.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Activate conda env
source /home/iowa/miniconda3/etc/profile.d/conda.sh
conda activate cerberus

LOG_DIR="$SCRIPT_DIR/build/_run_logs"
mkdir -p "$LOG_DIR"

ts() { date +'%Y-%m-%d %H:%M:%S'; }
banner() { printf '\n\033[1;36m=== [%s] %s ===\033[0m\n\n' "$(ts)" "$*"; }

banner "1/3 Curate auxiliary references (~10 GB transient, ~30 min)"
bash curate_aux_refs.sh 2>&1 | tee "$LOG_DIR/01_aux_refs.log"
df -h /home/iowa | tail -1

banner "2/3 Build Kraken2 GDPR compact DB (~15 GB transient, ~30-60 min)"
bash build_kraken2_gdpr.sh 2>&1 | tee "$LOG_DIR/02_kraken2_gdpr.log"
df -h /home/iowa | tail -1

banner "3/3 Build masked T2T+HLA indexes (~30 GB transient, 1-3 h)"
bash mask_t2t_hla.sh 2>&1 | tee "$LOG_DIR/03_mask_t2t_hla.log"
df -h /home/iowa | tail -1

banner "DONE. All three reference asset groups built."
ls -lh build/aux_refs/ build/kraken2_gdpr_compact.tar.zst build/masked_t2t_hla/ 2>/dev/null

banner "Manifest fragments (to be merged into cerberus/data/default_manifest.json after Zenodo upload):"
for f in build/*/manifest_fragment.json build/*manifest_fragment.json; do
    [ -f "$f" ] && echo "=== $f ===" && cat "$f"
done
