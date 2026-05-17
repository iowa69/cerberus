#!/usr/bin/env bash
# Build a Cerberus-compatible host reference directory from any FASTA.
#
# Produces a self-contained directory ready for ``cerberus --ref-dir <DIR>``:
#   - minimap2 short-read index   (sr preset)
#   - minimap2 long-read index    (default preset, suits ONT and PacBio)
#   - bowtie2 short-read index    (for the --profiling pass)
#   - manifest.json that satisfies RefManager
#   - optional copy of an existing Kraken2 DB for --gdpr
#
# Use cases:
#   - Mouse, pig, plant, or any non-human host
#   - A masked variant of T2T-CHM13 you built yourself
#   - A bespoke decoy reference (e.g. host + suspected contaminant)

set -euo pipefail

INPUT=""
OUT=""
NAME="masked_t2t_hla"
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 8)}"
DO_MASK=0
PLATFORM="both"
KRAKEN_DB=""
AUX_REF=""

usage() {
  cat <<USAGE
Usage: $0 -i FASTA -o OUT_DIR [options]

Build Cerberus-compatible host reference indexes from any FASTA file.

Required:
  -i, --input FASTA      Host genome (.fa, .fa.gz, .fna, .fna.gz)
  -o, --out DIR          Output directory (created if absent)

Options:
  -n, --name NAME        Index name prefix (default: masked_t2t_hla — keeps
                         the canonical Cerberus filenames so --ref-dir works).
  -t, --threads N        Threads (default: nproc)
  --mask                 Run bbmask repeat/entropy masking before indexing.
                         Recommended for any reference >100 Mb.
  --platform P           short | long | both  (default: both)
                         "short" builds minimap2 sr + bowtie2.
                         "long"  builds minimap2 default-preset only.
  --kraken-db PATH       Existing Kraken2 DB (a directory containing hash.k2d
                         or a .tar.zst archive) to copy/link for --gdpr.
                         If omitted, --gdpr will not be available.
  --aux-refs PATH        Optional auxiliary k-mer reference FASTA(.gz) used
                         by the profiling stage's bbduk k-mer pass.

Examples:
  # Mouse host removal
  $0 -i mouse_GRCm39.fa.gz -o ~/cerberus_mouse --mask -t 16

  # Just a long-read setup, no GDPR
  $0 -i ref.fa -o ./refs --platform long

  # Custom host but reuse the published Kraken2 GDPR DB (for human scrubbing)
  $0 -i custom.fa -o ./refs --kraken-db ~/.cerberus/refs/kraken2_gdpr_compact
USAGE
}

while [ $# -gt 0 ]; do
  case "$1" in
    -i|--input)     INPUT="$2"; shift 2 ;;
    -o|--out)       OUT="$2"; shift 2 ;;
    -n|--name)      NAME="$2"; shift 2 ;;
    -t|--threads)   THREADS="$2"; shift 2 ;;
    --mask)         DO_MASK=1; shift ;;
    --platform)     PLATFORM="$2"; shift 2 ;;
    --kraken-db)    KRAKEN_DB="$2"; shift 2 ;;
    --aux-refs)     AUX_REF="$2"; shift 2 ;;
    -h|--help)      usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[ -z "$INPUT" ] && { echo "ERROR: -i required" >&2; usage; exit 1; }
[ -z "$OUT" ]   && { echo "ERROR: -o required" >&2; usage; exit 1; }
[ -f "$INPUT" ] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }
case "$PLATFORM" in short|long|both) ;; *)
  echo "ERROR: --platform must be short|long|both" >&2; exit 1 ;;
esac

mkdir -p "$OUT"
WORK="$OUT/_work"
mkdir -p "$WORK"

log() { printf '\033[36m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*"; }

# ----- 1. Prepare uncompressed reference -----
if [[ "$INPUT" == *.gz ]]; then
  REF_FA="$WORK/host.fa"
  [ -f "$REF_FA" ] || { log "Decompressing $INPUT"; gunzip -c "$INPUT" > "$REF_FA"; }
else
  REF_FA="$INPUT"
fi

# ----- 2. Optional masking -----
if [ "$DO_MASK" -eq 1 ]; then
  log "bbmask repeats + entropy filtering"
  bbmask.sh in="$REF_FA" out="$WORK/host.masked.fa" \
    entropy=0.7 window=80 \
    maskrepeats=t kr=5 minlen=40 mincount=4 \
    threads="$THREADS" -Xmx12g
  REF_FA="$WORK/host.masked.fa"
fi
samtools faidx "$REF_FA"

# ----- 3. Short-read indexes -----
if [ "$PLATFORM" = "short" ] || [ "$PLATFORM" = "both" ]; then
  log "Building minimap2 short-read index ($NAME.mmi)"
  minimap2 -x sr -d "$OUT/$NAME.mmi" "$REF_FA"

  log "Building bowtie2 index ($NAME)"
  mkdir -p "$OUT/${NAME}_bt2"
  bowtie2-build --threads "$THREADS" "$REF_FA" "$OUT/${NAME}_bt2/$NAME"
fi

# ----- 4. Long-read index -----
if [ "$PLATFORM" = "long" ] || [ "$PLATFORM" = "both" ]; then
  log "Building minimap2 long-read index ($NAME.long.mmi)"
  minimap2 -d "$OUT/$NAME.long.mmi" "$REF_FA"
  if [ "$PLATFORM" = "long" ] && [ ! -f "$OUT/$NAME.mmi" ]; then
    log "Symlinking $NAME.mmi -> $NAME.long.mmi for asset compatibility"
    ln -sf "$NAME.long.mmi" "$OUT/$NAME.mmi"
  fi
fi

# ----- 5. Optional Kraken2 DB integration -----
if [ -n "$KRAKEN_DB" ]; then
  KDB_DST="$OUT/kraken2_gdpr_compact"
  if [ -d "$KRAKEN_DB" ]; then
    log "Linking Kraken2 DB from $KRAKEN_DB"
    ln -sfn "$(readlink -f "$KRAKEN_DB")" "$KDB_DST"
  elif [ -f "$KRAKEN_DB" ]; then
    log "Extracting Kraken2 DB archive $KRAKEN_DB"
    mkdir -p "$KDB_DST"
    case "$KRAKEN_DB" in
      *.tar.zst) zstd -dc "$KRAKEN_DB" | tar -x -C "$KDB_DST" ;;
      *.tar.gz)  tar -xzf "$KRAKEN_DB" -C "$KDB_DST" ;;
      *.tar)     tar -xf  "$KRAKEN_DB" -C "$KDB_DST" ;;
      *) echo "ERROR: unknown Kraken2 DB archive format: $KRAKEN_DB" >&2; exit 1 ;;
    esac
  else
    echo "WARNING: --kraken-db path not found: $KRAKEN_DB" >&2
  fi
fi

# ----- 6. Auxiliary refs -----
if [ -n "$AUX_REF" ] && [ -f "$AUX_REF" ]; then
  log "Copying auxiliary refs"
  cp "$AUX_REF" "$OUT/aux_refs.fa.gz"
else
  log "No --aux-refs provided; writing empty placeholder"
  : | gzip > "$OUT/aux_refs.fa.gz"
fi

# Placeholder for the GDPR k-mer set; not strictly needed since Cerberus
# uses minimap2 for the GDPR second-mechanism, but the manifest still lists it.
[ -f "$OUT/human_k27.fa.gz" ] || ln -sf "aux_refs.fa.gz" "$OUT/human_k27.fa.gz"

# ----- 7. Manifest -----
log "Writing manifest.json"
cat > "$OUT/manifest.json" <<EOF
{
  "schema_version": 1,
  "release": "custom-local",
  "zenodo_doi": "LOCAL",
  "notes": "Custom-built reference from $(basename "$INPUT"). Hashes empty; RefManager will skip verification.",
  "assets": {
    "masked_t2t_hla_minimap2": {"description":"","filename":"$NAME.mmi","url":"","sha256":"","size_bytes":null,"required_for":["meta","profiling-fast"]},
    "masked_t2t_hla_bowtie2":  {"description":"","filename":"${NAME}_bt2","url":"","sha256":"","size_bytes":null,"required_for":["profiling"]},
    "kraken2_gdpr_compact":    {"description":"","filename":"kraken2_gdpr_compact","url":"","sha256":"","size_bytes":null,"required_for":["gdpr"]},
    "aux_refs":                {"description":"","filename":"aux_refs.fa.gz","url":"","sha256":"","size_bytes":null,"required_for":["profiling"]},
    "human_kmer_set":          {"description":"","filename":"human_k27.fa.gz","url":"","sha256":"","size_bytes":null,"required_for":["gdpr"]}
  }
}
EOF

# ----- 8. Cleanup -----
rm -rf "$WORK"

log "Done."
echo
echo "Use:"
echo "    cerberus --ref-dir $OUT  -r1 R1.fq.gz -r2 R2.fq.gz -o out/ --all"
echo
ls -lh "$OUT" | head -20
