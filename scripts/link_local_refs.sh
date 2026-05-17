#!/usr/bin/env bash
# Place the freshly built reference artifacts in ~/.cerberus/refs/ so Cerberus
# can use them without re-downloading from Zenodo. Useful for local end-to-end
# testing before the Zenodo deposition is published.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$REPO_DIR/scripts/build_refs/build"
DST="${CERBERUS_REF_DIR:-$HOME/.cerberus/refs}"

if [ ! -d "$BUILD" ]; then
  echo "Build directory not found: $BUILD" >&2
  exit 1
fi

mkdir -p "$DST"

link_or_copy() {
  local src="$1" dst="$2"
  if [ ! -e "$src" ]; then
    echo "  ⚠ Source missing: $src" >&2
    return
  fi
  rm -f "$dst"
  ln -sf "$src" "$dst"
  echo "  ✓ $(basename "$dst")  →  $src"
}

extract_archive() {
  local archive="$1" extracted_dir="$2"
  if [ -d "$extracted_dir" ] && [ -n "$(ls -A "$extracted_dir" 2>/dev/null)" ]; then
    echo "  ✓ $(basename "$extracted_dir")/  already extracted"
    return
  fi
  mkdir -p "$extracted_dir"
  echo "  ⇪ Extracting $(basename "$archive")"
  zstd -dc "$archive" | tar -x -C "$extracted_dir"
}

echo "Linking refs into $DST"
link_or_copy "$BUILD/masked_t2t_hla/masked_t2t_hla.mmi" "$DST/masked_t2t_hla.mmi"
link_or_copy "$BUILD/aux_refs/aux_refs.fa.gz"            "$DST/aux_refs.fa.gz"
link_or_copy "$BUILD/aux_refs/human_k27.fa.gz"           "$DST/human_k27.fa.gz"

echo
echo "Extracting archives into $DST"
extract_archive "$BUILD/masked_t2t_hla/masked_t2t_hla_bt2.tar.zst" "$DST/masked_t2t_hla_bt2"
extract_archive "$BUILD/kraken2_gdpr_compact.tar.zst"              "$DST/kraken2_gdpr_compact"

echo
echo "Writing a placeholder manifest so RefManager treats files as satisfied"
cat > "$DST/manifest.json" <<'EOF'
{
  "schema_version": 1,
  "release": "0.1.0-local",
  "zenodo_doi": "LOCAL",
  "notes": "Local-build manifest. Hashes empty so RefManager skips verification.",
  "assets": {
    "masked_t2t_hla_minimap2": {"description":"","filename":"masked_t2t_hla.mmi","url":"","sha256":"","size_bytes":null,"required_for":["meta","profiling-fast"]},
    "masked_t2t_hla_bowtie2":  {"description":"","filename":"masked_t2t_hla_bt2.tar.zst","url":"","sha256":"","size_bytes":null,"required_for":["profiling"]},
    "kraken2_gdpr_compact":    {"description":"","filename":"kraken2_gdpr_compact.tar.zst","url":"","sha256":"","size_bytes":null,"required_for":["gdpr"]},
    "aux_refs":                {"description":"","filename":"aux_refs.fa.gz","url":"","sha256":"","size_bytes":null,"required_for":["profiling"]},
    "human_kmer_set":          {"description":"","filename":"human_k27.fa.gz","url":"","sha256":"","size_bytes":null,"required_for":["gdpr"]}
  }
}
EOF

echo
echo "Done. Set CERBERUS_REF_DIR=$DST or pass --ref-dir to cerberus."
ls -la "$DST"
