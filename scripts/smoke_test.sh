#!/usr/bin/env bash
# End-to-end smoke test.
#
# Builds miniature but REAL references from synthetic genomes (minimap2 index,
# bowtie2 index, Kraken2 database, auxiliary FASTA) and runs every Cerberus
# mode against reads of known provenance. Because each read's origin is encoded
# in its name, the test asserts that host reads were removed and microbial
# reads were kept — not merely that the pipeline exited 0.
#
# Runs in well under a minute and needs no network access.
#
# Usage:  bash scripts/smoke_test.sh [workdir]
set -euo pipefail

WORK="${1:-$(mktemp -d -t cerberus-smoke-XXXXXX)}"
REFS="$WORK/refs"
FIX="$WORK/fixture"
PASS=0
FAIL=0

say()  { printf '\n\033[1m== %s\033[0m\n' "$*"; }
ok()   { printf '  \033[32mPASS\033[0m %s\n' "$*"; PASS=$((PASS + 1)); }
bad()  { printf '  \033[31mFAIL\033[0m %s\n' "$*"; FAIL=$((FAIL + 1)); }
check() { if [ "$2" = "$3" ]; then ok "$1 ($2)"; else bad "$1: expected $3, got $2"; fi; }

for tool in cerberus minimap2 bowtie2-build kraken2-build bbduk.sh samtools fastp python3; do
  command -v "$tool" >/dev/null || { echo "missing required tool: $tool" >&2; exit 1; }
done

reads() { python3 -c "
import gzip,sys
n=0
with gzip.open(sys.argv[1],'rt') as f:
    for i,_ in enumerate(f): n=i+1
print(n//4)" "$1"; }

names_matching() { python3 -c "
import gzip,sys
pat=sys.argv[2]; n=0
with gzip.open(sys.argv[1],'rt') as f:
    for i,line in enumerate(f):
        if i%4==0 and pat in line: n+=1
print(n)" "$1" "$2"; }

say "Building fixture in $WORK"
mkdir -p "$FIX" "$REFS"
python3 - "$FIX" <<'PYEOF'
import gzip, random, sys
from pathlib import Path

out = Path(sys.argv[1])
rng = random.Random(20260803)

def genome(n, gc):
    at = (1 - gc) / 2
    return "".join(rng.choices("ACGT", weights=[at, gc/2, gc/2, at], k=n))

def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def mutate(s, rate):
    return "".join(rng.choice([c for c in "ACGT" if c != b]) if rng.random() < rate else b
                   for b in s)

def fasta(path, name, seq):
    with path.open("w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 70):
            f.write(seq[i:i+70] + "\n")

host, microbe = genome(600_000, 0.41), genome(400_000, 0.57)
fasta(out / "host.fa", "host_chr1", host)
fasta(out / "microbe.fa", "microbe_contig1", microbe)
fasta(out / "host_kraken.fa", "host_chr1|kraken:taxid|9606", host)
fasta(out / "aux_refs.fa", "host_aux_region", host[10_000:14_000])
with open(out / "aux_refs.fa", "rb") as fin, gzip.open(out / "aux_refs.fa.gz", "wb") as fo:
    fo.write(fin.read())

pairs = []
for src, label, n in ((host, "HOST", 2000), (microbe, "MICROBE", 2000)):
    for i in range(n):
        s = rng.randrange(0, len(src) - 350)
        frag = src[s:s+350]
        pairs.append((f"{label}_{i}",
                      mutate(frag[:150], 0.002), mutate(rc(frag[-150:]), 0.002)))
# pairs where exactly one mate is host: these separate the two drop strategies
for i in range(500):
    s = rng.randrange(0, len(host) - 150)
    t = rng.randrange(0, len(microbe) - 150)
    pairs.append((f"HALF_{i}", mutate(host[s:s+150], 0.002),
                  mutate(rc(microbe[t:t+150]), 0.002)))
rng.shuffle(pairs)

with gzip.open(out / "reads.R1.fq.gz", "wt") as f1, \
     gzip.open(out / "reads.R2.fq.gz", "wt") as f2:
    for name, r1, r2 in pairs:
        f1.write(f"@{name}/1\n{r1}\n+\n{'I'*len(r1)}\n")
        f2.write(f"@{name}/2\n{r2}\n+\n{'I'*len(r2)}\n")

with gzip.open(out / "reads.long.fq.gz", "wt") as f:
    for src, label in ((host, "HOST"), (microbe, "MICROBE")):
        for i in range(200):
            s = rng.randrange(0, len(src) - 6000)
            seq = mutate(src[s:s+6000], 0.05)
            f.write(f"@{label}long_{i}\n{seq}\n+\n{'5'*len(seq)}\n")
print("fixture: 2000 host + 2000 microbe + 500 half-host pairs, 400 long reads")
PYEOF

say "Building miniature references"
minimap2 -x sr -d "$REFS/masked_t2t_hla.mmi" "$FIX/host.fa" 2>/dev/null
mkdir -p "$REFS/masked_t2t_hla_bt2"
bowtie2-build --quiet --threads 2 "$FIX/host.fa" "$REFS/masked_t2t_hla_bt2/masked_t2t_hla"
cp "$FIX/aux_refs.fa.gz" "$REFS/aux_refs.fa.gz"
cp "$FIX/aux_refs.fa.gz" "$REFS/human_k27.fa.gz"

KDB="$REFS/kraken2_gdpr_compact"
mkdir -p "$KDB/taxonomy"
printf '1\t|\t1\t|\tno rank\t|\t-\t|\t\n9606\t|\t1\t|\tspecies\t|\tHS\t|\t\n' \
  > "$KDB/taxonomy/nodes.dmp"
printf '1\t|\troot\t|\t\t|\tscientific name\t|\n9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|\n' \
  > "$KDB/taxonomy/names.dmp"
kraken2-build --db "$KDB" --add-to-library "$FIX/host_kraken.fa" --no-masking >/dev/null 2>&1
kraken2-build --db "$KDB" --build --threads 2 --kmer-len 31 \
  --minimizer-len 25 --minimizer-spaces 6 >/dev/null 2>&1

python3 - "$REFS" <<'PYEOF'
import hashlib, json, sys
from pathlib import Path
refs = Path(sys.argv[1])

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def plain(name, desc, req):
    p = refs / name
    return {"description": desc, "filename": name, "url": "file://" + str(p),
            "sha256": sha(p), "size_bytes": p.stat().st_size, "required_for": req}

def extracted(name, desc, req):
    # already-extracted directories: no URL needed, RefManager finds them on disk
    return {"description": desc, "filename": name, "url": "PENDING",
            "sha256": "PENDING", "size_bytes": None, "required_for": req}

(refs / "manifest.json").write_text(json.dumps({
    "schema_version": 1, "release": "smoke-fixture",
    "notes": "Miniature synthetic references for the smoke test.",
    "assets": {
        "masked_t2t_hla_minimap2": plain("masked_t2t_hla.mmi", "mini minimap2 index",
                                         ["meta", "profiling-fast"]),
        "masked_t2t_hla_bowtie2": extracted("masked_t2t_hla_bt2.tar.zst",
                                            "mini bowtie2 index", ["profiling"]),
        "kraken2_gdpr_compact": extracted("kraken2_gdpr_compact.tar.zst",
                                          "mini kraken2 host DB", ["gdpr"]),
        "aux_refs": plain("aux_refs.fa.gz", "mini aux refs", ["profiling"]),
        "human_kmer_set": plain("human_k27.fa.gz", "mini human k-mers", ["gdpr"]),
    }}, indent=2))
PYEOF

say "cerberus doctor"
cerberus doctor --ref-dir "$REFS" >/dev/null && ok "doctor exits 0"

say "Short reads: --all"
OUT="$WORK/out_all"
cerberus -r1 "$FIX/reads.R1.fq.gz" -r2 "$FIX/reads.R2.fq.gz" \
  -o "$OUT" -t 2 --all --ref-dir "$REFS" -s SMOKE -q

check "meta R1 keeps microbial reads"   "$(names_matching "$OUT/SMOKE.meta.R1.fastq.gz" MICROBE)" 2000
check "meta R1 removes host reads"      "$(names_matching "$OUT/SMOKE.meta.R1.fastq.gz" HOST_)"      0
check "meta R1/R2 stay synchronised"    "$(reads "$OUT/SMOKE.meta.R1.fastq.gz")" \
                                        "$(reads "$OUT/SMOKE.meta.R2.fastq.gz")"
# The two drop strategies must genuinely differ. A pair with exactly one host
# mate is KEPT by meta (conservative: retain the microbial mate, accept the
# residual host) and DROPPED by profiling (aggressive). Under v0.1.1 both
# strategies compiled to the same filter and both produced 2000 here.
check "meta keeps half-host pairs"      "$(names_matching "$OUT/SMOKE.meta.R1.fastq.gz" HALF)"    500
check "meta pairs stay paired"          "$(names_matching "$OUT/SMOKE.meta.R2.fastq.gz" HALF)"    500
check "profiling drops half-host pairs" "$(names_matching "$OUT/SMOKE.profiling.fastq.gz" HALF)"    0
check "profiling removes host reads"    "$(names_matching "$OUT/SMOKE.profiling.fastq.gz" HOST_)"     0
check "GDPR R1 removes host reads"      "$(names_matching "$OUT/SMOKE.meta.R1_GDPR.fastq.gz" HOST_)"  0
# profiling's GDPR output is its deliverable, not a side stream, so it must not
# be named "orphans" the way meta's unpaired leftovers are
[ -s "$OUT/SMOKE.profiling_GDPR.fastq.gz" ] \
  && ok "profiling GDPR deliverable is named for what it is" \
  || bad "expected $OUT/SMOKE.profiling_GDPR.fastq.gz"
[ -e "$OUT/SMOKE.profiling.orphans_GDPR.fastq.gz" ] \
  && bad "profiling deliverable is still mislabelled as orphans" \
  || ok "no mislabelled profiling orphans file"

say "Output integrity"
for f in "$OUT"/*.fastq.gz; do
  if python3 -c "
import gzip,sys
with gzip.open(sys.argv[1],'rb') as fh:
    while fh.read(1<<20): pass" "$f"; then
    ok "valid gzip: $(basename "$f")"
  else
    bad "corrupt gzip: $(basename "$f")"
  fi
done

UNIQ=$(python3 -c "
import gzip
ids=[l for i,l in enumerate(gzip.open('$OUT/SMOKE.profiling.fastq.gz','rt')) if i%4==0]
print(1 if len(ids)==len(set(ids)) else 0)")
check "merged profiling read IDs are unique" "$UNIQ" 1

say "Reports"
[ -s "$OUT/reports/cerberus_report.html" ] && ok "HTML run report written" || bad "no HTML run report"
[ -s "$OUT/reports/run_record.json" ]      && ok "run_record.json written"  || bad "no run_record.json"
[ -s "$OUT/reports/accounting.tsv" ]       && ok "accounting.tsv written"   || bad "no accounting.tsv"
if grep -q "https://\|<script" "$OUT/reports/cerberus_report.html"; then
  bad "HTML report is not self-contained"
else
  ok "HTML report is self-contained"
fi
[ -d "$OUT/_work" ] && bad "_work/ was not cleaned up" || ok "_work/ cleaned up"

say "Long reads"
OUTL="$WORK/out_long"
cerberus --long -i "$FIX/reads.long.fq.gz" -o "$OUTL" -t 2 \
  --meta --profiling --ref-dir "$REFS" -s LONG -q
check "long meta removes host"      "$(names_matching "$OUTL/LONG.long_meta.fastq.gz" HOSTlong)"      0
check "long meta keeps microbial"   "$(names_matching "$OUTL/LONG.long_meta.fastq.gz" MICROBElong)" 200
check "long profiling removes host" "$(names_matching "$OUTL/LONG.long_profiling.fastq.gz" HOSTlong)" 0

say "Regression guards"
# A long-read preset on paired input used to yield empty output and exit 0.
OUTP="$WORK/out_ont"
cerberus -r1 "$FIX/reads.R1.fq.gz" -r2 "$FIX/reads.R2.fq.gz" -o "$OUTP" -t 2 \
  --meta --platform ont --ref-dir "$REFS" -s ONT -q
check "ont preset on paired reads still yields output" \
  "$(names_matching "$OUTP/ONT.meta.R1.fastq.gz" MICROBE)" 2000

# --dry-run used to crash, and to delete a previous run's outputs first.
BEFORE=$(reads "$OUT/SMOKE.meta.R1.fastq.gz")
cerberus -r1 "$FIX/reads.R1.fq.gz" -r2 "$FIX/reads.R2.fq.gz" -o "$OUT" -t 2 \
  --meta --ref-dir "$REFS" -s SMOKE --dry-run -q >/dev/null
ok "dry-run exits cleanly"
check "dry-run left existing outputs intact" "$(reads "$OUT/SMOKE.meta.R1.fastq.gz")" "$BEFORE"

say "Result: $PASS passed, $FAIL failed"
if [ "$FAIL" -gt 0 ]; then
  echo "Artefacts left in $WORK for inspection." >&2
  exit 1
fi
[ -n "${1:-}" ] || rm -rf "$WORK"
echo "Smoke test passed."
