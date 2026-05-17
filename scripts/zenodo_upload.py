#!/usr/bin/env python3
"""Upload the five Cerberus reference assets to Zenodo as a draft deposition.

Usage:
    CERBERUS_ZENODO_TOKEN=$(cat ~/.cerberus_zenodo) python scripts/zenodo_upload.py

Idempotent: if a file with matching name+size already exists in the deposition,
it is skipped. Deposition ID is cached at ~/.cerberus/zenodo_deposition.json
so re-runs continue where they left off.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ZENODO_BASE = "https://zenodo.org/api"
CACHE_PATH = Path.home() / ".cerberus" / "zenodo_deposition.json"
BUILD_DIR = Path("/home/iowa/Desktop/cerberus/scripts/build_refs/build")
REPO_URL = "https://github.com/iowa69/cerberus"

ASSETS = [
    ("aux_refs.fa.gz",                BUILD_DIR / "aux_refs/aux_refs.fa.gz"),
    ("human_k27.fa.gz",               BUILD_DIR / "aux_refs/human_k27.fa.gz"),
    ("kraken2_gdpr_compact.tar.zst",  BUILD_DIR / "kraken2_gdpr_compact.tar.zst"),
    ("masked_t2t_hla.mmi",            BUILD_DIR / "masked_t2t_hla/masked_t2t_hla.mmi"),
    ("masked_t2t_hla_bt2.tar.zst",    BUILD_DIR / "masked_t2t_hla/masked_t2t_hla_bt2.tar.zst"),
]

METADATA = {
    "metadata": {
        "title": "Cerberus reference data v0.1.0 — masked T2T-CHM13v2.0 + IPD-IMGT/HLA, Kraken2 GDPR DB, auxiliary k-mer references",
        "upload_type": "dataset",
        "description": (
            "<p>Reference data bundle for <a href='https://github.com/iowa69/cerberus'>Cerberus</a>, "
            "a three-headed host-removal pipeline for metagenomic data.</p>"
            "<p>Contents:</p>"
            "<ul>"
            "<li><code>masked_t2t_hla.mmi</code> — minimap2 index of T2T-CHM13v2.0 + IPD-IMGT/HLA, "
            "low-entropy and viral-homology masked.</li>"
            "<li><code>masked_t2t_hla_bt2.tar.zst</code> — bowtie2 index of the same reference.</li>"
            "<li><code>kraken2_gdpr_compact.tar.zst</code> — compact Kraken2 database covering human "
            "(T2T-CHM13v2.0), chimpanzee, gorilla, mouse, and rat for GDPR-pass host scrubbing.</li>"
            "<li><code>aux_refs.fa.gz</code> — Human host-decoy ncRNA (rRNA, snRNA, snoRNA, miRNA, "
            "Y_RNA, vaultRNA) + mitochondrion (NC_012920.1).</li>"
            "<li><code>human_k27.fa.gz</code> — Full T2T-CHM13v2.0 used as bbduk k-mer reference "
            "in Cerberus's GDPR pass (orthogonal k=31 mechanism alongside Kraken2).</li>"
            "</ul>"
            "<p>Source build scripts: <a href='https://github.com/iowa69/cerberus/tree/main/scripts/build_refs'>"
            "scripts/build_refs/</a>.</p>"
        ),
        "creators": [
            {"name": "Lorenzin, Giovanni", "affiliation": ""},
        ],
        "keywords": [
            "bioinformatics", "metagenomics", "host removal", "decontamination",
            "kraken2", "minimap2", "bowtie2", "T2T-CHM13", "IPD-IMGT/HLA", "GDPR",
        ],
        "license": "cc-by-4.0",
        "access_right": "open",
        "version": "0.1.0",
        "related_identifiers": [
            {"identifier": REPO_URL, "relation": "isSupplementTo",
             "resource_type": "software"},
        ],
    }
}


def _token() -> str:
    t = os.environ.get("CERBERUS_ZENODO_TOKEN")
    if not t:
        sys.exit("ERROR: set CERBERUS_ZENODO_TOKEN env var")
    return t.strip()


def _req(method: str, url: str, *, data: bytes | None = None,
         headers: dict | None = None) -> dict | bytes:
    headers = headers or {}
    headers.setdefault("Authorization", f"Bearer {_token()}")
    req = urllib.request.Request(url, data=data, method=method, headers=headers)
    try:
        with urllib.request.urlopen(req) as resp:
            body = resp.read()
            ct = resp.headers.get("Content-Type", "")
            if "json" in ct:
                return json.loads(body)
            return body
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", errors="replace")
        sys.exit(f"HTTP {e.code} on {method} {url}\n{body}")


def load_cache() -> dict:
    if CACHE_PATH.exists():
        return json.loads(CACHE_PATH.read_text())
    return {}


def save_cache(c: dict) -> None:
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    CACHE_PATH.write_text(json.dumps(c, indent=2))


def create_or_resume() -> dict:
    cache = load_cache()
    if dep_id := cache.get("deposition_id"):
        print(f"Resuming deposition {dep_id}", flush=True)
        info = _req("GET", f"{ZENODO_BASE}/deposit/depositions/{dep_id}")
        return info  # type: ignore[return-value]

    print("Creating new deposition", flush=True)
    info = _req(
        "POST", f"{ZENODO_BASE}/deposit/depositions",
        data=json.dumps({}).encode(),
        headers={"Content-Type": "application/json"},
    )
    cache["deposition_id"] = info["id"]                                # type: ignore[index]
    cache["bucket_url"]    = info["links"]["bucket"]                   # type: ignore[index]
    cache["created_at"]    = info.get("created")                       # type: ignore[union-attr]
    save_cache(cache)
    return info                                                        # type: ignore[return-value]


def set_metadata(dep_id: int) -> None:
    print("Setting deposition metadata", flush=True)
    _req(
        "PUT", f"{ZENODO_BASE}/deposit/depositions/{dep_id}",
        data=json.dumps(METADATA).encode(),
        headers={"Content-Type": "application/json"},
    )


def existing_files(dep_id: int) -> dict[str, int]:
    info = _req("GET", f"{ZENODO_BASE}/deposit/depositions/{dep_id}")
    out: dict[str, int] = {}
    for f in info.get("files", []):                                    # type: ignore[union-attr]
        out[f["filename"]] = f["filesize"]
    return out


def upload_one(bucket_url: str, filename: str, path: Path) -> dict:
    url = f"{bucket_url}/{filename}"
    size = path.stat().st_size
    t0 = time.time()
    print(f"⇧ {filename}  ({size/1024/1024:.0f} MB)", flush=True)

    chunk = 4 * 1024 * 1024
    sent = 0
    last_report = t0

    class StreamingFile:
        def __init__(self, fp):
            self.fp = fp
        def read(self, n=-1):
            nonlocal sent, last_report
            buf = self.fp.read(n if n > 0 else chunk)
            sent += len(buf)
            now = time.time()
            if now - last_report >= 10:
                pct = 100 * sent / size
                mbps = (sent / 1024 / 1024) / (now - t0)
                print(f"   {filename}: {pct:.1f}% ({mbps:.1f} MB/s)", flush=True)
                last_report = now
            return buf

    with path.open("rb") as f:
        req = urllib.request.Request(
            url,
            data=StreamingFile(f),                                     # type: ignore[arg-type]
            method="PUT",
            headers={
                "Authorization": f"Bearer {_token()}",
                "Content-Type": "application/octet-stream",
                "Content-Length": str(size),
            },
        )
        try:
            with urllib.request.urlopen(req) as resp:
                body = resp.read()
                return json.loads(body)
        except urllib.error.HTTPError as e:
            sys.exit(f"Upload failed for {filename}: HTTP {e.code}\n{e.read().decode(errors='replace')}")


def main() -> int:
    info = create_or_resume()
    dep_id = info["id"]
    bucket_url = info["links"]["bucket"]
    print(f"Deposition: id={dep_id}  draft URL: https://zenodo.org/deposit/{dep_id}", flush=True)

    set_metadata(dep_id)

    have = existing_files(dep_id)
    for name, path in ASSETS:
        if not path.exists():
            sys.exit(f"Missing asset: {path}")
        sz = path.stat().st_size
        if have.get(name) == sz:
            print(f"✓ {name}  already uploaded ({sz/1024/1024:.0f} MB)", flush=True)
            continue
        upload_one(bucket_url, name, path)

    print()
    print("=== Upload complete. Draft deposition has all 5 files. ===")
    print(f"View/edit:  https://zenodo.org/deposit/{dep_id}")
    print("Next: run scripts/zenodo_publish.py to publish and mint the DOI.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
