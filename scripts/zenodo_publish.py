#!/usr/bin/env python3
"""Publish the Zenodo draft created by zenodo_upload.py and update the in-repo manifest.

Usage:
    CERBERUS_ZENODO_TOKEN=$(cat ~/.cerberus_zenodo) python scripts/zenodo_publish.py
"""
from __future__ import annotations

import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

ZENODO_BASE = "https://zenodo.org/api"
CACHE_PATH = Path.home() / ".cerberus" / "zenodo_deposition.json"
MANIFEST_PATH = Path("/home/iowa/Desktop/cerberus/cerberus/data/default_manifest.json")

ASSET_KEY_BY_FILENAME = {
    "masked_t2t_hla.mmi":            "masked_t2t_hla_minimap2",
    "masked_t2t_hla_bt2.tar.zst":    "masked_t2t_hla_bowtie2",
    "kraken2_gdpr_compact.tar.zst":  "kraken2_gdpr_compact",
    "aux_refs.fa.gz":                "aux_refs",
    "human_k27.fa.gz":               "human_kmer_set",
}


def _token() -> str:
    t = os.environ.get("CERBERUS_ZENODO_TOKEN")
    if not t:
        sys.exit("ERROR: set CERBERUS_ZENODO_TOKEN env var")
    return t.strip()


def _req(method: str, url: str, *, data: bytes | None = None,
         headers: dict | None = None) -> dict:
    headers = headers or {}
    headers.setdefault("Authorization", f"Bearer {_token()}")
    req = urllib.request.Request(url, data=data, method=method, headers=headers)
    try:
        with urllib.request.urlopen(req) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as e:
        sys.exit(f"HTTP {e.code} on {method} {url}\n{e.read().decode(errors='replace')}")


def main() -> int:
    if not CACHE_PATH.exists():
        sys.exit("No deposition cache; run zenodo_upload.py first")
    cache = json.loads(CACHE_PATH.read_text())
    dep_id = cache["deposition_id"]

    print(f"Publishing deposition {dep_id}...", flush=True)
    published = _req("POST", f"{ZENODO_BASE}/deposit/depositions/{dep_id}/actions/publish")

    doi = published.get("doi") or published.get("metadata", {}).get("doi")
    concept_doi = published.get("conceptdoi")
    record_id = published.get("record_id") or dep_id
    record_url = f"https://zenodo.org/record/{record_id}"

    print(f"✓ Published. DOI: {doi}  conceptDOI: {concept_doi}")
    print(f"  Record: {record_url}")

    files_by_name: dict[str, dict] = {f["filename"]: f for f in published.get("files", [])}

    manifest = json.loads(MANIFEST_PATH.read_text())
    manifest["release"] = "0.1.0"
    manifest["zenodo_doi"] = doi
    manifest["zenodo_concept_doi"] = concept_doi
    manifest["zenodo_record_url"] = record_url

    updated_keys: list[str] = []
    for filename, asset_key in ASSET_KEY_BY_FILENAME.items():
        info = files_by_name.get(filename)
        if not info:
            print(f"⚠ Published deposition has no file {filename!r}; manifest entry unchanged.")
            continue
        download_url = info.get("links", {}).get("download") \
            or f"https://zenodo.org/record/{record_id}/files/{filename}?download=1"
        checksum = info.get("checksum", "")
        if checksum.startswith("md5:"):
            checksum_md5 = checksum[len("md5:"):]
        else:
            checksum_md5 = checksum
        entry = manifest["assets"][asset_key]
        entry["url"] = download_url
        entry["size_bytes"] = info.get("filesize")
        entry["md5"] = checksum_md5
        updated_keys.append(asset_key)

    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"✓ Updated manifest: {MANIFEST_PATH}  ({len(updated_keys)} entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
