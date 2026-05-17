"""Reference manager: lazy-download, hash-verify, cache.

Lifecycle:
  1. Cerberus is invoked; orchestrator asks RefManager for the assets each
     selected pipeline needs.
  2. RefManager checks ``<ref_dir>/manifest.json``. If missing, it copies the
     packaged default manifest. If a newer version is shipped in the package,
     the user is prompted (unless ``--update-refs``).
  3. For each required asset, RefManager verifies the on-disk file's SHA256
     against the manifest. Missing/mismatched ⇒ download (aria2c preferred,
     urllib fallback) into a ``.tmp`` next to the target, then atomic rename.
  4. Tarballed assets (``.tar.zst`` / ``.tar.gz``) are auto-extracted into a
     sibling directory.

The download step is the only network call Cerberus makes during a normal
pipeline run. Disable with ``--no-auto-download`` or pre-warm with
``cerberus fetch-refs``.
"""
from __future__ import annotations

import json
import shutil
import tarfile
import urllib.request
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Iterable

from cerberus.utils.hashing import sha256_file, verify_sha256
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import which

log = get_logger("refs")


_PIPELINE_TO_ASSETS = {
    "meta":           ["masked_t2t_hla_minimap2"],
    "profiling":      ["masked_t2t_hla_minimap2", "masked_t2t_hla_bowtie2", "aux_refs"],
    "profiling-fast": ["masked_t2t_hla_minimap2"],
    "gdpr":           ["kraken2_gdpr_compact", "masked_t2t_hla_minimap2"],
    "long-meta":      ["masked_t2t_hla_minimap2"],
    "long-profiling": ["masked_t2t_hla_minimap2", "aux_refs"],
    "long-gdpr":      ["kraken2_gdpr_compact", "masked_t2t_hla_minimap2"],
}


@dataclass
class Asset:
    key: str
    description: str
    filename: str
    url: str
    sha256: str
    size_bytes: int | None
    required_for: list[str]

    @property
    def is_archive(self) -> bool:
        return self.filename.endswith((".tar", ".tar.gz", ".tar.zst", ".tar.xz"))

    @property
    def extracted_dirname(self) -> str:
        return self.filename.split(".")[0]


class RefManagerError(RuntimeError):
    pass


class RefManager:
    def __init__(self, ref_dir: Path, *, auto_download: bool = True):
        self.ref_dir = ref_dir
        self.auto_download = auto_download
        self.ref_dir.mkdir(parents=True, exist_ok=True)
        self.manifest_path = ref_dir / "manifest.json"
        self.manifest: dict = self._load_manifest()

    def _load_manifest(self) -> dict:
        if not self.manifest_path.exists():
            log.info("No manifest at %s; seeding default", self.manifest_path)
            self._seed_default_manifest()
        with self.manifest_path.open() as f:
            return json.load(f)

    def _seed_default_manifest(self) -> None:
        with resources.files("cerberus.data").joinpath("default_manifest.json").open() as f:
            data = f.read()
        self.manifest_path.write_text(data)

    def asset(self, key: str) -> Asset:
        info = self.manifest["assets"].get(key)
        if not info:
            raise RefManagerError(f"Unknown asset key: {key!r}")
        return Asset(
            key=key,
            description=info["description"],
            filename=info["filename"],
            url=info["url"],
            sha256=info["sha256"],
            size_bytes=info.get("size_bytes"),
            required_for=info.get("required_for", []),
        )

    def required_assets_for(self, pipeline_keys: Iterable[str]) -> list[Asset]:
        seen: set[str] = set()
        out: list[Asset] = []
        for pkey in pipeline_keys:
            for akey in _PIPELINE_TO_ASSETS.get(pkey, []):
                if akey in seen:
                    continue
                seen.add(akey)
                out.append(self.asset(akey))
        return out

    def path_to(self, asset: Asset) -> Path:
        if asset.is_archive:
            return self.ref_dir / asset.extracted_dirname
        return self.ref_dir / asset.filename

    def is_satisfied(self, asset: Asset) -> bool:
        target = self.path_to(asset)
        if asset.is_archive:
            return target.is_dir() and any(target.iterdir())
        if not target.exists():
            return False
        if asset.sha256 and asset.sha256 != "PENDING":
            return verify_sha256(target, asset.sha256)
        return True

    def ensure(self, assets: Iterable[Asset]) -> None:
        for a in assets:
            if self.is_satisfied(a):
                log.info("✓ %s present and verified", a.key)
                continue
            if not self.auto_download:
                raise RefManagerError(
                    f"Asset {a.key!r} missing and --no-auto-download set. "
                    f"Run: cerberus fetch-refs"
                )
            if a.url in ("", "PENDING"):
                raise RefManagerError(
                    f"Asset {a.key!r} has no URL in manifest yet "
                    f"(awaiting first release). Provide --ref-dir with prebuilt assets."
                )
            self._download(a)

    def _download(self, asset: Asset) -> None:
        target = self.ref_dir / asset.filename
        tmp = target.with_suffix(target.suffix + ".tmp")
        log.info("⇩ Downloading %s (%s)", asset.key, asset.url)

        if which("aria2c"):
            self._aria2c(asset.url, tmp)
        else:
            self._urllib(asset.url, tmp)

        if asset.sha256 and asset.sha256 != "PENDING":
            log.info("Verifying SHA256 of %s", asset.filename)
            if not verify_sha256(tmp, asset.sha256, show_progress=True):
                tmp.unlink(missing_ok=True)
                raise RefManagerError(
                    f"SHA256 mismatch for {asset.key}. Expected {asset.sha256}. "
                    f"File deleted; please re-run."
                )
        else:
            actual = sha256_file(tmp, show_progress=True)
            log.warning("Manifest had no SHA256 for %s; observed %s", asset.key, actual)

        tmp.rename(target)
        if asset.is_archive:
            self._extract(target)

    def _aria2c(self, url: str, dst: Path) -> None:
        from cerberus.utils.shell import run
        log_path = dst.with_suffix(".aria2.log")
        run(
            [
                "aria2c", "-x16", "-s16", "-k1M", "--allow-overwrite=true",
                "--summary-interval=10", "--console-log-level=warn",
                "-d", str(dst.parent), "-o", dst.name, url,
            ],
            log_path=log_path,
        )

    def _urllib(self, url: str, dst: Path) -> None:
        try:
            from tqdm import tqdm
        except ImportError:
            tqdm = None  # type: ignore[assignment]

        log.info("aria2c not found; falling back to urllib (slower)")

        with urllib.request.urlopen(url) as resp:
            total = int(resp.headers.get("Content-Length", 0) or 0)
            bar = tqdm(total=total, unit="B", unit_scale=True, desc=dst.name) \
                if (tqdm and total) else None
            with dst.open("wb") as f:
                while chunk := resp.read(1024 * 1024):
                    f.write(chunk)
                    if bar:
                        bar.update(len(chunk))
            if bar:
                bar.close()

    def _extract(self, archive: Path) -> None:
        log.info("Extracting %s", archive.name)
        out_dir = self.ref_dir / archive.name.split(".")[0]
        out_dir.mkdir(parents=True, exist_ok=True)
        if archive.suffix in {".gz", ".xz"} or archive.name.endswith(".tar"):
            with tarfile.open(archive, "r:*") as tar:
                tar.extractall(out_dir, filter="data")
        elif archive.name.endswith(".tar.zst"):
            zstd = which("zstd")
            if not zstd:
                raise RefManagerError("zstd not found, required to extract .tar.zst")
            from cerberus.utils.shell import pipe
            pipe(
                [["zstd", "-dc", str(archive)], ["tar", "-x", "-C", str(out_dir)]],
                log_path=archive.with_suffix(".extract.log"),
            )
        else:
            raise RefManagerError(f"Unknown archive format: {archive.name}")

    def fetch_all(self) -> None:
        all_assets = [self.asset(k) for k in self.manifest["assets"]]
        self.ensure(all_assets)

    def doctor(self) -> list[str]:
        problems: list[str] = []
        for key in self.manifest["assets"]:
            a = self.asset(key)
            if not self.is_satisfied(a):
                problems.append(f"missing or corrupt: {key} ({a.filename})")
        return problems


def cleanup_partial(ref_dir: Path) -> int:
    """Remove leftover .tmp files from a previous aborted run. Returns count."""
    n = 0
    for f in ref_dir.glob("*.tmp"):
        f.unlink()
        n += 1
    return n
