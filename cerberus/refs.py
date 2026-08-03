"""Reference manager: lazy-download, hash-verify, cache.

Lifecycle:
  1. Cerberus is invoked; the orchestrator asks RefManager for the assets each
     selected pipeline needs.
  2. RefManager checks ``<ref_dir>/manifest.json``. If missing, it copies the
     packaged default manifest. If the packaged manifest is newer, that is
     reported and ``--update-refs`` adopts it.
  3. For each required asset, RefManager verifies the on-disk file against the
     manifest SHA256. Verification results are cached in ``.verified.json``
     keyed by (size, mtime), so a 7.7 GB index is not re-hashed on every run —
     but any change to the file invalidates the stamp and forces a re-hash.
  4. Tarballed assets (``.tar.zst`` / ``.tar.gz``) are extracted into a
     sibling directory *atomically*: extraction goes to a ``.partial``
     directory that is renamed into place only on success, so an interrupted
     extraction can never be mistaken for a complete one.

Concurrent runs sharing a ``--ref-dir`` coordinate through a lock file, so one
run's ``cleanup_partial`` cannot delete another run's in-flight download.

The download step is the only network call Cerberus makes during a normal
pipeline run. Disable with ``--no-auto-download`` or pre-warm with
``cerberus fetch-refs``.
"""
from __future__ import annotations

import json
import os
import shutil
import tarfile
import time
import urllib.error
import urllib.request
from collections.abc import Iterable
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

from cerberus.utils.hashing import sha256_file, verify_sha256
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import which

log = get_logger("refs")


_PIPELINE_TO_ASSETS = {
    "meta":           ["masked_t2t_hla_minimap2"],
    "profiling":      ["masked_t2t_hla_minimap2", "masked_t2t_hla_bowtie2", "aux_refs"],
    "profiling-fast": ["masked_t2t_hla_minimap2"],
    "gdpr":           ["kraken2_gdpr_compact", "masked_t2t_hla_minimap2", "human_kmer_set"],
    "long-meta":      ["masked_t2t_hla_minimap2"],
    "long-profiling": ["masked_t2t_hla_minimap2", "aux_refs"],
    "long-gdpr":      ["kraken2_gdpr_compact", "masked_t2t_hla_minimap2", "human_kmer_set"],
}

# Assets a run can proceed without, at reduced strength.
_OPTIONAL_ASSETS = {"human_kmer_set"}

_DOWNLOAD_TIMEOUT = 60
_DOWNLOAD_RETRIES = 3
_LOCK_STALE_SEC = 6 * 60 * 60


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
        """Directory an archive extracts into.

        Strips only the archive suffixes, so ``T2T-CHM13v2.0_bt2.tar.zst``
        becomes ``T2T-CHM13v2.0_bt2`` rather than being truncated at the
        first dot.
        """
        name = self.filename
        for suffix in (".tar.zst", ".tar.gz", ".tar.xz", ".tar"):
            if name.endswith(suffix):
                return name[: -len(suffix)]
        return name.split(".")[0]


class RefManagerError(RuntimeError):
    pass


class RefLock:
    """Advisory lock over a reference directory (atomic O_EXCL create)."""

    def __init__(self, ref_dir: Path, name: str = ".cerberus.lock"):
        self.path = ref_dir / name
        self.held = False

    def acquire(self, timeout: float = 1800.0, poll: float = 2.0) -> bool:
        deadline = time.time() + timeout
        while True:
            try:
                fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
                os.write(fd, f"{os.getpid()} {time.time():.0f}\n".encode())
                os.close(fd)
                self.held = True
                return True
            except FileExistsError:
                if self._is_stale():
                    log.warning("Removing stale reference lock %s", self.path)
                    self.path.unlink(missing_ok=True)
                    continue
                if time.time() >= deadline:
                    return False
                log.info("Waiting for another Cerberus run to finish with %s", self.path.parent)
                time.sleep(poll)
            except OSError as e:
                log.warning("Could not create reference lock (%s); proceeding unlocked", e)
                return True

    def _is_stale(self) -> bool:
        try:
            return (time.time() - self.path.stat().st_mtime) > _LOCK_STALE_SEC
        except OSError:
            return False

    def release(self) -> None:
        if self.held:
            self.path.unlink(missing_ok=True)
            self.held = False

    def __enter__(self):
        self.acquire()
        return self

    def __exit__(self, *exc) -> None:
        self.release()


class RefManager:
    def __init__(
        self,
        ref_dir: Path,
        *,
        auto_download: bool = True,
        kraken2_db_override: Path | None = None,
        aux_refs_override: Path | None = None,
    ):
        self.ref_dir = ref_dir
        self.auto_download = auto_download
        self.kraken2_db_override = kraken2_db_override
        self.aux_refs_override = aux_refs_override
        self.ref_dir.mkdir(parents=True, exist_ok=True)
        self.manifest_path = ref_dir / "manifest.json"
        self.manifest: dict = self._load_manifest()
        self._verify_cache_path = ref_dir / ".verified.json"
        self._verify_cache = self._load_verify_cache()
        self.skipped_optional: list[str] = []

    # ---------- manifest ----------

    def _load_manifest(self) -> dict:
        if not self.manifest_path.exists():
            log.info("No manifest at %s; seeding default", self.manifest_path)
            self._seed_default_manifest()
        try:
            with self.manifest_path.open() as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            raise RefManagerError(
                f"{self.manifest_path} is not valid JSON ({e}). "
                "Delete it to re-seed the packaged default."
            ) from e

    def _seed_default_manifest(self) -> None:
        data = resources.files("cerberus.data").joinpath("default_manifest.json").read_text()
        self.manifest_path.write_text(data)

    def _packaged_manifest(self) -> dict:
        return json.loads(
            resources.files("cerberus.data").joinpath("default_manifest.json").read_text()
        )

    def manifest_update_available(self) -> str | None:
        """Return the packaged release string when it differs from on-disk."""
        try:
            packaged = self._packaged_manifest()
        except (OSError, json.JSONDecodeError):
            return None
        mine = str(self.manifest.get("release", ""))
        theirs = str(packaged.get("release", ""))
        return theirs if theirs and theirs != mine else None

    def adopt_packaged_manifest(self) -> None:
        self._seed_default_manifest()
        self.manifest = self._load_manifest()
        log.info("Adopted packaged manifest (release %s)", self.manifest.get("release"))

    # ---------- verification cache ----------

    def _load_verify_cache(self) -> dict:
        try:
            with self._verify_cache_path.open() as f:
                data = json.load(f)
                return data if isinstance(data, dict) else {}
        except (OSError, json.JSONDecodeError):
            return {}

    def _save_verify_cache(self) -> None:
        try:
            self._verify_cache_path.write_text(json.dumps(self._verify_cache, indent=2))
        except OSError as e:
            log.debug("Could not persist verification cache: %s", e)

    @staticmethod
    def _stamp(path: Path) -> str:
        st = path.stat()
        return f"{st.st_size}:{int(st.st_mtime)}"

    def _cached_ok(self, path: Path, sha: str) -> bool:
        entry = self._verify_cache.get(str(path))
        return bool(entry and entry.get("stamp") == self._stamp(path)
                    and entry.get("sha256") == sha)

    def _remember(self, path: Path, sha: str) -> None:
        self._verify_cache[str(path)] = {"stamp": self._stamp(path), "sha256": sha}
        self._save_verify_cache()

    # ---------- assets ----------

    def asset(self, key: str) -> Asset:
        assets = self.manifest.get("assets")
        if not isinstance(assets, dict):
            raise RefManagerError(f"{self.manifest_path} has no 'assets' object")
        info = assets.get(key)
        if not info:
            raise RefManagerError(
                f"Unknown asset key: {key!r}. Known keys: {', '.join(sorted(assets))}"
            )
        return Asset(
            key=key,
            description=info.get("description", ""),
            filename=info["filename"],
            url=info.get("url", ""),
            sha256=info.get("sha256", ""),
            size_bytes=info.get("size_bytes"),
            required_for=info.get("required_for", []),
        )

    def has_asset(self, key: str) -> bool:
        assets = self.manifest.get("assets")
        return isinstance(assets, dict) and key in assets

    def required_assets_for(self, pipeline_keys: Iterable[str]) -> list[Asset]:
        seen: set[str] = set()
        out: list[Asset] = []
        for pkey in pipeline_keys:
            if pkey not in _PIPELINE_TO_ASSETS:
                raise RefManagerError(
                    f"No reference set defined for pipeline {pkey!r}. This is a Cerberus bug: "
                    "a mode combination was selected that the reference map does not cover."
                )
            for akey in _PIPELINE_TO_ASSETS[pkey]:
                if akey in seen or not self.has_asset(akey):
                    continue
                seen.add(akey)
                out.append(self.asset(akey))
        return out

    def path_to(self, asset: Asset) -> Path:
        if asset.is_archive:
            return self.ref_dir / asset.extracted_dirname
        return self.ref_dir / asset.filename

    # ---------- resolved paths used by pipelines ----------

    def kraken2_db_path(self) -> Path:
        """Directory holding hash.k2d, honouring --kraken2-db."""
        base = self.kraken2_db_override or self.path_to(self.asset("kraken2_gdpr_compact"))
        if (base / "hash.k2d").exists():
            return base
        if base.is_dir():
            for sub in sorted(base.iterdir()):
                if sub.is_dir() and (sub / "hash.k2d").exists():
                    return sub
        raise FileNotFoundError(
            f"No Kraken2 database (hash.k2d) found under {base}. "
            "Point --kraken2-db at a directory containing hash.k2d/opts.k2d/taxo.k2d."
        )

    def aux_refs_path(self) -> Path:
        if self.aux_refs_override:
            if not self.aux_refs_override.exists():
                raise FileNotFoundError(f"--aux-refs not found: {self.aux_refs_override}")
            return self.aux_refs_override
        return self.path_to(self.asset("aux_refs"))

    def human_kmer_path(self) -> Path | None:
        """Human k-mer reference for the GDPR belt-and-braces scrub, if present."""
        if not self.has_asset("human_kmer_set"):
            return None
        p = self.path_to(self.asset("human_kmer_set"))
        return p if p.exists() else None

    # ---------- satisfaction / download ----------

    def is_satisfied(self, asset: Asset) -> bool:
        target = self.path_to(asset)
        if asset.is_archive:
            if not (target.is_dir() and any(target.iterdir())):
                return False
            return (target / ".cerberus_extract_ok").exists() or _looks_complete(target)
        if not target.exists():
            return False
        if asset.sha256 and asset.sha256 != "PENDING":
            if self._cached_ok(target, asset.sha256):
                return True
            if verify_sha256(target, asset.sha256, show_progress=True):
                self._remember(target, asset.sha256)
                return True
            return False
        return True

    def ensure(self, assets: Iterable[Asset]) -> None:
        assets = list(assets)
        if not assets:
            return
        lock = RefLock(self.ref_dir)
        if not lock.acquire():
            raise RefManagerError(
                f"Timed out waiting for the reference lock at {lock.path}. "
                "Another Cerberus run may be downloading. Remove the file if it is stale."
            )
        try:
            for a in assets:
                if self.is_satisfied(a):
                    log.info("[ok] %s present and verified", a.key)
                    continue
                if not self.auto_download:
                    if a.key in _OPTIONAL_ASSETS:
                        self._skip_optional(a, "missing and --no-auto-download set")
                        continue
                    raise RefManagerError(
                        f"Asset {a.key!r} missing and --no-auto-download set. "
                        f"Run: cerberus fetch-refs"
                    )
                if a.url in ("", "PENDING"):
                    if a.key in _OPTIONAL_ASSETS:
                        self._skip_optional(a, "no URL in manifest")
                        continue
                    raise RefManagerError(
                        f"Asset {a.key!r} has no URL in the manifest. "
                        f"Provide --ref-dir with prebuilt assets."
                    )
                try:
                    self._download(a)
                except (RefManagerError, OSError, urllib.error.URLError) as e:
                    if a.key in _OPTIONAL_ASSETS:
                        self._skip_optional(a, str(e))
                        continue
                    raise
        finally:
            lock.release()

    def _skip_optional(self, asset: Asset, why: str) -> None:
        log.warning(
            "Optional asset %r unavailable (%s). Continuing with one fewer GDPR mechanism.",
            asset.key, why,
        )
        self.skipped_optional.append(asset.key)

    def _download(self, asset: Asset) -> None:
        target = self.ref_dir / asset.filename
        tmp = target.with_name(target.name + ".tmp")
        log.info("Downloading %s (%s)", asset.key, asset.url)

        last_err: Exception | None = None
        for attempt in range(1, _DOWNLOAD_RETRIES + 1):
            try:
                if which("aria2c"):
                    self._aria2c(asset.url, tmp)
                else:
                    self._urllib(asset.url, tmp)
                break
            except Exception as e:                      # noqa: BLE001 - retried below
                last_err = e
                log.warning("Download attempt %d/%d for %s failed: %s",
                            attempt, _DOWNLOAD_RETRIES, asset.key, e)
                if attempt == _DOWNLOAD_RETRIES:
                    tmp.unlink(missing_ok=True)
                    raise RefManagerError(
                        f"Could not download {asset.key} from {asset.url}: {e}"
                    ) from last_err
                time.sleep(2 ** attempt)

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

        tmp.replace(target)
        if asset.sha256 and asset.sha256 != "PENDING":
            self._remember(target, asset.sha256)
        if asset.is_archive:
            self._extract(target, asset)

    def _aria2c(self, url: str, dst: Path) -> None:
        from cerberus.utils.shell import run
        log_path = dst.with_name(dst.name + ".aria2.log")
        run(
            [
                "aria2c", "-x16", "-s16", "-k1M", "--allow-overwrite=true",
                "--continue=true", "--max-tries=3", "--retry-wait=5",
                "--summary-interval=10", "--console-log-level=warn",
                "-d", str(dst.parent), "-o", dst.name, url,
            ],
            log_path=log_path,
        )

    def _urllib(self, url: str, dst: Path) -> None:
        try:
            from tqdm import tqdm
        except ImportError:
            tqdm = None                                  # type: ignore[assignment]

        log.info("aria2c not found; falling back to urllib (slower)")

        with urllib.request.urlopen(url, timeout=_DOWNLOAD_TIMEOUT) as resp:
            if getattr(resp, "status", 200) >= 400:
                raise RefManagerError(f"HTTP {resp.status} for {url}")
            total = int(resp.headers.get("Content-Length", 0) or 0)
            bar = tqdm(total=total, unit="B", unit_scale=True, desc=dst.name) \
                if (tqdm and total) else None
            try:
                with dst.open("wb") as f:
                    while chunk := resp.read(1024 * 1024):
                        f.write(chunk)
                        if bar:
                            bar.update(len(chunk))
            finally:
                if bar:
                    bar.close()

    def _extract(self, archive: Path, asset: Asset) -> None:
        """Extract into a .partial dir, then rename — never a half-extracted tree."""
        log.info("Extracting %s", archive.name)
        final_dir = self.ref_dir / asset.extracted_dirname
        staging = self.ref_dir / f"{asset.extracted_dirname}.partial"
        if staging.exists():
            shutil.rmtree(staging, ignore_errors=True)
        staging.mkdir(parents=True)

        try:
            if archive.name.endswith(".tar.zst"):
                if not which("zstd"):
                    raise RefManagerError("zstd not found, required to extract .tar.zst")
                from cerberus.utils.shell import pipe
                pipe(
                    [["zstd", "-dc", str(archive)], ["tar", "-x", "-C", str(staging)]],
                    log_path=archive.with_name(archive.name + ".extract.log"),
                )
            elif archive.name.endswith((".tar", ".tar.gz", ".tar.xz")):
                with tarfile.open(archive, "r:*") as tar:
                    _safe_extract(tar, staging)
            else:
                raise RefManagerError(f"Unknown archive format: {archive.name}")

            (staging / ".cerberus_extract_ok").write_text(
                json.dumps({"asset": asset.key, "sha256": asset.sha256})
            )
            if final_dir.exists():
                shutil.rmtree(final_dir, ignore_errors=True)
            staging.replace(final_dir)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise

        if not self._keep_archives():
            archive.unlink(missing_ok=True)
            log.info("Removed %s after extraction", archive.name)

    @staticmethod
    def _keep_archives() -> bool:
        return os.environ.get("CERBERUS_KEEP_ARCHIVES", "").lower() in {"1", "true", "yes"}

    def fetch_all(self, *, force: bool = False) -> None:
        assets = [self.asset(k) for k in self.manifest.get("assets", {})]
        if force:
            for a in assets:
                self._verify_cache.pop(str(self.path_to(a)), None)
                if not a.is_archive:
                    self.path_to(a).unlink(missing_ok=True)
                else:
                    shutil.rmtree(self.path_to(a), ignore_errors=True)
            self._save_verify_cache()
        self.ensure(assets)

    def doctor(self) -> list[str]:
        problems: list[str] = []
        for key in self.manifest.get("assets", {}):
            a = self.asset(key)
            if self.is_satisfied(a):
                continue
            optional = " (optional)" if key in _OPTIONAL_ASSETS else ""
            problems.append(f"missing or corrupt{optional}: {key} ({a.filename})")
        return problems


def _looks_complete(extracted: Path) -> bool:
    """Heuristic for archives extracted by an older Cerberus without a marker."""
    known_markers = ("hash.k2d", ".1.bt2", ".1.bt2l")
    for p in extracted.rglob("*"):
        if p.name in known_markers or p.name.endswith(known_markers):
            return True
    return False


def _safe_extract(tar: tarfile.TarFile, dest: Path) -> None:
    """extractall with the data filter, working on Python 3.10-3.11 too."""
    try:
        tar.extractall(dest, filter="data")
    except TypeError:
        # filter= landed mid-3.11; fall back to a manual traversal check.
        base = dest.resolve()
        for member in tar.getmembers():
            target = (dest / member.name).resolve()
            if not str(target).startswith(str(base)):
                raise RefManagerError(f"Refusing unsafe tar path: {member.name}") from None
        tar.extractall(dest)


def cleanup_partial(ref_dir: Path) -> int:
    """Remove leftovers from a previous aborted run. Returns the count.

    Skipped entirely while another run holds the lock, so this can never
    delete a download that is currently in flight.
    """
    lock = RefLock(ref_dir)
    if lock.path.exists() and not lock._is_stale():
        log.info("Another Cerberus run holds %s; skipping partial cleanup", lock.path)
        return 0
    n = 0
    for pattern in ("*.tmp", "*.tmp.aria2", "*.aria2", "*.partial"):
        for f in ref_dir.glob(pattern):
            try:
                if f.is_dir():
                    shutil.rmtree(f, ignore_errors=True)
                else:
                    f.unlink()
                n += 1
            except OSError:
                pass
    return n
