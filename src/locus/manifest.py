"""Version manifest for Locus's external data sources.

``data/manifest.json`` is the source of truth for "what release of each database do
we currently have loaded" — it survives a DuckDB rebuild and turns today's
file-exists download skips into *version-aware* skips. The refresh engine
(``refresh.py``) reads it to decide what changed, and projects it into the DuckDB
``sources`` table for queryability.

Schema (per source name):
    {
      "version":      "<release id / md5 / date>",
      "url":          "<resolved download URL>",
      "checksum":     "<md5 or etag, if known>",
      "license":      "<spdx-ish label>",
      "last_checked": "<ISO8601, when we last probed>",
      "last_updated": "<ISO8601, when the version last changed>"
    }
"""

from __future__ import annotations

import datetime as _dt
import json
from pathlib import Path

from .config import settings

# Static license labels recorded per source (informational; all fine for personal use).
LICENSES = {
    "clinvar": "public-domain (NCBI ClinVar)",
    "pgs_catalog": "open (PGS Catalog, EBI)",
    "gwas_catalog": "open (NHGRI-EBI GWAS Catalog)",
    "pubmed": "public-domain (NCBI PubMed metadata)",
    "litvar": "public-domain (NCBI LitVar2 metadata)",
    "alphamissense": "CC-BY-NC 4.0",
    "gnomad": "open (gnomAD)",
}


def manifest_path() -> Path:
    return settings.data_dir / "manifest.json"


def now_iso() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


def load() -> dict:
    """Read the manifest (empty dict if absent or malformed)."""
    p = manifest_path()
    if not p.exists():
        return {}
    try:
        data = json.loads(p.read_text())
        return data if isinstance(data, dict) else {}
    except (json.JSONDecodeError, OSError):
        return {}


def save(manifest: dict) -> Path:
    """Write the manifest atomically (tmp + replace) so a crash can't truncate it."""
    p = manifest_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    tmp.replace(p)
    return p


def get(manifest: dict, source: str) -> dict:
    return manifest.get(source, {}) if isinstance(manifest.get(source), dict) else {}


def record(manifest: dict, source: str, *, version: str, url: str = "", checksum: str = "",
           changed: bool) -> dict:
    """Update one source's entry in-place and return it.

    ``last_checked`` always advances; ``last_updated`` advances only when ``changed``.
    """
    entry = get(manifest, source)
    ts = now_iso()
    entry.update({
        "version": version,
        "url": url or entry.get("url", ""),
        "checksum": checksum or entry.get("checksum", ""),
        "license": LICENSES.get(source, entry.get("license", "")),
        "last_checked": ts,
    })
    if changed or "last_updated" not in entry:
        entry["last_updated"] = ts
    manifest[source] = entry
    return entry
