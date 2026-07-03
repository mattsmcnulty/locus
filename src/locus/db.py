"""DuckDB access for Locus.

A single ``locus.duckdb`` file is the source of truth. The pipeline writes it
(``load`` phase); the MCP server and the FastAPI app open it **read-only** so
multiple readers can share it safely while you query.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import duckdb

from .config import settings

# Schema is versioned so we can detect/upgrade an old store.
# v2: added sources + watch_seen_ids + watch_findings (the "living" refresh spine).
# v3: traits + associations tables; ancestry_pca gained a `group` (continent) column.
# v4: watch_findings gained a `url` column (clickable citations for PubMed/GWAS findings).
SCHEMA_VERSION = 4


@contextmanager
def connect(read_only: bool = False, db_path: Path | None = None) -> Iterator[duckdb.DuckDBPyConnection]:
    """Yield a DuckDB connection to the Locus store.

    Use ``read_only=True`` from the MCP server and the API so several processes
    can query concurrently. The writer (load phase) opens with ``read_only=False``.
    """
    path = Path(db_path or settings.db_path)
    if read_only and not path.exists():
        raise FileNotFoundError(
            f"Locus DB not found at {path}. Run `locus pipeline` (or `locus load`) first."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = duckdb.connect(str(path), read_only=read_only)
    try:
        yield conn
    finally:
        conn.close()


def db_exists(db_path: Path | None = None) -> bool:
    return Path(db_path or settings.db_path).exists()
