"""FastAPI backend for the Locus debug SPA.

Thin HTTP layer over the same read-only ``queries`` module the MCP server uses.
Binds to localhost only (private data). Also serves the built React SPA from
``web/dist`` when present, so ``locus serve api`` is a single self-contained app.
"""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import queries
from .db import db_exists

app = FastAPI(title="Locus", description="Explore your genome locally", version="0.1.0")

# Allow the Vite dev server (localhost:5173) to call the API during development.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_methods=["*"],
    allow_headers=["*"],
)


def _guard_db() -> None:
    if not db_exists():
        raise HTTPException(503, "No Locus database. Build it with `locus pipeline`.")


@app.get("/api/overview")
def overview():
    _guard_db()
    return queries.overview()


@app.get("/api/variant/rsid/{rsid}")
def by_rsid(rsid: str, limit: int = 50, offset: int = 0):
    _guard_db()
    return queries.lookup_by_rsid(rsid, limit=limit, offset=offset)


@app.get("/api/gene/{gene}")
def by_gene(gene: str, limit: int = 100, offset: int = 0):
    _guard_db()
    return queries.lookup_by_gene(gene, limit=limit, offset=offset)


@app.get("/api/region")
def by_region(region: str, limit: int = 200, offset: int = 0):
    _guard_db()
    try:
        return queries.lookup_by_region(region, limit=limit, offset=offset)
    except ValueError as e:
        raise HTTPException(400, str(e)) from e


@app.get("/api/clinical")
def clinical(gene: str = "", significance: str = "", limit: int = 100, offset: int = 0):
    _guard_db()
    return queries.clinical_findings(
        gene=gene or None, significance=significance or None, limit=limit, offset=offset
    )


@app.get("/api/pgx")
def pgx(gene: str = "", drug: str = ""):
    _guard_db()
    return queries.pharmacogenomics(gene=gene or None, drug=drug or None)


@app.get("/api/structural")
def structural(region: str, limit: int = 100):
    _guard_db()
    try:
        return queries.structural_overlap(region, limit=limit)
    except ValueError as e:
        raise HTTPException(400, str(e)) from e


@app.get("/api/ancestry")
def ancestry():
    _guard_db()
    return queries.ancestry()


@app.get("/api/pgs")
def pgs():
    _guard_db()
    return queries.polygenic_risk()


@app.get("/api/whats_new")
def whats_new(since: str = "", tier: str = ""):
    _guard_db()
    return queries.whats_new(since=since or None, tier=tier or None)


@app.get("/api/secondary_findings")
def secondary_findings(limit: int = 100, offset: int = 0):
    _guard_db()
    return queries.secondary_findings(limit=limit, offset=offset)


@app.get("/api/traits")
def traits(category: str = ""):
    _guard_db()
    return queries.traits(category=category or None)


class SqlBody(BaseModel):
    query: str


@app.post("/api/sql")
def sql(body: SqlBody):
    _guard_db()
    try:
        return queries.run_sql(body.query)
    except ValueError as e:
        raise HTTPException(400, str(e)) from e


# Serve the built SPA if it exists (production); in dev you run Vite separately.
_DIST = Path(__file__).resolve().parents[2] / "web" / "dist"
if _DIST.exists():
    app.mount("/", StaticFiles(directory=str(_DIST), html=True), name="spa")
