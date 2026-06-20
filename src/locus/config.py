"""Central configuration for Locus.

All paths are env-overridable (prefix ``LOCUS_``) so the genome and the large
reference / annotation databases can live wherever you have room — e.g. an
external SSD — without touching code. Defaults keep everything under ``data/``,
which is fully ``.gitignore``d.

Example ``.env``::

    LOCUS_DATA_DIR=/Volumes/genome/locus-data
"""

from __future__ import annotations

from pathlib import Path

from pydantic import computed_field, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

_REPO_ROOT = Path(__file__).resolve().parents[2]


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="LOCUS_",
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )

    # Root for everything Locus stores locally (gitignored).
    data_dir: Path = _REPO_ROOT / "data"

    # The DuckDB store that both the MCP server and the SPA query. Defaults to
    # <data_dir>/locus.duckdb; override with LOCUS_DB_PATH to put it elsewhere
    # (e.g. a fast SSD) independent of the data dir.
    db_path: Path | None = None

    @model_validator(mode="after")
    def _default_db_path(self):
        if self.db_path is None:
            self.db_path = self.data_dir / "locus.duckdb"
        return self

    # Sample identity (used to label outputs / pick the sample column in VCFs).
    sample_id: str = "sample"

    # SPA backend bind address — localhost only by design (private data).
    api_host: str = "127.0.0.1"
    api_port: int = 8787

    @computed_field  # type: ignore[prop-decorator]
    @property
    def genome_dir(self) -> Path:
        """Raw inputs downloaded from sequencing.com (the VCFs, optional FASTQ)."""
        return self.data_dir / "genome"

    @computed_field  # type: ignore[prop-decorator]
    @property
    def reference_dir(self) -> Path:
        """GRCh38 reference FASTA + index."""
        return self.data_dir / "reference"

    @computed_field  # type: ignore[prop-decorator]
    @property
    def annotations_dir(self) -> Path:
        """Annotation databases: ClinVar, gnomAD, dbSNP, VEP cache, etc."""
        return self.data_dir / "annotations"

    @computed_field  # type: ignore[prop-decorator]
    @property
    def reports_dir(self) -> Path:
        """Tool reports (e.g. PharmCAT JSON/HTML) and intermediate VCFs."""
        return self.data_dir / "reports"

    @computed_field  # type: ignore[prop-decorator]
    @property
    def work_dir(self) -> Path:
        """Scratch space for normalized / intermediate VCFs."""
        return self.data_dir / "work"

    def ensure_dirs(self) -> None:
        for d in (
            self.data_dir,
            self.genome_dir,
            self.reference_dir,
            self.annotations_dir,
            self.reports_dir,
            self.work_dir,
        ):
            d.mkdir(parents=True, exist_ok=True)


settings = Settings()
