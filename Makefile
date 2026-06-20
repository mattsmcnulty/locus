# Locus — common tasks. Run `make help`.
.DEFAULT_GOAL := help
SHELL := /bin/bash

.PHONY: help setup tools doctor ingest annotate load pipeline serve-mcp serve-api web-dev web-build lint clean

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

tools: ## Install system tools (Homebrew)
	brew install uv bcftools samtools htslib openjdk

setup: ## Create the Python env and install deps (uv)
	uv sync

doctor: ## Check toolchain + data are ready
	uv run locus doctor

ingest: ## Index, QC, normalize the VCFs
	uv run locus ingest

annotate: ## Annotate against ClinVar / dbSNP / gnomAD / VEP / PharmCAT
	uv run locus annotate

load: ## Build the DuckDB store
	uv run locus load

pipeline: ## ingest -> annotate -> load
	uv run locus pipeline

serve-mcp: ## Start the MCP server (query with Claude)
	uv run locus serve mcp

serve-api: ## Start the FastAPI backend for the SPA
	uv run locus serve api

web-dev: ## Run the React debug SPA (Vite dev server)
	cd web && npm install && npm run dev

web-build: ## Build the SPA for production
	cd web && npm install && npm run build

lint: ## Lint Python
	uv run ruff check src

clean: ## Remove the DuckDB store and intermediate work (keeps raw inputs)
	rm -f data/locus.duckdb data/locus.duckdb.wal
	rm -rf data/work
