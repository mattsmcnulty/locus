#!/bin/bash
# Locus — launch the local genome explorer (API + SPA) and open it in your browser.
# Double-click this file (or the Desktop "Locus" launcher). Press Ctrl-C in this window to stop.
set -e
DIR="$(cd "$(dirname "$0")/.." && pwd)"   # repo root (this script lives in scripts/)
cd "$DIR"
UV="$(command -v uv || echo /opt/homebrew/bin/uv)"
echo "Starting Locus from $DIR …"
echo "Opening http://127.0.0.1:8787 — press Ctrl-C here to stop the server."
( sleep 3; open "http://127.0.0.1:8787" ) &   # open the browser once the server is up
exec "$UV" run locus serve api
