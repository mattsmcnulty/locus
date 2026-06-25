#!/bin/bash
# Locus — one-command setup for Apple Silicon Macs.
# Installs Homebrew + genomics tools + uv, then runs the guided `locus setup`.
# Safe to run repeatedly: every step is detect-then-skip.

cd "$(cd "$(dirname "$0")" && pwd)" || exit 1

say() { printf "\n\033[1m%s\033[0m\n" "$1"; }
die() { printf "\n\033[31m%s\033[0m\n" "$1"; exit 1; }

# 1. Apple Silicon only (the reference panels ship an arm64 PLINK2 binary).
if [ "$(uname -m)" != "arm64" ]; then
  die "Locus's setup targets Apple Silicon Macs (M-series). This looks like an Intel Mac."
fi

# 2. Xcode Command Line Tools (git + compilers Homebrew needs).
if ! xcode-select -p >/dev/null 2>&1; then
  say "Installing Xcode Command Line Tools — click 'Install' in the macOS dialog, wait for it to finish, then run setup again."
  xcode-select --install >/dev/null 2>&1 || true
  exit 0
fi

# 3. Homebrew (and put it on PATH for THIS shell — a fresh install isn't on PATH yet).
if [ -x /opt/homebrew/bin/brew ]; then
  eval "$(/opt/homebrew/bin/brew shellenv)"
elif command -v brew >/dev/null 2>&1; then
  eval "$(brew shellenv)"
else
  say "Installing Homebrew (you may be asked for your Mac password)…"
  NONINTERACTIVE=1 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" \
    || die "Homebrew install failed. Check your internet connection and run setup again."
  eval "$(/opt/homebrew/bin/brew shellenv)"
fi

# 4. System tools + uv (same set as `make tools`).
say "Installing genomics tools (bcftools, samtools, htslib), Java, and uv via Homebrew…"
brew install uv bcftools samtools htslib openjdk \
  || die "Homebrew couldn't install the tools. Check your connection and run setup again."
eval "$(/opt/homebrew/bin/brew shellenv)"   # ensure freshly-installed uv is on PATH
command -v uv >/dev/null 2>&1 || die "uv didn't install correctly. Run setup again."

# 5. Python environment.
say "Setting up the Python environment (uv sync)…"
uv sync || die "Could not set up the Python environment (uv sync failed). Run setup again."

# 6. Hand off to the guided installer (downloads DBs, builds your genome, registers with Claude).
say "Starting the guided Locus setup…"
exec uv run locus setup "$@"
