#!/bin/bash
# Build & install Locus.app — a native macOS Dock app that launches the Locus SPA/API
# server (uv run locus serve api), opens the browser, and stops the server when you quit.
#
# Usage:  scripts/build_macos_app.sh
# Installs to /Applications (falls back to ~/Applications). Requires: uv, Node/npm (to build the
# web UI the app serves), swiftc (Xcode command-line tools), and macOS iconutil. The .app itself is machine-local (not in git);
# this script reproduces it from the committed sources.
set -euo pipefail

DIR="$(cd "$(dirname "$0")/.." && pwd)"   # repo root
UV="$(command -v uv || echo /opt/homebrew/bin/uv)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "1/6  building the web UI (web/dist)…"
# The app opens the SPA that `locus serve api` serves from web/dist — without it the app opens a 404.
NPM="$(command -v npm || true)"
[ -z "$NPM" ] && [ -x /opt/homebrew/bin/npm ] && NPM=/opt/homebrew/bin/npm
[ -z "$NPM" ] && [ -x /usr/local/bin/npm ] && NPM=/usr/local/bin/npm
if [ -n "$NPM" ]; then
  (cd "$DIR/web" && PATH="$(dirname "$NPM"):$PATH" "$NPM" install --no-audit --no-fund --loglevel=error \
     && PATH="$(dirname "$NPM"):$PATH" "$NPM" run build --silent) >/dev/null
elif [ -f "$DIR/web/dist/index.html" ]; then
  echo "     npm not found — using the existing web/dist build"
else
  echo "✗ The Locus app serves a web UI that has to be built with Node.js, which isn't installed."
  echo "  Install it (brew install node), then re-run scripts/build_macos_app.sh."
  exit 1
fi

echo "2/6  rendering DNA icon…"
uv run --with pillow python "$DIR/scripts/make_app_icon.py" "$WORK" >/dev/null
iconutil -c icns "$WORK/Locus.iconset" -o "$WORK/Locus.icns"

echo "3/6  compiling native app (arm64/x86_64 per this Mac)…"
sed -e "s#@REPO@#$DIR#g" -e "s#@UV@#$UV#g" "$DIR/scripts/locus_app.swift" > "$WORK/main.swift"
compile() { swiftc -O "$WORK/main.swift" -o "$WORK/Locus" "$@"; }
if ! compile 2>"$WORK/swiftc.err"; then
  # Common right after a macOS update: the default SDK is newer than the Command Line Tools'
  # Swift compiler ("this SDK is not supported by the compiler"). Retry with older installed SDKs.
  echo "     default SDK failed; trying other installed macOS SDKs…"
  built=""
  for sdkdir in /Library/Developer/CommandLineTools/SDKs \
                "$(xcode-select -p 2>/dev/null)/Platforms/MacOSX.platform/Developer/SDKs"; do
    [ -d "$sdkdir" ] || continue
    for name in $(ls "$sdkdir" | grep -E '^MacOSX[0-9.]+\.sdk$' | sort -rV); do
      if compile -sdk "$sdkdir/$name" 2>/dev/null; then built="$name"; break 2; fi
    done
  done
  if [ -z "$built" ]; then
    cat "$WORK/swiftc.err" >&2
    echo "✗ swiftc couldn't build with any installed SDK. Usually the Command Line Tools are older than"
    echo "  the macOS SDK. Fix: update them (softwareupdate --list, or xcode-select --install), or point at"
    echo "  an SDK your compiler supports:  SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX<ver>.sdk \\"
    echo "    scripts/build_macos_app.sh"
    exit 1
  fi
  echo "     built with $built"
fi

echo "4/6  choosing install location…"
if [ -w /Applications ]; then APPDIR="/Applications"; else APPDIR="$HOME/Applications"; mkdir -p "$APPDIR"; fi
APP="$APPDIR/Locus.app"

echo "5/6  assembling bundle → $APP"
rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
cp "$WORK/Locus" "$APP/Contents/MacOS/Locus"; chmod +x "$APP/Contents/MacOS/Locus"
cp "$WORK/Locus.icns" "$APP/Contents/Resources/Locus.icns"
cat > "$APP/Contents/Info.plist" <<'PLIST'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>Locus</string>
  <key>CFBundleDisplayName</key><string>Locus</string>
  <key>CFBundleIdentifier</key><string>com.locus.app</string>
  <key>CFBundleVersion</key><string>1.0</string>
  <key>CFBundleShortVersionString</key><string>1.0</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleExecutable</key><string>Locus</string>
  <key>CFBundleIconFile</key><string>Locus</string>
  <key>LSMinimumSystemVersion</key><string>12.0</string>
  <key>NSHighResolutionCapable</key><true/>
</dict>
</plist>
PLIST

echo "6/6  signing + registering…"
codesign --force --deep --sign - "$APP"
LSREG="/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister"
[ -x "$LSREG" ] && "$LSREG" -f "$APP" 2>/dev/null || true
touch "$APP"; killall Dock 2>/dev/null || true

echo "✓ Installed $APP"
echo "  Launch it from Finder/Spotlight/Launchpad, or drag it to the Dock."
