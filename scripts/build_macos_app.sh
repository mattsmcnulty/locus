#!/bin/bash
# Build & install Locus.app — a native macOS Dock app that launches the Locus SPA/API
# server (uv run locus serve api), opens the browser, and stops the server when you quit.
#
# Usage:  scripts/build_macos_app.sh
# Installs to /Applications (falls back to ~/Applications). Requires: uv, swiftc (Xcode
# command-line tools), and macOS iconutil. The .app itself is machine-local (not in git);
# this script reproduces it from the committed sources.
set -euo pipefail

DIR="$(cd "$(dirname "$0")/.." && pwd)"   # repo root
UV="$(command -v uv || echo /opt/homebrew/bin/uv)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "1/5  rendering DNA icon…"
uv run --with pillow python "$DIR/scripts/make_app_icon.py" "$WORK" >/dev/null
iconutil -c icns "$WORK/Locus.iconset" -o "$WORK/Locus.icns"

echo "2/5  compiling native app (arm64/x86_64 per this Mac)…"
sed -e "s#@REPO@#$DIR#g" -e "s#@UV@#$UV#g" "$DIR/scripts/locus_app.swift" > "$WORK/main.swift"
swiftc -O "$WORK/main.swift" -o "$WORK/Locus"

echo "3/5  choosing install location…"
if [ -w /Applications ]; then APPDIR="/Applications"; else APPDIR="$HOME/Applications"; mkdir -p "$APPDIR"; fi
APP="$APPDIR/Locus.app"

echo "4/5  assembling bundle → $APP"
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
  <key>CFBundleIdentifier</key><string>com.mmcnulty.locus</string>
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

echo "5/5  signing + registering…"
codesign --force --deep --sign - "$APP"
LSREG="/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister"
[ -x "$LSREG" ] && "$LSREG" -f "$APP" 2>/dev/null || true
touch "$APP"; killall Dock 2>/dev/null || true

echo "✓ Installed $APP"
echo "  Launch it from Finder/Spotlight/Launchpad, or drag it to the Dock."
