"""Register the Locus MCP server with Claude (Desktop + Code).

Writes the small JSON that points Claude at this repo's ``locus-mcp`` server, using THIS
checkout's absolute path (so ``.mcp.json`` is generated per machine and gitignored, never
committed) and an absolute ``uv`` (Claude Desktop has no shell PATH). Existing config is backed
up and safe-merged — other MCP servers are preserved.

Claude Desktop must be QUIT while its config is written: a running Desktop rewrites
``claude_desktop_config.json`` from memory and silently drops an entry added underneath it.
"""

from __future__ import annotations

import datetime as _dt
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

from rich.console import Console

console = Console()

REPO_ROOT = Path(__file__).resolve().parents[2]


def _uv_path() -> str:
    return shutil.which("uv") or "/opt/homebrew/bin/uv"


def desktop_config_path() -> Path:
    return Path.home() / "Library" / "Application Support" / "Claude" / "claude_desktop_config.json"


def code_config_path() -> Path:
    return REPO_ROOT / ".mcp.json"


def _server_entry() -> dict:
    return {
        "type": "stdio",
        "command": _uv_path(),
        "args": ["run", "--directory", str(REPO_ROOT), "locus-mcp"],
    }


def _backup(path: Path) -> Path:
    # microseconds keep rapid re-runs from colliding on the same backup filename
    bak = path.with_name(path.name + f".bak-{_dt.datetime.now():%Y%m%d-%H%M%S-%f}")
    shutil.copy2(path, bak)
    return bak


def _install_into(path: Path, *, label: str) -> bool:
    """Safe-merge the ``locus`` server into a Claude JSON config. Returns True on success."""
    path.parent.mkdir(parents=True, exist_ok=True)
    data: dict = {}
    if path.exists():
        try:
            data = json.loads(path.read_text() or "{}")
            if not isinstance(data, dict):
                raise ValueError("config is not a JSON object")
        except (json.JSONDecodeError, ValueError) as e:
            bak = _backup(path)
            console.print(f"[yellow]{label}: existing config isn't valid JSON ({e}); backed it up to "
                          f"{bak.name} and left it unchanged. Fix or delete it, then re-run "
                          f"`locus mcp install`.[/]")
            return False
        _backup(path)  # preserve the prior good config before editing
    data.setdefault("mcpServers", {})["locus"] = _server_entry()
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    tmp.replace(path)
    console.print(f"[green]{label}: registered[/] → {path}")
    return True


def desktop_running() -> bool:
    """Is the Claude Desktop app open? (It owns its config while running.)"""
    try:
        return subprocess.run(["pgrep", "-f", "Claude.app/Contents/MacOS/Claude"],
                              capture_output=True).returncode == 0
    except OSError:
        return False


def _wait_for_desktop_quit(timeout_s: int = 300) -> bool:
    """Ask the user to quit Claude Desktop and wait for it. False if it's still open (or non-interactive)."""
    if not desktop_running():
        return True
    if not sys.stdin.isatty():
        return False
    console.print("[yellow]Claude Desktop is open. It rewrites its config while running and would silently "
                  "drop the Locus entry.[/]\n[bold]Quit Claude Desktop now (Cmd-Q, not just close the "
                  "window)[/] — waiting…")
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        time.sleep(2)
        if not desktop_running():
            time.sleep(1)  # let it finish its own last write before we touch the file
            return True
    return False


def run() -> bool:
    """Register the MCP server with Claude Code, and with Claude Desktop once it's quit.

    Returns True if Desktop was registered (False if it was left open and skipped).
    """
    console.rule("[bold]Register MCP server with Claude")
    desktop_ok = False
    if _wait_for_desktop_quit():
        desktop_ok = _install_into(desktop_config_path(), label="Claude Desktop")
    else:
        console.print("[yellow]Claude Desktop: NOT registered — it's still open, and it would overwrite the "
                      "entry. Quit it (Cmd-Q), then run `locus mcp install` again.[/]")
    _install_into(code_config_path(), label="Claude Code (.mcp.json)")
    if desktop_ok:
        console.print("[dim]Open Claude Desktop to load the server. "
                      "Claude Code picks up .mcp.json automatically in this folder.[/]")
    return desktop_ok


def _registered_dir(path: Path) -> str | None:
    """The ``--directory`` recorded for the locus server in a config, or None if absent."""
    if not path.exists():
        return None
    try:
        args = json.loads(path.read_text() or "{}").get("mcpServers", {}).get("locus", {}).get("args", [])
        if "--directory" in args:
            return args[args.index("--directory") + 1]
    except (json.JSONDecodeError, ValueError, IndexError, AttributeError):
        return None
    return None


def is_registered() -> dict:
    """For `doctor`: per-target registration + whether it points at THIS repo."""
    want = str(REPO_ROOT)
    out: dict = {}
    for key, path in (("desktop", desktop_config_path()), ("code", code_config_path())):
        d = _registered_dir(path)
        out[key] = d is not None
        out[f"{key}_matches"] = d == want
    return out


def status() -> None:
    info = is_registered()
    console.print(f"repo: {REPO_ROOT}")
    for label, key, path in (("Claude Desktop", "desktop", desktop_config_path()),
                             ("Claude Code", "code", code_config_path())):
        if info[key]:
            state = "[green]registered[/]" if info[f"{key}_matches"] else \
                    "[yellow]registered but points at a different folder[/]"
        else:
            state = "[yellow]not registered[/]"
        console.print(f"  {label}: {state}  ({path})")
    if not all(info[k] and info[f"{k}_matches"] for k in ("desktop", "code")):
        console.print("[dim]Run `locus mcp install` to (re)register for this checkout.[/]")
