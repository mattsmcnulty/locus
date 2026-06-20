"""Schedule ``locus refresh`` via a macOS launchd LaunchAgent (weekly).

launchd is the OS-correct scheduler on Apple Silicon: a ``StartCalendarInterval``
job runs on the next wake if the Mac was asleep (cron silently skips). The catch —
which CLAUDE.md already notes for the Claude Desktop config — is that launchd has
**no login PATH**, so the plist must invoke an absolute ``uv`` and set ``PATH`` /
``JAVA_HOME`` / ``LOCUS_DATA_DIR`` explicitly, or bcftools/SnpEff/uv won't resolve.
"""

from __future__ import annotations

import plistlib
import shutil
import subprocess
from pathlib import Path

from rich.console import Console

from . import shell
from .config import settings

console = Console()

LABEL = "com.locus.refresh"
_REPO_ROOT = Path(__file__).resolve().parents[2]


def _plist_path() -> Path:
    return Path.home() / "Library" / "LaunchAgents" / f"{LABEL}.plist"


def _uv_path() -> str:
    return shutil.which("uv") or "/opt/homebrew/bin/uv"


def _launchd_env() -> dict[str, str]:
    """The minimal environment a non-login launchd job needs to find our tools."""
    path_dirs = ["/opt/homebrew/bin", "/opt/homebrew/sbin", "/usr/local/bin",
                 "/usr/bin", "/bin", "/usr/sbin", "/sbin"]
    env = {"PATH": ":".join(path_dirs), "LOCUS_DATA_DIR": str(settings.data_dir)}
    java = shell.resolve_java()
    if java:  # JAVA_HOME = <jdk>/bin/java -> <jdk>
        env["JAVA_HOME"] = str(Path(java).resolve().parent.parent)
    return env


def build_plist(weekday: int, hour: int) -> dict:
    log_dir = settings.reports_dir
    return {
        "Label": LABEL,
        "ProgramArguments": [_uv_path(), "run", "--directory", str(_REPO_ROOT), "locus", "refresh"],
        "WorkingDirectory": str(_REPO_ROOT),
        "EnvironmentVariables": _launchd_env(),
        # Weekly: Weekday 0=Sunday. Runs at <hour>:00 local time.
        "StartCalendarInterval": {"Weekday": weekday, "Hour": hour, "Minute": 0},
        "RunAtLoad": False,
        "StandardOutPath": str(log_dir / "refresh.out.log"),
        "StandardErrorPath": str(log_dir / "refresh.err.log"),
    }


def install(weekday: int = 0, hour: int = 3) -> Path:
    """Write + load the LaunchAgent (weekly, default Sunday 03:00)."""
    settings.reports_dir.mkdir(parents=True, exist_ok=True)
    p = _plist_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("wb") as fh:
        plistlib.dump(build_plist(weekday, hour), fh)
    # Reload (ignore "not loaded" on first install).
    subprocess.run(["launchctl", "unload", str(p)], capture_output=True)
    res = subprocess.run(["launchctl", "load", str(p)], capture_output=True, text=True)
    if res.returncode != 0:
        console.print(f"[yellow]launchctl load warning:[/] {res.stderr.strip()}")
    days = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"]
    console.print(f"[green]Scheduled[/] weekly `locus refresh` — {days[weekday % 7]} {hour:02d}:00 → {p}")
    console.print(f"[dim]Logs: {settings.reports_dir}/refresh.{{out,err}}.log[/]")
    return p


def uninstall() -> None:
    p = _plist_path()
    if not p.exists():
        console.print("[yellow]No Locus schedule installed.[/]")
        return
    subprocess.run(["launchctl", "unload", str(p)], capture_output=True)
    p.unlink(missing_ok=True)
    console.print(f"[green]Removed[/] the Locus refresh schedule ({p}).")


def status() -> None:
    p = _plist_path()
    if not p.exists():
        console.print("[yellow]Not scheduled.[/] Install with `locus schedule install`.")
        return
    res = subprocess.run(["launchctl", "list", LABEL], capture_output=True, text=True)
    loaded = res.returncode == 0
    console.print(f"plist  : {p}")
    console.print(f"loaded : {'[green]yes[/]' if loaded else '[yellow]no (run `locus schedule install`)[/]'}")
    if loaded:
        console.print(res.stdout.strip())
