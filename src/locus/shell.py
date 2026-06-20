"""Thin subprocess helpers with consistent logging and error handling."""

from __future__ import annotations

import subprocess
from pathlib import Path

from rich.console import Console

console = Console()


class ToolError(RuntimeError):
    """A shelled-out command failed."""


def run(args: list[str], *, capture: bool = False, quiet: bool = False) -> subprocess.CompletedProcess:
    """Run a command from an argv list (no shell)."""
    if not quiet:
        console.log(f"[dim]$ {' '.join(str(a) for a in args)}[/]")
    proc = subprocess.run(
        [str(a) for a in args],
        capture_output=capture,
        text=True,
    )
    if proc.returncode != 0:
        err = proc.stderr if capture else ""
        raise ToolError(f"command failed ({proc.returncode}): {' '.join(str(a) for a in args)}\n{err}")
    return proc


def run_env(args: list[str], *, env: dict, quiet: bool = False) -> subprocess.CompletedProcess:
    """Like :func:`run` but with a custom environment (e.g. PATH/JAVA_HOME for PharmCAT)."""
    if not quiet:
        console.log(f"[dim]$ {' '.join(str(a) for a in args)}[/]")
    proc = subprocess.run([str(a) for a in args], env=env, text=True)
    if proc.returncode != 0:
        raise ToolError(f"command failed ({proc.returncode}): {' '.join(str(a) for a in args)}")
    return proc


def sh(cmd: str, *, quiet: bool = False) -> subprocess.CompletedProcess:
    """Run a shell pipeline string (for bcftools `|` chains). Uses bash, errexit on pipefail."""
    if not quiet:
        console.log(f"[dim]$ {cmd}[/]")
    proc = subprocess.run(
        ["bash", "-o", "pipefail", "-c", cmd],
        text=True,
    )
    if proc.returncode != 0:
        raise ToolError(f"pipeline failed ({proc.returncode}): {cmd}")
    return proc


def capture(args: list[str]) -> str:
    """Run a command and return stdout (raises on failure)."""
    proc = run(args, capture=True, quiet=True)
    return proc.stdout


_JAVA_CANDIDATES = (
    "/opt/homebrew/opt/openjdk/bin/java",
    "/usr/local/opt/openjdk/bin/java",
)


def resolve_java() -> str | None:
    """Return a *working* java binary.

    On macOS, ``/usr/bin/java`` is a stub that errors ("Unable to locate a Java
    Runtime") unless a JDK is installed, and the Homebrew openjdk is keg-only —
    so we test each candidate actually runs before returning it.
    """
    import shutil

    for cand in (shutil.which("java"), *_JAVA_CANDIDATES):
        if not cand or not Path(cand).exists():
            continue
        try:
            out = subprocess.run([cand, "-version"], capture_output=True, text=True, timeout=10)
            text = (out.stderr or out.stdout).lower()
            if "version" in text and "unable to locate" not in text:
                return cand
        except (OSError, subprocess.SubprocessError):
            continue
    return None


def java_cmd(args: list[str]) -> list[str]:
    """Build a ``java ...`` argv using a working java binary (falls back to 'java')."""
    return [resolve_java() or "java", *args]


def have(tool: str) -> bool:
    """Is a tool resolvable, accounting for Homebrew keg-only paths for java?"""
    import shutil

    if tool == "java":
        return resolve_java() is not None
    return shutil.which(tool) is not None


def docker_ready() -> bool:
    """Is the Docker CLI present AND the daemon actually running? (PharmCAT needs both.)"""
    import shutil

    if not shutil.which("docker"):
        return False
    try:
        return subprocess.run(["docker", "info"], capture_output=True, timeout=8).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False
