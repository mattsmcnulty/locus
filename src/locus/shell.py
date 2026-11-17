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
# Apple's launcher stub: it forwards to whatever JDK $JAVA_HOME / java_home names, so it is never
# itself a JDK home — and with a bogus JAVA_HOME (e.g. /usr) it can hang instead of failing.
_MACOS_JAVA_STUB = Path("/usr/bin/java")


def _java_candidates() -> list[str]:
    """Where to look for java, most-deliberate first.

    An explicit ``$JAVA_HOME`` wins, then macOS's registered default (``/usr/libexec/java_home``),
    then the Homebrew keg-only openjdk, and only then whatever ``java`` is on PATH (which on macOS
    is usually the stub).
    """
    import os
    import shutil

    cands: list[str] = []
    if home := os.environ.get("JAVA_HOME"):
        cands.append(str(Path(home) / "bin" / "java"))
    java_home_tool = Path("/usr/libexec/java_home")
    if java_home_tool.exists():
        try:
            out = subprocess.run([str(java_home_tool)], capture_output=True, text=True, timeout=10)
            if out.returncode == 0 and out.stdout.strip():
                cands.append(str(Path(out.stdout.strip()) / "bin" / "java"))
        except (OSError, subprocess.SubprocessError):
            pass
    cands.extend(_JAVA_CANDIDATES)
    if on_path := shutil.which("java"):
        cands.append(on_path)
    return cands


def _java_version(java: str) -> str | None:
    """First line of ``java -version`` if this binary actually runs a JVM, else None."""
    try:
        out = subprocess.run([java, "-version"], capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    text = (out.stderr or out.stdout).strip()
    if out.returncode != 0 or "version" not in text.lower() or "unable to locate" in text.lower():
        return None
    return text.splitlines()[0]


def resolve_java() -> str | None:
    """Return a *working* java binary.

    On macOS, ``/usr/bin/java`` is a stub that errors ("Unable to locate a Java
    Runtime") unless a JDK is installed, and the Homebrew openjdk is keg-only —
    so we test each candidate actually runs before returning it.
    """
    seen: set[str] = set()
    for cand in _java_candidates():
        if cand in seen or not Path(cand).exists():
            continue
        seen.add(cand)
        if _java_version(cand):
            return cand
    return None


def java_version() -> str | None:
    """``java -version`` banner of the java Locus will use (for `locus doctor`)."""
    java = resolve_java()
    return _java_version(java) if java else None


def java_home(java: str | None) -> str | None:
    """The JDK home a java binary belongs to, suitable for ``JAVA_HOME`` — or None.

    Resolves symlinks (Homebrew's ``opt/openjdk/bin/java`` points into the real
    ``…/Contents/Home``). Returns None for the macOS stub, whose ``parent.parent``
    is ``/usr`` — exporting that as JAVA_HOME makes the stub hang.
    """
    if not java:
        return None
    real = Path(java).resolve()
    if real == _MACOS_JAVA_STUB or Path(java) == _MACOS_JAVA_STUB:
        return None
    home = real.parent.parent
    return str(home) if (home / "bin" / "java").exists() else None


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
