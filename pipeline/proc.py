"""Running external tools as plain command lines.

Every external tool in this project is invoked through run_tool(); nothing is
imported as a library. A non-zero exit becomes TOOL_FAILED, a timeout kills the
whole process group and raises ToolTimeout (callers map it to a specific code).
"""
import logging
import os
import re
import signal
import subprocess
from pathlib import Path

from config import settings
from pipeline.errors import TOOL_FAILED, PipelineError

log = logging.getLogger("amr.proc")


class ToolTimeout(Exception):
    pass


def run_tool(cmd: list[str], *, timeout: int, log_path: Path | None = None, cwd: Path | None = None,
             name: str | None = None) -> str:
    """Run cmd, return stdout. stdout+stderr are also written to log_path (server-side only)."""
    name = name or Path(cmd[0]).name
    try:
        proc = subprocess.Popen(
            [str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            cwd=str(cwd) if cwd else None, env=settings.tool_env(), start_new_session=True,
        )
    except OSError as exc:
        log.error("cannot start %s: %s", name, exc)
        raise PipelineError(TOOL_FAILED, f"{name} could not be started") from None
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        _kill_group(proc)
        proc.communicate()
        raise ToolTimeout(name) from None
    except BaseException:
        _kill_group(proc)
        raise
    if log_path is not None:
        try:
            Path(log_path).parent.mkdir(parents=True, exist_ok=True)
            Path(log_path).write_text(f"$ {' '.join(map(str, cmd))}\n--- stdout\n{out}\n--- stderr\n{err}\n")
        except OSError:
            pass
    if proc.returncode != 0:
        tail = (err or out or "").strip().splitlines()[-5:]
        log.error("%s exited %s: %s", name, proc.returncode, " | ".join(tail))
        raise PipelineError(TOOL_FAILED, f"{name} exited with status {proc.returncode}")
    return out


def _kill_group(proc: subprocess.Popen) -> None:
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        proc.kill()


def tool_version(cmd: list[str], pattern: str = r"(\d+(?:\.\d+)+)") -> str:
    """Best-effort version string; never raises."""
    try:
        res = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=60,
                             env=settings.tool_env())
        m = re.search(pattern, (res.stdout or "") + "\n" + (res.stderr or ""))
        return m.group(1) if m else "unknown"
    except Exception:
        return "unknown"
