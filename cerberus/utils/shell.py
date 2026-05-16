"""Subprocess runner. Streams stdout/stderr to a per-step log file and to the JSONL logger."""
from __future__ import annotations

import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from cerberus.utils.logger import get_logger

log = get_logger("shell")


class ToolError(RuntimeError):
    def __init__(self, cmd: list[str], rc: int, log_path: Path | None):
        self.cmd = cmd
        self.rc = rc
        self.log_path = log_path
        super().__init__(
            f"Command failed (rc={rc}): {' '.join(shlex.quote(c) for c in cmd)}"
            + (f"\nSee log: {log_path}" if log_path else "")
        )


@dataclass
class ToolResult:
    cmd: list[str]
    returncode: int
    log_path: Path


def which(binary: str) -> str | None:
    from shutil import which as _which
    return _which(binary)


def require_tools(*binaries: str) -> None:
    missing = [b for b in binaries if which(b) is None]
    if missing:
        raise RuntimeError(
            f"Missing required tool(s): {', '.join(missing)}. "
            "Install via the cerberus conda environment."
        )


def run(
    cmd: list[str] | str,
    *,
    log_path: Path,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
    check: bool = True,
    timeout: int | None = None,
    dry_run: bool = False,
) -> ToolResult:
    if isinstance(cmd, str):
        cmd_list = shlex.split(cmd)
    else:
        cmd_list = list(cmd)

    pretty = " ".join(shlex.quote(c) for c in cmd_list)
    log.info("RUN  %s", pretty)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    if dry_run:
        log.info("DRY-RUN — skipped")
        log_path.write_text(f"# DRY-RUN\n{pretty}\n")
        return ToolResult(cmd=cmd_list, returncode=0, log_path=log_path)

    with log_path.open("w", encoding="utf-8") as logf:
        logf.write(f"# CMD: {pretty}\n")
        logf.flush()
        proc = subprocess.run(
            cmd_list,
            stdout=logf,
            stderr=subprocess.STDOUT,
            cwd=str(cwd) if cwd else None,
            env=dict(env) if env else None,
            timeout=timeout,
            check=False,
        )

    if check and proc.returncode != 0:
        raise ToolError(cmd_list, proc.returncode, log_path)

    return ToolResult(cmd=cmd_list, returncode=proc.returncode, log_path=log_path)


def pipe(
    cmds: list[list[str]],
    *,
    log_path: Path,
    final_stdout: Path | None = None,
    cwd: Path | None = None,
    dry_run: bool = False,
) -> ToolResult:
    """Run a chain of commands joined by pipes: cmds[0] | cmds[1] | ..."""
    pretty = " | ".join(" ".join(shlex.quote(c) for c in cmd) for cmd in cmds)
    log.info("PIPE %s", pretty)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    if dry_run:
        log_path.write_text(f"# DRY-RUN\n{pretty}\n")
        return ToolResult(cmd=cmds[-1], returncode=0, log_path=log_path)

    procs: list[subprocess.Popen] = []
    stderr_log = log_path.open("w", encoding="utf-8")
    stderr_log.write(f"# PIPE: {pretty}\n")
    stderr_log.flush()

    try:
        prev_stdout = None
        for i, cmd in enumerate(cmds):
            is_last = i == len(cmds) - 1
            stdout: int | object
            if is_last and final_stdout is not None:
                stdout = open(final_stdout, "wb")
            else:
                stdout = subprocess.PIPE if not is_last else subprocess.DEVNULL
            p = subprocess.Popen(
                cmd,
                stdin=prev_stdout,
                stdout=stdout,
                stderr=stderr_log,
                cwd=str(cwd) if cwd else None,
            )
            if prev_stdout is not None:
                prev_stdout.close()
            prev_stdout = p.stdout
            procs.append(p)
        for p in procs:
            p.wait()
    finally:
        stderr_log.close()

    rc = procs[-1].returncode
    if rc != 0:
        raise ToolError(cmds[-1], rc, log_path)
    return ToolResult(cmd=cmds[-1], returncode=rc, log_path=log_path)
