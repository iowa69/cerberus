"""Subprocess runner. Streams stdout/stderr to a per-step log file and to the JSONL logger.

Two entry points:

  ``run()``   one command. stderr (and stdout, unless ``stdout_path`` is
              given) goes to the step log.
  ``pipe()``  a chain joined by pipes. **Every** stage's exit status is
              checked, not just the last one — a dead aligner feeding a
              healthy ``samtools view`` used to look like success and produce
              a silently truncated dataset.

All spawned children are tracked so a Ctrl-C tears down the whole process
group instead of leaving a 32-thread aligner running detached.
"""
from __future__ import annotations

import shlex
import signal
import subprocess
import threading
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

from cerberus.utils.logger import get_logger

log = get_logger("shell")

# Children currently running, so an interrupt can terminate them.
_LIVE: set[subprocess.Popen] = set()
_LIVE_LOCK = threading.Lock()


def _track(proc: subprocess.Popen) -> None:
    with _LIVE_LOCK:
        _LIVE.add(proc)


def _untrack(proc: subprocess.Popen) -> None:
    with _LIVE_LOCK:
        _LIVE.discard(proc)


def terminate_all(sig: int = signal.SIGTERM) -> int:
    """Signal every tracked child. Returns how many were signalled."""
    with _LIVE_LOCK:
        procs = list(_LIVE)
    n = 0
    for p in procs:
        if p.poll() is None:
            try:
                p.send_signal(sig)
                n += 1
            except OSError:
                pass
    for p in procs:
        try:
            p.wait(timeout=10)
        except (subprocess.TimeoutExpired, OSError):
            try:
                p.kill()
            except OSError:
                pass
    return n


class ToolError(RuntimeError):
    def __init__(self, cmd: list[str], rc: int, log_path: Path | None, *, stage: str = ""):
        self.cmd = cmd
        self.rc = rc
        self.log_path = log_path
        where = f" [{stage}]" if stage else ""
        super().__init__(
            f"Command failed (rc={rc}){where}: {' '.join(shlex.quote(c) for c in cmd)}"
            + (f"\nSee log: {log_path}" if log_path else "")
        )


@dataclass
class ToolResult:
    cmd: list[str]
    returncode: int
    log_path: Path
    returncodes: list[int] = field(default_factory=list)


def which(binary: str) -> str | None:
    from shutil import which as _which
    return _which(binary)


def require_tools(*binaries: str) -> None:
    missing = [b for b in binaries if which(b) is None]
    if missing:
        raise RuntimeError(
            f"Missing required tool(s): {', '.join(missing)}. "
            "Install via the cerberus conda environment, or run 'cerberus doctor' "
            "to see the full picture."
        )


def run(
    cmd: list[str] | str,
    *,
    log_path: Path,
    stdout_path: Path | None = None,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
    check: bool = True,
    timeout: int | None = None,
    dry_run: bool = False,
) -> ToolResult:
    """Run one command.

    ``stdout_path`` captures stdout to its own file (used for deliverables
    such as ``samtools flagstat`` output, which must not be polluted by the
    log's ``# CMD:`` header).
    """
    cmd_list = shlex.split(cmd) if isinstance(cmd, str) else list(cmd)

    pretty = " ".join(shlex.quote(c) for c in cmd_list)
    log.info("RUN  %s", pretty)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if stdout_path is not None:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)

    if dry_run:
        log.info("DRY-RUN — skipped")
        log_path.write_text(f"# DRY-RUN\n{pretty}\n")
        if stdout_path is not None:
            stdout_path.write_text(f"# DRY-RUN\n{pretty}\n")
        return ToolResult(cmd=cmd_list, returncode=0, log_path=log_path, returncodes=[0])

    with log_path.open("w", encoding="utf-8") as logf:
        logf.write(f"# CMD: {pretty}\n")
        logf.flush()
        out_ctx = stdout_path.open("wb") if stdout_path is not None else None
        try:
            proc = subprocess.Popen(
                cmd_list,
                stdout=(out_ctx if out_ctx is not None else logf),
                stderr=logf if out_ctx is not None else subprocess.STDOUT,
                cwd=str(cwd) if cwd else None,
                env=dict(env) if env else None,
            )
            _track(proc)
            try:
                proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
                raise ToolError(cmd_list, -signal.SIGKILL, log_path, stage="timeout") from None
            finally:
                _untrack(proc)
        finally:
            if out_ctx is not None:
                out_ctx.close()

    if check and proc.returncode != 0:
        raise ToolError(cmd_list, proc.returncode, log_path)

    return ToolResult(cmd=cmd_list, returncode=proc.returncode, log_path=log_path,
                      returncodes=[proc.returncode])


def pipe(
    cmds: list[list[str]],
    *,
    log_path: Path,
    final_stdout: Path | None = None,
    cwd: Path | None = None,
    timeout: int | None = None,
    dry_run: bool = False,
) -> ToolResult:
    """Run a chain of commands joined by pipes: cmds[0] | cmds[1] | ...

    Every stage's exit status is checked. A non-zero status anywhere raises
    ``ToolError`` naming the stage that actually failed.
    """
    if not cmds:
        raise ValueError("pipe() needs at least one command")

    pretty = " | ".join(" ".join(shlex.quote(c) for c in cmd) for cmd in cmds)
    log.info("PIPE %s", pretty)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    if dry_run:
        log_path.write_text(f"# DRY-RUN\n{pretty}\n")
        if final_stdout is not None:
            final_stdout.parent.mkdir(parents=True, exist_ok=True)
            final_stdout.write_bytes(b"")
        return ToolResult(cmd=cmds[-1], returncode=0, log_path=log_path,
                          returncodes=[0] * len(cmds))

    if final_stdout is not None:
        final_stdout.parent.mkdir(parents=True, exist_ok=True)

    procs: list[subprocess.Popen] = []
    out_file = None
    with log_path.open("w", encoding="utf-8") as stderr_log:
        stderr_log.write(f"# PIPE: {pretty}\n")
        stderr_log.flush()
        try:
            if final_stdout is not None:
                out_file = final_stdout.open("wb")

            prev_stdout = None
            for i, cmd in enumerate(cmds):
                is_last = i == len(cmds) - 1
                if is_last:
                    stdout = out_file if out_file is not None else subprocess.DEVNULL
                else:
                    stdout = subprocess.PIPE
                p = subprocess.Popen(
                    cmd,
                    stdin=prev_stdout,
                    stdout=stdout,
                    stderr=stderr_log,
                    cwd=str(cwd) if cwd else None,
                )
                _track(p)
                # The parent must close its copy of the upstream pipe so the
                # downstream process sees EOF when the producer exits.
                if prev_stdout is not None:
                    prev_stdout.close()
                prev_stdout = p.stdout
                procs.append(p)

            for p in procs:
                try:
                    p.wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    for q in procs:
                        if q.poll() is None:
                            q.kill()
                    for q in procs:
                        q.wait()
                    raise ToolError(cmds[-1], -signal.SIGKILL, log_path,
                                    stage="timeout") from None
        finally:
            for p in procs:
                if p.poll() is None:
                    p.kill()
                    p.wait()
                _untrack(p)
                if p.stdout is not None and not p.stdout.closed:
                    p.stdout.close()
            if out_file is not None:
                out_file.close()

    codes = [p.returncode for p in procs]

    # SIGPIPE (rc -13 / 141) on a producer is normal when a consumer exits
    # early on purpose; it is not normal here because every consumer we build
    # reads to EOF. Report any non-zero status, naming the stage.
    for i, (cmd, rc) in enumerate(zip(cmds, codes)):
        if rc != 0:
            stage = f"stage {i + 1}/{len(cmds)}: {cmd[0]}"
            raise ToolError(cmd, rc, log_path, stage=stage)

    return ToolResult(cmd=cmds[-1], returncode=codes[-1], log_path=log_path,
                      returncodes=codes)
