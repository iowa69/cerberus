"""Logging setup. Two handlers: pretty stderr for humans, JSONL file for machines."""
from __future__ import annotations

import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path


class _JsonlFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        payload = {
            "ts": datetime.now(timezone.utc).isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "msg": record.getMessage(),
        }
        if record.exc_info:
            payload["exc"] = self.formatException(record.exc_info)
        return json.dumps(payload, ensure_ascii=False)


class _ConsoleFormatter(logging.Formatter):
    RESET = "\033[0m"
    COLORS = {
        "DEBUG": "\033[2;37m",
        "INFO": "\033[36m",
        "WARNING": "\033[33m",
        "ERROR": "\033[31m",
        "CRITICAL": "\033[1;31m",
    }

    def __init__(self, color: bool):
        super().__init__()
        self.color = color and sys.stderr.isatty()

    def format(self, record: logging.LogRecord) -> str:
        ts = datetime.now().strftime("%H:%M:%S")
        prefix = f"[{ts}] {record.levelname:<7s} {record.name}"
        msg = record.getMessage()
        if self.color:
            c = self.COLORS.get(record.levelname, "")
            return f"{c}{prefix}{self.RESET}  {msg}"
        return f"{prefix}  {msg}"


def setup_logging(log_dir: Path, verbose: bool = False, quiet: bool = False) -> Path:
    log_dir.mkdir(parents=True, exist_ok=True)
    logfile = log_dir / "cerberus.log.jsonl"

    root = logging.getLogger()
    root.handlers.clear()
    root.setLevel(logging.DEBUG)

    console = logging.StreamHandler(sys.stderr)
    console.setLevel(logging.WARNING if quiet else (logging.DEBUG if verbose else logging.INFO))
    console.setFormatter(_ConsoleFormatter(color=True))
    root.addHandler(console)

    fh = logging.FileHandler(logfile, mode="a", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(_JsonlFormatter())
    root.addHandler(fh)

    return logfile


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(f"cerberus.{name}")
