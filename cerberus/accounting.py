"""Per-stage read accounting.

A reviewer wants to see exactly how many reads survive each stage, so we keep
a structured TSV + JSON next to the outputs. The orchestrator builds this
incrementally and writes both formats at the end.

Counting convention
-------------------
Every ``reads`` value is a count of FASTQ **records**, not fragments. For
paired data that means a surviving pair contributes 2 to a merged file and 1
to each of R1/R2. The ``unit`` column states which is which so the numbers
can be reconciled rather than guessed at.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

from cerberus.utils.logger import get_logger

log = get_logger("accounting")


@dataclass
class StageCount:
    mode: str
    stage: str
    reads: int
    unit: str = "records"
    file: str = ""


@dataclass
class RunAccounting:
    sample_id: str
    input_r1_reads: int = 0
    input_r2_reads: int = 0
    input_long_reads: int = 0
    qc_paired: int = 0
    qc_orphans: int = 0
    qc_long: int = 0
    stages: list[StageCount] = field(default_factory=list)
    final_outputs: dict[str, dict] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def add_stage(
        self, mode: str, stage: str, reads: int,
        file: Path | None = None, unit: str = "records",
    ) -> None:
        self.stages.append(StageCount(
            mode=mode, stage=stage, reads=reads, unit=unit,
            file=str(file) if file else "",
        ))

    def warn(self, message: str) -> None:
        log.warning(message)
        self.warnings.append(message)

    def add_final(self, mode: str, paths: dict[str, Path | None]) -> None:
        from cerberus.utils.fastq import count_reads
        summary: dict[str, dict] = {}
        for key, path in paths.items():
            if path is None or not path.exists():
                continue
            summary[key] = {
                "path": str(path),
                "reads": count_reads(path),
                "size_bytes": path.stat().st_size,
            }
        self.final_outputs[mode] = summary

    def stages_for(self, mode: str) -> list[StageCount]:
        return [s for s in self.stages if s.mode == mode]

    def write(self, out_dir: Path) -> tuple[Path, Path]:
        out_dir.mkdir(parents=True, exist_ok=True)
        json_path = out_dir / "accounting.json"
        tsv_path = out_dir / "accounting.tsv"

        json_path.write_text(json.dumps(asdict(self), indent=2))

        lines = ["sample\tmode\tstage\treads\tunit\tfile"]

        def row(mode: str, stage: str, reads: int, unit: str, file: str = "") -> None:
            lines.append(f"{self.sample_id}\t{mode}\t{stage}\t{reads}\t{unit}\t{file}")

        if self.input_r1_reads or self.input_r2_reads:
            row("_input", "input_r1", self.input_r1_reads, "records")
            row("_input", "input_r2", self.input_r2_reads, "records")
        if self.input_long_reads:
            row("_input", "input_long", self.input_long_reads, "records")
        if self.qc_paired:
            row("_qc", "qc_paired_per_mate", self.qc_paired, "records")
        if self.qc_orphans:
            row("_qc", "qc_orphans", self.qc_orphans, "records")
        if self.qc_long:
            row("_qc", "qc_long", self.qc_long, "records")
        for s in self.stages:
            row(s.mode, s.stage, s.reads, s.unit, s.file)
        for mode, files in self.final_outputs.items():
            for key, info in files.items():
                row("_final", f"{mode}.{key}", info["reads"], "records", info["path"])

        tsv_path.write_text("\n".join(lines) + "\n")

        log.info("Wrote accounting: %s and %s", json_path.name, tsv_path.name)
        return json_path, tsv_path
