"""Per-stage read accounting.

A reviewer wants to see exactly how many reads survive each stage, so we keep
a structured TSV + JSON next to the outputs. The orchestrator builds this
incrementally and writes both formats at the end.
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
    file: str = ""


@dataclass
class RunAccounting:
    sample_id: str
    input_r1_reads: int = 0
    input_r2_reads: int = 0
    input_long_reads: int = 0
    qc_paired: int = 0
    qc_orphans: int = 0
    stages: list[StageCount] = field(default_factory=list)
    final_outputs: dict[str, dict] = field(default_factory=dict)

    def add_stage(self, mode: str, stage: str, reads: int, file: Path | None = None) -> None:
        self.stages.append(StageCount(
            mode=mode, stage=stage, reads=reads,
            file=str(file) if file else "",
        ))

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

    def write(self, out_dir: Path) -> tuple[Path, Path]:
        json_path = out_dir / "accounting.json"
        tsv_path = out_dir / "accounting.tsv"

        json_path.write_text(json.dumps(asdict(self), indent=2))

        lines = ["mode\tstage\treads\tfile"]
        lines.append(f"_input\tinput_r1\t{self.input_r1_reads}\t")
        if self.input_r2_reads:
            lines.append(f"_input\tinput_r2\t{self.input_r2_reads}\t")
        if self.input_long_reads:
            lines.append(f"_input\tinput_long\t{self.input_long_reads}\t")
        lines.append(f"_qc\tqc_paired\t{self.qc_paired}\t")
        if self.qc_orphans:
            lines.append(f"_qc\tqc_orphans\t{self.qc_orphans}\t")
        for s in self.stages:
            lines.append(f"{s.mode}\t{s.stage}\t{s.reads}\t{s.file}")
        for mode, files in self.final_outputs.items():
            for key, info in files.items():
                lines.append(f"_final\t{mode}.{key}\t{info['reads']}\t{info['path']}")
        tsv_path.write_text("\n".join(lines) + "\n")

        log.info("Wrote accounting: %s and %s", json_path.name, tsv_path.name)
        return json_path, tsv_path
