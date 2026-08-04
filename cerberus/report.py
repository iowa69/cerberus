"""Self-contained HTML run report, written on every Cerberus run.

The report answers the questions a reviewer actually asks after a run:

  * What did I run, on what, with which resolved parameters?
  * Which stages ran, and how many reads survived each one?
  * Which outputs exist, are they valid gzip, and are the pairs synchronised?
  * For the GDPR heads: what did each mechanism actually remove?
  * What went wrong, or looked suspicious?

Everything is inlined — no CDN, no external CSS — so the file can be attached
to an email or dropped into a supplementary bundle and still render.
"""
from __future__ import annotations

import gzip
import html
import json
import platform
import shutil
import subprocess
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from cerberus import __version__
from cerberus.accounting import RunAccounting
from cerberus.config import CerberusConfig
from cerberus.pipelines.base import PipelineResult
from cerberus.utils.fastq import count_reads, is_gzipped
from cerberus.utils.logger import get_logger
from cerberus.utils.shell import which

log = get_logger("report")

_VERSION_PROBES: dict[str, list[str]] = {
    "fastp": ["fastp", "--version"],
    "fastplong": ["fastplong", "--version"],
    "minimap2": ["minimap2", "--version"],
    "bowtie2": ["bowtie2", "--version"],
    "samtools": ["samtools", "--version"],
    "kraken2": ["kraken2", "--version"],
    "bbduk.sh": ["bbduk.sh", "--version"],
    "chopper": ["chopper", "--version"],
    "winnowmap": ["winnowmap", "--version"],
    "pigz": ["pigz", "--version"],
}


def tool_versions() -> dict[str, str]:
    """Best-effort version string for every tool on PATH."""
    out: dict[str, str] = {}
    for name, cmd in _VERSION_PROBES.items():
        path = which(name)
        if not path:
            continue
        try:
            proc = subprocess.run(cmd, capture_output=True, timeout=15, check=False)
            text = (proc.stdout or proc.stderr).decode(errors="replace")
            first = next((ln.strip() for ln in text.splitlines() if ln.strip()), "")
            out[name] = first[:120] or "(unknown)"
        except (OSError, subprocess.SubprocessError):
            out[name] = "(probe failed)"
    return out


def _gzip_ok(path: Path) -> bool | None:
    """True/False for gzip validity; None when the file is not gzip at all."""
    if not path.exists() or not is_gzipped(path):
        return None
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 20):
                pass
        return True
    except (OSError, EOFError):
        return False


def _fmt_int(n: Any) -> str:
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return str(n)


def _fmt_bytes(n: int) -> str:
    step = 1024.0
    val = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if val < step:
            return f"{val:.0f} {unit}" if unit == "B" else f"{val:.1f} {unit}"
        val /= step
    return f"{val:.1f} PB"


def _pct(part: float, whole: float) -> str:
    return f"{100.0 * part / whole:.2f}%" if whole else "—"


def collect_outputs(
    cfg: CerberusConfig,
    results: list[PipelineResult],
    gdpr_outputs: dict[str, dict[str, Path | None]],
) -> list[dict[str, Any]]:
    """Verify every published output: existence, size, records, gzip validity."""
    rows: list[dict[str, Any]] = []
    seen: set[Path] = set()

    def add(group: str, role: str, path: Path | None) -> None:
        if path is None or path in seen:
            return
        seen.add(path)
        exists = path.exists()
        rows.append({
            "group": group,
            "role": role,
            "path": str(path),
            "name": path.name,
            "exists": exists,
            "size": path.stat().st_size if exists else 0,
            "reads": count_reads(path) if exists else 0,
            "gzip_ok": _gzip_ok(path) if exists else None,
        })

    for r in results:
        add(r.mode, "R1", r.paired_r1)
        add(r.mode, "R2", r.paired_r2)
        add(r.mode, "unpaired", r.singletons)
        add(r.mode, "reads", r.long_reads)
    for mode, paths in gdpr_outputs.items():
        for role, p in paths.items():
            add(mode, role, p)

    # Pair-synchronisation check per group.
    by_group: dict[str, dict[str, dict]] = {}
    for row in rows:
        by_group.setdefault(row["group"], {})[row["role"]] = row
    for group, roles in by_group.items():
        r1, r2 = roles.get("R1"), roles.get("R2")
        if r1 and r2 and r1["exists"] and r2["exists"]:
            synced = r1["reads"] == r2["reads"]
            r1["paired_ok"] = synced
            r2["paired_ok"] = synced
    return rows


def build_context(
    cfg: CerberusConfig,
    *,
    accounting: RunAccounting,
    results: list[PipelineResult],
    gdpr_outputs: dict[str, dict[str, Path | None]],
    refs: Any,
    elapsed_sec: float,
) -> dict[str, Any]:
    prescan = cfg.prescan
    prescan_d = asdict(prescan) if is_dataclass(prescan) and not isinstance(prescan, type) else {}

    outputs = collect_outputs(cfg, results, gdpr_outputs)

    gdpr_mech: dict[str, dict[str, dict[str, float]]] = {}
    if cfg.gdpr:
        from cerberus.pipelines.gdpr import residual_host_estimate
        for r in results:
            key = f"{r.mode}_gdpr"
            stats = {s.stage: s.reads for s in accounting.stages_for(key)}
            est = residual_host_estimate(stats)
            if est:
                gdpr_mech[key] = est

    return {
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
        "version": __version__,
        "elapsed_sec": elapsed_sec,
        "cfg": cfg,
        "sample_id": cfg.sample_id,
        "command_line": cfg.command_line,
        "modes": (cfg.modes + (["gdpr"] if cfg.gdpr else [])) or ["(none)"],
        "tuned": cfg.tuned.as_dict(),
        "prescan": prescan_d,
        "accounting": accounting,
        "results": results,
        "outputs": outputs,
        "gdpr_mechanisms": gdpr_mech,
        "warnings": accounting.warnings,
        "tools": tool_versions(),
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "processor": platform.machine(),
            "threads_used": cfg.threads,
            "memory_budget_gb": cfg.memory_gb,
            "reference_dir": str(cfg.ref_dir),
            "disk_free": _fmt_bytes(shutil.disk_usage(cfg.out_dir).free)
            if cfg.out_dir.exists() else "unknown",
        },
        "skipped_optional": list(getattr(refs, "skipped_optional", []) or []),
        "dry_run": cfg.dry_run,
    }


# --------------------------------------------------------------------------
# rendering
# --------------------------------------------------------------------------

_CSS = """
:root{--bg:#fbfbfd;--fg:#16181d;--muted:#5d6470;--line:#e3e6ec;--card:#fff;
--ok:#0d7a4a;--warn:#a35a00;--bad:#b3261e;--accent:#2b4c9b;--accent-soft:#eef2fb}
@media(prefers-color-scheme:dark){:root{--bg:#14161a;--fg:#e8eaef;--muted:#9aa2b1;
--line:#2a2e37;--card:#1b1e24;--ok:#4ec98a;--warn:#e0a34a;--bad:#f2857c;
--accent:#8fb0ff;--accent-soft:#212836}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1120px;margin:0 auto;padding:32px 20px 72px}
header{border-bottom:2px solid var(--line);padding-bottom:20px;margin-bottom:28px}
h1{margin:0 0 6px;font-size:27px;letter-spacing:-.02em}
h2{margin:36px 0 12px;font-size:19px;letter-spacing:-.01em;
border-left:3px solid var(--accent);padding-left:10px}
h3{margin:22px 0 8px;font-size:15px;color:var(--muted);text-transform:uppercase;
letter-spacing:.06em}
.sub{color:var(--muted);font-size:14px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(168px,1fr));gap:12px;margin:20px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.card .k{font-size:12px;color:var(--muted);text-transform:uppercase;letter-spacing:.06em}
.card .v{font-size:21px;font-weight:600;margin-top:4px;letter-spacing:-.01em}
.tablewrap{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:var(--card)}
table{border-collapse:collapse;width:100%;font-size:14px}
th,td{text-align:left;padding:9px 13px;border-bottom:1px solid var(--line);
white-space:nowrap;vertical-align:top}
th{background:var(--accent-soft);font-weight:600;font-size:12.5px;
text-transform:uppercase;letter-spacing:.05em;color:var(--muted)}
tr:last-child td{border-bottom:none}
td.num{text-align:right;font-variant-numeric:tabular-nums}
td.wrap{white-space:normal;min-width:280px}
code,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:13px}
.pill{display:inline-block;padding:1px 9px;border-radius:999px;font-size:12px;font-weight:600}
.pill.ok{background:rgba(13,122,74,.13);color:var(--ok)}
.pill.warn{background:rgba(163,90,0,.14);color:var(--warn)}
.pill.bad{background:rgba(179,38,30,.13);color:var(--bad)}
.note{background:var(--card);border:1px solid var(--line);border-left:3px solid var(--warn);
border-radius:8px;padding:12px 16px;margin:12px 0}
.note.bad{border-left-color:var(--bad)}
.note.ok{border-left-color:var(--ok)}
.cmd{background:var(--card);border:1px solid var(--line);border-radius:8px;
padding:12px 16px;overflow-x:auto;margin:10px 0}
.bar{height:7px;border-radius:4px;background:var(--accent-soft);overflow:hidden;
min-width:90px;margin-top:5px}
.bar>span{display:block;height:100%;background:var(--accent)}
footer{margin-top:56px;padding-top:18px;border-top:1px solid var(--line);
color:var(--muted);font-size:13px}
"""


def _esc(x: Any) -> str:
    return html.escape(str(x), quote=True)


def _table(headers: list[str], rows: list[list[str]], aligns: str = "") -> str:
    if not rows:
        return '<p class="sub">Nothing recorded.</p>'
    head = "".join(f"<th>{_esc(h)}</th>" for h in headers)
    body = []
    for r in rows:
        tds = []
        for i, cell in enumerate(r):
            cls = ""
            if i < len(aligns):
                cls = {"n": ' class="num"', "w": ' class="wrap"'}.get(aligns[i], "")
            tds.append(f"<td{cls}>{cell}</td>")
        body.append("<tr>" + "".join(tds) + "</tr>")
    return (f'<div class="tablewrap"><table><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table></div>')


def _pill(text: str, kind: str) -> str:
    return f'<span class="pill {kind}">{_esc(text)}</span>'


def render_html(ctx: dict[str, Any]) -> str:
    cfg: CerberusConfig = ctx["cfg"]
    acct: RunAccounting = ctx["accounting"]
    parts: list[str] = []

    # ---- header ----
    parts.append(f"""<div class="wrap"><header>
<h1>Cerberus run report — {_esc(ctx['sample_id'])}</h1>
<div class="sub">Cerberus {_esc(ctx['version'])} &middot; {_esc(ctx['generated'])} &middot;
modes: {_esc(', '.join(ctx['modes']))}
{' &middot; <strong>DRY RUN</strong>' if ctx['dry_run'] else ''}</div>
</header>""")

    if ctx["command_line"]:
        parts.append(f'<div class="cmd mono">{_esc(ctx["command_line"])}</div>')

    # ---- headline numbers ----
    total_in = acct.input_r1_reads + acct.input_r2_reads + acct.input_long_reads
    total_out = sum(r["reads"] for r in ctx["outputs"] if r["exists"])
    n_ok = sum(1 for r in ctx["outputs"] if r["exists"] and r["gzip_ok"] is not False)
    cards = [
        ("Input records", _fmt_int(total_in)),
        ("Output records", _fmt_int(total_out)),
        ("Retained", _pct(total_out, total_in) if total_in else "—"),
        ("Outputs written", f"{n_ok}/{len(ctx['outputs'])}"),
        ("Runtime", f"{ctx['elapsed_sec']:.1f}s"),
        ("Warnings", str(len(ctx["warnings"]))),
    ]
    parts.append('<div class="cards">' + "".join(
        f'<div class="card"><div class="k">{_esc(k)}</div><div class="v">{_esc(v)}</div></div>'
        for k, v in cards
    ) + "</div>")

    # ---- warnings ----
    if ctx["warnings"]:
        parts.append("<h2>Warnings</h2>")
        for w in ctx["warnings"]:
            parts.append(f'<div class="note bad">{_esc(w)}</div>')
    else:
        parts.append('<div class="note ok">No warnings were raised during this run.</div>')

    if ctx["skipped_optional"]:
        parts.append(
            '<div class="note">Optional reference(s) unavailable, so the run used fewer '
            f'mechanisms than usual: <code>{_esc(", ".join(ctx["skipped_optional"]))}</code>.'
            "</div>"
        )

    # ---- outputs ----
    parts.append("<h2>Outputs</h2>")
    rows = []
    for r in ctx["outputs"]:
        if not r["exists"]:
            status = _pill("missing", "bad")
        elif r["gzip_ok"] is False:
            status = _pill("corrupt gzip", "bad")
        elif r.get("paired_ok") is False:
            status = _pill("pairs desynced", "bad")
        elif r["reads"] == 0:
            status = _pill("empty", "warn")
        else:
            status = _pill("ok", "ok")
        rows.append([
            _esc(r["group"]), _esc(r["role"]),
            f'<code>{_esc(r["name"])}</code>',
            _fmt_int(r["reads"]), _fmt_bytes(r["size"]), status,
        ])
    parts.append(_table(
        ["Head", "Role", "File", "Records", "Size", "Status"], rows, aligns="   nn ",
    ))
    parts.append(
        '<p class="sub">"Records" counts FASTQ entries. A synchronised pair therefore '
        "shows the same number in R1 and R2, and a merged file shows the sum.</p>"
    )

    # ---- per-stage accounting ----
    parts.append("<h2>Stage-by-stage read accounting</h2>")
    modes = []
    for s in acct.stages:
        if s.mode not in modes:
            modes.append(s.mode)
    if not modes:
        parts.append('<p class="sub">No per-stage counts (dry run).</p>')
    for mode in modes:
        stages = acct.stages_for(mode)
        rows = []
        prev = None
        baseline = 0
        for s in stages:
            # A head interleaves several streams (paired, orphans, merged...).
            # A count that goes *up* means a new stream started, so the delta
            # against the previous row would be meaningless — restart there.
            new_stream = prev is None or s.reads > prev
            if new_stream:
                baseline = s.reads
                removed = share = "—"
            else:
                removed = _fmt_int(prev - s.reads)
                share = _pct(prev - s.reads, prev) if prev else "—"
            width = (100.0 * s.reads / baseline) if baseline else 0
            rows.append([
                f'<code>{_esc(s.stage)}</code>',
                _fmt_int(s.reads), removed, share,
                f'<div class="bar"><span style="width:{width:.1f}%"></span></div>',
            ])
            prev = s.reads
        parts.append(f"<h3>{_esc(mode)}</h3>")
        parts.append(_table(
            ["Stage", "Records surviving", "Removed here", "% removed", "Retained"],
            rows, aligns=" nnn ",
        ))

    # ---- GDPR mechanisms ----
    if ctx["gdpr_mechanisms"]:
        parts.append("<h2>GDPR scrub — contribution of each mechanism</h2>")
        rows = []
        for mode, streams in ctx["gdpr_mechanisms"].items():
            for stream, mech in streams.items():
                for name, pct in mech.items():
                    flag = _pill("no effect", "warn") if pct == 0 else _pill("active", "ok")
                    rows.append([_esc(mode), _esc(stream),
                                 f"<code>{_esc(name)}</code>", f"{pct:.4f}%", flag])
        parts.append(_table(["Head", "Stream", "Mechanism", "% removed at this step", "Status"],
                            rows, aligns="   n "))
        parts.append(
            '<div class="note">A mechanism marked <em>no effect</em> removed nothing. That '
            "has two very different readings: either the upstream heads had already taken the "
            "host out (expected when GDPR post-processes a cleaned output), or the mechanism "
            "could not act at all because of a preset or database mismatch. Compare it against "
            "the stage table above before trusting the result.<br><br>"
            "Cerberus reports measured removal rather than asserting a zero-host guarantee: "
            "all three mechanisms derive from the same reference assemblies, so host sequence "
            "absent from those assemblies — population-specific insertions, V(D)J junctions, "
            "novel structural variants — is invisible to every one of them.</div>"
        )

    # ---- resolved parameters ----
    parts.append("<h2>Resolved parameters</h2>")
    parts.append("<h3>Autotuned</h3>")
    parts.append(_table(
        ["Parameter", "Value"],
        [[f"<code>{_esc(k)}</code>", _esc(v)] for k, v in ctx["tuned"].items()],
    ))
    if ctx["prescan"]:
        p = ctx["prescan"]
        parts.append("<h3>Input prescan</h3>")
        parts.append(_table(
            ["Metric", "Value"],
            [
                ["Reads sampled", _fmt_int(p.get("reads_sampled", 0))],
                ["Mean read length", f'{p.get("mean_length", 0):.1f} bp'],
                ["Q20 rate", f'{100 * p.get("q20_rate", 0):.2f}%'],
                ["Q30 rate", f'{100 * p.get("q30_rate", 0):.2f}%'],
                ["Source", f'<code>{_esc(p.get("source", ""))}</code>'],
            ],
        ))

    parts.append("<h3>Run configuration</h3>")
    cfg_rows = [
        ["Sample ID", _esc(cfg.sample_id)],
        ["Output directory", f"<code>{_esc(cfg.out_dir)}</code>"],
        ["Reference directory", f"<code>{_esc(cfg.ref_dir)}</code>"],
        ["Threads", str(cfg.threads)],
        ["Memory budget", f"{cfg.memory_gb} G"],
        ["Platform", _esc(cfg.platform.value)],
        ["Fast mode", "yes" if cfg.fast else "no"],
        ["Double pass", "yes" if cfg.double_pass else "no"],
        ["Keep intermediates", "yes" if cfg.keep_intermediates else "no"],
    ]
    if cfg.gdpr:
        cfg_rows.append(["Kraken2 confidence (GDPR)", str(cfg.gdpr_confidence)])
        cfg_rows.append(["Human k-mer scrub", "yes" if cfg.gdpr_kmer_scrub else "no"])
    parts.append(_table(["Setting", "Value"], cfg_rows))

    # ---- environment ----
    parts.append("<h2>Environment</h2>")
    parts.append(_table(
        ["Property", "Value"],
        [[_esc(k.replace("_", " ").title()), _esc(v)] for k, v in ctx["environment"].items()],
    ))
    parts.append("<h3>Tool versions</h3>")
    parts.append(_table(
        ["Tool", "Version"],
        [[f"<code>{_esc(k)}</code>", _esc(v)] for k, v in sorted(ctx["tools"].items())],
    ))

    parts.append(
        f"<footer>Generated by Cerberus {_esc(ctx['version'])}. "
        "Machine-readable equivalents of this page are in "
        "<code>accounting.json</code> and <code>run_record.json</code> "
        "in the same directory.</footer></div>"
    )

    body = "\n".join(parts)
    return (
        "<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
        "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
        f"<title>Cerberus run report — {_esc(ctx['sample_id'])}</title>"
        f"<style>{_CSS}</style></head><body>{body}</body></html>\n"
    )


def write_run_report(
    cfg: CerberusConfig,
    *,
    accounting: RunAccounting,
    results: list[PipelineResult],
    gdpr_outputs: dict[str, dict[str, Path | None]],
    refs: Any,
    elapsed_sec: float,
) -> Path:
    """Render the HTML report and the machine-readable run record."""
    ctx = build_context(
        cfg, accounting=accounting, results=results,
        gdpr_outputs=gdpr_outputs, refs=refs, elapsed_sec=elapsed_sec,
    )
    cfg.reports_dir.mkdir(parents=True, exist_ok=True)

    html_path = cfg.reports_dir / "cerberus_report.html"
    html_path.write_text(render_html(ctx), encoding="utf-8")

    record = {
        "cerberus_version": ctx["version"],
        "generated_utc": ctx["generated"],
        "sample_id": ctx["sample_id"],
        "command_line": ctx["command_line"],
        "modes": ctx["modes"],
        "elapsed_sec": round(elapsed_sec, 3),
        "dry_run": ctx["dry_run"],
        "resolved_parameters": ctx["tuned"],
        "prescan": ctx["prescan"],
        "environment": ctx["environment"],
        "tool_versions": ctx["tools"],
        "outputs": [
            {k: v for k, v in row.items() if k != "name"} for row in ctx["outputs"]
        ],
        "gdpr_mechanisms": ctx["gdpr_mechanisms"],
        "warnings": ctx["warnings"],
        "skipped_optional_assets": ctx["skipped_optional"],
    }
    (cfg.reports_dir / "run_record.json").write_text(json.dumps(record, indent=2))

    log.info("Wrote run report: %s", html_path)
    return html_path
