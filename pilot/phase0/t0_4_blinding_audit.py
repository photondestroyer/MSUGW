"""Task 0.4  Blinding audit.

"List every earlier analysis in the repository that examined ET by aquifer condition, with the periods it
covered. The repository shows SSEBop ET used up to 21 May 2022."
Physical idea: the pilot's test is only clean if nobody has already looked at drought-year ET by aquifer
condition; where someone has, the periods they saw must be on record before the freeze.
Decisions: none.   Deliverables: blinding_audit.csv (one row per file), blinding_audit.md.
"""
import os
import re
import subprocess

SKIP_DIRS = {".git", "__pycache__", "node_modules", "worktrees", ".kilo", ".claude", "pilot", "pilot test",
             "merged_datasets", "zarr", "plan_v3_build", "data audit"}
CODE_EXT = (".py", ".ipynb", ".r", ".rmd", ".md", ".txt", ".json", ".csv", ".log")
ET_RE = re.compile(r"(?i)\b(evapotranspiration|ssebop|openet|et_ensemble|ecostress|echo-?et|residual[ _-]?et|"
                   r"buffering residual|\bET\b|Bi_|et_anom|eta\b)")
AQ_RE = re.compile(r"(?i)(water[ _-]?table|\bwtd\b|depth[ _-]?to[ _-]?water|saturated[ _-]?thickness|aquifer[ _-]?"
                   r"(thickness|condition|depletion)|hres[ _-]?wtd|dwl_|thickness|well[ _-]?capacity|dtw\b)")
DATE_RE = re.compile(r"\b(19[89]\d|20[0-2]\d)-(0[1-9]|1[0-2])-(0[1-9]|[12]\d|3[01])\b")
YEAR_RANGE_RE = re.compile(r"\b(20[0-2]\d)\s*(?:-|–|to)\s*(20[0-2]\d)\b")
KNOWN = {   # outputs whose period is documented elsewhere in the repository
    "gw_buffering_xgb": "SSEBop dekadal ET 2015-04 to 2022-05-21 (RESIDUAL_ET_MODELS.md; catalogue item 3); "
                        "residual ET by water-table depth (buffering_by_wtd.png, effects.csv)",
    "residual_et_store": "feature store of the corrected residual-ET models; dekads 2015-04-01 to 2022-05-21",
    "XGB model": "XGBoost residual-ET script; water-table depth among the summaries",
    "temporal buffering": "temporal buffering index (TBI) scripts over USGS water-level change",
    "tbi_maps_v3_usgs": "TBI maps v3 (USGS water-level data)",
}


def _files(root):
    for dp, dirs, files in os.walk(root):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS and not d.startswith(".")]
        for f in files:
            if f.lower().endswith(CODE_EXT) or f.lower().endswith(".png"):
                yield os.path.join(dp, f)


def _git(path):
    try:
        out = subprocess.run(["git", "-C", os.path.dirname(path), "log", "--format=%ad", "--date=short", "--", path],
                             capture_output=True, text=True, timeout=60).stdout.split()
        return (out[-1], out[0]) if out else ("untracked", "untracked")
    except (OSError, subprocess.SubprocessError):
        return ("", "")


def run(ctx):
    import pandas as pd
    from .. import config as C
    from ..common import report
    rows = []
    # ══════ STEP 0.4-1  every code, notebook, text and figure file outside the pilot and the data folders ══════
    for p in _files(C.REPO_ROOT):
        rel = os.path.relpath(p, C.REPO_ROOT).replace("\\", "/")
        top = rel.split("/")[0]
        # ══════ STEP 0.4-2  figures: file names that pair ET with an aquifer variable ══════
        if p.lower().endswith(".png"):
            if ET_RE.search(rel) or re.search(r"(?i)buffer|Bi_", rel):
                if AQ_RE.search(rel):
                    rows.append(dict(file=rel, kind="figure", et_hits="name", aquifer_hits="name", lines="",
                                     dates_min="", dates_max="", year_ranges="", known=KNOWN.get(top, "")))
            continue
        if os.path.getsize(p) > 30e6:
            continue
        # ══════ STEP 0.4-3  text: lines naming ET and lines naming aquifer condition, with dates and year ranges ══════
        try:
            txt = open(p, encoding="utf-8", errors="ignore").read()
        except OSError:
            continue
        et_l = [i + 1 for i, l in enumerate(txt.splitlines()) if ET_RE.search(l)]
        aq_l = [i + 1 for i, l in enumerate(txt.splitlines()) if AQ_RE.search(l)]
        if not et_l or not aq_l:
            continue
        dates = sorted(m.group(0) for m in DATE_RE.finditer(txt))
        ranges = sorted({f"{a}-{b}" for a, b in YEAR_RANGE_RE.findall(txt) if a < b})
        rows.append(dict(file=rel, kind=os.path.splitext(p)[1].lstrip("."), et_hits=len(et_l), aquifer_hits=len(aq_l),
                         lines=f"ET {et_l[:5]} / aquifer {aq_l[:5]}", dates_min=dates[0] if dates else "",
                         dates_max=dates[-1] if dates else "", year_ranges=", ".join(ranges[:12]),
                         known=KNOWN.get(top, "")))
    df = pd.DataFrame(rows)
    if df.empty:
        df = pd.DataFrame(columns=["file", "kind", "et_hits", "aquifer_hits", "lines", "dates_min", "dates_max",
                                   "year_ranges", "known"])
    # ══════ STEP 0.4-4  git history of each file (when it was written and last changed) ══════
    g = [_git(os.path.join(C.REPO_ROOT, f)) for f in df["file"]]
    df["git_first"], df["git_last"] = [a for a, _ in g], [b for _, b in g]
    df["overlaps_panel"] = df["year_ranges"].astype(str).str.contains(r"201[6-9]|202[0-4]") | \
        df["dates_max"].astype(str).str[:4].isin([str(y) for y in C.PANEL_YEARS])
    df = df.sort_values(["known", "file"], ascending=[False, True])
    ctx.log.info("STEP 0.4-3  %d files pair ET with an aquifer variable", len(df))

    # ══════ STEP 0.4-5  write the audit ══════
    report.table(ctx, df, "blinding_audit")
    lines = ["# Blinding audit (Task 0.4)", "",
             "Earlier work in this repository that put ET next to an aquifer variable, found by scanning code, "
             "notebooks, documents and figure names. A row is a lead to read, not proof that an analysis was run.", "",
             "## Known analyses", ""]
    for k, v in KNOWN.items():
        lines.append(f"- **{k}/**: {v}")
    lines += ["", f"## Files ({len(df)})", "", "| file | kind | dates in file | year ranges | git first / last | "
              "overlaps 2016-2024 |", "|---|---|---|---|---|---|"]
    for r in df.itertuples():
        lines.append(f"| `{r.file}` | {r.kind} | {r.dates_min} .. {r.dates_max} | {r.year_ranges} | "
                     f"{r.git_first} / {r.git_last} | {'yes' if r.overlaps_panel else ''} |")
    lines += ["", "Rule for the freeze: results from these analyses for 2016-2024 dry seasons are treated as "
              "already seen; the decision log records any choice that could have been influenced by them."]
    report.text(ctx, "blinding_audit.md", "\n".join(lines))
    ctx.complete(dict(files=int(len(df)), overlapping_panel=int(df["overlaps_panel"].sum())))
