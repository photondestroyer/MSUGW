"""The blinding rule (pilot Section 3, rule 1) in one place.

"Do not view, tabulate or plot cohort ET or surface temperature by aquifer thickness for dry conditions.
 Allowed before the freeze: support masks and observation counts; ET of rainfed cells; cohort ET from non-dry
 conditions, used only to estimate residual variance. Replication region: metadata only."

How the code keeps it:
  STEP B-1  dry flags (task `dry`) come from weather (Task 1.1) and S (Task 1.3) only, via definitions.py.
  STEP B-2  gate(et, ...) is the ONLY way an ET table reaches disk: it keeps rainfed rows and cohort rows of
            non-dry months, drops everything else, refuses tables that carry thickness columns, and refuses
            cells outside Kansas (the replication region). Dry cohort ET is computed in memory and dropped
            before any write, so it never exists on disk.
  STEP B-3  audit(root) re-reads every table written under the output root and fails when any ET table holds
            a cohort row in dry conditions, a thickness column, or a row without a dry flag.
No surface temperature is extracted anywhere in the pilot.
"""
import json
import os

import numpy as np

from . import config as C
from . import definitions as DEF
from .core.errors import PipelineStop

GATE_VERSION = "blinding-gate-1"
ET_PREFIX = ("et_", "ET_")
THICKNESS_COLUMNS = {"A_m", "A_se_m", "A_third", "A_poor", "wte_m", "base_m", "thickness", "satthick"}
# ... and any column whose NAME says thickness, water table or aquifer base, whatever table it came from
# (the cell table holds kgs_satthick_2022_2024, kgs_predev_satthick_updated, kgs_hpa_bedrock, kgs_wte_2022_2024).
THICKNESS_PATTERN = r"thick|wte|bedrock|^A_"


def thickness_columns(cols):
    """The columns of a table that carry aquifer thickness or what it is made of."""
    import re
    pat = re.compile(THICKNESS_PATTERN)
    return sorted(c for c in map(str, cols) if c in THICKNESS_COLUMNS or pat.search(c))
ROLE_ALLOWED = {"rainfed": "all months", "cohort": "non-dry months only"}


class BlindingError(PipelineStop):
    pass


# ------------------------------------------------------------------------------------------- STEP B-1 dry flags
def dry_table(weather_cy, s_cm, d20=None):
    """County-month dry flags for the season months of every panel year.
    weather_cy: county, year, deficit_season_mm (Task 1.1); s_cm: county, year, month, S (Task 1.3)."""
    import pandas as pd
    w = weather_cy[["county", "year", "deficit_season_mm"]].copy()
    w["dry_season"], season_thr = DEF.season_dry(w["deficit_season_mm"].to_numpy(), d20)
    grid = pd.MultiIndex.from_product([sorted(w.county.unique()), C.PANEL_YEARS, C.SEASON_MONTHS],
                                      names=["county", "year", "month"]).to_frame(index=False)
    s = grid.merge(s_cm[["county", "year", "month", "S"]], on=["county", "year", "month"], how="left")
    s["dry_month"], s_thr = DEF.month_dry(s["S"].to_numpy(), d20)
    out = s.merge(w[["county", "year", "dry_season"]], on=["county", "year"], how="left")
    out["dry_season"] = out["dry_season"].astype("boolean").fillna(True).astype(bool)   # no weather -> dry
    out["dry"] = out["dry_season"] | out["dry_month"]
    return out, dict(season_threshold_mm=season_thr, s_threshold=s_thr)


# ------------------------------------------------------------------------------------------- STEP B-2 the gate
def gate(et, cells, dry):
    """et: cell, year, month, et_* columns (any). cells: cell_id, county, role. dry: county, year, month, dry.
    Returns the rows that may be written, with role and dry attached, and a count of what was dropped."""
    bad = thickness_columns(et.columns)
    if bad:
        raise BlindingError(f"an ET table must not carry thickness columns: {bad}")
    c = cells[["cell_id", "county", "role"]]
    if not c["county"].astype(str).str.startswith(C.STATE_FIPS).all():
        raise BlindingError("cells outside Kansas reached the ET gate (the replication region is metadata only)")
    x = et.merge(c, left_on="cell", right_on="cell_id", how="left").drop(columns="cell_id")
    if x["role"].isna().any():
        raise BlindingError("ET rows for cells that are not in the cell table")
    x = x.merge(dry[["county", "year", "month", "dry"]], on=["county", "year", "month"], how="left")
    x["dry"] = x["dry"].astype("boolean").fillna(True).astype(bool)     # unknown -> treated as dry
    keep = (x["role"] == "rainfed") | ((x["role"] == "cohort") & ~x["dry"])
    dropped = dict(rows_in=int(len(x)), rows_kept=int(keep.sum()),
                   cohort_dry_dropped=int(((x["role"] == "cohort") & x["dry"]).sum()),
                   other_dropped=int((~x["role"].isin(["rainfed", "cohort"])).sum()))
    out = x[keep].copy()
    out["blinding_gate"] = GATE_VERSION
    return out.reset_index(drop=True), dropped


def write_et(ctx, et, cells, dry, name, units=None, notes=""):
    """The one function that writes ET values. Gate first, write second; the gate result goes in the sidecar."""
    from .common import report
    out, dropped = gate(et, cells, dry)
    report.table(ctx, out, name, units=units or {"et_*": "mm per month"}, parquet=True, csv=False,
                 mask="blinded: rainfed all months; cohort non-dry months only",
                 notes=f"{notes} gate={json.dumps(dropped)}")
    return out, dropped


# ------------------------------------------------------------------------------------------- STEP B-3 audit
def _et_columns(cols):
    """Columns that hold ET values: the ensemble, the six models, and their per-cell sums."""
    import re
    names = "|".join(["ensemble_mad", "ensemble_sam"] + C.ET_MODELS)
    pat = re.compile(rf"^et_({names})(_sum)?$")
    return [c for c in cols if pat.match(str(c))]


def audit(root=None):
    """Scan every CSV/parquet under the output root. Returns (passes, problems, files_checked)."""
    import pandas as pd
    root = root or C.out_root()
    problems, checked = [], 0
    for dp, _, files in os.walk(root):
        for f in files:
            if not f.endswith((".csv", ".parquet")):
                continue
            p = os.path.join(dp, f)
            try:                                        # column names only: no table is loaded to be skipped
                if f.endswith(".parquet"):
                    import pyarrow.parquet as pq
                    cols = pq.read_schema(p).names
                else:
                    cols = pd.read_csv(p, nrows=0).columns
            except Exception as e:                      # unreadable file: report, never ignore
                problems.append(f"{p}: unreadable ({e})")
                continue
            ets = _et_columns(cols)
            if not ets:
                continue
            checked += 1
            df = pd.read_parquet(p) if f.endswith(".parquet") else pd.read_csv(p, low_memory=False)
            if thickness_columns(df.columns):
                problems.append(f"{p}: ET table with thickness columns {thickness_columns(df.columns)}")
            if "role" not in df.columns or "dry" not in df.columns:
                problems.append(f"{p}: ET table without role/dry columns (did not pass the gate)")
                continue
            leak = (df["role"] == "cohort") & df["dry"].astype(bool)
            if leak.any():
                problems.append(f"{p}: {int(leak.sum())} cohort rows in dry conditions")
            if (~df["role"].isin(["cohort", "rainfed"])).any():
                problems.append(f"{p}: ET rows for cells that are neither cohort nor rainfed")
    return not problems, problems, checked


def run(ctx):
    """Task `dry`: STEP B-1 for every Kansas county-month of the panel."""
    from .common import report
    w = report.read(ctx.dep("t1_1"), "weather_county_year")
    s = report.read(ctx.dep("t1_3"), "county_month_S")
    out, thr = dry_table(w, s)
    report.table(ctx, out, "dry_flags", parquet=True, units={"S": "z-score"},
                 notes=f"thresholds={json.dumps(thr)}")
    ctx.complete(dict(thresholds=thr, county_months=int(len(out)), dry_share=float(out["dry"].mean()),
                      dry_season_share=float(out.drop_duplicates(["county", "year"])["dry_season"].mean())))
