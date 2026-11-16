"""Task 1.7  Complete support table (Stage 3), and the preliminary county-year inventory.

"Join Tasks 1.1 to 1.5 for all eligible cells. Keep the reason for every exclusion. Build four sample masks: the
soil-moisture contrast; the peak-demand contrast; the sample common to both; and the monthly timing analysis.
Recount the dry episodes and the range of deficit and S inside each mask."
"First deliverable. A preliminary county-year inventory ... from Tasks 1.1 to 1.5. Mark items that are not yet
available as pending."
Physical idea: each contrast needs the same things at the same place and month - a cohort cell, its ET seen by
at least four models, the drought signal S of its county, and wells close enough to know its aquifer. The masks
are where all of these hold; what they keep (episodes, ranges of deficit and S) is what the tests can learn from.
Decisions: D21 (ET support), D25 (well distance for thickness support).
Deliverables: support_county_year.csv, cell_exclusions.csv, masks (parquet), mask_summary.csv (Paper 1 Table 2),
support_map.png, support_summary.md; inventory: county_year_inventory.csv.
"""
import numpy as np


def _cell_year(ctx, cells, years):
    import pandas as pd
    from ..common import report
    w = report.read(ctx.dep("t1_1"), "weather_county_year")
    w["county"] = w["county"].astype(str)
    cy = cells[["cell_id", "county", "role", "crop_class", "area_ha"]].merge(pd.DataFrame(dict(year=years)),
                                                                              how="cross")
    return cy.merge(w, on=["county", "year"], how="left")


def inventory(ctx):
    """The first deliverable: county-year inventory from Tasks 1.1 and 1.5 (+ 0.2); S and ET pending."""
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from ..common import report
    from ..phase0 import t0_2_cells as T02
    cells = T02.load(ctx.dep("t0_2"))
    w = report.read(ctx.dep("t1_1"), "weather_county_year")
    w["county"] = w["county"].astype(str)
    dist = report.read(ctx.dep("t1_5"), "cell_well_distance")
    from .. import definitions as DEF
    near = (dist.assign(ok=DEF.well_support(dist.well_km)).groupby(["county", "year"]).ok.mean()
            .rename("cohort_share_well_8km").reset_index())
    near["county"] = near["county"].astype(str)
    n = cells.groupby(["county", "role"]).size().unstack(fill_value=0).reindex(columns=["cohort", "rainfed"],
                                                                               fill_value=0).reset_index()
    inv = w.merge(n, on="county", how="left").merge(near, on=["county", "year"], how="left")
    for col in ("S_peak", "S_valid_peak", "n_reference_footprints", "cohort_share_et_supported_peak"):
        inv[col] = "pending"
    report.table(ctx, inv, "county_year_inventory",
                 units=dict(deficit_season_mm="mm", cohort="cells", rainfed="cells", cohort_share_well_8km="fraction"),
                 notes="preliminary: S and ET items are pending until Tasks 1.3 and 1.4 run")
    ctx.complete(dict(county_years=int(len(inv)), counties=int(inv.county.nunique())))


def run(ctx):
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from .. import definitions as DEF
    from ..common import report, cells as CL
    from ..phase0 import t0_2_cells as T02, counties as CTY

    # ══════ STEP 1.7-1  cells in scope; cell-months of the season with weather, S, dry, ET support, wells ══════
    cells = T02.load(ctx.dep("t0_2"))
    sample = report.read(ctx.dep("sample"), "sample_cells").cell_id
    cells = cells[CL.scope(cells, sample)].copy()
    work = cells[cells.role.isin(["cohort", "rainfed"])]
    cy = _cell_year(ctx, work, C.PANEL_YEARS)
    cm = cy.merge(pd.DataFrame(dict(month=C.SEASON_MONTHS)), how="cross")
    cm["is_peak"] = cm.month == cm.peak_month
    S = report.read(ctx.dep("t1_3"), "county_month_S")
    S["county"] = S["county"].astype(str)
    cm = cm.merge(S[["county", "year", "month", "S", "S_valid", "n_footprints"]], on=["county", "year", "month"],
                  how="left")
    cm["S_valid"] = cm["S_valid"].astype("boolean").fillna(False).astype(bool)
    dry = report.read(ctx.dep("dry"), "dry_flags")
    dry["county"] = dry["county"].astype(str)
    cm = cm.merge(dry[["county", "year", "month", "dry"]], on=["county", "year", "month"], how="left")
    cov = report.read(ctx.dep("t1_4"), "et_coverage")[["cell", "year", "month", "n_models", "supported"]]
    cm = cm.merge(cov.rename(columns={"cell": "cell_id", "supported": "et_supported"}),
                  on=["cell_id", "year", "month"], how="left")
    cm["et_supported"] = cm["et_supported"].astype("boolean").fillna(False).astype(bool)
    wk = report.read(ctx.dep("t1_5"), "cell_well_distance")[["cell_id", "year", "well_km"]]
    cm = cm.merge(wk, on=["cell_id", "year"], how="left")
    cm["a_support"] = DEF.well_support(cm.well_km)
    ctx.log.info("STEP 1.7-1  %d cell-months", len(cm))

    # ══════ STEP 1.7-2  the four sample masks (definitions.sample_masks) ══════
    mk = DEF.sample_masks(cm.role, cm.et_supported, cm.S_valid, cm.a_support.fillna(False), cm.is_peak)
    for k, v in mk.items():
        cm[k] = v

    # ══════ STEP 1.7-3  exclusions: every cell of the study, with the first reason it drops out ══════
    allc = T02.load(ctx.dep("t0_2"), columns=["cell_id", "county", "role", "reason", "area_ha"])
    allc = allc[allc.cell_id.isin(cells.cell_id)]
    used = cm[cm.mask_peak | cm.mask_sm].groupby("cell_id").size()
    any_et = cm.groupby("cell_id").et_supported.any()
    any_a = cm.groupby("cell_id").a_support.any()
    ex = allc.copy()
    ex["final_reason"] = ex["reason"]
    w = ex.role.isin(["cohort", "rainfed"])
    ex.loc[w & ~ex.cell_id.map(any_et).astype("boolean").fillna(False).astype(bool), "final_reason"] = "no month with >= 4 ET models"
    ex.loc[w & (ex.role == "cohort") & ~ex.cell_id.map(any_a).astype("boolean").fillna(False).astype(bool) &
           ex.final_reason.eq(""), "final_reason"] = "no measured well within 8 km in any year"
    ex.loc[w & ex.final_reason.eq("") & ~ex.cell_id.isin(used.index), "final_reason"] = \
        "in no peak-month mask (S or ET support never coincide)"
    ex["retained"] = ex.final_reason.eq("")

    # ══════ STEP 1.7-4  inside each mask: counties, cells, county-years, dry episodes, ranges of deficit and S ══════
    rows = []
    for k in DEF.MASKS:
        sub = cm[cm[k]]
        co = sub[sub.role == "cohort"]
        dry_cy = co[co.dry_season.astype(bool)].drop_duplicates(["county", "year"])
        rows.append(dict(mask=k, counties=int(co.county.nunique()), cohort_cells=int(co.cell_id.nunique()),
                         rainfed_cells=int(sub[sub.role == "rainfed"].cell_id.nunique()),
                         cohort_cell_months=int(len(co)), county_years=int(co.drop_duplicates(["county", "year"]).shape[0]),
                         dry_county_years=int(len(dry_cy)), dry_episodes=int(dry_cy.episode_id.nunique()),
                         dry_episodes_multiyear=int(dry_cy.episode_id_multiyear.nunique()),
                         deficit_peak_min=float(co.deficit_peak_mm.min()) if len(co) else np.nan,
                         deficit_peak_max=float(co.deficit_peak_mm.max()) if len(co) else np.nan,
                         S_min=float(co.S.min()) if len(co) else np.nan, S_max=float(co.S.max()) if len(co) else np.nan,
                         area_ha=float(cells.set_index("cell_id").loc[co.cell_id.unique(), "area_ha"].sum())
                         if len(co) else 0.0))
    ms = pd.DataFrame(rows)

    # ══════ STEP 1.7-5  county-year support table ══════
    pk = cm[cm.is_peak]
    cyt = (pk.groupby(["county", "year"])
           .agg(cohort_cells=("role", lambda r: int((r == "cohort").sum())),
                rainfed_cells=("role", lambda r: int((r == "rainfed").sum())),
                peak_month=("peak_month", "first"), deficit_season_mm=("deficit_season_mm", "first"),
                deficit_peak_mm=("deficit_peak_mm", "first"), dry_season=("dry_season", "first"),
                episode_id=("episode_id", "first"), S_peak=("S", "first"), S_valid_peak=("S_valid", "first"),
                n_footprints_peak=("n_footprints", "first"),
                cohort_share_et_supported_peak=("et_supported", lambda s: float(s[pk.loc[s.index, "role"] == "cohort"].mean())
                                                if (pk.loc[s.index, "role"] == "cohort").any() else np.nan),
                cohort_share_well_8km=("a_support", lambda s: float(s[pk.loc[s.index, "role"] == "cohort"].mean())
                                       if (pk.loc[s.index, "role"] == "cohort").any() else np.nan))
           .reset_index())

    # ══════ STEP 1.7-6  write: tables, map, one-page summary ══════
    keep = ["cell_id", "county", "role", "crop_class", "year", "month", "is_peak", "peak_month", "dry", "dry_season",
            "episode_id", "deficit_peak_mm", "S", "S_valid", "et_supported", "n_models", "well_km", "a_support"] + DEF.MASKS
    report.table(ctx, cm[keep], "masks", parquet=True, csv=False, units=dict(S="z", deficit_peak_mm="mm", well_km="km"))
    report.table(ctx, ms, "mask_summary", units=dict(area_ha="ha", deficit_peak_min="mm", S_min="z"),
                 notes="Paper 1 Table 2 (sample) draft")
    report.table(ctx, cyt, "support_county_year", units=dict(deficit_season_mm="mm", S_peak="z"))
    report.table(ctx, ex, "cell_exclusions", units=dict(area_ha="ha"))
    area = ex.groupby(["role", "retained"]).area_ha.sum().reset_index()
    report.table(ctx, area, "area_retained", units=dict(area_ha="ha"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    geo = T02.load_geometry(ctx.dep("t0_2")).merge(ex[["cell_id", "retained", "role"]], on=["cell_id", "role"])
    fig, ax = plt.subplots(figsize=(8, 5))
    CTY.load(ctx.dep("counties")).boundary.plot(ax=ax, color="#888888", linewidth=0.4)
    geo[~geo.retained].plot(ax=ax, color=C.COLORS["excluded"], linewidth=0)
    geo[geo.retained & (geo.role == "cohort")].plot(ax=ax, color=C.COLORS["cohort"], linewidth=0)
    geo[geo.retained & (geo.role == "rainfed")].plot(ax=ax, color=C.COLORS["rainfed"], linewidth=0)
    ax.set_axis_off()
    ax.set_title("Retained cohort (purple) and rainfed (green) cells; excluded (grey)", fontsize=9)
    report.figure(ctx, fig, "support_map", mask="cells in a peak-month mask")
    L = ["# Support table: one-page summary (Task 1.7)", "", f"Mode {C.RUN_MODE}. Cells in scope: {len(cells)}.", "",
         "## Sample masks", "", report.md(ms), "", "## Exclusions", "",
         report.md(ex.final_reason.replace("", "retained").value_counts().rename_axis("reason")
                   .reset_index(name="cells")), "", "## Area", "", report.md(area)]
    report.text(ctx, "support_summary.md", "\n".join(L))
    ctx.complete(dict(cell_months=int(len(cm)), masks=ms.set_index("mask")[["cohort_cells", "county_years",
                                                                             "dry_episodes"]].to_dict("index"),
                      retained_cells=int(ex.retained.sum())))
