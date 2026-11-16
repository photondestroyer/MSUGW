"""Check 2.4  Aquifer thickness.

"Build the pre-season water-level panel from recovered levels. Interpolate the water table and validate it by
holding out whole well locations and whole years. Include the uncertainty of the base elevation. Compare the
estimates with existing products in the years they share."
Targets: each year at least 60% of cohort cells within 8 km of a measured well; at least 15 counties with 50 or
more cells in both the thin and the thick third; a thin-to-thick thickness difference of at least three times its
standard error. For a county without both thirds, report its own range of A. Exclude cells with poor estimates.
Physical idea: A is the saturated aquifer under a field - what its wells can draw on. It must be known well enough
that "thin" and "thick" are real differences, not interpolation noise.
Decisions: D22-D26.   Deliverables: thickness_cell_year (parquet; Paper 2 data product), thirds.json,
validation_location.csv, validation_year.csv, validation_stats.csv (Paper 2 Table 2), comparison.csv,
targets.csv, thickness_map.png (Paper 2 Fig 2), validation.png (Paper 2 Fig 3).
"""
import json

import numpy as np

FT = 0.3048


def _sample_nc(path, var, x, y):
    """Nearest-cell value of a 2-D variable on an EPSG:5070 x/y grid (USGS 500 m change file)."""
    import netCDF4 as nc
    with nc.Dataset(path) as d:
        xs, ys = np.asarray(d["x"][:], float), np.asarray(d["y"][:], float)
        a = np.ma.filled(d[var][:].astype(float), np.nan)
    dx, dy = xs[1] - xs[0], ys[1] - ys[0]                 # regular axes of pixel centres (dy < 0 here)
    ci = np.rint((np.asarray(x) - xs[0]) / dx).astype(int)
    ri = np.rint((np.asarray(y) - ys[0]) / dy).astype(int)
    inside = (ci >= 0) & (ci < len(xs)) & (ri >= 0) & (ri < len(ys))
    out = np.full(len(ci), np.nan)
    out[inside] = a[ri[inside], ci[inside]]
    return out


def _err_stats(df, by):
    g = df.groupby(by)
    return g.apply(lambda d: dict(n=len(d), bias_m=float(d.err.mean()), mae_m=float(d.err.abs().mean()),
                                  rmse_m=float(np.sqrt((d.err ** 2).mean())),
                                  r2=float(1 - (d.err ** 2).sum() / ((d.obs - d.obs.mean()) ** 2).sum())),
                   include_groups=False).apply(__import__("pandas").Series).reset_index()


def run(ctx):
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from .. import definitions as DEF
    from ..common import io, report, stats, cells as CL
    from ..phase0 import t0_2_cells as T02
    d23, d24 = D.value("D23"), D.value("D24")

    # ══════ STEP 2.4-1  inputs: usable pre-season wells (Task 1.5) and the cells (Task 0.2) ══════
    wy = report.read(ctx.dep("t1_5"), "well_year")
    wy = wy[wy.use.astype(bool) & np.isfinite(wy.wte_m) & np.isfinite(wy.land_m)]
    cells = T02.load(ctx.dep("t0_2"))
    sample = report.read(ctx.dep("sample"), "sample_cells").cell_id
    cells = cells[CL.scope(cells, sample) & np.isfinite(cells.land_m.to_numpy())].reset_index(drop=True)
    cxy = np.c_[cells.x5070, cells.y5070]
    ctx.log.info("STEP 2.4-1  %d usable well-years; %d cells", len(wy), len(cells))

    # ══════ STEP 2.4-2  per year: regression kriging of water-table elevation on land surface (D23) ══════
    models, pred, loc_cv = {}, [], []
    rs = np.random.default_rng(0)
    for y in C.WELL_YEARS:
        w = wy[wy.year == y]
        if len(w) < 30:
            ctx.log.warning("   %d: only %d wells, skipped", y, len(w))
            continue
        xy = np.c_[w.x5070, w.y5070]
        rk = stats.RegressionKriging(xy, w.wte_m.to_numpy(), w.land_m.to_numpy(), k=d23["neighbours"],
                                     max_km=d23["max_km"])
        models[y] = (rk, w)
        p, se = rk.predict(cxy, cells.land_m.to_numpy())
        pred.append(pd.DataFrame(dict(cell_id=cells.cell_id, year=y, wte_m=p, wte_se_m=se,
                                      well_km=stats.nearest_km(cxy, xy))))
        # hold out whole well locations: 10 folds by well
        wells = w.well_id.unique()
        fold = dict(zip(wells, rs.integers(0, 10, len(wells))))
        f = w.well_id.map(fold).to_numpy()
        for k in range(10):
            tr, te = f != k, f == k
            if te.sum() == 0 or tr.sum() < 30:
                continue
            m = stats.RegressionKriging(xy[tr], w.wte_m.to_numpy()[tr], w.land_m.to_numpy()[tr], k=d23["neighbours"],
                                        max_km=d23["max_km"])
            pp, ss = m.predict(xy[te], w.land_m.to_numpy()[te])
            loc_cv.append(pd.DataFrame(dict(year=y, well_id=w.well_id.to_numpy()[te], obs=w.wte_m.to_numpy()[te],
                                            pred=pp, se=ss)))
        ctx.log.info("STEP 2.4-2  %d: %d wells, variogram %s", y, len(w), {k: round(v, 1) for k, v in rk.vario.items()})
        ctx.check_budget()
    pred = pd.concat(pred, ignore_index=True)
    loc_cv = pd.concat(loc_cv, ignore_index=True)
    loc_cv["err"] = loc_cv.pred - loc_cv.obs

    # ══════ STEP 2.4-3  hold out whole years: predict year y from the surfaces of y-1 and y+1 ══════
    yr_cv = []
    for y in sorted(models):
        if y - 1 in models and y + 1 in models:
            _, w = models[y]
            xy = np.c_[w.x5070, w.y5070]
            pa, _ = models[y - 1][0].predict(xy, w.land_m.to_numpy())
            pb, _ = models[y + 1][0].predict(xy, w.land_m.to_numpy())
            yr_cv.append(pd.DataFrame(dict(year=y, well_id=w.well_id, obs=w.wte_m, pred=(pa + pb) / 2)))
    yr_cv = pd.concat(yr_cv, ignore_index=True)
    yr_cv["err"] = yr_cv.pred - yr_cv.obs

    # ══════ STEP 2.4-4  base (KGS) and its uncertainty from the USGS-KGS spread by GMD (D24); thickness A ══════
    diff = (cells.usgs_bottom_ft - cells.kgs_hpa_bedrock) * FT
    gm = cells.gmd.fillna(0).astype(int)
    mad = pd.Series(diff.to_numpy()).groupby(gm.to_numpy()).apply(lambda s: 1.4826 * np.nanmedian(np.abs(s - np.nanmedian(s))))
    overall = 1.4826 * np.nanmedian(np.abs(diff - np.nanmedian(diff)))
    cells["base_m"] = cells.kgs_hpa_bedrock * FT
    cells["base_se_m"] = gm.map(mad).fillna(overall).to_numpy()
    t = pred.merge(cells[["cell_id", "county", "role", "base_m", "base_se_m", "land_m"]], on="cell_id")
    t["A_m"], t["A_se_m"] = DEF.thickness(t.wte_m, t.wte_se_m, t.base_m, t.base_se_m)
    t["A_poor"] = DEF.thickness_poor(t.A_m, t.A_se_m, t.well_km)

    # ══════ STEP 2.4-5  thirds, cut once on the cohort's 2016-2024 mean A (D26) ══════
    panel = t[t.year.isin(C.PANEL_YEARS)]
    cm = panel[panel.role == "cohort"].groupby("cell_id").A_m.mean()
    cuts = DEF.thickness_thirds(cm.to_numpy())
    third = pd.Series(DEF.assign_third(cm.to_numpy(), cuts), index=cm.index)
    t["A_third"] = t.cell_id.map(third).fillna("")

    # ══════ STEP 2.4-6  comparison with existing products in shared years ══════
    comp = []
    hdir = io.path("haacker_dir")
    for y in range(2011, 2017):
        p = f"{hdir}/WTelev_Clip_{y}.tif"
        if not __import__("os").path.exists(p) or y not in models:
            continue
        hv = io.sample_raster(p, cells.x5070, cells.y5070, "EPSG:5070")
        ours = pred[pred.year == y].set_index("cell_id").reindex(cells.cell_id).wte_m.to_numpy()
        d = ours - hv
        comp.append(dict(product="Haacker et al. water table", year=y, n=int(np.isfinite(d).sum()),
                         median_diff_m=float(np.nanmedian(d)), mad_m=float(np.nanmedian(np.abs(d - np.nanmedian(d))))))
    sat = io.table("kgs_satthick_1921")[["TRS", "satthick_2019_2021"]]
    a1921 = t[t.year.isin([2019, 2020, 2021])].groupby("cell_id").A_m.mean()
    k = cells[["cell_id", "kgs_TRS"]].merge(sat, left_on="kgs_TRS", right_on="TRS", how="left")
    d = a1921.reindex(k.cell_id).to_numpy() - k.satthick_2019_2021.to_numpy() * FT
    comp.append(dict(product="KGS section saturated thickness 2019-2021", year="2019-2021", n=int(np.isfinite(d).sum()),
                     median_diff_m=float(np.nanmedian(d)), mad_m=float(np.nanmedian(np.abs(d - np.nanmedian(d))))))
    w2224 = t[t.year.isin([2022, 2023, 2024])].groupby("cell_id").wte_m.mean().reindex(cells.cell_id).to_numpy()
    d = w2224 - cells.kgs_wte_2022_2024.to_numpy() * FT
    comp.append(dict(product="KGS section water table 2022-2024", year="2022-2024", n=int(np.isfinite(d).sum()),
                     median_diff_m=float(np.nanmedian(d)), mad_m=float(np.nanmedian(np.abs(d - np.nanmedian(d))))))
    if 2017 in models and 2019 in models:
        usg = _sample_nc(io.path("usgs_change"), "dWL_2017_to_2019", cells.x5070.to_numpy(), cells.y5070.to_numpy()) * FT
        ch = (pred[pred.year == 2019].set_index("cell_id").wte_m - pred[pred.year == 2017].set_index("cell_id").wte_m)
        d = ch.reindex(cells.cell_id).to_numpy() - usg
        comp.append(dict(product="USGS water-level change 2017-2019 (McGuire & Strauch 2024)", year="2017-2019",
                         n=int(np.isfinite(d).sum()), median_diff_m=float(np.nanmedian(d)),
                         mad_m=float(np.nanmedian(np.abs(d - np.nanmedian(d))))))
    comp = pd.DataFrame(comp)

    # ══════ STEP 2.4-7  targets (reported, never enforced) ══════
    co = panel[panel.role == "cohort"]
    share8 = co.groupby("year").well_km.apply(lambda s: float((s <= C.TARGETS["well_km"]).mean()))
    good = co[~co.A_poor.astype(bool)].drop_duplicates("cell_id")
    good = good.assign(third=good.cell_id.map(third))
    per_cty = good.groupby("county").third.value_counts().unstack(fill_value=0).reindex(columns=["thin", "middle",
                                                                                                 "thick"], fill_value=0)
    both = per_cty[(per_cty.thin >= C.TARGETS["both_thirds_cells"]) & (per_cty.thick >= C.TARGETS["both_thirds_cells"])]
    cell_mean = co[~co.A_poor.astype(bool)].groupby("cell_id").agg(A=("A_m", "mean"), se=("A_se_m", "mean"))
    cell_mean["third"] = cell_mean.index.map(third)
    thin, thick = cell_mean[cell_mean.third == "thin"], cell_mean[cell_mean.third == "thick"]
    gap = float(thick.A.median() - thin.A.median()) if len(thin) and len(thick) else np.nan
    se_gap = float(np.sqrt(thin.se.median() ** 2 + thick.se.median() ** 2)) if len(thin) and len(thick) else np.nan
    ranges = co.groupby("county").A_m.agg(["min", "median", "max"]).reset_index()
    ranges["has_both_thirds"] = ranges.county.isin(both.index)
    tg = pd.DataFrame([
        dict(target="cohort cells within 8 km of a measured well, worst year", target_value=C.TARGETS["well_share"],
             observed=float(share8.min()) if len(share8) else np.nan),
        dict(target="counties with >= 50 cells in both thin and thick thirds", target_value=C.TARGETS["both_thirds_counties"],
             observed=int(len(both))),
        dict(target="thin-to-thick difference / its SE (cell-level errors)", target_value=C.TARGETS["thin_thick_se_ratio"],
             observed=gap / se_gap if se_gap else np.nan),
        dict(target="cohort cell-years with poor thickness (excluded)", target_value=None,
             observed=float(co.A_poor.astype(bool).mean())),
    ])
    vstats = pd.concat([_err_stats(loc_cv, "year").assign(holdout="well locations"),
                        _err_stats(yr_cv, "year").assign(holdout="whole years")], ignore_index=True)

    # ══════ STEP 2.4-8  write: product, validation, comparison, targets, figures ══════
    out = t[["cell_id", "county", "role", "year", "wte_m", "wte_se_m", "base_m", "base_se_m", "A_m", "A_se_m",
             "well_km", "A_poor", "A_third"]]
    report.table(ctx, out, "thickness_cell_year", parquet=True, csv=False,
                 units=dict(wte_m="m", base_m="m", A_m="m", well_km="km"), notes="Paper 2 data product (Kansas)")
    report.text(ctx, "thirds.json", json.dumps(dict(cuts_m=cuts, rule=D.value("D26"), n_cohort_cells=int(len(cm))),
                                               indent=1))
    report.table(ctx, loc_cv, "validation_location", units=dict(err="m"))
    report.table(ctx, yr_cv, "validation_year", units=dict(err="m"))
    report.table(ctx, vstats, "validation_stats", units=dict(bias_m="m", rmse_m="m"), notes="Paper 2 Table 2 draft")
    report.table(ctx, comp, "comparison", units=dict(median_diff_m="m (ours - product)"))
    report.table(ctx, tg, "targets")
    report.table(ctx, ranges, "county_A_range", units=dict(min="m", median="m", max="m"), mask="cohort cell-years")
    report.table(ctx, share8.rename("share").reset_index(), "well_support_by_year", units=dict(share="fraction"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    yr = 2020 if 2020 in models else max(models)
    g = T02.load_geometry(ctx.dep("t0_2")).merge(out[out.year == yr], on=["cell_id", "county", "role"])
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    g.plot(ax=ax[0], column="A_m", cmap="viridis", legend=True, linewidth=0, legend_kwds=dict(label="A (m)"))
    g.plot(ax=ax[1], column="A_se_m", cmap="magma_r", legend=True, linewidth=0, legend_kwds=dict(label="SE of A (m)"))
    for a, ttl in zip(ax, [f"Saturated thickness, {yr} pre-season", "Its standard error"]):
        a.set_axis_off()
        a.set_title(ttl, fontsize=9)
    report.figure(ctx, fig, "thickness_map", units={"A": "m"}, mask="cells in scope", caption="Paper 2 Figure 2 draft")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].scatter(loc_cv.obs, loc_cv.pred, s=2, alpha=0.3)
    lim = [np.nanmin(loc_cv.obs), np.nanmax(loc_cv.obs)]
    ax[0].plot(lim, lim, "k-", lw=0.8)
    ax[0].set_xlabel("measured water table (m)")
    ax[0].set_ylabel("kriged, location held out (m)")
    for h, s in vstats.groupby("holdout"):
        ax[1].plot(s.year, s.rmse_m, marker="o", label=h)
    ax[1].set_ylabel("RMSE (m)")
    ax[1].legend(fontsize=8)
    report.figure(ctx, fig, "validation", units={"water table": "m"}, mask="usable pre-season wells",
                  caption="Paper 2 Figure 3 draft")
    ctx.complete(dict(years=sorted(int(y) for y in models), cuts_m=cuts, targets=tg.to_dict("records"),
                      rmse_location_m=float(np.sqrt((loc_cv.err ** 2).mean())),
                      rmse_year_m=float(np.sqrt((yr_cv.err ** 2).mean())) if len(yr_cv) else None,
                      base_se_by_gmd={str(k): float(v) for k, v in mad.items()}))
