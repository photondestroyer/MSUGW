"""Task 1.1  Weather: deficit, peak month, dry seasons, dry episodes.

"Compute the monthly and seasonal deficit and the peak month for every county and year. Group dry seasons
into dry episodes."
Physical idea: a well's capacity matters most when the crop's demand for water is highest. The peak month is
that moment, and it is chosen from weather only (reference ET minus rain), so it cannot be tuned to any ET
result. Dry seasons, and the episodes that link neighbouring dry counties, set how many independent droughts
the pilot can learn from.
Decisions: D01 (reference ET), D10 (episodes), D20 (dry season percentile).
Deliverables: weather_county_month.csv, weather_county_year.csv, dry_episodes.csv, deficit_heatmap.png.
"""
import numpy as np

QC_NO_DATA = 4


def county_weights(counties4326, lat, lon):
    """For each county: gridMET cell indices (row, col) and weights = area of the cell inside the county."""
    import geopandas as gpd
    from ..common import grids
    out = {}
    le, oe = grids.cell_edges(lat), grids.cell_edges(lon)
    for cty, geom in zip(counties4326["county"], counties4326.geometry):
        x0, y0, x1, y1 = geom.bounds
        rows = np.where((np.maximum(le[:-1], le[1:]) >= y0) & (np.minimum(le[:-1], le[1:]) <= y1))[0]
        cols = np.where((oe[1:] >= x0) & (oe[:-1] <= x1))[0]
        R, Cc = np.meshgrid(rows, cols, indexing="ij")
        R, Cc = R.ravel(), Cc.ravel()
        boxes = grids.latlon_cell_polygons(lat, lon, R, Cc)
        inter = gpd.GeoSeries(boxes, crs=4326).intersection(geom)
        area = inter.to_crs(5070).area.to_numpy()
        keep = area > 0
        out[cty] = (R[keep], Cc[keep], area[keep] / area[keep].sum())
    return out


def run(ctx):
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from .. import definitions as DEF
    from ..common import io, report
    from ..phase0 import counties as CTY

    var = D.value("D01")["variable"]
    cty = CTY.load(ctx.dep("counties")).to_crs(4326)

    # ══════ STEP 1.1-1  county weights on the gridMET lab grid (area of each 1/24 deg cell inside the county) ══════
    ds = io.gridmet()
    lat, lon, dates = io.gridmet_axes(ds)
    W = county_weights(cty, lat, lon)
    rr = np.concatenate([w[0] for w in W.values()])
    cc = np.concatenate([w[1] for w in W.values()])
    r0, r1, c0, c1 = rr.min(), rr.max() + 1, cc.min(), cc.max() + 1
    ctx.log.info("STEP 1.1-1  %d counties on %d gridMET cells", len(W), len(set(zip(rr, cc))))

    # ══════ STEP 1.1-2  daily reference ET (D01) and precipitation, May-September of every panel year ══════
    status = np.asarray(ds["status"][:])
    qc = {v: np.asarray(ds[f"{v}_qc"][:]) for v in (var, "pr")}
    rows = []
    for y in C.PANEL_YEARS:
        sel = np.where((dates.year == y) & dates.month.isin(C.SEASON_MONTHS))[0]
        i0, i1 = sel[0], sel[-1] + 1
        et = io.gridmet_read(ds, var, i0, i1, r0, r1, c0, c1)
        pr = io.gridmet_read(ds, "pr", i0, i1, r0, r1, c0, c1)
        for v, a in ((var, et), ("pr", pr)):
            bad = (qc[v][i0:i1] & QC_NO_DATA) > 0
            a[bad] = np.nan
        d_month = dates.month[i0:i1].to_numpy()
        for c, (R, Cc, w) in W.items():
            e = (et[:, R - r0, Cc - c0] * w).sum(axis=1)               # NaN when any cell is NaN that day
            p = (pr[:, R - r0, Cc - c0] * w).sum(axis=1)
            for m in C.SEASON_MONTHS:
                k = d_month == m
                rows.append(dict(county=c, year=y, month=m, eto_mm=e[k].sum() if np.isfinite(e[k]).all() else np.nan,
                                 p_mm=p[k].sum() if np.isfinite(p[k]).all() else np.nan, days=int(k.sum()),
                                 days_not_permanent=int((status[i0:i1][k] > 0).sum()),
                                 days_qc_flagged=int(((qc[var][i0:i1][k] | qc["pr"][i0:i1][k]) > 0).sum())))
        ctx.log.info("STEP 1.1-2  %d done", y)
    ds.close()
    cm = pd.DataFrame(rows)
    cm["deficit_mm"] = DEF.deficit(cm.eto_mm, cm.p_mm)

    # ══════ STEP 1.1-3  season sums and the peak month (June-August, largest deficit) ══════
    piv = cm.pivot_table(index=["county", "year"], columns="month", values="deficit_mm")
    cy = piv.reset_index()[["county", "year"]].copy()
    cy["deficit_season_mm"] = piv[C.SEASON_MONTHS].sum(axis=1, min_count=len(C.SEASON_MONTHS)).to_numpy()
    cy["peak_month"] = DEF.peak_month(piv[6].to_numpy(), piv[7].to_numpy(), piv[8].to_numpy())
    cy["deficit_peak_mm"] = [piv.loc[(c, y), m] if m else np.nan for c, y, m in zip(cy.county, cy.year, cy.peak_month)]
    s = cm.groupby(["county", "year"])[["eto_mm", "p_mm"]].sum(min_count=5).reset_index()
    cy = cy.merge(s.rename(columns={"eto_mm": "eto_season_mm", "p_mm": "p_season_mm"}), on=["county", "year"])

    # ══════ STEP 1.1-4  dry seasons: above the 70th percentile over all Kansas county-years (D20) ══════
    cy["dry_season"], thr = DEF.season_dry(cy["deficit_season_mm"].to_numpy())
    cy["deficit_season_pct"] = cy["deficit_season_mm"].rank(pct=True) * 100

    # ══════ STEP 1.1-5  dry episodes: dry counties sharing a border in the same season form one (D10) ══════
    nb = CTY.neighbours(ctx.dep("counties"))
    cy["episode_id"], cy["episode_id_multiyear"] = DEF.dry_episodes(cy.year.to_numpy(), cy.county.to_numpy(),
                                                                    cy.dry_season.to_numpy(), nb)
    ep = (cy[cy.episode_id >= 0].groupby("episode_id")
          .agg(year=("year", "first"), counties=("county", lambda x: " ".join(sorted(x))), n_counties=("county", "size"),
               deficit_mean_mm=("deficit_season_mm", "mean"), multiyear_id=("episode_id_multiyear", "first"))
          .reset_index())
    ctx.log.info("STEP 1.1-5  threshold %.0f mm; %d dry county-years; %d episodes (%d multi-year)", thr,
                 int(cy.dry_season.sum()), len(ep), ep.multiyear_id.nunique() if len(ep) else 0)

    # ══════ STEP 1.1-6  write tables and the deficit heat map ══════
    u = {"eto_mm": "mm/month", "p_mm": "mm/month", "deficit_mm": "mm/month"}
    report.table(ctx, cm, "weather_county_month", units=u, parquet=True, notes=f"reference ET = gridMET {var} (D01)")
    report.table(ctx, cy, "weather_county_year", parquet=True,
                 units={"deficit_season_mm": "mm (May-Sep)", "deficit_peak_mm": "mm/month"},
                 notes=f"dry season = deficit above {thr:.1f} mm (70th percentile, D20)")
    report.table(ctx, ep, "dry_episodes", units={"deficit_mean_mm": "mm"})
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    hm = cy.pivot(index="county", columns="year", values="deficit_season_mm")
    dm = cy.pivot(index="county", columns="year", values="dry_season")
    fig, ax = plt.subplots(figsize=(7, 0.18 * len(hm) + 1.5))
    im = ax.imshow(hm.to_numpy(), aspect="auto", cmap="YlOrBr")
    yy, xx = np.where(dm.to_numpy().astype(bool))
    ax.scatter(xx, yy, marker="x", color="k", s=12, label="dry season (D20)")
    ax.set_xticks(range(hm.shape[1]), hm.columns)
    ax.set_yticks(range(hm.shape[0]), hm.index, fontsize=6)
    fig.colorbar(im, ax=ax, label="May-Sep deficit (mm)")
    ax.legend(loc="upper right", fontsize=7)
    report.figure(ctx, fig, "deficit_heatmap", units={"deficit": "mm"}, mask="all study county-years")
    ctx.complete(dict(threshold_mm=thr, county_years=int(len(cy)), dry_county_years=int(cy.dry_season.sum()),
                      episodes=int(len(ep)), episodes_multiyear=int(ep.multiyear_id.nunique()) if len(ep) else 0,
                      peak_month_share={int(k): float(v) for k, v in cy.peak_month.value_counts(normalize=True).items()},
                      missing_months=int(cm.deficit_mm.isna().sum())))
