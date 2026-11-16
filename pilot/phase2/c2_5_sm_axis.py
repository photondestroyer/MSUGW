"""Check 2.5  Soil-moisture axis.

"Rainfed cells only. Relate peak-month ET in mm to S, with cell effects, and account for which month is the peak
and for crop and soil. Compare S with the weather deficit as a consistency check."
Report the slope and its 95% interval. A positive slope is expected. If the interval includes zero, find out
whether the cause is low precision or unsuitable reference footprints.
Physical idea: a rainfed crop has no well; its ET should fall as the soil dries. If S does not move rainfed ET,
it cannot be the drought axis for the irrigated cohort either. Crops are never pooled (pilot Section 4).
Deliverables: c2_5_slopes.csv, c2_5_consistency.csv, rainfed_curve.csv + rainfed_curve.png (Paper 1 Figure 1,
rainfed curve), c2_5.md.
"""
import numpy as np


def run(ctx):
    import pandas as pd
    from .. import config as C
    from ..common import report, stats
    from ..phase0 import t0_2_cells as T02
    # ══════ STEP 2.5-1  rainfed peak-month ET (blinded table holds every rainfed month) with S and weather ══════
    et = report.read(ctx.dep("et_cells"), "et_cell_month")
    et = et[et.role == "rainfed"].copy()
    w = report.read(ctx.dep("t1_1"), "weather_county_year")
    wm = report.read(ctx.dep("t1_1"), "weather_county_month")
    S = report.read(ctx.dep("t1_3"), "county_month_S")
    for d in (w, wm, S, et):
        d["county"] = d["county"].astype(str)
    cells = T02.load(ctx.dep("t0_2"), columns=["cell_id", "crop_class", "rzws_mm"])
    d = (et.merge(w[["county", "year", "peak_month", "deficit_peak_mm"]], on=["county", "year"])
         .query("month == peak_month")
         .merge(S[["county", "year", "month", "S", "S_valid"]], on=["county", "year", "month"], how="left")
         .merge(wm[["county", "year", "month", "eto_mm"]], on=["county", "year", "month"], how="left")
         .merge(cells, left_on="cell", right_on="cell_id", how="left"))
    d = d[d.S_valid.astype("boolean").fillna(False).astype(bool) & np.isfinite(d[C.ET_ENSEMBLE])]
    d["et_ratio"] = d[C.ET_ENSEMBLE] / d.eto_mm
    ctx.log.info("STEP 2.5-1  %d rainfed peak-month cell-years with S", len(d))

    # ══════ STEP 2.5-2  ET ~ S + peak-month dummies | cell FE, clustered by county; per crop (never pooled) ══════
    rows = []
    for crop, g in d.groupby("crop_class"):
        if g.cell.nunique() < 5 or g.year.nunique() < 3:
            rows.append(dict(crop=crop, n=len(g), slope_mm_per_sd=np.nan, note="too few cells or years"))
            continue
        X = np.column_stack([g.S] + [(g.peak_month == m).astype(float) for m in (7, 8)])
        cl, cl_label = stats.cluster_ids(g.county, year=g.year, fallback="county_year")
        for y_name in (C.ET_ENSEMBLE, "et_ratio"):
            r = stats.fe_ols(g[y_name].to_numpy(), X, fe=[g.cell.to_numpy()], cluster=cl,
                             names=["S", "peak_jul", "peak_aug"])
            lo, hi = stats.ci95(r["beta"][0], r["se"][0], r["n_clusters"])
            rows.append(dict(crop=crop, outcome=y_name, n=r["n"], cells=int(g.cell.nunique()),
                             counties=int(g.county.nunique()), clusters=r["n_clusters"], cluster_level=cl_label,
                             slope=float(r["beta"][0]), se=float(r["se"][0]), ci_low=float(lo), ci_high=float(hi),
                             includes_zero=bool(lo <= 0 <= hi),
                             mean_footprints=float(g.drop_duplicates(["county", "year"]).merge(
                                 S[["county", "year", "month", "n_footprints"]], on=["county", "year", "month"],
                                 how="left").n_footprints.mean())))
    sl = pd.DataFrame(rows)

    # ══════ STEP 2.5-3  consistency of S with the weather deficit (county-years, peak month) ══════
    cyp = w.merge(S, left_on=["county", "year", "peak_month"], right_on=["county", "year", "month"])
    cyp = cyp[cyp.S_valid.astype("boolean").fillna(False).astype(bool)]
    cons = pd.DataFrame([dict(pairs=len(cyp), corr_S_peak_deficit=float(np.corrcoef(cyp.S, cyp.deficit_peak_mm)[0, 1])
                              if len(cyp) > 2 else np.nan,
                              corr_S_season_deficit=float(np.corrcoef(cyp.S, cyp.deficit_season_mm)[0, 1])
                              if len(cyp) > 2 else np.nan)])

    # ══════ STEP 2.5-4  the rainfed curve: ET (mm and ratio to ETo) by S bin, per crop ══════
    d["S_bin"] = pd.cut(d.S, np.arange(-3, 3.01, 0.5))
    curve = (d.groupby(["crop_class", "S_bin"], observed=True)
             .agg(n=("cell", "size"), et_mm=(C.ET_ENSEMBLE, "mean"), et_ratio=("et_ratio", "mean")).reset_index())
    curve["S_mid"] = curve.S_bin.apply(lambda b: b.mid).astype(float)
    # ══════ STEP 2.5-5  write ══════
    report.table(ctx, sl, "c2_5_slopes", units=dict(slope="mm per unit S (ET) or ratio per unit S"),
                 mask="rainfed cells, peak month, S valid")
    report.table(ctx, cons, "c2_5_consistency")
    report.table(ctx, curve.drop(columns="S_bin"), "rainfed_curve", units=dict(et_mm="mm/month", S_mid="z"),
                 mask="rainfed cells, peak month, S valid")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for crop, g in curve.groupby("crop_class"):
        g = g[g.n >= 20]
        ax[0].plot(g.S_mid, g.et_mm, marker="o", label=crop)
        ax[1].plot(g.S_mid, g.et_ratio, marker="o", label=crop)
    ax[0].set_ylabel("peak-month ET (mm)")
    ax[1].set_ylabel("ET / reference ET")
    for a in ax:
        a.set_xlabel("soil-moisture state S (z)")
    ax[0].legend(fontsize=7, title="rainfed crop")
    report.figure(ctx, fig, "rainfed_curve", units={"ET": "mm/month", "S": "z"},
                  mask="rainfed cells, peak month, S valid", caption="Paper 1 Figure 1 (rainfed curve) draft")
    main = sl[(sl.get("outcome") == C.ET_ENSEMBLE)] if "outcome" in sl else sl
    L = ["# Check 2.5 Soil-moisture axis", "", "Rainfed cells only; slope of peak-month ET on S with cell effects and "
         "peak-month dummies, clustered by county; crops never pooled.", "", report.md(sl), "",
         "## S against the weather deficit", "", report.md(cons), "",
         "If an interval includes zero: compare its SE with the number of cells (precision) and the mean number of "
         "footprints behind S (reference quality). The reference definition may be revised once, tested on "
         "county-years not used to make it (pilot)."]
    report.text(ctx, "c2_5.md", "\n".join(L))
    ctx.complete(dict(slopes=main.to_dict("records"), consistency=cons.to_dict("records")[0]))
