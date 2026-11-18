"""Extensive margin, long window (pilot Section 11 fallback): a cohort labelled from 2006-2010, outcome 2011-2020.

SPECIFIED 2026-11-18, BEFORE ANY RESULT OF IT WAS SEEN. Secondary analysis; not in the frozen decision log; no ET.
Why: the 2016-2020 window has 5 map years and no severe drought. LANID covers 1997-2020 and thickness is kriged
for 2011-2024, so 2011-2020 gives 10 years including the 2011-2012 drought, with no new data.
Cohort F: regular cells inside the aquifer, not served by surface water (D09), LANID irrigated (D04 rule, D06
coverage) in at least COHORT_MIN_IRR_YEARS (4) of 2006-2010, LANID known in all five years, and LANID and
AIM-HPA agreeing in at least D06 min_agree_years of 2006-2010 (definitions.maps_agree; AIM-HPA 2006-2010 read from
"pilot test/AIM HPA"; this rule added 2026-11-18 before the second run, replacing the no-agreement shortcut).
  shortcut: no corn rule (the outcome is irrigation of any crop); add it if this becomes a primary analysis.
S4 also draws AIM-HPA (2011-2017) for the same cells, dashed, as a check on single-year map artefacts.
Outcome: irrigated per LANID (D04), 2011-2020; unknown cell-years left out. Thickness: cell-year A of Check 2.4
(frozen main run), poor cell-years (D25) left out; thirds re-cut for cohort F by the D26 rule on its 2011-2020 means.
Deficit: seasonal (May-Sep) reference ET minus rain, D01, county-area-weighted gridMET, 2011-2020 (the Task 1.1
computation, kept here so the frozen module stays untouched).
Control: root-zone storage. Corn years 2013-2015 are NOT a control here: they fall inside the outcome window and are
partly an outcome of irrigation.
  S1  share of baseline cropland (the 2016-2020 analysis's population) irrigated, 2011-2020, by A bin.
  S2  cohort F still irrigated, by A bin, and within-county slope per 10 m (A error corrected).
  S3  cohort F: thin minus thick change in irrigated (pp) in a dry season (q90 vs median deficit), cell and
      county x year effects, A error corrected (D30 form).
  S4  (display) share of cohort F irrigated each year 2011-2020, by thickness third, with the deficit.
Expected signs: S1, S2 slopes positive; S3 contrast negative. Evidence: 95% interval excluding zero.
Run:  PILOT_RUN_MODE=main python -m pilot.phase3.extensive_margin_2011_2020
"""
import os

import numpy as np

COHORT_YEARS = list(range(2006, 2011))
YEARS = list(range(2011, 2021))


def lanid_shares(ctx, years):
    """LANID irrigated share per cell for the given years (D06 coverage rule), as t0_2 computes 2011-2020."""
    from .. import definitions as DEF
    from ..common import io
    from ..phase0 import t0_2_cells as T02
    geo = T02.load_geometry(ctx.dep("t0_2"))
    ln = io.lanid()
    ly = [int(v) for v in ln["time"][:]]
    lvar = ln["irrigated"]

    def win(name, r0, r1, c0, c1):
        a = np.asarray(lvar[ly.index(int(name[1:])), r0:r1, c0:c1])
        return a == 1, (a == 0) | (a == 1)
    acc = T02._accumulate_5070(ctx, geo, io.lanid_transform(ln), lvar.shape[1:], win, [f"l{y}" for y in years],
                               "LANID 2006-2010")
    ln.close()
    out = geo[["cell_id"]].copy()
    for y in years:
        out[f"lanid_{y}"] = DEF.map_share(acc.cnt[f"l{y}_hit"], acc.cnt[f"l{y}_map"], acc.cnt["px"])[0]
    return out


def aim_shares(ctx, years):
    """AIM-HPA irrigated share per cell (D06 coverage rule), from the full 2002-2017 AIM-HPA folder."""
    import rasterio
    from .. import config as C
    from .. import definitions as DEF
    from ..common import io
    from ..phase0 import t0_2_cells as T02
    geo = T02.load_geometry(ctx.dep("t0_2"))
    src = {y: rasterio.open(f"{C.PILOT_TEST}/AIM HPA/{y}_AIM-HPA_rse_finalMapBinary_20190909.tif") for y in years}
    first = src[years[0]]

    def win(name, r0, r1, c0, c1):
        a = src[int(name[1:])].read(1, window=((r0, r1), (c0, c1)))
        return a == io.AIM_IRRIGATED, (a == io.AIM_IRRIGATED) | (a == io.AIM_NOT_IRRIGATED)
    acc = T02._accumulate_5070(ctx, geo, first.transform, first.shape, win, [f"a{y}" for y in years],
                               "AIM-HPA 2006-2010")
    for s in src.values():
        s.close()
    out = geo[["cell_id"]].copy()
    for y in years:
        out[f"aim_{y}"] = DEF.map_share(acc.cnt[f"a{y}_hit"], acc.cnt[f"a{y}_map"], acc.cnt["px"])[0]
    return out


def season_deficit(ctx, years):
    """County seasonal deficit (mm), D01: area-weighted gridMET reference ET minus rain, May-September."""
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from .. import definitions as DEF
    from ..common import io
    from ..phase0 import counties as CTY
    from ..phase1.t1_1_weather import county_weights, QC_NO_DATA
    var = D.value("D01")["variable"]
    ds = io.gridmet()
    lat, lon, dates = io.gridmet_axes(ds)
    W = county_weights(CTY.load(ctx.dep("counties")).to_crs(4326), lat, lon)
    rr = np.concatenate([w[0] for w in W.values()])
    cc = np.concatenate([w[1] for w in W.values()])
    r0, r1, c0, c1 = rr.min(), rr.max() + 1, cc.min(), cc.max() + 1
    qc = {v: np.asarray(ds[f"{v}_qc"][:]) for v in (var, "pr")}
    rows = []
    for y in years:
        sel = np.where((dates.year == y) & dates.month.isin(C.SEASON_MONTHS))[0]
        i0, i1 = sel[0], sel[-1] + 1
        a = {v: io.gridmet_read(ds, v, i0, i1, r0, r1, c0, c1) for v in (var, "pr")}
        for v in a:
            a[v][(qc[v][i0:i1] & QC_NO_DATA) > 0] = np.nan
        for c, (R, Cc, w) in W.items():
            e = (a[var][:, R - r0, Cc - c0] * w).sum(axis=1)
            p = (a["pr"][:, R - r0, Cc - c0] * w).sum(axis=1)
            rows.append(dict(county=str(c), year=y, deficit_season_mm=float(DEF.deficit(e, p).sum())
                             if np.isfinite(e).all() and np.isfinite(p).all() else np.nan))
    ds.close()
    return pd.DataFrame(rows)


def main():
    import json
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from .. import definitions as DEF
    from .. import freeze
    from ..common import report
    from ..core.errors import PipelineStop
    from ..phase0 import t0_2_cells as T02
    from . import extensive_margin as EM
    from .explore_t1_by_month import _Ctx
    if C.RUN_MODE != "main":
        raise PipelineStop("reads the frozen main-mode thickness: run with PILOT_RUN_MODE=main")
    rec = freeze.check()

    class Ctx(_Ctx):                                  # what t0_2's accumulator needs from a task context
        log = __import__("logging").getLogger("pilot.extensive_margin_2011_2020")
        check_budget = staticmethod(lambda: None)
    ctx = Ctx()
    out = f"{C.out_root()}/explore/extensive_margin_2011_2020"
    os.makedirs(out, exist_ok=True)

    # ══════ STEP L-1  cells, cohort F from 2006-2010 labels, outcome 2011-2020 ══════
    crop = [f"cdl_cropland_{y}" for y in range(2011, 2016)]
    lan = [f"lanid_{y}" for y in YEARS]
    aim_late = [f"aim_{y}" for y in range(2011, 2018)]
    cells = T02.load(ctx.dep("t0_2"), columns=["cell_id", "county", "regular", "in_aquifer", "surface_water",
                                               "rzws_mm", *crop, *lan, *aim_late])
    cells = cells.merge(lanid_shares(ctx, COHORT_YEARS), on="cell_id", how="left") \
                 .merge(aim_shares(ctx, COHORT_YEARS), on="cell_id", how="left")
    e = cells[[f"lanid_{y}" for y in COHORT_YEARS]].to_numpy(float)
    ea = cells[[f"aim_{y}" for y in COHORT_YEARS]].to_numpy(float)
    eligible = cells.regular.astype(bool) & cells.in_aquifer.astype(bool) & ~cells.surface_water.astype(bool)
    cells["cohort_f"] = eligible & np.isfinite(e).all(axis=1) & \
        (DEF.is_irrigated_year(e).sum(axis=1) >= C.COHORT_MIN_IRR_YEARS) & DEF.maps_agree(e, ea)[1]
    z = np.zeros((len(cells), len(crop)))           # the D07 cropland test without its irrigation test
    cells["base_cropland"] = DEF.is_rainfed(z, z, cells[crop].to_numpy(), cells.regular)[0] & eligible
    long = cells.melt(id_vars=["cell_id", "county", "cohort_f", "base_cropland", "rzws_mm"], value_vars=lan,
                      var_name="y", value_name="share")
    long["year"] = long.y.str[-4:].astype(int)
    long = long[np.isfinite(long.share)]
    long["irr"] = 100.0 * DEF.is_irrigated_year(long.share.to_numpy())
    long["county"] = long.county.astype(str)
    long["corn_years"] = 0.0                          # not a control here (see the docstring); constant drops out

    th = report.read(ctx.dep("c2_4"), "thickness_cell_year")[["cell_id", "year", "A_m", "wte_se_m", "base_se_m",
                                                              "A_poor"]]
    wy = season_deficit(ctx, YEARS)
    d = long.merge(th, on=["cell_id", "year"]).merge(wy, on=["county", "year"])
    d = d[~d.A_poor.astype(bool) & np.isfinite(d.A_m) & np.isfinite(d.rzws_mm) &
          np.isfinite(d.deficit_season_mm)].reset_index(drop=True)
    d["A_bin"] = EM.bins(d.A_m)
    co = d[d.cohort_f].copy()
    cm = co.groupby("cell_id").A_m.mean()
    cuts = DEF.thickness_thirds(cm.to_numpy())
    co["A_third"] = co.cell_id.map(pd.Series(DEF.assign_third(cm.to_numpy(), cuts), index=cm.index))
    rs = np.random.default_rng(D.value("D28")["seed"])

    # ══════ STEP L-2  S1-S3 with the estimators of the 2016-2020 analysis ══════
    s1 = d[d.base_cropland]
    s1b, s1w = EM.between(s1), EM.within_bins(s1, ["rzws_mm"])
    s2b, s2w = EM.between(co), EM.within_bins(co, ["rzws_mm"])
    s2s, _ = EM.slope_per_10m(co, rs)
    s3 = EM.drought(co, rs)
    s4 = co.groupby(["year", "A_third"]).irr.mean().unstack()[["thin", "middle", "thick"]]
    am = cells[cells.cohort_f].melt(id_vars="cell_id", value_vars=aim_late, var_name="y", value_name="share")
    am["year"] = am.y.str[-4:].astype(int)
    am = am[np.isfinite(am.share)].merge(co[["cell_id", "year", "A_third"]], on=["cell_id", "year"])
    am["irr"] = 100.0 * DEF.is_irrigated_year(am.share.to_numpy())
    s4_aim = am.groupby(["year", "A_third"]).irr.mean().unstack()[["thin", "middle", "thick"]]
    defy = wy.groupby("year").deficit_season_mm.mean()
    tables = {"S1_between": s1b, "S1_within": s1w[0], "S2_between": s2b, "S2_within": s2w[0],
              "S4_by_year": s4.reset_index(), "S4_by_year_AIM": s4_aim.reset_index()}
    for k, t in tables.items():
        t.to_csv(f"{out}/{k}.csv", index=False)
    summ = dict(freeze=rec["tag"], years=f"{YEARS[0]}-{YEARS[-1]}", cohort_from=f"{COHORT_YEARS[0]}-{COHORT_YEARS[-1]}",
                cohort_f_cells=int(cells.cohort_f.sum()), cohort_f_cell_years=len(co), thirds_cuts_m=cuts,
                S1_cell_years=len(s1), S1_counties=s1w[1], S2_slope_per_10m=s2s, S3_drought_contrast_pp=s3)
    with open(f"{out}/summary.json", "w") as fh:
        json.dump(summ, fh, indent=1, default=float)

    # ══════ STEP L-3  figure: S1, S2 by thickness; S4 by year with the deficit ══════
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    for a, (b, w, title) in zip(ax[:2], ((s1b, s1w[0], "S1: baseline cropland irrigated, 2011-2020"),
                                         (s2b, s2w[0], "S2: 2006-2010 cohort still irrigated, 2011-2020"))):
        x = np.arange(len(b))
        a.errorbar(x - 0.12, b["mean"], yerr=b.half95, fmt="o", color="k", capsize=2, label="between places")
        w = w.set_index("A_bin").reindex(b.A_bin)
        anchor = float(b["mean"].iloc[0])
        a.errorbar(x + 0.12, anchor + w.effect_pp, yerr=[w.effect_pp - w.lo_pp, w.hi_pp - w.effect_pp], fmt="s",
                   color="#1d4ed8", mfc="none", capsize=2, label="within county-year (anchored at first bin)")
        a.set_xticks(x, b.A_bin.astype(str), rotation=45, fontsize=8)
        a.set_xlabel("Aquifer thickness A (m)")
        a.set_ylabel("Irrigated (% of cell-years)")
        a.set_title(title, fontsize=10)
    ax[0].legend(frameon=False, fontsize=8)
    a = ax[2]
    for g in ("thin", "middle", "thick"):
        a.plot(s4.index, s4[g], marker="o", color=C.COLORS[g], label=f"{g} (LANID)")
        a.plot(s4_aim.index, s4_aim[g], linestyle="--", color=C.COLORS[g], alpha=0.7, label=f"{g} (AIM-HPA)")
    a.set_ylabel("Cohort F irrigated (%)")
    a.set_title("S4: cohort F irrigated by year and thickness third", fontsize=10)
    a.legend(frameon=False, fontsize=8, loc="lower left")
    a2 = a.twinx()
    a2.bar(defy.index, defy.values, color="#d4d4d4", alpha=0.6, zorder=0)
    a2.set_ylabel("Mean seasonal deficit (mm)", color="#6b7280")
    a.set_zorder(a2.get_zorder() + 1)
    a.patch.set_visible(False)
    fig.text(0.01, -0.03, "LANID (D04, D06); cohort F irrigated >= 4 of 2006-2010, LANID and AIM-HPA agreeing >= 3 "
             "years; A from Check 2.4, poor cell-years "
             f"out; 95% across counties (between) or county-clustered (within); specified 2026-11-18 before running; "
             f"{rec['tag']}", fontsize=6)
    fig.tight_layout()
    fig.savefig(f"{out}/extensive_margin_2011_2020.png", dpi=150, bbox_inches="tight")
    with open(f"{out}/README.md", "w", encoding="utf-8") as fh:
        fh.write("# Extensive margin 2011-2020 (cohort from 2006-2010; secondary, specified before running)\n\n" +
                 __doc__ + "\n\n" + "\n\n".join(f"## {k}\n\n{report.md(t.round(2))}" for k, t in tables.items()) +
                 "\n\n## Summary\n\n```\n" + json.dumps(summ, indent=1, default=float) + "\n```\n")
    print(json.dumps(summ, indent=1, default=float))
    for k, t in tables.items():
        print(k)
        print(t.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
