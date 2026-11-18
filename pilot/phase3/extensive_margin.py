"""Extensive margin: irrigated area against aquifer thickness (the analogue of Mieno et al. 2024, Fig. 3).

SPECIFIED 2026-11-18, BEFORE ANY RESULT OF IT WAS SEEN. Not in the frozen decision log: a secondary analysis.
It reads no ET (the blinding rule covers ET and surface temperature only).
Outcome: a cell-year is irrigated when its LANID irrigated share meets D04 (definitions.is_irrigated_year),
2016-2020 (LANID ends in 2020); unknown where LANID covers too little of the cell (D06), and then left out.
Thickness: cell-year A of Check 2.4 (frozen main run); poor-thickness cell-years (D25) left out.
Baseline cropland: regular cells inside the aquifer, not served by surface water (D09), cropland in every year
2011-2015 (the D07 cropland rule).
  S1  share of baseline cropland irrigated, by 10 m bins of A: (a) mean of county means (between places, the
      Mieno Fig. 3 analogue); (b) bin effects with county x year effects and root-zone storage (within county).
  S2  cohort cells (irrigated in 2011-2015): share still irrigated, by A bin as in S1, and the within-county slope
      per 10 m of A (county x year effects, root-zone storage, corn years; A error corrected as in main_models).
  S3  cohort cells: does irrigation drop more on thin aquifer in dry seasons? Cell and county x year effects,
      A, A x seasonal deficit, controls x deficit (the D30 form); contrast = beta[A x deficit] x (A_thin - A_thick)
      x (deficit q90 - median), in percentage points.
Expected signs: irrigated share rises with A (S1, S2 positive slopes); S3 contrast negative.
Evidence: a 95% interval excluding zero (county-clustered; t with G-1 df below 50 counties).
Run:  PILOT_RUN_MODE=main python -m pilot.phase3.extensive_margin
"""
import os

import numpy as np

EDGES = [0, 10, 20, 30, 40, 50, 60, 70, 80, 100, 1000]
LAN = [f"lanid_{y}" for y in range(2016, 2021)]


def bins(a):
    import pandas as pd
    return pd.cut(a, EDGES, right=False, labels=[f"{lo}-{hi}" if hi < 1000 else f">={lo}"
                                                for lo, hi in zip(EDGES[:-1], EDGES[1:])])


def between(d):
    """Mean of county means per A bin (percent irrigated) and 1.96 SE across counties."""
    cm = d.groupby(["A_bin", "county"], observed=True).irr.mean().reset_index()
    out = cm.groupby("A_bin", observed=True).irr.agg(["mean", "std", "count"]).reset_index()
    out["half95"] = 1.96 * out["std"].fillna(0) / np.sqrt(out["count"].clip(lower=1))
    return out.rename(columns={"count": "counties"})


def within_bins(d, controls):
    """Bin effects (percentage points, against the first populated bin) with county x year effects."""
    import pandas as pd
    from ..common import stats
    cats = [c for c in d.A_bin.cat.categories if (d.A_bin == c).any()]
    X = np.column_stack([(d.A_bin == c).to_numpy(float) for c in cats[1:]] + [d[c].to_numpy(float) for c in controls])
    fe = [(d.county.astype(str) + "_" + d.year.astype(str)).to_numpy()]
    f = stats.fe_ols(d.irr.to_numpy(float), X, fe=fe, cluster=d.county.to_numpy())
    rows = [dict(A_bin=cats[0], effect_pp=0.0, lo_pp=0.0, hi_pp=0.0, n=int((d.A_bin == cats[0]).sum()))]
    for k, c in enumerate(cats[1:]):
        lo, hi = stats.ci95(f["beta"][k], f["se"][k], f["n_clusters"])
        rows.append(dict(A_bin=c, effect_pp=float(f["beta"][k]), lo_pp=float(lo), hi_pp=float(hi),
                         n=int((d.A_bin == c).sum())))
    return pd.DataFrame(rows), int(f["n_clusters"])


def slope_per_10m(d, rs):
    """S2: within-county slope of irrigated (pp) on A, per 10 m, A error corrected (main_models.fit)."""
    import pandas as pd
    from ..common import stats
    from . import main_models as MM
    X = np.column_stack([d.A_m.to_numpy(float), d.rzws_mm.to_numpy(float), d.corn_years.to_numpy(float)])
    fe = [(d.county.astype(str) + "_" + d.year.astype(str)).to_numpy()]
    T = np.ones((len(d), 1))
    G = MM.expected_error_gram(T, fe, d.wte_se_m.to_numpy(), d.base_se_m.to_numpy(), pd.factorize(d.cell_id)[0], 2, rs)
    f = MM.fit(d.irr.to_numpy(float), X, fe, d.county.to_numpy(), G)
    out = {}
    for tag, b, V in (("corrected", f["beta_c"], f["cov_c"]), ("uncorrected", f["beta"], f["cov"])):
        s = float(np.sqrt(max(V[0, 0], 0))) * 10
        lo, hi = stats.ci95(float(b[0]) * 10, s, f["n_clusters"])
        out[tag] = dict(pp_per_10m=float(b[0]) * 10, lo=float(lo), hi=float(hi))
    return out, int(f["n_clusters"])


def drought(d, rs):
    """S3: the D30 form on the irrigated indicator; contrast thin - thick at a dry season, in pp."""
    import pandas as pd
    from . import main_models as MM
    X, T, names = MM.design(d, [])
    fe = [d.cell_id.to_numpy(), (d.county.astype(str) + "_" + d.year.astype(str)).to_numpy()]
    G = MM.expected_error_gram(T, fe, d.wte_se_m.to_numpy(), d.base_se_m.to_numpy(), pd.factorize(d.cell_id)[0],
                               X.shape[1] - T.shape[1], rs)
    f = MM.fit(d.irr.to_numpy(float), X, fe, d.county.to_numpy(), G)
    q90, q50 = np.quantile(d.deficit_season_mm, [0.9, 0.5])
    a_thin, a_thick = MM.thirds_mean_A(d)
    return dict(MM.contrast(f, names.index("A:D_season"), a_thin, a_thick, q90 - q50), A_thin_m=a_thin,
                A_thick_m=a_thick, deficit_q90_mm=float(q90), deficit_median_mm=float(q50),
                clusters=int(f["n_clusters"]), n=len(d), cells=int(d.cell_id.nunique()))


def main():
    import json
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .. import config as C
    from .. import decisions as D
    from .. import definitions as DEF
    from .. import freeze
    from ..common import report
    from ..core.errors import PipelineStop
    from ..phase0 import t0_2_cells as T02
    from .explore_t1_by_month import _Ctx
    if C.RUN_MODE != "main":
        raise PipelineStop("reads the frozen main-mode thickness: run with PILOT_RUN_MODE=main")
    rec = freeze.check()
    ctx = _Ctx()
    out = f"{C.out_root()}/explore/extensive_margin"
    os.makedirs(out, exist_ok=True)
    crop = [f"cdl_cropland_{y}" for y in range(2011, 2016)]
    cells = T02.load(ctx.dep("t0_2"), columns=["cell_id", "county", "regular", "in_aquifer", "surface_water", "role",
                                               "rzws_mm", "corn_years", *crop, *LAN])
    z = np.zeros((len(cells), len(crop)))           # the D07 cropland test without its irrigation test
    cells["base_cropland"] = DEF.is_rainfed(z, z, cells[crop].to_numpy(), cells.regular)[0] & \
        cells.in_aquifer.astype(bool) & ~cells.surface_water.astype(bool)
    long = cells.melt(id_vars=["cell_id", "county", "role", "base_cropland", "rzws_mm", "corn_years"],
                      value_vars=LAN, var_name="y", value_name="share")
    long["year"] = long.y.str[-4:].astype(int)
    long = long[np.isfinite(long.share)]
    long["irr"] = 100.0 * DEF.is_irrigated_year(long.share.to_numpy())          # percent
    th = report.read(ctx.dep("c2_4"), "thickness_cell_year")[["cell_id", "year", "A_m", "wte_se_m", "base_se_m",
                                                              "A_poor", "A_third"]]
    wy = report.read(ctx.dep("t1_1"), "weather_county_year")[["county", "year", "deficit_season_mm"]]
    wy["county"] = wy.county.astype(str)
    long["county"] = long.county.astype(str)
    d = long.merge(th, on=["cell_id", "year"]).merge(wy, on=["county", "year"])
    d = d[~d.A_poor.astype(bool) & np.isfinite(d.A_m) & np.isfinite(d.rzws_mm)].reset_index(drop=True)
    d["A_bin"] = bins(d.A_m)
    rs = np.random.default_rng(D.value("D28")["seed"])

    s1 = d[d.base_cropland]
    co = d[d.role == "cohort"]
    s1b, s1w = between(s1), within_bins(s1, ["rzws_mm"])
    s2b, s2w = between(co), within_bins(co, ["rzws_mm", "corn_years"])
    s2s, _ = slope_per_10m(co, rs)
    s3 = drought(co, rs)
    tables = {"S1_between": s1b, "S1_within": s1w[0], "S2_between": s2b, "S2_within": s2w[0]}
    for k, t in tables.items():
        t.to_csv(f"{out}/{k}.csv", index=False)
    summ = dict(freeze=rec["tag"], S1_cell_years=len(s1), S1_cells=int(s1.cell_id.nunique()),
                S1_counties=s1w[1], S2_cell_years=len(co), S2_cells=int(co.cell_id.nunique()),
                S2_slope_per_10m=s2s, S3_drought_contrast_pp=s3)
    with open(f"{out}/summary.json", "w") as fh:
        json.dump(summ, fh, indent=1, default=float)

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
    for a, (b, w, title) in zip(ax, ((s1b, s1w[0], "S1: share of baseline cropland irrigated"),
                                     (s2b, s2w[0], "S2: cohort cells still irrigated"))):
        x = np.arange(len(b))
        a.errorbar(x - 0.12, b["mean"], yerr=b.half95, fmt="o", color="k", capsize=2, label="between places")
        w = w.set_index("A_bin").reindex(b.A_bin)
        anchor = float(b["mean"].iloc[0])
        a.errorbar(x + 0.12, anchor + w.effect_pp, yerr=[w.effect_pp - w.lo_pp, w.hi_pp - w.effect_pp], fmt="s",
                   color="#1d4ed8", mfc="none", capsize=2, label="within county-year (anchored at first bin)")
        a.set_xticks(x, b.A_bin.astype(str), rotation=45, fontsize=8)
        a.set_xlabel("Aquifer thickness A (m)")
        a.set_ylabel("Irrigated (% of cell-years), 2016-2020")
        a.set_title(title, fontsize=10)
    ax[0].legend(frameon=False, fontsize=8)
    fig.text(0.01, -0.02, "LANID 2016-2020, D04 rule; A from Check 2.4, poor cell-years out; 95% across counties "
             f"(between) or county-clustered (within); specified 2026-11-18 before running; {rec['tag']}", fontsize=6)
    fig.tight_layout()
    fig.savefig(f"{out}/extensive_margin.png", dpi=150, bbox_inches="tight")
    with open(f"{out}/README.md", "w", encoding="utf-8") as fh:
        fh.write("# Extensive margin (secondary; specified before running)\n\n" + __doc__ + "\n\n" +
                 "\n\n".join(f"## {k}\n\n{report.md(t.round(2))}" for k, t in tables.items()) +
                 "\n\n## Summary\n\n```\n" + json.dumps(summ, indent=1, default=float) + "\n```\n")
    print(json.dumps(summ, indent=1, default=float))
    for k, t in tables.items():
        print(k)
        print(t.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
