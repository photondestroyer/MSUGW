"""EXPLORATORY (not pre-registered): the T1 timing model with its thickness terms split by calendar month.

Frozen T1 (D30) has A x peak-month indicator and A x exposure x peak indicator. Here the same monthly data, fixed
effects (cell x month, county x month x year), controls and thickness-error correction (main_models) are used,
with one A term and one A x monthly-deficit term PER MONTH (May-September):
    A x month_m            within a cell, how month m's ET moves with that cell's own A from year to year
    A x deficit x month_m  in month m, how the ET response to that month's deficit differs with A
For each month: contrast_m = beta[A x deficit x month_m] * (A_thin - A_thick) * (deficit q90_m - median_m),
the thin-minus-thick ET difference (mm) at a high deficit of that month against its median.
Reads the frozen main-mode outputs; writes only to <out>/explore/t1_by_month (never into a frozen task folder).
Run:  PILOT_RUN_MODE=main python -m pilot.phase3.explore_t1_by_month
"""
import os

import numpy as np

MONTHS = [5, 6, 7, 8, 9]


class _Ctx:
    """Read-only stand-in for a task context: dep() resolves the frozen outputs that m_models reads."""
    def __init__(self):
        from .. import pipeline
        self.runner = pipeline.runner()
        self.stage = next(s for s in pipeline.STAGES if s["name"] == "m_models")

    def dep(self, name):
        from .. import config as C
        from ..core.errors import PipelineStop
        d = f"{C.out_root()}/{name}/{self.runner.sig(name)}"
        if not os.path.exists(f"{d}/_complete.json"):
            raise PipelineStop(f"frozen output {name} is not complete ({d})")
        return d


def fit_by_month(d, outcome, rs):
    """One regression on monthly rows; returns one row per month with the corrected estimates and contrast."""
    import pandas as pd
    from ..common import stats
    from . import main_models as MM
    need = ["A_m", "deficit_mm", outcome, *MM.CONTROLS, "wte_se_m", "base_se_m"]
    d = d[np.isfinite(d[need].to_numpy(float)).all(axis=1) & d.month.isin(MONTHS)].reset_index(drop=True)
    ex = d.deficit_mm.to_numpy(float)
    M = np.column_stack([(d.month == m).to_numpy(float) for m in MONTHS])
    T = np.column_stack([M, M * ex[:, None]])                       # month_m, deficit x month_m
    Z = d[list(MM.CONTROLS)].to_numpy(float)
    X = np.column_stack([d.A_m.to_numpy(float)[:, None] * T, Z * ex[:, None]])
    m = d.month.astype(str)
    fe = [(d.cell_id + "_" + m).to_numpy(), (d.county.astype(str) + "_" + m + "_" + d.year.astype(str)).to_numpy()]
    cl, cl_label = stats.cluster_ids(d.county, d.x5070.to_numpy(), d.y5070.to_numpy(), fallback="block")
    G = MM.expected_error_gram(T, fe, d.wte_se_m.to_numpy(), d.base_se_m.to_numpy(), pd.factorize(d.cell_id)[0],
                               Z.shape[1], rs)
    f = MM.fit(d[outcome].to_numpy(float), X, fe, cl, G)
    a_thin, a_thick = MM.thirds_mean_A(d)
    rows = []
    for i, mo in enumerate(MONTHS):
        j = len(MONTHS) + i
        x = d.loc[d.month == mo, "deficit_mm"].to_numpy(float)
        q90, q50 = np.quantile(x, [0.9, 0.5])
        b_a, s_a = float(f["beta_c"][i]), float(np.sqrt(max(f["cov_c"][i, i], 0)))
        lo_a, hi_a = stats.ci95(b_a, s_a, f["n_clusters"])
        c = MM.contrast(f, j, a_thin, a_thick, q90 - q50)
        rows.append(dict(month=mo, outcome=outcome, n=int((d.month == mo).sum()), clusters=int(f["n_clusters"]),
                         cluster_level=cl_label, beta_A=b_a, lo_A=float(lo_a), hi_A=float(hi_a),
                         beta_A_deficit=float(f["beta_c"][j]), se_A_deficit=float(np.sqrt(max(f["cov_c"][j, j], 0))),
                         deficit_q90_mm=float(q90), deficit_median_mm=float(q50), A_thin_m=a_thin, A_thick_m=a_thick,
                         **{f"contrast_{k}": v for k, v in c.items()}))
    return rows


def write_readme(out, e):
    from ..common import report
    with open(f"{out}/README.md", "w", encoding="utf-8") as fh:
        fh.write("# EXPLORATORY: T1 by month (not pre-registered)\n\n" + __doc__ + "\n\n" + report.md(
            e[["month", "n", "beta_A", "lo_A", "hi_A", "contrast_value_mm", "contrast_lo_mm", "contrast_hi_mm"]]
            .round(3)))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from .. import freeze
    from ..core.errors import PipelineStop
    from . import main_models as MM
    if C.RUN_MODE != "main":
        raise PipelineStop("this exploratory analysis reads the main-mode outputs: run with PILOT_RUN_MODE=main")
    rec = freeze.check()
    out = f"{C.out_root()}/explore/t1_by_month"
    os.makedirs(out, exist_ok=True)
    d, ets = MM.assemble(_Ctx())
    d = d[d.mask_timing.astype(bool)]
    rs = np.random.default_rng(D.value("D28")["seed"])
    res = pd.DataFrame([r for e in ets for r in fit_by_month(d, e, rs)])
    res.to_csv(f"{out}/t1_by_month.csv", index=False)
    e = res[res.outcome == C.ET_ENSEMBLE]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    lab = ["May", "Jun", "Jul", "Aug", "Sep"]
    ax[0].errorbar(e.month, e.beta_A, yerr=[e.beta_A - e.lo_A, e.hi_A - e.beta_A], fmt="o-", color="k", capsize=3)
    ax[0].set_ylabel("A x month: mm ET per m of thickness\n(within a cell, year to year)")
    ax[1].errorbar(e.month, e.contrast_value_mm, yerr=[e.contrast_value_mm - e.contrast_lo_mm,
                   e.contrast_hi_mm - e.contrast_value_mm], fmt="o-", color="#c2410c", capsize=3)
    for m, s in res[res.outcome != C.ET_ENSEMBLE].groupby("outcome"):
        ax[1].plot(s.month, s.contrast_value_mm, color="#9ca3af", linewidth=0.8, alpha=0.8)
    ax[1].set_ylabel("Thin minus thick ET (mm) at the month's\n90th-percentile deficit vs its median")
    for a in ax:
        a.axhline(0, color="k", linewidth=0.8)
        a.set_xticks(MONTHS, lab)
    fig.suptitle("EXPLORATORY (not pre-registered): T1 thickness terms by month; ensemble ET, thickness error "
                 "corrected, 95% intervals; grey = the six ET models", fontsize=9)
    fig.text(0.01, 0.01, f"{rec['tag']} | cohort cell-months of the timing mask | {len(d)} rows", fontsize=6)
    fig.savefig(f"{out}/t1_by_month.png", dpi=150, bbox_inches="tight")
    write_readme(out, e)
    print(e[["month", "n", "beta_A", "lo_A", "hi_A", "beta_A_deficit", "deficit_q90_mm", "deficit_median_mm",
             "contrast_value_mm", "contrast_lo_mm", "contrast_hi_mm"]].round(3).to_string(index=False))
    print(res.pivot(index="outcome", columns="month", values="contrast_value_mm").round(2).to_string())
    print("written to", out)


if __name__ == "__main__":
    main()
