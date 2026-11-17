"""Main analysis figures (after the freeze only).

Figure 1 (lead, pilot Section 10): "Peak-month ET, in mm and as a ratio to reference ET, against the regional
soil-moisture state, for rainfed cells and for cohort cells over thin and thick aquifer."
Figure 3: "The two contrasts with 95% intervals ... each ET model, with the timing and joint-model results."
ET by thickness: the ET version of Mieno et al. (2024) Figures 3-5 (phase2.fig_thickness.draw) on real ET.
Bins of S: deciles over all plotted rows; a bin's value is the mean of its county means, its interval 1.96 SE
across counties (counties are the clusters of the models, D31).
Decisions: D29, D31, D32.   Deliverables: Fig1_ET_vs_S, Fig3_contrasts, ET_by_thickness_Fig3/4/5 (+ data).
"""
import numpy as np


def binned(d, value, groups, bins):
    """Mean of county means of `value` per S bin and group, and 1.96 x SE across counties."""
    import pandas as pd
    d = d.assign(bin=pd.cut(d.S, bins, labels=False, include_lowest=True))
    cm = d.groupby(["group", "bin", "county"], observed=True)[value].mean().reset_index()
    out = cm.groupby(["group", "bin"], observed=True)[value].agg(["mean", "std", "count"]).reset_index()
    out["half95"] = 1.96 * out["std"].fillna(0) / np.sqrt(out["count"].clip(lower=1))
    mid = (bins[:-1] + bins[1:]) / 2
    out["S_mid"] = mid[out["bin"].astype(int)]
    return out[out.group.isin(groups)]


def run(ctx):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    from .. import config as C
    from .. import freeze
    from ..common import report
    from ..core.errors import PipelineStop
    from ..phase2 import fig_thickness as F
    from . import main_models as MM
    if C.RUN_MODE not in ("main", "smoke", "test"):
        raise PipelineStop(f"the main figures are drawn in main mode only (after the freeze), not {C.RUN_MODE}")
    rec = freeze.check() if C.RUN_MODE == "main" else dict(tag=f"{C.RUN_MODE} CODE TEST: cohort ET SIMULATED")
    E = C.ET_ENSEMBLE
    d, ets = MM.assemble(ctx)
    pk = d[d.is_peak.astype(bool) & (d.mask_sm | d.mask_peak).astype(bool)].copy()

    # ══════ STEP MF-1  Figure 1: peak-month ET against S; rainfed, thin, thick ══════
    et = report.read(ctx.dep("et_cells"), "et_cell_month")
    mk = report.read(ctx.dep("t1_7"), "masks")
    wm = report.read(ctx.dep("t1_1"), "weather_county_month")[["county", "year", "month", "eto_mm"]]
    for t in (et, mk, wm):
        t["county"] = t["county"].astype(str)
    rf = (mk[(mk.role == "rainfed") & mk.is_peak.astype(bool)]
          .merge(et.rename(columns={"cell": "cell_id"})[["cell_id", "year", "month", E]], on=["cell_id", "year", "month"])
          .merge(wm, on=["county", "year", "month"]))
    rf["group"] = "rainfed"
    co = pk[pk.mask_sm.astype(bool)].assign(group=pk.A_third)
    f1 = pd.concat([rf[["group", "county", "S", E, "eto_mm"]], co[["group", "county", "S", E, "eto_mm"]]])
    f1 = f1[np.isfinite(f1.S) & np.isfinite(f1[E])]
    f1["ratio"] = f1[E] / f1.eto_mm
    bins = np.unique(np.quantile(f1.S, np.linspace(0, 1, 11)))
    groups = ["rainfed", "thin", "thick"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    tabs = []
    for ax, (v, lab) in zip(axs, ((E, "Peak-month ET (mm)"), ("ratio", "Peak-month ET / reference ET"))):
        b = binned(f1, v, groups, bins)
        tabs.append(b.assign(value=v))
        for g in groups:
            s = b[b.group == g]
            ax.errorbar(s.S_mid, s["mean"], yerr=s.half95, color=C.COLORS[g], marker="o", ms=3, capsize=2, label=g)
        ax.set_xlabel("Regional soil-moisture state S (SD from the monthly mean)")
        ax.set_ylabel(lab)
    axs[0].legend(frameon=False, fontsize=8)
    fig.subplots_adjust(bottom=0.2)
    mask = "peak months; rainfed cells, and cohort cells of the soil-moisture mask (thin and thick thirds, D26)"
    report.figure(ctx, fig, "Fig1_ET_vs_S", units={"ET": "mm", "S": "SD"}, mask=mask,
                  caption=f"bin = decile of S; point = mean of county means; bar = 95% across counties; {rec['tag']}")
    report.table(ctx, pd.concat(tabs), "Fig1_data", units={"mean": "mm or ratio"}, mask=mask)

    # ══════ STEP MF-2  Figure 3: the two contrasts for each ET model ══════
    con = report.read(ctx.dep("m_models"), "main_contrasts")
    if "corrected" not in con:                        # no model could be estimated (e.g. a one-county code test)
        con = con.assign(corrected=False, value_mm=np.nan, lo_mm=np.nan, hi_mm=np.nan, contrast="")
    con = con[con.corrected.astype("boolean").fillna(False).astype(bool) & (con.base_error == "cell")
              & con.value_mm.notna()]
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, (name, prim) in zip(axs, (("soil moisture", "sm"), ("peak demand", "peak"))):
        for k, (smp, mod, mk_) in enumerate(((prim, None, "o"), ("common", "J1", "s"))):
            s = con[(con.contrast == name) & (con["sample"] == smp) & ((con.model == mod) if mod else True)]
            s = s.set_index("outcome").reindex(ets)
            yy = np.arange(len(ets)) + (0.18 if k else -0.18)
            ax.errorbar(s.value_mm, yy, xerr=[s.value_mm - s.lo_mm, s.hi_mm - s.value_mm], fmt=mk_, ms=4,
                        color="k" if k == 0 else "#6b7280", mfc="k" if k == 0 else "none", capsize=2,
                        label="primary sample" if k == 0 else "common sample, J1")
        ax.axvline(0, color="k", linewidth=0.8)
        ax.set_title(f"{name} contrast", fontsize=10)
        ax.set_xlabel("Thin minus thick, dry minus median (mm peak-month ET)")
    axs[0].set_yticks(range(len(ets)), [e.replace("et_", "") for e in ets])
    axs[0].legend(frameon=False, fontsize=8)
    fig.subplots_adjust(bottom=0.2, left=0.15)
    report.figure(ctx, fig, "Fig3_contrasts", units={"contrast": "mm"}, mask="D31 samples; thickness error corrected",
                  caption="95% intervals; the range across ET models is not a confidence interval")

    # ══════ STEP MF-3  ET by aquifer thickness: Mieno Figs 3-5 on real ET ══════
    dz = pk[np.isfinite(pk[E]) & np.isfinite(pk.deficit_season_mm)].reset_index(drop=True)
    note = f"REAL peak-month ET (ensemble): solid = fitted surface, band = 95% interval; {rec['tag']}"
    summ = F.draw(ctx, dz, dz[E].to_numpy(), dz.x5070.to_numpy(), dz.y5070.to_numpy(),
                  "cohort peak-month cell-years of the contrast masks (D32 intent to treat)", note, simulated=False)
    ctx.complete(dict(freeze=rec["tag"], fig_thickness=summ))
