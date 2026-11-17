"""ET by aquifer thickness: the ET version of Figures 3-5 of Mieno et al. (2024), drawn on SIMULATED ET.

Mieno et al.: Fig. 3 irrigated share against aquifer thickness; Fig. 4 yield against seasonal water deficit at
thickness 10, 40, 70 and 100 m; Fig. 5 the yield difference from the 10 m level with its 95% interval.
Here peak-month ET (mm) takes the place of yield and share.
Blinding (pilot Section 3): cohort ET by aquifer thickness in dry conditions may not be viewed, tabulated or
plotted before the freeze. So this task reads NO ET value. The design is real (cohort peak-month cell-years of
the contrast masks, A of Check 2.4, seasonal deficit of Task 1.1); the outcome is planted_effect(A, D) plus
noise with the variance components of Check 2.6. Each panel shows the fitted curve of one simulated draw, its
95% interval (the width the real figure would have) and the planted effect (dashed).
After the freeze, phase3.main_figs calls draw() with the real peak-month ET in place of the simulated one.
Decisions: D28 (seed).   Deliverables: ET_by_thickness_Fig3/4/5 (.png, _data.csv), fig_thickness.md.
"""
import numpy as np

A_LEVELS = [10.0, 40.0, 70.0, 100.0]        # Mieno Figs 4-5; levels outside the observed range of A are dropped
MASK = "cohort peak-month cell-years of the contrast masks (design only); ET SIMULATED, no ET value read"


def planted_effect(A, D):
    """Hypothetical peak-month ET (mm, up to a constant) at aquifer thickness A (m) and seasonal deficit D (mm).
    It is added to the simulated noise and drawn as the dashed line. Zero = no thickness effect."""
    # TODO(human): the response the project expects (see the Learn by Doing note in the chat)
    return np.zeros(np.broadcast(A, D).shape)


def basis(A, D, a0, sa, d0, sd):
    """Response surface: quadratic in thickness and in deficit, with thickness x deficit interactions.
    shortcut: quadratics, not the splines of Mieno et al.; replace if the real fit shows a threshold they miss."""
    a, d = (np.asarray(A, float) - a0) / sa, (np.asarray(D, float) - d0) / sd
    return np.column_stack([a, a * a, d, d * d, a * d, a * d * d, a * a * d])


def fit(A, D, y, fe=(), cluster=None):
    """Fit the surface to y (mm) with fixed effects swept out and cluster-robust covariance."""
    from ..common import stats
    c = dict(a0=float(np.mean(A)), sa=float(np.std(A)), d0=float(np.mean(D)), sd=float(np.std(D)))
    X = basis(A, D, **c)
    f = stats.fe_ols(y, X, fe=fe, cluster=cluster)
    return dict(beta=f["beta"], cov=f["cov"], xbar=X.mean(axis=0), n_clusters=f["n_clusters"], scale=c)


def predict(f, A, D, ref_A=None):
    """ET (mm) at (A, D) relative to the sample mean, or to thickness ref_A at the same deficit; and its SE."""
    A, D = (v.ravel() for v in np.broadcast_arrays(np.asarray(A, float), np.asarray(D, float)))
    x = basis(A, D, **f["scale"]) - (f["xbar"] if ref_A is None else basis(np.full(A.size, ref_A), D, **f["scale"]))
    return x @ f["beta"], np.sqrt(np.clip(np.einsum("ij,jk,ik->i", x, f["cov"], x), 0, None))


def draw(ctx, dz, yv, x, y, mask, note, simulated, planted=None):
    """Fit the surface to yv on the design dz (A_m, deficit_season_mm, county, year) and write Figs 3-5 with their
    data tables. x, y: cell coordinates (m) for clustering. planted(A, D): the dashed line (simulation only).
    Returns the summary dict."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.colors import LinearSegmentedColormap
    from .. import config as C
    from ..common import report, stats
    A, Dm = dz.A_m.to_numpy(), dz.deficit_season_mm.to_numpy()
    # ══════ STEP F-2  fit the response surface ══════
    # one county: the deficit differs only between years, so year effects would absorb it
    fe = [dz.year.to_numpy()] if dz.county.nunique() > 1 else []
    cl, cl_label = stats.cluster_ids(dz.county, x, y, fallback="block")
    f = fit(A, Dm, yv, fe, cl)
    q = stats.ci95(0.0, 1.0, f["n_clusters"])[1]
    levels = [a for a in A_LEVELS if A.min() <= a <= A.max()]
    if len(levels) < 2:                               # Mieno's levels lie outside this sample: use its own range
        levels = sorted({float(round(v)) for v in np.quantile(A, [0.05, 0.35, 0.65, 0.95])})
    base = levels[0]
    ag = np.linspace(*np.quantile(A, [0.01, 0.99]), 60)
    dg = np.linspace(Dm.min(), Dm.max(), 60)          # inside the observed range: no extrapolation
    p_mean = float(np.mean(planted(A, Dm))) if planted else 0.0
    col = "sim_et_mm" if simulated else "et_mm"

    def curve(a, d, ref=None):
        a, d = np.broadcast_arrays(np.asarray(a, float), np.asarray(d, float))
        m, s = predict(f, a, d, ref)
        out = pd.DataFrame({"A_m": a, "deficit_mm": d, "ref_A_m": np.nan if ref is None else ref, col: m,
                            "se_mm": s, "lo_mm": m - q * s, "hi_mm": m + q * s})
        if planted:
            out["planted_mm"] = planted(a, d) - (p_mean if ref is None else planted(ref, d))
        return out

    def panel(values, xlabel, ylabel):                # main axes with the histogram of the x variable on top
        fig, (h, ax) = plt.subplots(2, 1, figsize=(7, 4.8), sharex=True,
                                    gridspec_kw=dict(height_ratios=[1, 4], hspace=0.04))
        h.hist(values, bins=30, facecolor="none", edgecolor="#166534")
        h.set_axis_off()
        fig.subplots_adjust(bottom=0.2)               # room for the units and sample line report.figure adds
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        return fig, ax

    units = {col: "mm peak-month ET" + (" (simulated)" if simulated else ""), "A_m": "m",
             "deficit_mm": "mm May-Sep"}
    # ══════ STEP F-3  Fig. 3: ET against aquifer thickness at the mean deficit ══════
    c3 = curve(ag, Dm.mean())
    fig, ax = panel(A, "Aquifer thickness (m)", "Peak-month ET relative to the sample mean (mm)")
    ax.fill_between(c3.A_m, c3.lo_mm, c3.hi_mm, color="#9ca3af", alpha=0.6, linewidth=0)
    ax.plot(c3.A_m, c3[col], color="k")
    if planted:
        ax.plot(c3.A_m, c3.planted_mm, color="k", linestyle="--", linewidth=0.8)
    report.figure(ctx, fig, "ET_by_thickness_Fig3", units=units, mask=mask, caption=note)
    # ══════ STEP F-4  Fig. 4: ET against seasonal deficit at fixed thickness levels ══════
    c4 = pd.concat([curve(a, dg) for a in levels])
    cmap = LinearSegmentedColormap.from_list("thin_thick", [C.COLORS["thin"], C.COLORS["middle"], C.COLORS["thick"]])
    fig, ax = panel(Dm, "Seasonal water deficit, May-September (mm)", "Peak-month ET relative to the sample mean (mm)")
    for i, a in enumerate(levels):
        g = c4[c4.A_m == a]
        ax.plot(g.deficit_mm, g[col], color=cmap(i / (len(levels) - 1)), label=f"{a:g} m")
    ax.legend(title="Aquifer thickness", fontsize=8, frameon=False)
    report.figure(ctx, fig, "ET_by_thickness_Fig4", units=units, mask=mask, caption=note)
    # ══════ STEP F-5  Fig. 5: ET difference from the thinnest level, with its 95% interval ══════
    c5 = pd.concat([curve(a, dg, base) for a in levels[1:]])
    fig, axs = plt.subplots(1, len(levels) - 1, figsize=(3.6 * (len(levels) - 1), 3.6), sharey=True, squeeze=False)
    for ax, a in zip(axs[0], levels[1:]):
        g = c5[c5.A_m == a]
        ax.fill_between(g.deficit_mm, g.lo_mm, g.hi_mm, color="#a5b4fc", alpha=0.7, linewidth=0)
        ax.plot(g.deficit_mm, g[col], color="#dc2626")
        if planted:
            ax.plot(g.deficit_mm, g.planted_mm, color="k", linestyle="--", linewidth=0.8)
        ax.axhline(0, color="k", linewidth=0.8)
        ax.set_title(f"{a:g} m minus {base:g} m", fontsize=9)
        ax.set_xlabel("Seasonal water deficit (mm)")
    axs[0][0].set_ylabel("Difference in peak-month ET (mm)")
    fig.subplots_adjust(bottom=0.32)
    report.figure(ctx, fig, "ET_by_thickness_Fig5", units=units, mask=mask, caption=note)
    # ══════ STEP F-6  tables ══════
    for name, c in (("Fig3", c3), ("Fig4", c4), ("Fig5", c5)):
        report.table(ctx, c, f"ET_by_thickness_{name}_data", units=units, mask=mask, notes=note)
    top = c5[c5.A_m == levels[-1]].iloc[-1]           # thickest level against the thinnest, at the largest deficit
    return dict(rows=int(len(dz)), cells=int(dz.cell_id.nunique()), counties=int(dz.county.nunique()),
                clusters=int(f["n_clusters"]), cluster_level=cl_label, year_effects=bool(fe), levels_m=levels,
                A_range_m=[float(A.min()), float(A.max())], deficit_range_mm=[float(Dm.min()), float(Dm.max())],
                fig5_diff_mm_at_max_deficit=float(top[col]), fig5_halfwidth_mm_at_max_deficit=float(q * top.se_mm))


def run(ctx):
    import pandas as pd
    from .. import decisions as DEC
    from ..common import report
    from ..phase0 import t0_2_cells as T02
    from . import c2_6_precision as C26
    # ══════ STEP F-1  the real design and one simulated outcome (no ET value is read) ══════
    w = report.read(ctx.dep("t1_1"), "weather_county_year")
    w["county"] = w["county"].astype(str)
    dz = C26.design(ctx, w)
    vc = report.read(ctx.dep("c2_6"), "c2_6_variance").iloc[0].to_dict()
    xy = T02.load(ctx.dep("t0_2"), columns=["cell_id", "x5070", "y5070"]).set_index("cell_id").loc[dz.cell_id]
    x, y = xy.x5070.to_numpy(), xy.y5070.to_numpy()
    block = pd.factorize(pd.Series((x // 10000).astype(int).astype(str)) + "_" +
                         pd.Series((y // 10000).astype(int).astype(str)) + "_" + dz.year.astype(str))[0]
    rs = np.random.default_rng(DEC.value("D28")["seed"])
    ysim = planted_effect(dz.A_m.to_numpy(), dz.deficit_season_mm.to_numpy()) + \
        C26.noise(rs, vc, block, pd.factorize(dz.cell_id)[0])
    note = "SIMULATED ET on the real design: solid = fitted, band = 95% interval, dashed = planted effect"
    summ = draw(ctx, dz, ysim, x, y, MASK, note, simulated=True, planted=planted_effect)
    report.text(ctx, "fig_thickness.md", "\n".join([
        "# ET by aquifer thickness (ET version of Mieno et al. 2024, Figures 3-5): SIMULATED", "",
        "No ET value was read. The blinding rule forbids viewing cohort ET by aquifer thickness in dry conditions "
        "before the freeze, so the outcome is simulated: the planted effect plus noise with the variance "
        "components of Check 2.6, on the real design. The bands show the precision the real figures would have; "
        "they leave out the thickness error (Check 2.6 includes it).", "",
        report.md(pd.DataFrame([{k: str(v) if isinstance(v, list) else v for k, v in summ.items()}]))]))
    ctx.complete(summ)
