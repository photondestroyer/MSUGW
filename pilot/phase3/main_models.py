"""Main analysis: models C1, P1, J1, T1 (D30) and the two contrasts (D31). Runs only after the freeze.

Pilot Section 8: "Terms in every model. Peak-month ET in mm for each cell and year; cell effects; county-by-year
effects; thickness A as a continuous variable; A x seasonal deficit; and interactions of root-zone water storage
and baseline crop history with the drought terms. C1 adds A x S in the peak month. P1 adds A x peak-month
deficit. J1 contains both added terms. T1 uses monthly ET for May-September ..."
Contrasts (D31): soil-moisture = (thin - thick) peak-month ET at S = 10th percentile minus the same at median S
(from C1); peak-demand = the same at the 90th percentile of the peak-month deficit against its median (from P1).
In these models the contrast is beta[A x term] * (A_thin - A_thick) * (x_eval - x_median).
Thickness error (D23, D24): A is measured with error shared by every A term; the estimates are corrected with the
expected error Gram matrix E[U'U] (Monte Carlo over the error model of Check 2.4), as Check 2.6 projected.
Estimand (D32): every cohort cell in every year (intent to treat); the same models on cell-years known to be
irrigated are the secondary result.
Decisions: D28 (seed), D29, D30, D31, D32.   Deliverables: main_estimates.csv, main_contrasts.csv,
main_support.csv, main_results.md.
"""
import numpy as np

N_ERR_DRAWS = 50          # Monte Carlo draws of the thickness error for E[U'U]; the correction is a mean, not a test
SAMPLES = (("C1", "sm", ["S"]), ("P1", "peak", ["deficit_peak_mm"]), ("C1", "common", ["S"]),
           ("P1", "common", ["deficit_peak_mm"]), ("J1", "common", ["S", "deficit_peak_mm"]))
TERM = {"S": "A:S", "deficit_peak_mm": "A:D_peak"}
CONTRAST = {"S": ("soil moisture", 0.10), "deficit_peak_mm": ("peak demand", 0.90)}
CONTROLS = ("rzws_mm", "corn_years")       # D30: root-zone storage and corn years 2013-2015, x each drought term


def design(d, drought, exposure="deficit_season_mm"):
    """X = [A x (1, exposure, drought...), controls x (exposure, drought...)], T = (1, exposure, drought...).
    Returns X, T (the factor multiplying A, for the error), and the column names."""
    n = len(d)
    T = np.column_stack([np.ones(n), d[exposure].to_numpy(float)] + [d[c].to_numpy(float) for c in drought])
    A = d["A_m"].to_numpy(float)
    Z = d[list(CONTROLS)].to_numpy(float)
    ctrl = (Z[:, :, None] * T[:, None, 1:]).reshape(n, -1)
    lab = ["D_season" if exposure == "deficit_season_mm" else "exposure"] + [TERM.get(c, c)[2:] for c in drought]
    names = ["A"] + [f"A:{x}" for x in lab] + [f"{z}:{x}" for z in ("rzws", "corn_years") for x in lab]
    return np.column_stack([A[:, None] * T, ctrl]), T, names


def expected_error_gram(T, fe, wte_se, base_se, base_group, k_ctrl, rs, n_draws=N_ERR_DRAWS):
    """E[Uw'Uw]: U = a_err x T for the A columns (zeros for the controls), a_err = base error (one draw per
    base_group: cell, or district) + water-table error (per row); swept of the fixed effects like X."""
    from ..common import stats
    p = T.shape[1] + k_ctrl
    G = np.zeros((p, p))
    for _ in range(n_draws):
        a = rs.normal(size=base_group.max() + 1)[base_group] * base_se + rs.normal(size=len(T)) * wte_se
        U = np.column_stack([a[:, None] * T, np.zeros((len(T), k_ctrl))])
        Uw = stats.within(U, fe)
        G += Uw.T @ Uw
    return G / n_draws


def fit(y, X, fe, cluster, G=None):
    """fe_ols, then (G given) the errors-in-variables correction: L = (Xw'Xw)^-1 (Xw'Xw - G),
    beta_c = L^-1 beta, cov_c = L^-1 cov L^-T. Rows must be finite (the caller filters)."""
    from ..common import stats
    f = stats.fe_ols(y, X, fe=fe, cluster=cluster)
    out = dict(beta=f["beta"], cov=f["cov"], n_clusters=f["n_clusters"], beta_c=f["beta"], cov_c=f["cov"])
    if G is not None:
        Xw = stats.within(X, fe)
        S = Xw.T @ Xw
        Li = np.linalg.pinv(np.linalg.pinv(S) @ (S - G))
        out.update(beta_c=Li @ f["beta"], cov_c=Li @ f["cov"] @ Li.T)
    return out


def contrast(f, j, a_thin, a_thick, dx, corrected=True):
    """D31: beta[j] * (A_thin - A_thick) * dx, its SE and 95% interval (t with G-1 df below 50 clusters)."""
    from ..common import stats
    b, V = (f["beta_c"], f["cov_c"]) if corrected else (f["beta"], f["cov"])
    k = (a_thin - a_thick) * dx
    v, se = float(b[j] * k), float(np.sqrt(max(V[j, j], 0.0)) * abs(k))
    lo, hi = stats.ci95(v, se, f["n_clusters"])
    return dict(value_mm=v, se_mm=se, lo_mm=float(lo), hi_mm=float(hi))


def thirds_mean_A(d):
    m = d.groupby("cell_id").agg(A=("A_m", "mean"), third=("A_third", "first"))
    return float(m.A[m.third == "thin"].mean()), float(m.A[m.third == "thick"].mean())


def estimate(d, outcome, model, sample, drought, rs, base_error="cell"):
    """One model on one sample and outcome: estimates (corrected and not) and the contrast(s)."""
    import pandas as pd
    from ..common import stats
    need = ["A_m", "deficit_season_mm", outcome, *drought, *CONTROLS, "wte_se_m", "base_se_m"]
    d = d[np.isfinite(d[need].to_numpy(float)).all(axis=1)].reset_index(drop=True)
    if len(d) < 50 or d.county.nunique() < 2:
        return [], [dict(model=model, sample=sample, outcome=outcome, base_error=base_error, n=len(d),
                         note="too few rows")]
    X, T, names = design(d, drought)
    fe = [d.cell_id.to_numpy(), (d.county.astype(str) + "_" + d.year.astype(str)).to_numpy()]
    cl, cl_label = stats.cluster_ids(d.county, d.x5070.to_numpy(), d.y5070.to_numpy(), fallback="block")
    grp = pd.factorize(d.cell_id if base_error == "cell" else d.gmd.fillna(0).astype(str))[0]
    G = expected_error_gram(T, fe, d.wte_se_m.to_numpy(), d.base_se_m.to_numpy(), grp, X.shape[1] - T.shape[1], rs)
    f = fit(d[outcome].to_numpy(float), X, fe, cl, G)
    meta = dict(model=model, sample=sample, outcome=outcome, base_error=base_error, n=len(d),
                cells=int(d.cell_id.nunique()), counties=int(d.county.nunique()), clusters=int(f["n_clusters"]),
                cluster_level=cl_label)
    est = [dict(meta, term=t, beta=float(f["beta"][i]), se=float(np.sqrt(max(f["cov"][i, i], 0))),
                beta_corrected=float(f["beta_c"][i]), se_corrected=float(np.sqrt(max(f["cov_c"][i, i], 0))))
           for i, t in enumerate(names)]
    a_thin, a_thick = thirds_mean_A(d)
    con = []
    for c in drought:
        name, q = CONTRAST[c]
        x = d[c].to_numpy(float)
        ev, med = np.quantile(x, [q, 0.5])
        j = names.index(TERM[c])
        for corr in (True, False):
            con.append(dict(meta, contrast=name, corrected=corr, eval_point=float(ev), median=float(med),
                            A_thin_m=a_thin, A_thick_m=a_thick, **contrast(f, j, a_thin, a_thick, ev - med, corr)))
    return est, con


def support(d, drought):
    """D31 support: is each evaluation point inside the 5-95% range of both the thin and the thick third?"""
    rows = []
    for c in drought:
        _, q = CONTRAST[c]
        x = d[c].to_numpy(float)
        ev = float(np.nanquantile(x, q))
        for g in ("thin", "thick"):
            s = d.loc[d.A_third == g, c].dropna()
            lo, hi = (float(s.quantile(0.05)), float(s.quantile(0.95))) if len(s) else (np.nan, np.nan)
            rows.append(dict(variable=c, third=g, eval_point=ev, p05=lo, p95=hi, inside=bool(lo <= ev <= hi)))
    return rows


def timing(dm, outcome, rs):
    """T1 (pilot Section 8): monthly ET May-September; cell x month and county x month x year effects; A x monthly
    exposure (monthly deficit, D30), A x peak-month indicator, A x exposure x peak indicator; controls x exposure."""
    import pandas as pd
    from ..common import stats
    need = ["A_m", "deficit_mm", outcome, *CONTROLS, "wte_se_m", "base_se_m"]
    d = dm[np.isfinite(dm[need].to_numpy(float)).all(axis=1)].reset_index(drop=True)
    if len(d) < 50:
        return []
    pk = d.is_peak.astype(float).to_numpy()
    ex = d.deficit_mm.to_numpy(float)
    T = np.column_stack([np.ones(len(d)), ex, pk, ex * pk])
    Z = d[list(CONTROLS)].to_numpy(float)
    X = np.column_stack([d.A_m.to_numpy(float)[:, None] * T, Z * ex[:, None]])
    names = ["A", "A:exposure", "A:peak", "A:exposure:peak", "rzws:exposure", "corn_years:exposure"]
    m = d.month.astype(str)
    fe = [(d.cell_id + "_" + m).to_numpy(), (d.county.astype(str) + "_" + m + "_" + d.year.astype(str)).to_numpy()]
    cl, cl_label = stats.cluster_ids(d.county, d.x5070.to_numpy(), d.y5070.to_numpy(), fallback="block")
    G = expected_error_gram(T, fe, d.wte_se_m.to_numpy(), d.base_se_m.to_numpy(), pd.factorize(d.cell_id)[0],
                            Z.shape[1], rs)
    f = fit(d[outcome].to_numpy(float), X, fe, cl, G)
    out = []
    for i, t in enumerate(names):
        b, s = float(f["beta_c"][i]), float(np.sqrt(max(f["cov_c"][i, i], 0)))
        lo, hi = stats.ci95(b, s, f["n_clusters"])
        out.append(dict(model="T1", sample="timing", outcome=outcome, base_error="cell", n=len(d),
                        cells=int(d.cell_id.nunique()), counties=int(d.county.nunique()),
                        clusters=int(f["n_clusters"]), cluster_level=cl_label, term=t, beta=float(f["beta"][i]),
                        se=float(np.sqrt(max(f["cov"][i, i], 0))), beta_corrected=b, se_corrected=s,
                        lo_corrected=float(lo), hi_corrected=float(hi)))
    return out


def assemble(ctx):
    """Peak-month (and monthly, for T1) cohort rows: ET (unblinded in main mode), masks, A, deficits, controls,
    and the D32 'known irrigated' flag. Poor-thickness cell-years (D25) are left out (D30)."""
    import pandas as pd
    from .. import config as C
    from .. import definitions as DEF
    from ..common import report
    from ..phase0 import t0_2_cells as T02
    lan = [f"lanid_{y}" for y in C.PANEL_YEARS if y <= 2020]
    cells = T02.load(ctx.dep("t0_2"), columns=["cell_id", "gmd", "x5070", "y5070", *CONTROLS, *lan])
    et = report.read(ctx.dep("et_cells"), "et_cell_month")
    et = et[et.role == "cohort"].rename(columns={"cell": "cell_id"})
    ets = [C.ET_ENSEMBLE] + [f"et_{m}" for m in C.ET_MODELS]
    mk = report.read(ctx.dep("t1_7"), "masks")
    mk = mk[mk.role == "cohort"]
    th = report.read(ctx.dep("c2_4"), "thickness_cell_year")[["cell_id", "year", "A_m", "wte_se_m", "base_se_m",
                                                              "A_poor", "A_third"]]
    wy = report.read(ctx.dep("t1_1"), "weather_county_year")[["county", "year", "deficit_season_mm"]]
    wm = report.read(ctx.dep("t1_1"), "weather_county_month")[["county", "year", "month", "deficit_mm", "eto_mm"]]
    for t in (et, mk, wy, wm):
        t["county"] = t["county"].astype(str)
    d = (mk.merge(et[["cell_id", "year", "month", *ets]], on=["cell_id", "year", "month"], how="left")
           .merge(th, on=["cell_id", "year"], how="left")
           .merge(wy, on=["county", "year"], how="left").merge(wm, on=["county", "year", "month"], how="left")
           .merge(cells.drop(columns=lan), on="cell_id", how="left"))
    d = d[~d.A_poor.astype("boolean").fillna(True).astype(bool) & np.isfinite(d.A_m)].reset_index(drop=True)
    if C.RUN_MODE != "main":              # code test (smoke/test): cohort ET replaced by noise, never a real value
        rs = np.random.default_rng(0)
        for e in ets:
            d[e] = np.where(np.isfinite(d[e]), 100.0 + rs.normal(0, 20, len(d)), np.nan)
    # D32 secondary: cell-years known to be irrigated (LANID 2016-2020; Check 2.8 clear links 2021-2023)
    irr = cells.melt(id_vars="cell_id", value_vars=lan, var_name="y", value_name="share")
    irr["year"] = irr.y.str[-4:].astype(int)
    irr = irr[DEF.is_irrigated_year(irr.share.to_numpy())][["cell_id", "year"]]
    st = report.read(ctx.dep("c2_8"), "c2_8_status")
    irr = pd.concat([irr, st[st.status == "irrigated"][["cell_id", "year"]]]).drop_duplicates()
    d = d.merge(irr.assign(known_irrigated=True), on=["cell_id", "year"], how="left")
    d["known_irrigated"] = d.known_irrigated.astype("boolean").fillna(False).astype(bool)
    return d, ets


def run(ctx):
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from .. import freeze
    from ..common import report
    from ..core.errors import PipelineStop
    if C.RUN_MODE not in ("main", "smoke", "test"):
        raise PipelineStop(f"the main analysis runs in main mode only (after the freeze), not {C.RUN_MODE}")
    rec = freeze.check() if C.RUN_MODE == "main" else dict(tag=f"{C.RUN_MODE} CODE TEST: cohort ET SIMULATED")
    rs = np.random.default_rng(D.value("D28")["seed"])

    # ══════ STEP M-1  cohort rows: ET, masks, A, deficits, controls ══════
    d, ets = assemble(ctx)
    pk = d[d.is_peak.astype(bool)]
    ctx.log.info("STEP M-1  %d peak-month cohort cell-years (%d cells); freeze %s", len(pk), pk.cell_id.nunique(),
                 rec["tag"])

    # ══════ STEP M-2  C1, P1, J1 per sample and ET outcome; the contrasts (D31) ══════
    est, con = [], []
    mask = dict(sm=pk.mask_sm, peak=pk.mask_peak, common=pk.mask_common)
    for model, sample, drought in SAMPLES:
        s = pk[mask[sample].astype(bool)]
        for outcome in ets:
            for be in (("cell", "district") if outcome == C.ET_ENSEMBLE and sample != "common" else ("cell",)):
                e, c = estimate(s, outcome, model, sample, drought, rs, be)
                est += e
                con += c
            ctx.check_budget()
        ctx.log.info("STEP M-2  %s on %s sample done", model, sample)
    # D32 secondary: cell-years known to be irrigated, ensemble only
    for model, sample, drought in SAMPLES[:2]:
        s = pk[mask[sample].astype(bool) & pk.known_irrigated]
        e, c = estimate(s, C.ET_ENSEMBLE, model, f"{sample}, known irrigated", drought, rs)
        est += e
        con += c

    # ══════ STEP M-3  T1 timing model (monthly) ══════
    tm = d[d.mask_timing.astype(bool)]
    est += timing(tm, C.ET_ENSEMBLE, rs)
    ctx.log.info("STEP M-3  T1 done on %d cell-months", len(tm))

    # ══════ STEP M-4  support of the evaluation points; write ══════
    sup = pd.DataFrame(support(pk[pk.mask_peak.astype(bool)], ["deficit_peak_mm"]) +
                       support(pk[pk.mask_sm.astype(bool)], ["S"]))
    est, con = pd.DataFrame(est), pd.DataFrame(con)
    mask_txt = "cohort cells, all years (D32 intent to treat), poor-thickness cell-years left out (D25)"
    report.table(ctx, est, "main_estimates", units={"beta*": "mm per (m x term unit)"}, mask=mask_txt)
    report.table(ctx, con, "main_contrasts", units={"value_mm": "mm peak-month ET", "*_mm": "mm"}, mask=mask_txt)
    report.table(ctx, sup, "main_support", mask=mask_txt)
    main = con[(con.outcome == C.ET_ENSEMBLE) & con.corrected.astype("boolean").fillna(False).astype(bool)
               & (con.base_error == "cell") & con["sample"].isin(["sm", "peak"])] if "corrected" in con else con.iloc[0:0]
    report.text(ctx, "main_results.md", "\n".join([
        "# Main analysis: the two contrasts", "", f"Freeze: {rec['tag']}. Expected sign (D31): negative. Evidence: a "
        "95% interval that excludes zero. Ensemble ET, thickness error corrected, primary samples:", "",
        report.md(main[["contrast", "model", "sample", "value_mm", "lo_mm", "hi_mm", "n", "cells", "clusters"]])
        if len(main) else "(no estimate)", "", "Every ET model, the common sample (J1), the district base error, the "
        "uncorrected estimates and the known-irrigated subset are in main_contrasts.csv; T1 in main_estimates.csv. "
        "The range across ET models is not a confidence interval.", "", "## Support of the evaluation points", "",
        report.md(sup)]))
    ctx.complete(dict(freeze=rec["tag"], rows=int(len(pk)),
                      primary=main[["contrast", "value_mm", "lo_mm", "hi_mm"]].to_dict("records") if len(main) else []))
