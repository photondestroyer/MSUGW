"""Check 2.6  Precision of the two planned contrasts (simulation).

"Simulate the two planned contrasts on the full sample design. Use the design matrix, the residual variance from
non-dry conditions, spatial and episode dependence, and the thickness error. Add one scenario with higher variance
in drought."   Targets: half-width of the 95% interval at or under 15 mm of peak-month ET; at least three
independent dry episodes. Report how uncertain the projected precision is.
Physical idea: before looking at any drought result, ask whether the design could tell a real thin-versus-thick
difference from noise. The answer depends on how many independent droughts there are, how much ET varies on its
own, and how well A is known - all of which are measurable without the outcome.
Blinding: cohort ET from NON-DRY months only, used only to estimate residual variance; only variances are written.
Decisions: D28; the fitted model has the D30 terms (not listed as a task decision, so a proposed D30 does not block
the check; the terms are in the code below).   Deliverables: c2_6_variance.csv, c2_6_precision.csv, c2_6.md.
"""
import numpy as np


def _variance_components(r, cell, block_year):
    """Residual SD split into a shared 10 km block-year part, a persistent cell part, and the rest."""
    import pandas as pd
    d = pd.DataFrame(dict(r=r, cell=cell, by=block_year))
    v_tot = float(np.var(d.r))
    bm = d.groupby("by").r.transform("mean")
    v_block = float(np.var(bm))
    cm = (d.r - bm).groupby(d.cell).transform("mean")
    v_cell = float(np.var(cm))
    v_rest = max(v_tot - v_block - v_cell, 1e-9)
    return dict(sd_total=np.sqrt(v_tot), sd_block_year=np.sqrt(v_block), sd_cell=np.sqrt(v_cell), sd_rest=np.sqrt(v_rest))


def reliability(Xw, U, fe):
    """Attenuation matrix L of regressors measured with error: OLS on X + U estimates L @ beta, where
    L = (X'X + U'U)^-1 X'X, both swept of the same fixed effects (Xw already is, U is swept here). A covariance
    in true units is L^-1 V L^-T. Regressors that share one error (A, A x deficit, A x S) must be corrected
    together; for one regressor L is the usual reliability var(real) / (var(real) + var(error))."""
    from ..common import stats
    Xw = np.asarray(Xw, float).reshape(len(Xw), -1)
    Uw = stats.within(np.asarray(U, float).reshape(len(Xw), -1), fe)
    sx = Xw.T @ Xw
    return np.linalg.pinv(sx + Uw.T @ Uw) @ sx


def design(ctx, w):
    """The design: cohort peak-month cell-years of the contrast masks with usable A, S and the deficits (no ET).
    w: weather_county_year of Task 1.1 with county as text. Read by this check and by fig_thickness."""
    from ..common import report
    from ..core.errors import PipelineStop
    mk = report.read(ctx.dep("t1_7"), "masks")
    mk = mk[(mk.role == "cohort") & mk.is_peak & (mk.mask_sm | mk.mask_peak)]
    th = report.read(ctx.dep("c2_4"), "thickness_cell_year")[["cell_id", "year", "A_m", "A_se_m", "wte_se_m",
                                                              "base_se_m", "A_poor", "A_third"]]
    dz = mk.merge(th, on=["cell_id", "year"]).merge(w[["county", "year", "deficit_season_mm"]], on=["county", "year"])
    dz = dz[~dz.A_poor.astype(bool) & np.isfinite(dz.A_m)].reset_index(drop=True)
    if len(dz) < 50:
        raise PipelineStop(f"design has only {len(dz)} cohort cell-years with usable A")
    return dz


def noise(rs, vc, block, cellc):
    """One draw of peak-month ET noise (mm) with the measured structure: 10 km block-year, cell, and the rest."""
    return (vc["sd_block_year"] * rs.normal(size=block.max() + 1)[block] +
            vc["sd_cell"] * rs.normal(size=cellc.max() + 1)[cellc] + vc["sd_rest"] * rs.normal(size=len(block)))


def run(ctx):
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from .. import definitions as DEF
    from ..common import report, stats
    from ..phase0 import t0_2_cells as T02
    d28 = dict(D.value("D28"))
    if C.RUN_MODE in ("test", "smoke"):
        d28["n_sim"] = min(d28["n_sim"], 20)            # code tests only; stage2/stage3 use the decided count
    rs = np.random.default_rng(d28["seed"])
    cells = T02.load(ctx.dep("t0_2"), columns=["cell_id", "county", "gmd", "x5070", "y5070", "rzws_mm", "corn_years"])
    xy = cells.set_index("cell_id")[["x5070", "y5070"]]

    # ══════ STEP 2.6-1  residual variance from NON-DRY cohort peak months (in memory; variances only leave) ══════
    et = report.read(ctx.dep("et_cells"), "et_cell_month")
    w = report.read(ctx.dep("t1_1"), "weather_county_year")
    w["county"] = w["county"].astype(str)
    et["county"] = et["county"].astype(str)
    nd = et[(et.role == "cohort") & ~et.dry.astype(bool)].merge(w[["county", "year", "peak_month"]], on=["county", "year"])
    nd = nd[(nd.month == nd.peak_month) & np.isfinite(nd[C.ET_ENSEMBLE])]
    if len(nd) < 50:
        from ..core.errors import PipelineStop
        raise PipelineStop(f"only {len(nd)} non-dry cohort peak-month values: too few to estimate residual variance")
    r = stats.within(nd[C.ET_ENSEMBLE].to_numpy(),                 # residual after cell and county-year effects
                     [nd.cell.to_numpy(), (nd.county + "_" + nd.year.astype(str)).to_numpy()])
    blk = (pd.Series((xy.loc[nd.cell, "x5070"].to_numpy() // 10000).astype(int)).astype(str) + "_" +
           pd.Series((xy.loc[nd.cell, "y5070"].to_numpy() // 10000).astype(int)).astype(str) + "_" +
           nd.year.astype(str).reset_index(drop=True))                  # 10 km block x year
    vc = _variance_components(r, nd.cell.to_numpy(), blk.to_numpy())
    del et, nd
    ctx.log.info("STEP 2.6-1  residual SD components (mm): %s", {k: round(v, 2) for k, v in vc.items()})

    # ══════ STEP 2.6-2  the design: cohort cell-years of the contrast masks, with A, S, deficits (no ET) ══════
    dz = design(ctx, w)
    # thin and thick are the thirds cut ONCE over the whole cohort by Check 2.4 (D26), never cut again here
    by_cell = dz.groupby("cell_id").agg(A=("A_m", "mean"), third=("A_third", "first"))
    A_thin = float(by_cell.A[by_cell.third == "thin"].mean())
    A_thick = float(by_cell.A[by_cell.third == "thick"].mean())
    if not (np.isfinite(A_thin) and np.isfinite(A_thick)):
        from ..core.errors import PipelineStop
        raise PipelineStop("the design holds no thin or no thick cohort cells (D26 thirds of Check 2.4)")
    n_ep = int(dz[dz.dry_season.astype(bool)].drop_duplicates(["county", "year"]).episode_id.nunique())

    # ══════ STEP 2.6-3  simulate: y = noise with the measured structure; fit C1 and P1; CI half-width ══════
    # The fit uses A itself. Thickness error enters through the attenuation matrix of the three A terms, which
    # share it: the water-table error is drawn per cell-year, the base error once per cell (the base does not
    # change between years) or, in the "district" scenario, once per groundwater management district.
    bx = (xy.loc[dz.cell_id, "x5070"].to_numpy() // 10000).astype(int)
    byy = (xy.loc[dz.cell_id, "y5070"].to_numpy() // 10000).astype(int)
    block = pd.factorize(pd.Series(bx.astype(str)) + "_" + pd.Series(byy.astype(str)) + "_" + dz.year.astype(str))[0]
    cellc = pd.factorize(dz.cell_id)[0]
    gmdc = pd.factorize(cells.set_index("cell_id").gmd.loc[dz.cell_id].astype(str))[0]     # no district = its own group
    cy = (dz.county + "_" + dz.year.astype(str)).to_numpy()
    dry = dz.dry.astype(bool).to_numpy()
    sm = dz.mask_sm.to_numpy()
    cl, cl_label = stats.cluster_ids(dz.county, xy.loc[dz.cell_id, "x5070"].to_numpy(),
                                     xy.loc[dz.cell_id, "y5070"].to_numpy(), fallback="block")
    # D30 controls in every model: root-zone storage and corn years 2013-2015 (static cell traits) x each drought
    # term. Measured without error, so they carry zero columns in the attenuation matrix.
    Z = cells.set_index("cell_id").loc[dz.cell_id, ["rzws_mm", "corn_years"]].to_numpy(float)
    out = []
    scenarios =[(1.0, "cell"), (1.0, "district")] + [(f, "cell") for f in d28["dry_variance_factors"]]
    for factor, base_error in scenarios:
        bcode = cellc if base_error == "cell" else gmdc
        for model, var, mask in (("C1", "S", sm), ("P1", "deficit_peak_mm", np.ones(len(dz), bool))):
            hw = []
            sub = np.where(mask & np.isfinite(dz[var].to_numpy()) & np.isfinite(Z).all(axis=1))[0]
            if len(sub) < 50:
                out.append(dict(model=model, dry_variance_factor=factor, base_error=base_error, n=len(sub),
                                note="too few rows"))
                continue
            X0 = dz.iloc[sub]
            q = np.nanquantile(X0[var], [0.1, 0.5, 0.9])
            dX = (q[0] - q[1]) if model == "C1" else (q[2] - q[1])
            A = dz.A_m.to_numpy()
            T = np.column_stack([np.ones(len(dz)), dz.deficit_season_mm.to_numpy(), dz[var].to_numpy()])
            drought = T[:, 1:]                              # seasonal deficit, S or peak deficit
            ctrl = (Z[:, :, None] * drought[:, None, :]).reshape(len(dz), -1)
            Xs = np.column_stack([A[:, None] * T, ctrl])[sub]   # A, A x seas. deficit, A x (S or peak), controls
            zero = np.zeros((len(dz), ctrl.shape[1]))
            fe = [dz.cell_id.to_numpy()[sub], cy[sub]]
            Xw = stats.within(Xs, fe)
            lam = []
            for s in range(d28["n_sim"]):
                e = noise(rs, vc, block, cellc)
                e = np.where(dry, e * np.sqrt(factor), e)
                a_err = (rs.normal(size=bcode.max() + 1)[bcode] * dz.base_se_m.to_numpy()
                         + rs.normal(size=len(dz)) * dz.wte_se_m.to_numpy())
                Li = np.linalg.pinv(reliability(Xw, np.column_stack([a_err[:, None] * T, zero])[sub], fe))
                fit = stats.fe_ols(e[sub], Xs, fe=fe, cluster=cl[sub])
                se_true = float(np.sqrt(max((Li @ fit["cov"] @ Li.T)[2, 2], 0.0)))
                lam.append(fit["se"][2] / se_true if se_true > 0 else np.nan)
                hw.append(1.96 * se_true * abs((A_thin - A_thick) * dX))
            hw, lam = np.array(hw), np.array(lam)
            out.append(dict(model=model, contrast="soil moisture" if model == "C1" else "peak demand",
                            dry_variance_factor=factor, base_error=base_error, n=len(sub),
                            cells=int(X0.cell_id.nunique()),
                            counties=int(X0.county.nunique()), cluster_level=cl_label,
                            clusters=int(len(np.unique(cl[sub]))), halfwidth_median_mm=float(np.median(hw)),
                            halfwidth_p10_mm=float(np.percentile(hw, 10)), halfwidth_p90_mm=float(np.percentile(hw, 90)),
                            target_mm=C.TARGETS["ci_halfwidth_mm"], A_thin_m=A_thin, A_thick_m=A_thick, delta_x=float(dX),
                            reliability_median=float(np.median(lam))))
            ctx.log.info("STEP 2.6-3  %s factor %.1f base error per %s: median half-width %.1f mm", model, factor,
                         base_error, np.median(hw))
            ctx.check_budget()
    pr = pd.DataFrame(out)
    # ══════ STEP 2.6-4  write (variances and half-widths only) ══════
    report.table(ctx, pd.DataFrame([vc]), "c2_6_variance", units={"sd_*": "mm (peak-month ET)"},
                 mask="cohort, non-dry peak months (variance only)")
    report.table(ctx, pr, "c2_6_precision", units={"halfwidth_*": "mm peak-month ET"},
                 mask="cohort cell-years in the contrast masks (design only, no outcome)")
    L = ["# Check 2.6 Precision (simulated)", "", f"Independent dry episodes in the design: {n_ep} "
         f"(minimum {C.TARGETS['min_dry_episodes']}).", "", report.md(pr), "",
         "Half-widths are 1.96 x SE of the interaction times |A_thin - A_thick| x (evaluation point - median). "
         "The SE is corrected for thickness error jointly over the three thickness terms, which share that "
         "error (covariance L^-1 V L^-T, L the attenuation matrix after the fixed effects); reliability_median "
         "is the uncorrected SE divided by the corrected one. Thickness error: water table per cell-year; base "
         "once per cell (base_error = cell, the main case) or once per groundwater management district "
         "(base_error = district). p10-p90 across simulations show how uncertain the projection is."]
    report.text(ctx, "c2_6.md", "\n".join(L))
    ctx.complete(dict(dry_episodes=n_ep, variance=vc, precision=pr.to_dict("records")))
