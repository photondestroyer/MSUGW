"""Check 2.9  Model support.

"For the models in Section 8, check that every term can be estimated; that thin and thick cells cover similar
ranges of deficit and S; how strongly S and the peak-month deficit are correlated; and whether the peak month
changes between years within a county."   Report any contrast that would need extrapolation or that has little
variation left after the fixed effects.
Physical idea: with cell and county-year effects, only variation WITHIN a county-year across thickness, and within
a cell across years, identifies the interactions. If thin and thick cells never meet the same drought, the
contrast is an extrapolation.
Deliverables: c2_9_terms.csv, c2_9_overlap.csv, c2_9_peak_changes.csv, c2_9.md. Uses no ET.
"""
import numpy as np


def run(ctx):
    import pandas as pd
    from .. import config as C
    from .. import definitions as DEF
    from ..common import report, stats
    # ══════ STEP 2.9-1  design rows: cohort peak months in the contrast masks, with A (Check 2.4) ══════
    mk = report.read(ctx.dep("t1_7"), "masks")
    mk = mk[(mk.role == "cohort") & mk.is_peak & (mk.mask_sm | mk.mask_peak)]
    th = report.read(ctx.dep("c2_4"), "thickness_cell_year")[["cell_id", "year", "A_m", "A_poor", "A_third"]]
    w = report.read(ctx.dep("t1_1"), "weather_county_year")
    w["county"] = w["county"].astype(str)
    d = mk.merge(th, on=["cell_id", "year"]).merge(w[["county", "year", "deficit_season_mm"]], on=["county", "year"])
    d = d[~d.A_poor.astype(bool) & np.isfinite(d.A_m)]
    fe = [d.cell_id.to_numpy(), (d.county + "_" + d.year.astype(str)).to_numpy()]

    # ══════ STEP 2.9-2  every term: share of its variance left after cell and county-year effects ══════
    terms = {"A": d.A_m, "A:D_season": d.A_m * d.deficit_season_mm, "A:S_peak": d.A_m * d.S,
             "A:D_peak": d.A_m * d.deficit_peak_mm}
    rows = []
    for name, x in terms.items():
        ok = np.isfinite(x.to_numpy())
        xx = x.to_numpy()[ok]
        left = stats.within(xx, [f[ok] for f in fe])
        rows.append(dict(term=name, rows=int(ok.sum()), variance_share_left=float(np.var(left) / np.var(xx))
                         if np.var(xx) > 0 else np.nan, models=", ".join(k for k, v in DEF.MODELS.items()
                                                                          if name in v["adds"]) or "all"))
    terms_t = pd.DataFrame(rows)

    # ══════ STEP 2.9-3  overlap of thin and thick on deficit and S; extrapolation of the evaluation points ══════
    ov = []
    for var, at in (("deficit_peak_mm", 0.9), ("S", 0.1)):
        x = d[np.isfinite(d[var])]
        ev = float(np.nanquantile(x[var], at))
        for g in ("thin", "thick"):
            s = x[x.A_third == g][var]
            ov.append(dict(variable=var, third=g, n=len(s), p05=float(s.quantile(0.05)) if len(s) else np.nan,
                           p50=float(s.median()) if len(s) else np.nan, p95=float(s.quantile(0.95)) if len(s) else np.nan,
                           evaluation_point=ev, needs_extrapolation=bool(len(s) == 0 or not (s.min() <= ev <= s.max()))))
    ov = pd.DataFrame(ov)
    cyp = d.drop_duplicates(["county", "year"])
    corr = float(np.corrcoef(cyp.S.dropna(), cyp.loc[cyp.S.notna(), "deficit_peak_mm"])[0, 1]) if cyp.S.notna().sum() > 2 \
        else np.nan

    # ══════ STEP 2.9-4  does the peak month change between years within a county? ══════
    pc = w.groupby("county").peak_month.agg(distinct="nunique", modal=lambda s: s.mode().iloc[0]).reset_index()
    pc["share_years_not_modal"] = w.merge(pc[["county", "modal"]], on="county").assign(
        x=lambda q: q.peak_month != q.modal).groupby("county").x.mean().to_numpy()
    # ══════ STEP 2.9-5  write ══════
    report.table(ctx, terms_t, "c2_9_terms", units=dict(variance_share_left="fraction"))
    report.table(ctx, ov, "c2_9_overlap", units=dict(p05="mm or z"))
    report.table(ctx, pc, "c2_9_peak_changes")
    L = ["# Check 2.9 Model support", "", report.md(terms_t), "", report.md(ov), "",
         f"Correlation of S and peak-month deficit across county-years: {corr:.2f}", "",
         f"Counties whose peak month changes between years: {int((pc.distinct > 1).sum())} of {len(pc)}"]
    report.text(ctx, "c2_9.md", "\n".join(L))
    ctx.complete(dict(terms=terms_t.to_dict("records"), corr_S_peak_deficit=corr,
                      needs_extrapolation=ov[ov.needs_extrapolation][["variable", "third"]].to_dict("records"),
                      counties_peak_changes=int((pc.distinct > 1).sum())))
