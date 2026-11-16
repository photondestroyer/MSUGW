"""Task 1.6  Processing test (Stage 2): the test counties and cells, and the test log.

"Run the whole pipeline on the test counties. Apply the blinding mask before any ET is exported. Export only
rainfed ET and cohort ET from non-dry conditions. Fix errors, then run the support extraction for all eligible
Kansas cells." ... "Three to five counties that span thin to thick aquifer and different soils and rainfall.
About 50 to 100 cells per county, drawn with a recorded random seed. ... No feasibility conclusion is drawn from
this sample."
Physical idea: the code is tested where the contrast lives - from thin to thick aquifer - so every join and mask
is exercised on the kinds of cells the analysis will need.
Decisions: D27.
select() (task `sample`) -> sample_counties.csv, sample_cells.csv. In stage2 every ET task reads only these cells.
run()    (task `t1_6`)   -> processing_test_log.md: counts, joins, masks, gate, audit, timings. No conclusions.
"""
import numpy as np


def select(ctx):
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from ..common import io, report
    from ..phase0 import t0_2_cells as T02
    d27 = D.value("D27")
    rs = np.random.default_rng(d27["seed"])
    cells = T02.load(ctx.dep("t0_2"))

    # ══════ STEP 1.6-1  candidate counties: enough cohort cells (D27) ══════
    n_co = cells[cells.role == "cohort"].groupby("county").size()
    cand = n_co[n_co >= d27["min_cohort_cells"]].index.tolist()
    if C.RUN_MODE == "smoke":
        cand = [C.SMOKE_COUNTY]

    # ══════ STEP 1.6-2  county thickness (KGS 2019-2021 sections), rainfall (2016-2024 May-Sep), soils ══════
    st = io.table("kgs_satthick_1921")[["TRS", "satthick_2019_2021"]]
    co = cells[cells.role == "cohort"].merge(st, left_on="kgs_TRS", right_on="TRS", how="left")
    thick = co.groupby("county").satthick_2019_2021.median()
    rain = report.read(ctx.dep("t1_1"), "weather_county_year").assign(county=lambda d: d.county.astype(str)) \
        .groupby("county").p_season_mm.mean()
    soil = cells.groupby("county").rzws_mm.mean()
    tab = pd.DataFrame(dict(cohort_cells=n_co, thickness_ft=thick, rain_mm=rain, rzws_mm=soil)).loc[cand].dropna()

    # ══════ STEP 1.6-3  one county per thickness quartile, the most different in rainfall and soil ══════
    chosen = []
    if len(tab) <= d27["n_counties"]:
        chosen = tab.index.tolist()
    else:
        tab["quartile"] = pd.qcut(tab.thickness_ft, d27["n_counties"], labels=False)
        z = (tab[["rain_mm", "rzws_mm"]] - tab[["rain_mm", "rzws_mm"]].mean()) / tab[["rain_mm", "rzws_mm"]].std()
        for q in range(d27["n_counties"]):
            pool = tab[tab.quartile == q].index
            if not chosen:
                med = tab.loc[pool, "thickness_ft"].median()
                pick = (tab.loc[pool, "thickness_ft"] - med).abs().idxmin()
            else:
                dist = [np.min(np.hypot(*(z.loc[c] - z.loc[chosen]).to_numpy().T)) for c in pool]
                pick = pool[int(np.argmax(dist))]
            chosen.append(pick)
    tab["chosen"] = tab.index.isin(chosen)

    # ══════ STEP 1.6-4  cells: cohort and rainfed per county, drawn with the recorded seed ══════
    rows = []
    for c in chosen:
        for role, k in (("cohort", d27["cohort_per_county"]), ("rainfed", d27["rainfed_per_county"])):
            pool = cells[(cells.county == c) & (cells.role == role)].cell_id.sort_values().to_numpy()
            take = rs.choice(pool, size=min(k, len(pool)), replace=False) if len(pool) else []
            rows += [dict(cell_id=x, county=c, role=role) for x in sorted(take)]
    sample = pd.DataFrame(rows)
    report.table(ctx, tab.reset_index(names="county"), "sample_counties",
                 units=dict(thickness_ft="ft (KGS 2019-2021 median of cohort cells)", rain_mm="mm (May-Sep mean)"))
    report.table(ctx, sample, "sample_cells", notes=f"seed {d27['seed']}")
    ctx.log.info("STEP 1.6-4  counties %s; %d cells", chosen, len(sample))
    ctx.complete(dict(counties=chosen, cells=int(len(sample)), by_role=sample.role.value_counts().to_dict(),
                      seed=d27["seed"]))


def run(ctx):
    import json
    import pandas as pd
    from .. import config as C
    from .. import blinding
    from ..common import report
    from ..phase0 import t0_2_cells as T02

    # ══════ STEP 1.6-5  what the run touched, and the checks of every join and mask ══════
    sample = report.read(ctx.dep("sample"), "sample_cells")
    cells = T02.load(ctx.dep("t0_2"), columns=["cell_id", "county", "role"])
    cov = report.read(ctx.dep("t1_4"), "et_coverage")
    et = report.read(ctx.dep("et_cells"), "et_cell_month")
    dry = report.read(ctx.dep("dry"), "dry_flags")
    checks = []

    def check(name, ok, detail=""):
        checks.append(dict(check=name, passes=bool(ok), detail=str(detail)))

    in_scope = set(sample.cell_id) if C.RUN_MODE == "stage2" else set(cells.cell_id)
    check("coverage rows only for cells in scope", set(cov.cell) <= in_scope, f"{cov.cell.nunique()} cells")
    check("ET rows only for cohort or rainfed cells", et.role.isin(["cohort", "rainfed"]).all())
    check("no cohort ET in dry months (gate)", not ((et.role == "cohort") & et.dry).any())
    check("every ET row has a county and a dry flag", et.county.notna().all() and et.dry.notna().all())
    check("one row per cell-month", not et.duplicated(["cell", "year", "month"]).any())
    vals = et[C.ET_ENSEMBLE].dropna()
    check("ensemble ET within 0-400 mm/month", vals.between(0, 400).all(),
          f"min {vals.min() if len(vals) else 'n/a'}, max {vals.max() if len(vals) else 'n/a'}")
    check("dry flags cover every study county-month", dry.groupby(["year", "month"]).county.nunique().min() > 0)
    ok, problems, n = blinding.audit()
    check("blinding audit of the whole output root", ok, f"{n} ET tables; {problems[:3]}")
    ev = ctx.ledger.read()
    took = {}
    for e in ev:
        if e.get("unit") == "_stage" and e.get("status") == "done":
            took[e["stage"]] = e.get("seconds")

    # ══════ STEP 1.6-6  the test log (no feasibility conclusion is drawn from this sample) ══════
    chk = pd.DataFrame(checks)
    report.table(ctx, chk, "processing_test_checks")
    L = ["# Processing test log (Task 1.6)", "", f"Mode: {C.RUN_MODE}. Counties: {sorted(sample.county.unique())}. "
         f"Cells: {len(sample)} ({sample.role.value_counts().to_dict()}).", "",
         "No feasibility conclusion is drawn from this sample (pilot Section 1).", "", "## Checks", "",
         report.md(chk), "", "## Seconds per task (last completed run)", "",
         "```", json.dumps(took, indent=1), "```"]
    report.text(ctx, "processing_test_log.md", "\n".join(L))
    ctx.complete(dict(checks_passed=int(chk.passes.sum()), checks=int(len(chk)),
                      failed=chk.loc[~chk.passes, "check"].tolist()))
