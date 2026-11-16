"""Feasibility memo (pilot Section 9, deliverable 4).

"Feasibility memo that answers four questions: Can each contrast be estimated with useful precision? Is the
thickness estimate reliable? Is well capacity available for a subset of cells? Are the satellite inputs
supported? For each answer give the observed value, the target, the projected precision and what remains open."
The memo only collects what the checks wrote; it draws no conclusion the checks did not, and it never treats a
missed target as a stop (pilot rule 4: targets are planning targets).
"""
import json


def _s(ctx, task):
    with open(f"{ctx.dep(task)}/_complete.json") as f:
        return json.load(f)["summary"]


def run(ctx):
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from ..common import report
    s = {t: _s(ctx, t) for t in ("c2_1", "c2_2", "c2_3", "c2_4", "c2_5", "c2_6", "c2_7", "c2_8", "c2_9", "t1_7")}
    T = C.TARGETS
    prec = pd.DataFrame(s["c2_6"]["precision"])
    base = prec[prec.dry_variance_factor == 1.0] if "dry_variance_factor" in prec else prec
    tg4 = {r["target"]: r["observed"] for r in s["c2_4"]["targets"]}
    rows = [
        ("1 Precision", "soil-moisture contrast: median 95% half-width (mm)",
         base[base.model == "C1"].halfwidth_median_mm.iloc[0] if (base.model == "C1").any() else None, T["ci_halfwidth_mm"]),
        ("1 Precision", "peak-demand contrast: median 95% half-width (mm)",
         base[base.model == "P1"].halfwidth_median_mm.iloc[0] if (base.model == "P1").any() else None, T["ci_halfwidth_mm"]),
        ("1 Precision", "independent dry episodes in the design", s["c2_6"]["dry_episodes"], T["min_dry_episodes"]),
        ("2 Thickness", "cohort cells within 8 km of a well, worst year",
         tg4.get("cohort cells within 8 km of a measured well, worst year"), T["well_share"]),
        ("2 Thickness", "counties with >= 50 cells in both thin and thick thirds",
         tg4.get("counties with >= 50 cells in both thin and thick thirds"), T["both_thirds_counties"]),
        ("2 Thickness", "thin-to-thick difference / SE", tg4.get("thin-to-thick difference / its SE (cell-level errors)"),
         T["thin_thick_se_ratio"]),
        ("2 Thickness", "water-table RMSE, locations held out (m)", s["c2_4"]["rmse_location_m"], None),
        ("3 Capacity", "cells with a clear link", s["c2_7"]["share_clear_link"], None),
        ("3 Capacity", "cells with a usable rate (driller yield)", s["c2_7"]["share_usable_rate"], None),
        ("3 Capacity", "cohort cell-years 2021-2024 with known irrigation status", s["c2_8"]["share_known"], None),
        ("4 Satellite", "dry county-years kept by the soil-moisture mask", s["c2_1"]["share_dry_county_years_kept"],
         T["smap_dry_share"]),
        ("4 Satellite", "counties with >= 3 usable reference footprints", s["c2_2"]["share_counties_meeting"],
         T["fp_county_share"]),
        ("4 Satellite", "peak-month cohort cells with >= 4 ET models", s["c2_3"]["peak_share_supported"], T["et_cell_share"]),
    ]
    tab = pd.DataFrame(rows, columns=["question", "item", "observed", "target"])

    def meets(r):
        obs, tgt = r["observed"], r["target"]
        if tgt is None or obs is None or pd.isna(tgt) or pd.isna(obs):
            return ""
        lower_is_better = "half-width" in r["item"]
        return "yes" if (obs <= tgt if lower_is_better else obs >= tgt) else "below target"
    tab["meets"] = tab.apply(meets, axis=1)
    report.table(ctx, tab, "feasibility_table")
    open_items = [
        f"Reference footprints use a {s['c2_2'].get('radius_km')} km radius (D13); the SMAP file used is named in "
        "Task 1.2 (smap_source)",
        "Authorised rate, pumping hours and verified capacity: wimas_lookup_sheet.csv (Check 2.7)",
        "Footprint check against 2018-2020 labels (LANID): see check_2018_2020 in the Task 1.2 summary and "
        "the irr_max_2018_2020 column of footprints.csv",
        f"Decisions still proposed: {[k for k, e in D.entries().items() if e['status'] != 'fixed']}",
    ]
    L = ["# Feasibility memo", "", f"Mode: {C.RUN_MODE}. Targets are planning targets: a value below target is "
         "reported and discussed and does not stop the work (pilot Section 3).", "", report.md(tab), "",
         "## Projected precision across scenarios", "", report.md(prec) if len(prec) else "", "",
         "## Open", ""] + [f"- {x}" for x in open_items] + [
        "", "## Next step (pilot Section 11)", "",
        "Read the table against Section 11 of the pilot: which row of that table describes these findings decides "
        "the next step, with the supervisor."]
    report.text(ctx, "feasibility_memo.md", "\n".join(L))
    ctx.complete(dict(rows=tab.to_dict("records")))
