"""Check 2.3  ET coverage.

"Check peak months and the other months separately, June included."   Target: at least four models in at least
80% of cohort cells. Mask unsupported cell-months. A missing non-peak month does not remove a valid peak-month
observation.
Physical idea: the peak month carries the test; June matters because irrigation starts there and clouds are
common in early summer.
Deliverables: c2_3_by_month.csv, c2_3_peak_vs_other.csv, c2_3.md.
"""


def run(ctx):
    import pandas as pd
    from .. import config as C
    from ..common import report
    from ..phase0 import t0_2_cells as T02
    # ══════ STEP 2.3-1  cohort cell-months with their peak-month flag ══════
    cov = report.read(ctx.dep("t1_4"), "et_coverage")
    cov = cov[cov.role == "cohort"].copy()
    w = report.read(ctx.dep("t1_1"), "weather_county_year")
    w["county"] = w["county"].astype(str)
    cov["county"] = cov["county"].astype(str)
    cov = cov.merge(w[["county", "year", "peak_month"]], on=["county", "year"], how="left")
    cov["is_peak"] = cov.month == cov.peak_month
    # ══════ STEP 2.3-2  share with >= 4 models: by month (June included), peak vs other ══════
    by_month = cov.groupby("month").agg(cell_months=("cell", "size"), share_supported=("supported", "mean"),
                                        **{f"share_{m}": (f"ok_{m}", "mean") for m in C.ET_MODELS}).reset_index()
    pv = cov.groupby("is_peak").agg(cell_months=("cell", "size"), share_supported=("supported", "mean")).reset_index()
    pv["months"] = pv.is_peak.map({True: "peak month", False: "other season months"})
    by_year_peak = cov[cov.is_peak].groupby("year").supported.mean().rename("share_supported").reset_index()
    # ══════ STEP 2.3-3  write ══════
    report.table(ctx, by_month, "c2_3_by_month", units={"share_*": "fraction of cohort cell-months"})
    report.table(ctx, pv, "c2_3_peak_vs_other", units={"share_supported": "fraction"})
    report.table(ctx, by_year_peak, "c2_3_peak_by_year", units={"share_supported": "fraction"})
    peak_share = float(cov[cov.is_peak].supported.mean()) if cov.is_peak.any() else float("nan")
    L = ["# Check 2.3 ET coverage", "", f"Peak-month cohort cells with >= 4 models: {peak_share:.1%} "
         f"(target {C.TARGETS['et_cell_share']:.0%}).", "", report.md(by_month), "",
         "Unsupported cell-months are masked (Task 1.7); a missing non-peak month never removes a valid peak month."]
    report.text(ctx, "c2_3.md", "\n".join(L))
    ctx.complete(dict(peak_share_supported=peak_share, target=C.TARGETS["et_cell_share"],
                      june_share=float(by_month.set_index("month").share_supported.get(6, float("nan")))))
