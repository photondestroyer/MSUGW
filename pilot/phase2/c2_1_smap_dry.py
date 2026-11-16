"""Check 2.1  SMAP coverage in dry conditions.

"Count the dry county-years and dry episodes that remain after the soil-moisture mask. Show what remains of 2019
and 2022."   Target: at least 70% of dry county-years. If lower, report the range of drought conditions that
remain and whether they represent the region.
Physical idea: the soil-moisture contrast can only be tested in droughts the satellite saw; if SMAP misses the
worst months, the test would rest on mild droughts only.
Deliverables: c2_1_dry_coverage.csv, c2_1_by_year.csv, c2_1_ranges.csv, c2_1.md.
"""
import numpy as np


def run(ctx):
    import pandas as pd
    from .. import config as C
    from ..common import report
    # ══════ STEP 2.1-1  dry county-years and whether S is valid in their peak month ══════
    w = report.read(ctx.dep("t1_1"), "weather_county_year")
    S = report.read(ctx.dep("t1_3"), "county_month_S")
    for d in (w, S):
        d["county"] = d["county"].astype(str)
    j = w.merge(S[["county", "year", "month", "S", "S_valid"]], left_on=["county", "year", "peak_month"],
                right_on=["county", "year", "month"], how="left")
    j["S_valid"] = j["S_valid"].astype("boolean").fillna(False).astype(bool)
    dry = j[j.dry_season.astype(bool)]
    kept = dry[dry.S_valid]
    share = float(len(kept) / len(dry)) if len(dry) else np.nan

    # ══════ STEP 2.1-2  episodes kept; what remains of 2019 and 2022; ranges of drought kept vs all ══════
    ep_all, ep_kept = dry.episode_id.nunique(), kept.episode_id.nunique()
    by_year = dry.groupby("year").agg(dry_county_years=("county", "size"),
                                      kept=("S_valid", "sum")).reset_index()
    by_year["share_kept"] = by_year.kept / by_year.dry_county_years
    q = [0.1, 0.5, 0.9]
    ranges = pd.DataFrame(dict(quantile=q, deficit_all_dry_mm=dry.deficit_season_mm.quantile(q).to_numpy(),
                               deficit_kept_mm=kept.deficit_season_mm.quantile(q).to_numpy() if len(kept) else np.nan))
    tab = pd.DataFrame([dict(item="dry county-years", all=len(dry), kept=len(kept), share=share,
                             target=C.TARGETS["smap_dry_share"]),
                        dict(item="dry episodes", all=ep_all, kept=ep_kept, share=ep_kept / ep_all if ep_all else np.nan,
                             target=None)])
    # ══════ STEP 2.1-3  write ══════
    report.table(ctx, tab, "c2_1_dry_coverage", mask="dry seasons (D20), peak month S valid (D19)")
    report.table(ctx, by_year, "c2_1_by_year")
    report.table(ctx, ranges, "c2_1_ranges", units=dict(deficit_all_dry_mm="mm (May-Sep)"))
    focus = by_year[by_year.year.isin([2019, 2022])]
    L = ["# Check 2.1 SMAP coverage in dry conditions", "", report.md(tab), "",
         "## What remains of 2019 and 2022", "", report.md(focus) if len(focus) else "no dry county-years",
         "", "## Range of drought kept", "", report.md(ranges)]
    report.text(ctx, "c2_1.md", "\n".join(L))
    ctx.complete(dict(share_dry_county_years_kept=share, target=C.TARGETS["smap_dry_share"], dry_county_years=len(dry),
                      episodes_all=int(ep_all), episodes_kept=int(ep_kept),
                      y2019_2022=focus.to_dict("records")))
