"""Check 2.2  Reference footprints.

"Count the usable reference footprints for each county and year, as actual footprints and not as 9 km grid
cells. Report their overlap."   Target: at least three within 100 km for at least 80% of counties. One widening
to 150 km is allowed if it is decided beforehand (D13).
Physical idea: several independent footprints average out one footprint's local storm or retrieval error, so S
stands for the region's drought.
Usable = a reference footprint whose month is supported (Task 1.3) in the county's peak month of that year.
Deliverables: c2_2_county_year.csv, c2_2_overlap.csv, c2_2.md.
"""
import numpy as np


def run(ctx):
    import pandas as pd
    from .. import config as C
    from ..common import report
    from ..phase0 import counties as CTY
    # ══════ STEP 2.2-1  usable reference footprints per county-year (supported in the peak month) ══════
    ref = report.read(ctx.dep("t1_2"), "county_reference_weights")
    fm = report.read(ctx.dep("t1_3"), "footprint_month")
    w = report.read(ctx.dep("t1_1"), "weather_county_year")
    for d in (ref, w):
        d["county"] = d["county"].astype(str)
    cy = w[["county", "year", "peak_month"]]
    j = cy.merge(ref[["county", "fp"]], on="county").merge(
        fm[["fp", "year", "month", "supported"]], left_on=["fp", "year", "peak_month"],
        right_on=["fp", "year", "month"], how="left")
    j["supported"] = j["supported"].astype("boolean").fillna(False).astype(bool)
    n = j.groupby(["county", "year"]).supported.sum().rename("usable").reset_index()
    n = cy.merge(n, on=["county", "year"], how="left").fillna({"usable": 0})
    n["meets"] = n.usable >= C.TARGETS["fp_min"]
    per_year = n.groupby("year").meets.mean().rename("share_counties").reset_index()
    per_county = n.groupby("county").usable.median()
    share_overall = float((per_county >= C.TARGETS["fp_min"]).mean())

    # ══════ STEP 2.2-2  overlap between neighbouring counties (shared footprints, Jaccard index) ══════
    nb = CTY.neighbours(ctx.dep("counties"))
    sets = ref.groupby("county").fp.apply(set).to_dict()
    ov = []
    for a, bs in nb.items():
        for b in bs:
            if a < b and a in sets and b in sets:
                inter, union = len(sets[a] & sets[b]), len(sets[a] | sets[b])
                ov.append(dict(county_a=a, county_b=b, shared=inter, jaccard=inter / union if union else np.nan))
    ov = pd.DataFrame(ov)
    radius = __import__("json").load(open(f"{ctx.dep('t1_2')}/_complete.json"))["summary"]["radius_km"]
    # ══════ STEP 2.2-3  write ══════
    report.table(ctx, n, "c2_2_county_year", units=dict(usable="footprints"))
    report.table(ctx, per_year, "c2_2_by_year", units=dict(share_counties="fraction"))
    report.table(ctx, ov, "c2_2_overlap", units=dict(jaccard="0-1"))
    L = ["# Check 2.2 Reference footprints", "", f"Radius used: {radius} km (decision D13).",
         f"Counties whose median usable count is at least {C.TARGETS['fp_min']}: {share_overall:.0%} "
         f"(target {C.TARGETS['fp_county_share']:.0%}).", "", report.md(per_year), "",
         f"Neighbouring-county overlap: median Jaccard {ov.jaccard.median():.2f}" if len(ov) else ""]
    report.text(ctx, "c2_2.md", "\n".join(L))
    ctx.complete(dict(radius_km=radius, share_counties_meeting=share_overall, target=C.TARGETS["fp_county_share"],
                      by_year=per_year.set_index("year").share_counties.round(3).to_dict(),
                      median_jaccard=float(ov.jaccard.median()) if len(ov) else None))
