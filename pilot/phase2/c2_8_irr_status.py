"""Check 2.8  Irrigation status after 2020.

"Check whether water-right and field records identify which cells were irrigated in 2021-2024. Reported irrigated
acres alone do not locate irrigation inside a quarter-section."   No target. Report coverage. Keep status as
unknown unless the link is clear.
Physical idea: a cohort cell that stopped irrigating after 2020 would show low ET for reasons unrelated to well
capacity in a drought; knowing which did keeps that from being misread.
Clear link: the cell holds at least one irrigation point of diversion and the acres they report fit inside the cell
(at most 110% of its area). Status: irrigated (acres > 0) / not irrigated (acres and use 0) / unknown otherwise.
WIMAS 2024 is not in the downloaded table, so 2024 is unknown by construction.
Deliverables: c2_8_status.csv (cohort cell-years 2021-2024), c2_8_coverage.csv, c2_8.md.
"""
import numpy as np

YEARS = [2021, 2022, 2023, 2024]
HA_PER_ACRE = 0.404686


def run(ctx):
    import geopandas as gpd
    import pandas as pd
    from ..common import io, report
    from ..phase0 import t0_2_cells as T02
    # ══════ STEP 2.8-1  cohort cells and the irrigation points of diversion inside them ══════
    cells = T02.load(ctx.dep("t0_2"), columns=["cell_id", "county", "role", "area_ha"])
    geo = T02.load_geometry(ctx.dep("t0_2"))
    geo = geo[geo.cell_id.isin(cells[cells.role == "cohort"].cell_id)].to_crs(5070)
    pod = io.table("wimas_pods")
    pts = gpd.GeoDataFrame(pod, geometry=gpd.points_from_xy(pod.long_nad83, pod.lat_nad83), crs=4269).to_crs(5070)
    j = gpd.sjoin(pts, geo[["cell_id", "geometry"]], predicate="within")
    # ══════ STEP 2.8-2  status per cohort cell-year, 2021-2024 ══════
    rows = []
    area = cells.set_index("cell_id").area_ha
    by = {k: g for k, g in j.groupby("cell_id")}
    for cid in geo.cell_id:
        g = by.get(cid)
        for y in YEARS:
            ac, af = f"ACRES_{y}", f"AF_USED_IRR_{y}"
            if g is None:
                rows.append(dict(cell_id=cid, year=y, pods=0, status="unknown", why="no point of diversion in the cell"))
                continue
            if ac not in g.columns:
                rows.append(dict(cell_id=cid, year=y, pods=len(g), status="unknown", why="year not in WIMAS table"))
                continue
            acres = float(g[ac].fillna(0).sum())
            use = float(g[af].fillna(0).sum()) if af in g.columns else np.nan
            fits = acres * HA_PER_ACRE <= 1.1 * area[cid]
            if not fits:
                st, why = "unknown", "reported acres exceed the cell (right serves more land)"
            elif acres > 0:
                st, why = "irrigated", "clear link"
            elif use == 0:
                st, why = "not irrigated", "clear link, no acres and no use reported"
            else:
                st, why = "unknown", "use without acres"
            rows.append(dict(cell_id=cid, year=y, pods=len(g), acres=acres, af_used=use, status=st, why=why))
    s = pd.DataFrame(rows)
    cov = s.groupby("year").status.value_counts(normalize=True).unstack(fill_value=0).reset_index()
    # ══════ STEP 2.8-3  write ══════
    report.table(ctx, s, "c2_8_status", units=dict(acres="acres", af_used="acre-feet"), mask="cohort cells")
    report.table(ctx, cov, "c2_8_coverage", units={"*": "share of cohort cells"})
    report.text(ctx, "c2_8.md", "# Check 2.8 Irrigation status after 2020\n\n" + report.md(cov) +
                "\n\nStatus stays unknown unless the link is clear (see the module docstring).")
    ctx.complete(dict(share_known=float((s.status != "unknown").mean()) if len(s) else None,
                      by_year=cov.set_index("year").round(3).to_dict("index")))
