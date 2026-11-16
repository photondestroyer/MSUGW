"""Check 2.7  Well capacity.

"Link about 20 cohort cells to WIMAS water rights and points of diversion. If the linking works, extend toward 200
cells across the thickness range. Keep three quantities apart: authorised rate; operating rate (reported use
divided by pumping hours); and verified capacity. Normalise by the area served."   No target. Report the share of
cells with a clear link and a usable rate.
Physical idea: thickness is a proxy; the mechanism is the rate a well can deliver. Where rates can be tied to
cells, the mechanism can be stated directly for a subset.
Links: cell -> irrigation points of diversion inside it (WIMAS) -> wells of the KGS Master Well Inventory with
that PDIV_ID -> WWC5 driller records (estimated yield at completion). Rates not in the downloaded tables
(authorised rate, pumping hours, verified capacity) go to a lookup sheet to be filled from WIMAS.
Deliverables: c2_7_links.csv, wimas_lookup_sheet.csv, c2_7.md.
"""
import numpy as np

HA_PER_ACRE = 0.404686


def _link(cells_geo, pod, mwi, wwc5):
    import geopandas as gpd
    import pandas as pd
    pts = gpd.GeoDataFrame(pod, geometry=gpd.points_from_xy(pod.long_nad83, pod.lat_nad83), crs=4269).to_crs(5070)
    j = gpd.sjoin(pts, cells_geo[["cell_id", "geometry"]], predicate="within")
    irr_cols = [c for c in pod.columns if c.startswith("AF_USED_IRR_")]
    j = j[j[irr_cols].fillna(0).sum(axis=1) > 0]
    m = mwi.dropna(subset=["WRIS_PDIV_ID"]).copy()
    m["WRIS_PDIV_ID"] = pd.to_numeric(m.WRIS_PDIV_ID, errors="coerce")
    jw = j.merge(m[["WRIS_PDIV_ID", "WELL_KID"]], left_on="PDIV_ID", right_on="WRIS_PDIV_ID", how="left")
    y = wwc5[["WELL_KID", "EST_YIELD", "COMPLE_DATE", "WELL_DEPTH", "STATUS"]].copy()
    y["EST_YIELD"] = pd.to_numeric(y.EST_YIELD, errors="coerce")
    jw = jw.merge(y, on="WELL_KID", how="left")
    return jw


def run(ctx):
    import pandas as pd
    from .. import config as C
    from ..common import io, report
    from ..phase0 import t0_2_cells as T02
    rs = np.random.default_rng(C.CAPACITY_SEED)
    # ══════ STEP 2.7-1  cohort cells across the thickness range (stratified by third) ══════
    th = report.read(ctx.dep("c2_4"), "thickness_cell_year")
    co = th[(th.role == "cohort") & th.year.isin(C.PANEL_YEARS) & ~th.A_poor.astype(bool)]
    third = co.groupby("cell_id").A_third.first()
    pod = io.table("wimas_pods")
    mwi = io.table("master_wells")
    wwc5 = io.table("wwc5")
    geo = T02.load_geometry(ctx.dep("t0_2")).to_crs(5070)
    cells = T02.load(ctx.dep("t0_2"), columns=["cell_id", "area_ha"])
    results = {}
    for n_target in (C.CAPACITY_N_START, C.CAPACITY_N_EXTEND):
        pick = []
        for g in ("thin", "middle", "thick"):
            pool = third[third == g].index.to_numpy()
            k = min(len(pool), int(np.ceil(n_target / 3)))
            pick += list(rs.choice(np.sort(pool), size=k, replace=False)) if k else []
        # ══════ STEP 2.7-2  links: cell -> points of diversion -> wells -> WWC5 yields ══════
        jw = _link(geo[geo.cell_id.isin(pick)], pod, mwi, wwc5)
        per = (jw.groupby("cell_id").agg(pods=("PDIV_ID", "nunique"), wells=("WELL_KID", "nunique"),
                                         driller_yield_gpm=("EST_YIELD", "sum"),
                                         acres_2016=("ACRES_2016", "sum") if "ACRES_2016" in jw else ("PDIV_ID", "size"))
               .reindex(pick).reset_index())
        per["third"] = per.cell_id.map(third)
        per["area_ha"] = per.cell_id.map(cells.set_index("cell_id").area_ha)
        per["clear_link"] = (per.pods >= 1) & (per.wells >= 1)
        per["usable_rate"] = per.clear_link & (per.driller_yield_gpm > 0)
        per["driller_yield_gpm_per_acre"] = per.driller_yield_gpm / (per.area_ha / HA_PER_ACRE)
        results[n_target] = (per, jw)
        share = float(per.usable_rate.mean()) if len(per) else 0.0
        ctx.log.info("STEP 2.7-2  %d cells: clear link %.0f%%, usable driller yield %.0f%%", len(per),
                     100 * per.clear_link.mean(), 100 * share)
        if share < 0.5:                       # "If the linking works, extend toward 200 cells"
            break
    n_used = max(results)
    per, jw = results[n_used]

    # ══════ STEP 2.7-3  the lookup sheet for rates not in the tables (authorised, operating, verified) ══════
    sheet = jw[["cell_id", "PDIV_ID", "WELL_KID"]].drop_duplicates().assign(
        authorised_rate_gpm="", authorised_quantity_af="", pumping_hours_2016="", pumping_hours_2020="",
        verified_capacity_gpm="", source="", notes="")
    # ══════ STEP 2.7-4  write ══════
    report.table(ctx, per, "c2_7_links", units=dict(driller_yield_gpm="gal/min (initial, WWC5)", area_ha="ha"),
                 notes="driller yield is a fourth quantity: not authorised, operating or verified rate")
    report.table(ctx, sheet, "wimas_lookup_sheet", notes="fill from WIMAS by PDIV_ID; keep the three rates apart")
    L = ["# Check 2.7 Well capacity", "", f"Cells linked: {len(per)} (extended to {n_used} when at least half the "
         f"first {C.CAPACITY_N_START} had a usable rate).", "",
         f"Clear link (point of diversion and well): {per.clear_link.mean():.0%}. Usable driller yield: "
         f"{per.usable_rate.mean():.0%}.", "", "Authorised rate, operating rate (reported use / pumping hours) and "
         "verified capacity are not in the downloaded tables: wimas_lookup_sheet.csv lists the PDIV_IDs to look up."]
    report.text(ctx, "c2_7.md", "\n".join(L))
    ctx.complete(dict(cells=int(len(per)), share_clear_link=float(per.clear_link.mean()),
                      share_usable_rate=float(per.usable_rate.mean()), extended=n_used > C.CAPACITY_N_START))
