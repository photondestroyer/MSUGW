"""Task 0.2  Cells and masks.

"Build quarter-section cells for all Kansas High Plains aquifer counties. Attach county, SMAP footprint
identifiers, root-zone water storage (gSSURGO) and a surface-water flag. Write down the rules that turn annual
crop and irrigation pixels into cell labels."
Physical idea: the cohort is the set of fields whose irrigation history makes them depend on wells; rainfed
cells beside them give the baseline response of ET to soil moisture. Both are fixed from 2011-2015 labels,
before the panel, so no panel-year outcome can move a cell between groups.
Decisions: D03 (cells), D04-D06, D09 (cohort), D07-D08 (rainfed). D06 names the primary irrigation map (LANID)
and the check map (AIM-HPA); the code reads both names from the decision, never from here.
Deliverables: cells.gpkg (geometry), cells.parquet (attributes and yearly shares), label_rules.md,
role_summary.csv.
"""
import numpy as np

AIM_LAST_YEAR, LANID_LAST_YEAR = 2017, 2020   # what each product holds (AIM-HPA 2017; LANID file through 2020)
CDL_YEARS = list(range(2011, 2016))


def _map_years(C, last):
    """Label years 2011-2015 plus the later-irrigation years a product holds."""
    return [y for y in C.LABEL_YEARS + C.LATER_IRR_YEARS if y <= last]


def _sections_in_scope(ctx, counties):
    import geopandas as gpd
    from .. import config as C
    from ..common import io
    sec = io.vector("plss_sections").to_crs(5070)
    sec = sec[sec.geometry.notna() & ~sec.geometry.is_empty]
    if C.RUN_MODE == "smoke":
        counties = counties[counties.county == C.SMOKE_COUNTY]
    cen = gpd.GeoDataFrame(dict(i=np.arange(len(sec))), geometry=sec.geometry.centroid.values, crs=5070)
    j = gpd.sjoin(cen, counties[["county", "geometry"]], predicate="within", how="inner")
    out = sec.iloc[j["i"].to_numpy()].copy()
    out["county"] = j["county"].to_numpy()
    return out.reset_index(drop=True)


def _accumulate_5070(ctx, cells, grid_transform, grid_shape, read_window, names, label):
    """Strips over a 30 m EPSG:5070 grid: ids once per strip, then every layer read_window(name, window)."""
    from ..common import cells as CL, zonal
    acc = zonal.Accumulator(len(cells), ["px"] + [f"{n}_{k}" for n in names for k in ("hit", "map")])
    n_strips = 0
    for r0, r1, c0, c1, ids in CL.strips(cells, grid_transform, grid_shape):
        acc.add("px", ids)
        for n in names:
            hit, mapped = read_window(n, r0, r1, c0, c1)
            acc.add(f"{n}_hit", ids, mask=hit)
            acc.add(f"{n}_map", ids, mask=mapped)
        n_strips += 1
        ctx.check_budget()
    ctx.log.info("   %s: %d strips", label, n_strips)
    return acc


def run(ctx):
    import geopandas as gpd
    import pandas as pd
    import rasterio
    from .. import config as C
    from .. import decisions as D
    from .. import definitions as DEF
    from ..common import io, report, grids, cells as CL, zonal, stats
    from . import counties as CTY

    counties = CTY.load(ctx.dep("counties"))

    # ══════ STEP 0.2-1  PLSS sections whose centre lies in a study county ══════
    sec = _sections_in_scope(ctx, counties)
    ctx.log.info("STEP 0.2-1  %d sections in %d counties", len(sec), sec.county.nunique())

    # ══════ STEP 0.2-2  cells = quarter sections (D03); county, GMD, aquifer, EASE ids, centroids ══════
    cells = DEF.build_cells(sec[["FRSTDIVID", "geometry"]])
    cells["county"] = cells["section_id"].map(dict(zip(sec.FRSTDIVID.astype(str), sec.county)))
    cen = cells.geometry.centroid
    cells["x5070"], cells["y5070"] = cen.x.to_numpy(), cen.y.to_numpy()
    ll = gpd.GeoSeries(cen, crs=5070).to_crs(4326)
    cells["lon"], cells["lat"] = ll.x.to_numpy(), ll.y.to_numpy()
    ea = gpd.GeoSeries(cen, crs=5070).to_crs(6933)
    cells["r9"], cells["c9"] = grids.ease_index(ea.x, ea.y, C.EASE_M09)
    cells["r36"], cells["c36"] = grids.ease_index(ea.x, ea.y, C.EASE_M36)
    cells["fp"] = grids.fp_id(cells.r36, cells.c36)
    pts = gpd.GeoDataFrame(dict(i=np.arange(len(cells))), geometry=cen.values, crs=5070)
    gmd = io.vector("gmd").to_crs(5070)[["GMD_ID", "geometry"]]
    j = gpd.sjoin(pts, gmd, predicate="within", how="left").drop_duplicates("i").sort_values("i")
    cells["gmd"] = j["GMD_ID"].to_numpy()
    aq = io.vector("aquifer_shp").to_crs(5070).union_all()
    cells["in_aquifer"] = pts.within(aq).to_numpy()
    ctx.log.info("STEP 0.2-2  %d cells; %.1f%% regular (D03); %.1f%% inside the aquifer polygon", len(cells),
                 100 * cells.regular.mean(), 100 * cells.in_aquifer.mean())

    # ══════ STEP 0.2-3  static attributes: root-zone storage, surface water, aquifer base, land surface ══════
    # (a) gSSURGO rootznaws (mm; the file once labelled it cm), mean over the cell's 120 m pixels
    gs = io.gssurgo()
    gtr = io.gssurgo_transform(gs)
    var = gs["rootznaws"]
    acc = zonal.Accumulator(len(cells), ["rz"])
    for r0, r1, c0, c1, ids in CL.strips(cells, gtr, var.shape, strip_rows=2048):
        a = np.ma.filled(var[r0:r1, c0:c1].astype("float64"), np.nan)
        acc.add("rz", ids, a)                     # gSSURGO rootznaws is mm (Valu1)
    cells["rzws_mm"] = acc.mean("rz")
    cells["rzws_px"] = acc.cnt["rz"]
    gs.close()
    # (b) surface-water points of diversion with irrigation use (D09)
    pod = io.table("wimas_pods")
    irr_cols = [c for c in pod.columns if c.startswith(("AF_USED_IRR_", "ACRES_"))]
    sw = pod[(pod["source"] == "S") & (pod[irr_cols].fillna(0).to_numpy() > 0).any(axis=1)]
    swp = gpd.GeoSeries(gpd.points_from_xy(sw.long_nad83, sw.lat_nad83), crs=4269).to_crs(5070)
    cells["sw_km"] = stats.nearest_km(np.c_[cells.x5070, cells.y5070], np.c_[swp.x, swp.y])
    cells["surface_water"] = DEF.served_by_surface_water(cells["sw_km"].to_numpy())
    # (c) KGS section properties (ft): base, land surface, thickness; the section holding the cell centre
    ks = io.vector("kgs_sections").to_crs(5070)
    keep = ["TRS", "hpa_bedrock", "ned_elev", "predev_wte_updated", "predev_satthick_updated", "satthick_2022_2024",
            "wte_2022_2024", "eul_3yr_2013_2023_200gpm"]
    j = gpd.sjoin(pts, ks[keep + ["geometry"]], predicate="within", how="left").drop_duplicates("i").sort_values("i")
    for k in keep:
        cells[f"kgs_{k}"] = j[k].to_numpy()
    # (d) USGS bottom surface (ft NAVD88, 100 m; an independent base estimate) and land surface (GLO30, m)
    cells["usgs_bottom_ft"] = io.sample_raster("usgs_bottom", cells.x5070, cells.y5070, "EPSG:5070")
    cells["land_m"] = io.sample_raster("dem", cells.x5070, cells.y5070, "EPSG:5070")
    # (e) KGS extent classes ("little" saturated thickness, "outcrop") at the cell centre
    ext = io.vector("kgs_extent").to_crs(5070)[["DESCRIPTIO", "geometry"]]
    j = gpd.sjoin(pts, ext, predicate="within", how="left").drop_duplicates("i").sort_values("i")
    cells["kgs_extent_class"] = j["DESCRIPTIO"].fillna("").to_numpy()
    ctx.log.info("STEP 0.2-3  rzws known %.1f%%; surface-water cells %d; KGS base known %.1f%%",
                 100 * np.isfinite(cells.rzws_mm).mean(), int(cells.surface_water.sum()),
                 100 * np.isfinite(cells.kgs_hpa_bedrock).mean())

    # ══════ STEP 0.2-4  irrigation shares per cell-year from both maps (D06 says which is primary) ══════
    # A share is known only where the map labels enough of the cell (DEF.map_share, D06 min_map_coverage).
    primary, check = DEF.map_roles()
    AIM_YEARS, LANID_YEARS = _map_years(C, AIM_LAST_YEAR), _map_years(C, LANID_LAST_YEAR)
    with rasterio.open(io.aim_path(AIM_YEARS[0])) as r0_:
        aim_tr, aim_shape = r0_.transform, r0_.shape
    aim_src = {y: rasterio.open(io.aim_path(y)) for y in AIM_YEARS}

    def aim_window(name, r0, r1, c0, c1):
        a = aim_src[int(name[1:])].read(1, window=((r0, r1), (c0, c1)))
        return a == io.AIM_IRRIGATED, (a == io.AIM_IRRIGATED) | (a == io.AIM_NOT_IRRIGATED)
    acc = _accumulate_5070(ctx, cells, aim_tr, aim_shape, aim_window, [f"a{y}" for y in AIM_YEARS], "AIM-HPA")
    for y in AIM_YEARS:
        cells[f"aim_{y}"], cells[f"aim_cov_{y}"] = DEF.map_share(acc.cnt[f"a{y}_hit"], acc.cnt[f"a{y}_map"],
                                                                  acc.cnt["px"])
    for s in aim_src.values():
        s.close()
    ln = io.lanid()
    ln_years = [int(v) for v in ln["time"][:]]
    missing = [y for y in LANID_YEARS if y not in ln_years]
    if missing:
        from ..core.errors import SourceError
        raise SourceError(f"LANID file {io.path('lanid')} lacks the years {missing}: build the 1997-2020 file "
                          f"(local data processing/build_lanid_1997_2020.py)")
    lvar = ln["irrigated"]

    def lanid_window(name, r0, r1, c0, c1):
        a = np.asarray(lvar[ln_years.index(int(name[1:])), r0:r1, c0:c1])
        return a == 1, (a == 0) | (a == 1)
    acc = _accumulate_5070(ctx, cells, io.lanid_transform(ln), lvar.shape[1:], lanid_window,
                           [f"l{y}" for y in LANID_YEARS], "LANID")
    for y in LANID_YEARS:
        cells[f"lanid_{y}"], cells[f"lanid_cov_{y}"] = DEF.map_share(acc.cnt[f"l{y}_hit"], acc.cnt[f"l{y}_map"],
                                                                      acc.cnt["px"])
    ln.close()
    ctx.log.info("STEP 0.2-4  primary map %s known for %.1f%% of cells in 2015; check map %s for %.1f%%",
                 primary, 100 * np.isfinite(cells[f"{primary}_2015"]).mean(),
                 check, 100 * np.isfinite(cells[f"{check}_2015"]).mean())

    # ══════ STEP 0.2-5  crop shares per cell-year from CDL tiles, 2011-2015 ══════
    classes = list(C.CDL_CLASSES) + ["fallow", "noncrop"]
    lut = np.full(256, classes.index("noncrop"), "int64")
    for k, codes in C.CDL_CLASSES.items():
        lut[codes] = classes.index(k)
    lut[61] = classes.index("fallow")
    K = len(classes)
    counts = {y: np.zeros((len(cells), K)) for y in CDL_YEARS}
    frames = io.tile_frames("cdl")
    sig = CL.cells_signature(cells)
    cache = f"{C.shared_root()}/cellids"
    for tile in C.KS_TILES:
        F = frames.get(tile)
        if F is None:
            continue
        ids = CL.tile_ids(cells, F, cache, sig)
        m = ids >= 0
        if not m.any():
            continue
        idv = ids[m].astype("int64")
        for y in CDL_YEARS:
            p = io.cdl_file(tile, y)
            ds = io.open_tile(p)
            io.tile_frame_of(ds, frames, tile)
            c = np.asarray(ds["cropland"][:])[m]
            ds.close()
            ok = c != 255
            counts[y] += np.bincount(idv[ok] * K + lut[c[ok]], minlength=len(cells) * K).reshape(len(cells), K)
        ctx.log.info("   CDL tile %s: %d cells", tile, len(np.unique(idv)))
        ctx.check_budget()
    for y in CDL_YEARS:
        tot = counts[y].sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            sh = counts[y] / tot[:, None]
        for k, name in enumerate(classes):
            cells[f"cdl_{name}_{y}"] = np.where(tot > 0, sh[:, k], np.nan)
        cells[f"cdl_cropland_{y}"] = np.where(tot > 0, 1.0 - sh[:, classes.index("noncrop")], np.nan)
        cells[f"cdl_cov_{y}"] = tot * C.PIXEL_M ** 2 / 1e4 / cells["area_ha"].to_numpy()
    ctx.log.info("STEP 0.2-5  CDL known for %.1f%% of cells in 2015", 100 * np.isfinite(cells["cdl_corn_2015"]).mean())

    # ══════ STEP 0.2-6  labels: cohort, rainfed, crop class, matching, later irrigation (definitions.py) ══════
    L = C.LABEL_YEARS
    prim = cells[[f"{primary}_{y}" for y in L]].to_numpy()          # decides "irrigated" (D06: LANID)
    chk = cells[[f"{check}_{y}" for y in L]].to_numpy()             # the secondary check (D06: AIM-HPA)
    corn = cells[[f"cdl_corn_{y}" for y in C.CORN_YEARS]].to_numpy()
    crop = cells[[f"cdl_cropland_{y}" for y in L]].to_numpy()
    cells["primary_map"], cells["check_map"] = primary, check
    cells["irr_years"] = DEF.is_irrigated_year(prim).sum(axis=1)
    cells["irr_years_check"] = DEF.is_irrigated_year(chk).sum(axis=1)
    cells["corn_years"] = DEF.is_corn_year(corn).sum(axis=1)
    cells["maps_agree_years"], _ = DEF.maps_agree(prim, chk)
    cohort, why_c = DEF.is_cohort(prim, chk, corn, cells["surface_water"].to_numpy(), cells["regular"].to_numpy())
    rainfed, why_r = DEF.is_rainfed(prim, chk, crop, cells["regular"].to_numpy())
    names = [k for k in C.CDL_CLASSES]
    shares = np.stack([cells[[f"cdl_{k}_{y}" for k in names]].to_numpy() for y in C.CORN_YEARS], axis=1)
    dom = DEF.dominant_class(shares, names)
    for i, y in enumerate(C.CORN_YEARS):
        cells[f"dominant_{y}"] = dom[:, i]
    cells["crop_class"] = DEF.crop_class(dom)
    nan = np.full(len(cells), np.nan)                     # a map without that year (AIM-HPA after 2017)
    later = np.stack([np.column_stack([cells[f"{m}_{y}"].to_numpy() if f"{m}_{y}" in cells else nan
                                       for y in C.LATER_IRR_YEARS]) for m in (primary, check)], axis=-1)
    cells["irrigated_later"] = DEF.irrigated_later(later)
    role = np.where(cohort, "cohort", np.where(rainfed, "rainfed", "excluded")).astype(object)
    reason = np.where(cohort | rainfed, "", "").astype(object)
    irr_any = DEF.is_irrigated_year(prim).any(axis=1)
    reason[~cohort & ~rainfed & irr_any] = "irrigated, not cohort: " + why_c[~cohort & ~rainfed & irr_any]
    reason[~cohort & ~rainfed & ~irr_any] = "not rainfed: " + why_r[~cohort & ~rainfed & ~irr_any]
    out_of_aq = ~cells["in_aquifer"].to_numpy()
    reason[out_of_aq & (role == "excluded")] = "outside the aquifer polygon"
    cells["role"], cells["reason"] = role, reason
    matched, n_match = DEF.match_rainfed(cells.county, cells.crop_class, cells.rzws_mm, cohort, rainfed)
    cells["rainfed_matched"], cells["cohort_n_rainfed_matches"] = matched, n_match
    ctx.log.info("STEP 0.2-6  cohort %d, rainfed %d (matched %d), excluded %d", int(cohort.sum()), int(rainfed.sum()),
                 int(matched.sum()), int((role == "excluded").sum()))

    # ══════ STEP 0.2-7  write: cells (geometry + table), label rules, role summary ══════
    import os
    from ..core.drive import deliver
    tmp = f"{ctx.scratch}/cells.gpkg"
    if os.path.exists(tmp):
        os.remove(tmp)
    cells[["cell_id", "county", "role", "geometry"]].to_file(tmp, layer="cells", driver="GPKG")
    deliver(tmp, f"{ctx.out}/cells.gpkg")
    tab = pd.DataFrame(cells.drop(columns="geometry"))
    report.table(ctx, tab, "cells", parquet=True, csv=False,
                 units=dict(area_ha="ha", rzws_mm="mm", sw_km="km", kgs_hpa_bedrock="ft", usgs_bottom_ft="ft",
                            land_m="m", aim_YYYY="share", lanid_YYYY="share", cdl_class_YYYY="share"),
                 notes=f"primary irrigation map {primary}, check map {check} (D06); a yearly share is NaN where "
                       f"the map labels too little of the cell (D06 min_map_coverage)")
    summ = tab.groupby(["county", "role"]).size().unstack(fill_value=0).reset_index()
    report.table(ctx, summ, "role_summary", units={"cells": "count"})
    reasons = tab.loc[tab.role == "excluded", "reason"].value_counts().rename_axis("reason").reset_index(name="cells")
    report.table(ctx, reasons, "exclusion_reasons")
    from .. import docs_gen
    report.text(ctx, "label_rules.md", docs_gen.definitions_md())
    # Pilot: "If rainfed corn is scarce, report results by crop and do not pool different crops."
    by_crop = (tab[tab.role != "excluded"].groupby(["county", "crop_class", "role"]).size().unstack(fill_value=0)
               .reindex(columns=["cohort", "rainfed"], fill_value=0).reset_index())
    report.table(ctx, by_crop, "cells_by_county_crop", units={"cohort": "cells", "rainfed": "cells"},
                 notes="matching is within county and crop class only (D08); crops are never pooled")
    corn_rf = int(((tab.role == "rainfed") & (tab.crop_class == "corn")).sum())
    ctx.complete(dict(cells=int(len(tab)), cohort=int(cohort.sum()), rainfed=int(rainfed.sum()),
                      rainfed_matched=int(matched.sum()), rainfed_corn=corn_rf,
                      rainfed_corn_scarce=bool(corn_rf < max(30, 0.05 * int(cohort.sum()))),
                      regular_share=float(tab.regular.mean()), primary_map=primary, check_map=check,
                      cohort_by_check_map_alone=int(DEF.is_cohort(chk, prim, corn, cells["surface_water"].to_numpy(),
                                                                  cells["regular"].to_numpy())[0].sum()),
                      rainfed_irrigated_later=int((cells["irrigated_later"].to_numpy() & rainfed).sum()),
                      crop_classes_cohort=tab[tab.role == "cohort"].crop_class.value_counts().to_dict(),
                      crop_classes_rainfed=tab[tab.role == "rainfed"].crop_class.value_counts().to_dict(),
                      signature_cells=sig))


def load(stage_dir, columns=None):
    import pandas as pd
    return pd.read_parquet(f"{stage_dir}/cells.parquet", columns=columns)


def load_geometry(stage_dir):
    import geopandas as gpd
    return gpd.read_file(f"{stage_dir}/cells.gpkg", layer="cells")
