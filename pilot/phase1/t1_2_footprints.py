"""Task 1.2  Reference footprints.

"Define the reference footprints before counting SMAP dates. Apply the rule with 2011-2015 labels and check it
against 2018-2020 labels. Record which footprints are shared between counties."
Physical idea: irrigation wets the soil SMAP sees, so soil moisture over irrigated land would echo the very
water the cohort pumps. Footprints with almost no irrigation measure the drought the region is having, not the
response to it - the irrigation-free drought signal the design rests on.
Decisions: D11 (footprint grid), D12 (irrigation share), D13 (radius and its one widening), D14 (similarity),
D15 (weights).
Deliverables: footprints.csv, county_footprints.csv, shared_footprints.csv, footprints_map.png.
"""
import numpy as np

DEC = 4                     # irrigation maps read at every 4th pixel (120 m) for 36 km footprint shares


def _footprint_grid(cty5070, max_km):
    """36 km EASE cells whose centre lies within max_km + 36 km of any county centroid."""
    import geopandas as gpd
    import shapely
    from .. import config as C
    from ..common import grids
    c = gpd.GeoSeries(gpd.points_from_xy(cty5070.x5070, cty5070.y5070), crs=5070).to_crs(6933)
    pad = (max_km + 40) * 1000
    r_a, c_a = grids.ease_index(c.x.min() - pad, c.y.max() + pad, C.EASE_M36)
    r_b, c_b = grids.ease_index(c.x.max() + pad, c.y.min() - pad, C.EASE_M36)
    R, Cc = np.meshgrid(np.arange(r_a, r_b + 1), np.arange(c_a, c_b + 1), indexing="ij")
    R, Cc = R.ravel(), Cc.ravel()
    boxes = [shapely.box(*grids.ease_box(r, cc, C.EASE_M36)) for r, cc in zip(R, Cc)]
    fp = gpd.GeoDataFrame(dict(fp=grids.fp_id(R, Cc), r36=R, c36=Cc), geometry=boxes, crs=6933).to_crs(5070)
    cen = fp.geometry.centroid
    fp["x5070"], fp["y5070"] = cen.x, cen.y
    d = np.hypot(fp.x5070.to_numpy()[:, None] - cty5070.x5070.to_numpy()[None, :],
                 fp.y5070.to_numpy()[:, None] - cty5070.y5070.to_numpy()[None, :]) / 1000
    return fp[d.min(axis=1) <= max_km + 36].reset_index(drop=True)


def _map_shares(ctx, fp, paths, irrigated, labelled, label):
    """Per footprint: largest yearly irrigated share (of the labelled pixels) and smallest yearly label coverage.
    Coverage = labelled pixels / pixels of the WHOLE footprint (from its area), so the part of a footprint that
    lies beyond the edge of the raster counts as not labelled.
    paths: {year: raster file}, all on one EPSG:5070 grid; irrigated(a), labelled(a): boolean masks of raster
    values. Decimated read (every DEC-th pixel)."""
    import rasterio
    from rasterio import features, windows
    from ..core.errors import SourceError
    years = sorted(paths)
    with rasterio.open(paths[years[0]]) as r:
        grid = (r.crs, r.transform, r.shape)
        g = fp.geometry if fp.crs == r.crs else fp.geometry.to_crs(r.crs)
        x0, y0, x1, y1 = g.total_bounds
        w = windows.from_bounds(x0, y0, x1, y1, r.transform).round_offsets().round_lengths()
        c_a, r_a = max(int(w.col_off), 0), max(int(w.row_off), 0)
        c_b, r_b = min(int(w.col_off + w.width), r.width), min(int(w.row_off + w.height), r.height)
        if c_b <= c_a or r_b <= r_a:
            raise SourceError(f"{label}: {paths[years[0]]} does not reach the footprints")
        win = windows.Window(c_a, r_a, c_b - c_a, r_b - r_a)
        shape = (int(win.height // DEC), int(win.width // DEC))
        tr = windows.transform(win, r.transform) * rasterio.Affine.scale(win.width / shape[1], win.height / shape[0])
    ids = features.rasterize(((geom, i) for i, geom in enumerate(g)), out_shape=shape, transform=tr, fill=-1,
                             dtype="int32")
    m = ids >= 0
    idv = ids[m]
    n = len(fp)
    whole = np.array([geom.area for geom in g]) / abs(tr.a * tr.e)       # EPSG:5070 is equal-area
    shares, cov = [], []
    for y in years:
        with rasterio.open(paths[y]) as r:
            if (r.crs, r.transform, r.shape) != grid:
                raise SourceError(f"{label} {y}: grid differs from {years[0]} ({paths[y]})")
            a = r.read(1, window=win, out_shape=shape)[m]
        lab = labelled(a)
        hit = np.bincount(idv[lab & irrigated(a)], minlength=n)
        mapped = np.bincount(idv[lab], minlength=n)
        with np.errstate(invalid="ignore", divide="ignore"):
            shares.append(np.where(mapped > 0, hit / mapped, np.nan))
            cov.append(np.minimum(mapped / whole, 1.0))
        ctx.log.info("   %s %d read", label, y)
    shares = np.vstack(shares)
    known = np.isfinite(shares).any(axis=0)
    largest = np.max(np.where(np.isfinite(shares), shares, -1.0), axis=0)
    return np.where(known, largest, np.nan), np.min(np.vstack(cov), axis=0)


def _aim_shares(ctx, fp, years):
    """Largest yearly irrigated share and label coverage per footprint from AIM-HPA (D12's map)."""
    from ..common import io
    return _map_shares(ctx, fp, {y: io.aim_path(y) for y in years},
                       lambda a: a == io.AIM_IRRIGATED,
                       lambda a: (a == io.AIM_IRRIGATED) | (a == io.AIM_NOT_IRRIGATED), "AIM-HPA")


def _lanid_later_shares(ctx, fp, years):
    """The same from LANID 2018-2020 as released (CONUS rasters: 1 = irrigated, 0 = not irrigated or no data, so
    every pixel inside the raster counts as labelled). None when the files are not on disk."""
    from ..common import io
    paths = {y: p for y, p in io.lanid_1820_paths().items() if y in years}
    if sorted(paths) != sorted(years):
        return None
    return _map_shares(ctx, fp, paths, lambda a: a == 1, lambda a: np.ones(a.shape, bool), "LANID")


def _rain_normal(ctx, polys4326, years):
    """Mean May-September precipitation (mm) per polygon over the normal years, from gridMET cells whose centre
    lies inside it; NaN when fewer than half of the polygon's area is on the gridMET lab grid."""
    import shapely
    from .. import config as C
    from ..common import io
    ds = io.gridmet()
    lat, lon, dates = io.gridmet_axes(ds)
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    cell_area = np.cos(np.radians(LA))
    owner = np.full(LA.shape, -1, int)
    for i, g in enumerate(polys4326):
        inside = shapely.contains_xy(g, LO, LA)
        owner[inside & (owner < 0)] = i
    rows, cols = np.where(owner >= 0)
    r0, r1, c0, c1 = rows.min(), rows.max() + 1, cols.min(), cols.max() + 1
    acc = np.zeros(LA.shape)
    for y in range(years[0], years[1] + 1):
        sel = np.where((dates.year == y) & dates.month.isin(C.SEASON_MONTHS))[0]
        a = io.gridmet_read(ds, "pr", sel[0], sel[-1] + 1, r0, r1, c0, c1)
        acc[r0:r1, c0:c1] += np.nansum(a, axis=0)
    ds.close()
    acc /= (years[1] - years[0] + 1)
    n = len(polys4326)
    o = owner.ravel()
    k = o >= 0
    w = cell_area.ravel()[k]
    s = np.bincount(o[k], weights=acc.ravel()[k] * w, minlength=n)
    sw = np.bincount(o[k], weights=w, minlength=n)
    from pyproj import Geod
    geod = Geod(ellps="WGS84")
    full = np.array([abs(geod.geometry_area_perimeter(g)[0]) for g in polys4326])
    covered = sw * (1 / 24) ** 2 * (111320.0 ** 2)                 # approx area of the cells found (m2)
    with np.errstate(invalid="ignore", divide="ignore"):
        val = s / sw
    return np.where(covered >= 0.5 * full, val, np.nan)


def _rzws(polys5070):
    """Mean root-zone storage (mm) per polygon and the share of its area gSSURGO covers."""
    from ..common import io, cells as CL, zonal
    import geopandas as gpd
    ds = io.gssurgo()
    tr = io.gssurgo_transform(ds)
    v = ds["rootznaws"]
    g = gpd.GeoDataFrame(geometry=list(polys5070), crs=5070)
    acc = zonal.Accumulator(len(g), ["rz", "px"])
    for r0, r1, c0, c1, ids in CL.strips(g, tr, v.shape, strip_rows=2048):
        a = np.ma.filled(v[r0:r1, c0:c1].astype("float64"), np.nan)
        acc.add("rz", ids, a)                     # gSSURGO rootznaws is mm (Valu1)
        acc.add("px", ids)
    ds.close()
    area_px = np.array([geom.area for geom in g.geometry]) / (tr.a ** 2)
    return acc.mean("rz"), acc.cnt["rz"] / area_px


def run(ctx):
    import geopandas as gpd
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from .. import definitions as DEF
    from ..common import io, report, grids
    from ..phase0 import counties as CTY
    d12, d13, d14 = D.value("D12"), D.value("D13"), D.value("D14")
    cty = CTY.load(ctx.dep("counties"))

    # ══════ STEP 1.2-1  36 km EASE footprints around the study counties (D11) ══════
    fp = _footprint_grid(cty, d13["fallback_km"])
    ctx.log.info("STEP 1.2-1  %d candidate footprints", len(fp))

    # ══════ STEP 1.2-2  irrigated share 2011-2015 (largest yearly share) and label coverage (D12) ══════
    fp["irr_max"], fp["label_cov"] = _aim_shares(ctx, fp, list(range(d12["years"][0], d12["years"][1] + 1)))

    # ══════ STEP 1.2-3  how much of each footprint the SMAP file holds (nested 9 km cells with data) ══════
    p, smap_label = io.smap_source()
    ds = io.smap_open(p)
    r9, c9 = grids.check_ease_axes(ds["x"][:], ds["y"][:])
    roi = np.asarray(ds["roi_mask"][:]) > 0 if "roi_mask" in ds.variables else np.ones((len(r9), len(c9)), bool)
    ds.close()
    R9, C9 = np.meshgrid(r9, c9, indexing="ij")
    key = pd.Series(1, index=pd.MultiIndex.from_arrays([R9[roi] // 4, C9[roi] // 4])).groupby(level=[0, 1]).size()
    fp["n9_in_file"] = [int(key.get((r, c), 0)) for r, c in zip(fp.r36, fp.c36)]
    fp["smap_ok"] = DEF.footprint_observed(fp["n9_in_file"])

    # ══════ STEP 1.2-4  rainfall normal (May-Sep, D14 years) and root-zone storage, footprints and counties ══════
    ny = d14["normal_years"]
    fp4326 = fp.to_crs(4326)
    fp["rain_normal_mm"] = _rain_normal(ctx, list(fp4326.geometry), ny)
    fp["rzws_mm"], fp["rzws_cov"] = _rzws(fp.geometry)
    cty["rain_normal_mm"] = _rain_normal(ctx, list(cty.to_crs(4326).geometry), ny)
    cty["rzws_mm"], cty["rzws_cov"] = _rzws(cty.geometry)
    ctx.log.info("STEP 1.2-4  rainfall known for %d of %d footprints; soils for %d", int(np.isfinite(fp.rain_normal_mm).sum()),
                 len(fp), int((fp.rzws_cov >= 0.5).sum()))

    # ══════ STEP 1.2-5  reference test per county at 100 and 150 km; the radius rule (D13) ══════
    rows = []
    for c in cty.itertuples():
        dist = np.hypot(fp.x5070 - c.x5070, fp.y5070 - c.y5070).to_numpy() / 1000
        for radius in (d13["radius_km"], d13["fallback_km"]):
            ok, why = DEF.is_reference_footprint(fp.irr_max, fp.label_cov, dist, fp.rain_normal_mm / c.rain_normal_mm,
                                                 fp.rzws_mm / c.rzws_mm, (fp.rzws_cov >= 0.5) & (c.rzws_cov >= 0.5),
                                                 radius)
            sel = dist <= radius
            rows.append(pd.DataFrame(dict(county=c.county, radius_km=radius, fp=fp.fp[sel], dist_km=dist[sel],
                                          reference=ok[sel], reason=why[sel], smap_ok=fp.smap_ok[sel])))
    cf = pd.concat(rows, ignore_index=True)
    n_at = {r: cf[(cf.radius_km == r) & cf.reference].groupby("county").size().reindex(cty.county, fill_value=0).to_numpy()
            for r in (d13["radius_km"], d13["fallback_km"])}
    radius, share100 = DEF.footprint_radius(n_at)
    cf = cf[cf.radius_km == radius].copy()
    ctx.log.info("STEP 1.2-5  %.0f%% of counties reach %d footprints at %d km -> radius %d km for every county",
                 100 * share100, d13["min_footprints"], d13["radius_km"], radius)

    # ══════ STEP 1.2-6  equal weights over each county's reference footprints (D15); shared footprints ══════
    ref = cf[cf.reference].copy()
    ref["weight"] = 1.0 / ref.groupby("county").fp.transform("size")
    shared = ref.groupby("fp").county.agg(lambda x: " ".join(sorted(x))).reset_index()
    shared["n_counties"] = shared.county.str.split().str.len()
    shared = shared[shared.n_counties > 1]

    # ══════ STEP 1.2-7  check against 2018-2020 labels: is a reference footprint still irrigation-free? ══════
    # The rule itself stays defined from 2011-2015 (D12); the check is reported and moves no footprint.
    later = _lanid_later_shares(ctx, fp, C.FOOTPRINT_CHECK_YEARS)
    if later is None:
        fp["irr_max_2018_2020"], fp["label_cov_2018_2020"] = np.nan, np.nan
        fp["check_2018_2020"] = "pending: LANID 2018-2020 not on disk"
        check_summary = dict(status="pending: LANID 2018-2020 not on disk")
    else:
        fp["irr_max_2018_2020"], fp["label_cov_2018_2020"] = later
        still = DEF.irrigation_free(fp.irr_max_2018_2020, fp.label_cov_2018_2020)
        fp["check_2018_2020"] = np.where(still, "irrigation-free in 2018-2020",
                                         np.where(np.isfinite(fp.irr_max_2018_2020),
                                                  "irrigated share at or above the limit",
                                                  "not covered by LANID 2018-2020"))
        is_ref = fp.fp.isin(ref.fp).to_numpy()
        lost = set(fp.fp[is_ref & ~still])
        per_now = ref.groupby("county").size().reindex(cty.county, fill_value=0)
        per_then = ref[~ref.fp.isin(lost)].groupby("county").size().reindex(cty.county, fill_value=0)
        check_summary = dict(status="done", map="LANID 2018-2020 (CONUS release)",
                             reference_footprints=int(is_ref.sum()), fail_2018_2020=int(len(lost)),
                             counties_with_min_footprints=int((per_now >= d13["min_footprints"]).sum()),
                             counties_with_min_footprints_after_check=int((per_then >= d13["min_footprints"]).sum()))
        ctx.log.info("STEP 1.2-7  %d of %d reference footprints fail the 2018-2020 check", len(lost), int(is_ref.sum()))

    # ══════ STEP 1.2-8  write ══════
    fp["smap_source"] = smap_label
    fp["is_reference_any"] = fp.fp.isin(ref.fp)
    report.table(ctx, pd.DataFrame(fp.drop(columns="geometry")), "footprints", parquet=True,
                 units=dict(irr_max="share", irr_max_2018_2020="share", rain_normal_mm="mm (May-Sep)", rzws_mm="mm"),
                 notes=f"SMAP file: {smap_label}; irrigation maps sampled every {DEC} pixels; irr_max from "
                       f"{d12['map']} 2011-2015 (the rule), irr_max_2018_2020 from LANID (the check: a different "
                       f"product, so a difference can be change over time or a difference between the maps)")
    report.table(ctx, cf, "county_footprints", parquet=True, units={"dist_km": "km"},
                 notes=f"radius {radius} km (D13; share of counties with >=3 at 100 km = {share100:.2f})")
    report.table(ctx, ref[["county", "fp", "dist_km", "weight", "smap_ok"]], "county_reference_weights", parquet=True)
    report.table(ctx, shared, "shared_footprints")
    report.table(ctx, pd.DataFrame(cty.drop(columns="geometry"))[["county", "name", "rain_normal_mm", "rzws_mm",
                                                                    "rzws_cov"]], "county_normals",
                 units=dict(rain_normal_mm="mm (May-Sep)", rzws_mm="mm"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 6))
    fp.plot(ax=ax, column="irr_max", cmap="Blues", vmin=0, vmax=0.3, edgecolor="#bbbbbb", linewidth=0.3, legend=True,
            legend_kwds=dict(label="largest irrigated share 2011-2015"))
    fp[fp.fp.isin(ref.fp)].boundary.plot(ax=ax, color=C.COLORS["footprint"], linewidth=1.2)
    cty.boundary.plot(ax=ax, color="k", linewidth=0.4)
    ax.set_axis_off()
    ax.set_title(f"Reference footprints (outlined), radius {radius} km; SMAP: {smap_label}", fontsize=9)
    report.figure(ctx, fig, "footprints_map", units={"share": "fraction"}, mask="study counties")
    per = ref.groupby("county").size().reindex(cty.county, fill_value=0)
    ctx.complete(dict(radius_km=radius, share_counties_3_at_100km=share100, footprints=int(len(fp)),
                      reference=int(ref.fp.nunique()), shared=int(len(shared)), smap_source=smap_label,
                      counties_with_3=float((per >= d13["min_footprints"]).mean()),
                      check_2018_2020=check_summary))
