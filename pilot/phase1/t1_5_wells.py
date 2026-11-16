"""Task 1.5  Well inventory.

"Combine KGS WIZARD water levels with the existing USGS records. Check aquifer code, units, datum, static or
pumping status, and duplicates. Map the wells measured in each pre-season and their distance from cohort cells.
Interpolation comes later, in Check 2.4."
Physical idea: aquifer thickness - how much saturated aquifer a well can draw on - is known only where someone
measured the water level after the pumps rested. Pre-season (January-April) recovered levels are those
measurements; their density near the cohort decides how well thickness can be mapped.
Decisions: D22.   Deliverables: well_year.csv (Paper 2 Table 1 source), wells_by_source_year.csv,
cell_well_distance (parquet), wells_support.png.
"""
import numpy as np

FT = 0.3048


def _kgs(C):
    """KGS winter (pre-season) network -> one row per well-year."""
    import pandas as pd
    from ..common import io
    w = io.table("kgs_winter")
    years = [y for y in C.WELL_YEARS if f"winter_1yr_{y}" in w.columns]
    base = pd.DataFrame(dict(well_id="KGS:" + w["USGS_ID"].astype(str), usgs_site=w["USGS_ID"].astype(str),
                             lon=w["LONGITUDE"], lat=w["LATITUDE"], aquifer=w["AQUIFER_CODE"].fillna("").astype(str),
                             gmd=w["GROUNDWATER_MGMT_DISTRICT"], well_depth_ft=w["DEPTH_OF_WELL"],
                             reported_land_ft=w["LAND_SURFACE_ALTITUDE"], datum=w["ALTITUDE_DATUM"].fillna(""),
                             depth_to_bedrock_ft=w["DEPTH_TO_BEDROCK"]))
    rows = []
    for y in years:
        d = base.copy()
        d["year"] = y
        d["dtw_ft"] = pd.to_numeric(w[f"winter_1yr_{y}"], errors="coerce")
        rows.append(d[d.dtw_ft.notna()])
    out = pd.concat(rows, ignore_index=True)
    out["source"], out["month"], out["status"] = "KGS winter", np.nan, "winter (recovered)"
    return out


def _usgs(C, d22, kgs_sites):
    """USGS/NWIS records of Kansas wells: January-April, status Static (D22), one median per well-year."""
    import pandas as pd
    import netCDF4 as nc
    from ..common import io
    ds = nc.Dataset(io.path("wells_ts"))
    st = np.array([str(x) for x in ds["state"][:]])
    ks = np.where(st == "Kansas")[0]
    sidx = np.asarray(ds["station_index"][:])
    m = np.isin(sidx, ks)
    year, month, ok_t = io.well_dates(np.ma.filled(ds["time"][:].astype(float), np.nan)[m])
    sel_idx = np.where(m)[0]
    stext = [str(x) for x in ds["status_text"][:]]
    status = np.asarray(ds["status"][:])[m]
    dtw = np.ma.filled(ds["depth_to_water"][:][m].astype(float), np.nan)
    df = pd.DataFrame(dict(si=sidx[sel_idx], year=year, month=month, dtw_ft=dtw,
                           status=[stext[k] if 0 <= k < len(stext) else "" for k in status]))
    n_all = len(df)
    df = df[ok_t & np.isin(month, d22["usgs_months"]) & np.isfinite(dtw)]
    n_season = len(df)
    df = df[df.status.str.split(";").str[0].str.strip().isin(d22["usgs_status"])]
    st_vars = dict(usgs_site=np.array([str(x) for x in ds["usgs_site_no"][:]]), lon=np.asarray(ds["lon"][:], float),
                   lat=np.asarray(ds["lat"][:], float), aquifer=np.array([str(x) for x in ds["aquifer"][:]]),
                   well_depth_ft=np.ma.filled(ds["well_depth"][:].astype(float), np.nan))
    ds.close()
    for k, v in st_vars.items():
        df[k] = v[df.si.to_numpy()]
    df = df[df.year.isin(C.WELL_YEARS)]
    dup = df.usgs_site.isin(kgs_sites)
    n_dup = int(dup.sum())
    df = df[~dup]
    g = (df.groupby(["usgs_site", "year"])
         .agg(dtw_ft=("dtw_ft", "median"), month=("month", "min"), lon=("lon", "first"), lat=("lat", "first"),
              aquifer=("aquifer", "first"), well_depth_ft=("well_depth_ft", "first"), status=("status", "first"),
              n_meas=("dtw_ft", "size")).reset_index())
    g["well_id"] = "USGS:" + g.usgs_site
    g["source"] = "USGS Jan-Apr static"
    return g, dict(usgs_records_ks=n_all, jan_apr=n_season, static=int(len(df) + n_dup), duplicates_of_kgs=n_dup)


def run(ctx):
    import geopandas as gpd
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from ..common import io, report, stats
    from ..phase0 import counties as CTY, t0_2_cells as T02
    d22 = D.value("D22")

    # ══════ STEP 1.5-1  KGS winter levels (recovered, pre-season) ══════
    kg = _kgs(C)
    ctx.log.info("STEP 1.5-1  KGS: %d well-years, %d wells", len(kg), kg.well_id.nunique())

    # ══════ STEP 1.5-2  USGS January-April static records; KGS wins duplicates ══════
    us, ucount = _usgs(C, d22, set(kg.usgs_site))
    ctx.log.info("STEP 1.5-2  USGS: %d well-years after the D22 rules (%s)", len(us), ucount)
    wy = pd.concat([kg, us], ignore_index=True)

    # ══════ STEP 1.5-3  checks: aquifer code (D22), location, units, sign; one value per well-year ══════
    pts = gpd.GeoDataFrame(wy, geometry=gpd.points_from_xy(wy.lon, wy.lat), crs=4269).to_crs(5070)
    ext = io.vector("kgs_extent").to_crs(5070)
    in_ext = gpd.sjoin(pts[["geometry"]].reset_index(), ext[["DESCRIPTIO", "geometry"]], predicate="within",
                       how="left").drop_duplicates("index").set_index("index").sort_index()
    wy["extent_class"] = in_ext["DESCRIPTIO"].reindex(wy.index).fillna("outside").to_numpy()
    code = wy.aquifer.str.strip()
    hpa_code = code.str.startswith(tuple(d22["keep_code_prefixes"]))
    blank = code.isin(["", "nan", "None"])
    wy["code_ok"] = hpa_code | (blank & ~wy.extent_class.isin(["outside", "outcrop"]))
    wy["x5070"], wy["y5070"] = pts.geometry.x.to_numpy(), pts.geometry.y.to_numpy()
    wy["dtw_m"] = wy.dtw_ft * FT
    wy["flag"] = np.where(~wy.code_ok, "aquifer code not HPA", np.where(wy.dtw_m < 0, "above land surface",
                          np.where(~np.isfinite(wy.lon), "no location", "")))
    # ══════ STEP 1.5-4  land surface from one source (GLO30 DEM, D22); water-table elevation ══════
    sites = wy.drop_duplicates("well_id")[["well_id", "x5070", "y5070"]]
    sites["land_m"] = io.sample_raster("dem", sites.x5070, sites.y5070, "EPSG:5070")
    wy = wy.merge(sites[["well_id", "land_m"]], on="well_id", how="left")
    wy["wte_m"] = wy.land_m - wy.dtw_m
    wy["reported_land_m"] = wy.reported_land_ft * FT
    wy.loc[wy.flag.eq("") & ~np.isfinite(wy.land_m), "flag"] = "no DEM"
    wy["use"] = wy.flag.eq("")
    datum_check = (wy[wy.use & np.isfinite(wy.reported_land_m)].assign(d=lambda x: x.land_m - x.reported_land_m)
                   .groupby("datum").d.agg(["median", "count"]).reset_index())
    ctx.log.info("STEP 1.5-4  %d usable well-years (%d wells)", int(wy.use.sum()), wy[wy.use].well_id.nunique())

    # ══════ STEP 1.5-5  support: nearest measured well of each pre-season from every cohort cell ══════
    cells = T02.load(ctx.dep("t0_2"), columns=["cell_id", "county", "role", "x5070", "y5070"])
    co = cells[cells.role == "cohort"]
    dist = []
    for y in C.PANEL_YEARS:
        w = wy[wy.use & (wy.year == y)]
        dist.append(pd.DataFrame(dict(cell_id=co.cell_id, county=co.county, year=y,
                                      well_km=stats.nearest_km(np.c_[co.x5070, co.y5070], np.c_[w.x5070, w.y5070]))))
    dist = pd.concat(dist, ignore_index=True)
    share = dist.groupby("year").well_km.apply(lambda s: float((s <= C.TARGETS["well_km"]).mean()))
    ctx.log.info("STEP 1.5-5  share of cohort cells within %g km of a measured well: %s", C.TARGETS["well_km"],
                 share.round(2).to_dict())

    # ══════ STEP 1.5-6  write ══════
    keep = ["well_id", "source", "year", "month", "lon", "lat", "x5070", "y5070", "aquifer", "extent_class",
            "well_depth_ft", "status", "dtw_ft", "dtw_m", "land_m", "wte_m", "reported_land_m", "datum", "gmd",
            "depth_to_bedrock_ft", "flag", "use"]
    report.table(ctx, wy[keep], "well_year", parquet=True,
                 units=dict(dtw_m="m below land", wte_m="m (GLO30 datum)", dtw_ft="ft"),
                 notes=f"D22 rules; USGS counts {ucount}")
    t1 = wy.groupby(["source", "year"]).agg(well_years=("well_id", "size"), usable=("use", "sum")).reset_index()
    report.table(ctx, t1, "wells_by_source_year", units={"well_years": "count"})
    report.table(ctx, datum_check, "land_surface_vs_reported", units={"median": "m (DEM - reported)"})
    report.table(ctx, dist, "cell_well_distance", parquet=True, csv=False, units={"well_km": "km"})
    report.table(ctx, share.rename("share_within_8km").reset_index(), "cohort_well_support",
                 units={"share_within_8km": "fraction"}, mask="cohort cells")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cty = CTY.load(ctx.dep("counties"))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), gridspec_kw=dict(width_ratios=[1.6, 1]))
    cty.boundary.plot(ax=axes[0], color="#888888", linewidth=0.4)
    yr = C.PANEL_YEARS[len(C.PANEL_YEARS) // 2]
    w = wy[wy.use & (wy.year == yr)]
    axes[0].scatter(w.x5070, w.y5070, s=3, c=np.where(w.source.str.startswith("KGS"), "#1d4ed8", "#c2410c"))
    axes[0].set_title(f"Wells measured in the {yr} pre-season (blue KGS, orange USGS)", fontsize=9)
    axes[0].set_axis_off()
    t1p = t1.pivot(index="year", columns="source", values="usable").fillna(0)
    t1p.plot(kind="bar", stacked=True, ax=axes[1], color=["#1d4ed8", "#c2410c"][:t1p.shape[1]])
    axes[1].set_ylabel("usable well-years")
    report.figure(ctx, fig, "wells_support", units={"wells": "count"}, mask="usable wells (D22)")
    ctx.complete(dict(well_years=int(len(wy)), usable=int(wy.use.sum()), wells=int(wy[wy.use].well_id.nunique()),
                      by_flag=wy.flag.value_counts().to_dict(), usgs=ucount,
                      cohort_share_within_8km={int(k): round(v, 3) for k, v in share.items()}))
