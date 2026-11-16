"""Replication region screen (Stage 4): metadata only.

"One region outside Kansas. Candidates: Texas High Plains, southwest Nebraska. Screen metadata only, to choose
the region. Its ET and surface temperature are not extracted."
Physical idea: the replication must offer the same ingredients - droughts, irrigated cropland on wells, a range
of aquifer thickness, measured wells, irrigation-free footprints - or it cannot confirm anything.
Decisions: D01, D12, D20 (applied as in Kansas).   Deliverables: rep_screen_county.csv, rep_screen.csv, rep_screen.md.
No ET, no surface temperature: this module never opens an OpenET or Landsat file.
"""
import numpy as np


def run(ctx):
    import pandas as pd
    import netCDF4 as nc
    from .. import config as C
    from .. import decisions as D
    from .. import definitions as DEF
    from ..common import io, report
    from ..phase0 import counties as CTY
    from .t1_1_weather import county_weights
    from .t1_2_footprints import _footprint_grid, _aim_shares

    allc = CTY.load(ctx.dep("counties"), study_only=False)
    rep = allc[allc.region.isin(C.REP_REGIONS)].copy()

    # ══════ STEP R-1  weather: May-September deficit 2016-2024 and dry seasons (D01, D20 within the region) ══════
    var = D.value("D01")["variable"]
    ds = io.gridmet()
    lat, lon, dates = io.gridmet_axes(ds)
    W = county_weights(rep.to_crs(4326), lat, lon)
    rows = []
    for y in C.PANEL_YEARS:
        sel = np.where((dates.year == y) & dates.month.isin(C.SEASON_MONTHS))[0]
        for c, (R, Cc, w) in W.items():
            r0, r1, c0, c1 = R.min(), R.max() + 1, Cc.min(), Cc.max() + 1
            e = io.gridmet_read(ds, var, sel[0], sel[-1] + 1, r0, r1, c0, c1)[:, R - r0, Cc - c0] @ w
            p = io.gridmet_read(ds, "pr", sel[0], sel[-1] + 1, r0, r1, c0, c1)[:, R - r0, Cc - c0] @ w
            rows.append(dict(county=c, year=y, deficit_season_mm=float(np.sum(e - p))))
    ds.close()
    wx = pd.DataFrame(rows).merge(rep[["county", "region"]], on="county")
    wx["dry_season"] = False
    for reg, g in wx.groupby("region"):
        wx.loc[g.index, "dry_season"] = DEF.season_dry(g.deficit_season_mm.to_numpy())[0]

    # ══════ STEP R-2  irrigation (AIM-HPA 2015) and irrigation-free footprints within 100 km (D12) ══════
    fp = _footprint_grid(rep, 100)
    fp["irr_max"], fp["label_cov"] = _aim_shares(ctx, fp, list(range(2011, 2016)))
    clean = fp[DEF.irrigation_free(fp.irr_max, fp.label_cov)]
    rep["fp_clean_100km"] = [int((np.hypot(clean.x5070 - c.x5070, clean.y5070 - c.y5070) / 1000 <= 100).sum())
                             for c in rep.itertuples()]

    # ══════ STEP R-3  wells measured January-April per year (combined well file) ══════
    d = nc.Dataset(io.path("wells_ts"))
    sidx = np.asarray(d["station_index"][:])
    yr, mo, known = io.well_dates(np.ma.filled(d["time"][:].astype(float), np.nan))
    lonw, latw = np.asarray(d["lon"][:], float), np.asarray(d["lat"][:], float)
    d.close()
    ok = known & np.isin(mo, C.PRESEASON_MONTHS) & np.isin(yr, C.PANEL_YEARS)
    import geopandas as gpd
    sp = gpd.GeoDataFrame(dict(si=np.arange(len(lonw))), geometry=gpd.points_from_xy(lonw, latw), crs=4269).to_crs(5070)
    j = gpd.sjoin(sp, rep[["county", "geometry"]], predicate="within")
    st2c = dict(zip(j.si, j.county))
    obs = pd.DataFrame(dict(si=sidx[ok], year=yr[ok]))
    obs["county"] = obs.si.map(st2c)
    wells = obs.dropna().drop_duplicates(["si", "year"]).groupby(["county", "year"]).size().rename("wells").reset_index()

    # ══════ STEP R-4  aquifer thickness available (Haacker 2015 water table and base at county centres) ══════
    hb = io.sample_raster(f"{io.path('haacker_dir')}/AqBase_m.tif", rep.x5070, rep.y5070, "EPSG:5070")
    hw = io.sample_raster(f"{io.path('haacker_dir')}/WTelev_Clip_2015.tif", rep.x5070, rep.y5070, "EPSG:5070")
    rep["haacker_thickness_2015_m"] = hw - hb

    # ══════ STEP R-5  write: per county and per region ══════
    per = (rep[["county", "name", "region", "fp_clean_100km", "haacker_thickness_2015_m"]]
           .merge(wx.groupby("county").agg(dry_seasons=("dry_season", "sum"),
                                           deficit_mean_mm=("deficit_season_mm", "mean")).reset_index(), on="county")
           .merge(wells.groupby("county").wells.mean().rename("wells_per_preseason").reset_index(), on="county",
                  how="left"))
    reg = (per.groupby("region").agg(counties=("county", "size"), dry_county_years=("dry_seasons", "sum"),
                                      deficit_mean_mm=("deficit_mean_mm", "mean"),
                                      counties_with_3_clean_footprints=("fp_clean_100km", lambda s: int((s >= 3).sum())),
                                      wells_per_preseason=("wells_per_preseason", "sum"),
                                      thickness_p10_m=("haacker_thickness_2015_m", lambda s: float(np.nanpercentile(s, 10))),
                                      thickness_p90_m=("haacker_thickness_2015_m", lambda s: float(np.nanpercentile(s, 90))))
           .reset_index())
    report.table(ctx, per, "rep_screen_county", units=dict(deficit_mean_mm="mm", haacker_thickness_2015_m="m"),
                 notes="metadata only; no ET or surface temperature")
    report.table(ctx, reg, "rep_screen", units=dict(deficit_mean_mm="mm", thickness_p10_m="m"))
    report.text(ctx, "rep_screen.md", "# Replication region screen (metadata only)\n\n" + report.md(reg)
                + "\n\nThe region stays reserved; this table is for choosing it, recorded in the decision log.")
    ctx.complete(dict(regions=reg.set_index("region").to_dict("index")))
