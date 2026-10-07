"""
Build one CF NetCDF with every depth-to-water measurement found for wells of the High Plains aquifer.

Sources (downloaded unchanged beforehand):
  USGS   Water Data API, national aquifer code N100HGHPLN, parameter 72019   (download_usgs_hpa_water_levels.py)
         This already carries the measurements that state agencies share with USGS.
  TWDB   Texas Water Development Board Groundwater Database, GWDBDownload.zip (wells whose aquifer names Ogallala)
  NECSD  University of Nebraska Conservation and Survey Division, WLDB.zip    (wells inside the aquifer polygon;
         the database has no aquifer field)

Nothing is averaged or filtered by quality here. Every record keeps its source's status text. Two flags help
later use:
  well_group                 wells of different sources within 25 m of each other share a group number
  duplicate_of_other_source  a record that repeats (same group, same day, depth within 0.05 ft) a record of a
                             source listed earlier (USGS, then TWDB, then NECSD)
"""
import gzip
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import netCDF4 as nc4
import numpy as np
import pandas as pd
import shapely
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")

RAW_USGS = Path(r"G:/MSU_GWB/_usgs_wl_raw")
RAW_STATE = Path(r"G:/MSU_GWB/_state_gw_raw")
OUT = Path(r"G:/MSU_GWB/datasets/merged_datasets/HPA_wells_water_level_timeseries.nc")
BOUND = "zip://G:/MSU_GWB/datasets/high_plains_quifer.zip!high_plains_quifer/hp_bound2010.shp"
MATCH_M = 25.0
EPOCH = pd.Timestamp("1900-01-01")
COMP = dict(zlib=True, complevel=4, shuffle=True)


def log(m):
    print(f"[hpa-timeseries] {m}", flush=True)


roi = gpd.read_file(BOUND).to_crs(4326).geometry.union_all()
shapely.prepare(roi)


def read_pages(stem):
    out = []
    for p in sorted(RAW_USGS.glob(f"{stem}_p*.json.gz")):
        with gzip.open(p, "rt", encoding="utf-8") as f:
            out.extend(json.load(f))
    return out


# ───────────────────────── USGS ─────────────────────────
feats = read_pages("wells")
uw = pd.DataFrame([dict(f["properties"], lon=(f.get("geometry") or {}).get("coordinates", [np.nan, np.nan])[0],
                        lat=(f.get("geometry") or {}).get("coordinates", [np.nan, np.nan])[1]) for f in feats])
uw = uw.rename(columns={"id": "site"})
rows = []
for p in sorted(RAW_USGS.glob("fm_*_p*.json.gz")):
    with gzip.open(p, "rt", encoding="utf-8") as f:
        for ft in json.load(f):
            q = ft["properties"]
            rows.append((q["monitoring_location_id"], q["time"], q["value"], q.get("unit_of_measure"),
                         "; ".join(q["qualifier"]) if isinstance(q.get("qualifier"), list) else (q.get("qualifier") or ""),
                         q.get("approval_status") or "", q.get("measuring_agency") or "",
                         q.get("observing_procedure") or ""))
um = pd.DataFrame(rows, columns=["site", "time", "value", "unit", "status", "approval", "agency", "method"])
um["depth_ft"] = pd.to_numeric(um["value"], errors="coerce")
um["date"] = pd.to_datetime(um["time"], errors="coerce", utc=True).dt.tz_localize(None).dt.floor("D")
log(f"USGS: {len(uw)} wells listed, {len(um)} measurements; units {um.unit.value_counts().to_dict()}; "
    f"unparseable value {int(um.depth_ft.isna().sum())}, unparseable date {int(um.date.isna().sum())}")
um = um[um.depth_ft.notna() & um.date.notna() & (um.unit == "ft")]
uw = uw[uw.site.isin(um.site.unique()) & uw.lat.notna()]
usgs_w = pd.DataFrame(dict(
    well_id="USGS:" + uw.site.astype(str), source="USGS", lat=uw.lat.astype(float), lon=uw.lon.astype(float),
    well_depth_ft=pd.to_numeric(uw.well_constructed_depth, errors="coerce"), state=uw.state_name.astype(str),
    county=uw.county_name.astype(str), name=uw.monitoring_location_name.astype(str),
    agency=uw.agency_code.astype(str), aquifer=uw.aquifer_code.fillna("").astype(str),
    usgs_site_no=uw.monitoring_location_number.astype(str), key=uw.site.astype(str)))
usgs_m = pd.DataFrame(dict(key=um.site, date=um.date, year=um.date.dt.year, depth_ft=um.depth_ft, status=um.status,
                           approval=um.approval, agency=um.agency, method=um.method, source="USGS"))

# ───────────────────────── Texas (TWDB) ─────────────────────────
tx = RAW_STATE / "TX" / "GWDBDownload"
wm = pd.read_csv(tx / "WellMain.txt", sep="|", dtype=str, encoding="latin-1", on_bad_lines="skip")
lv = pd.concat([pd.read_csv(tx / f, sep="|", dtype=str, encoding="latin-1", on_bad_lines="skip")
                for f in ("WaterLevelsMajor.txt", "WaterLevelsMinor.txt", "WaterLevelsCombination.txt",
                          "WaterLevelsOtherUnassigned.txt")], ignore_index=True)
wm = wm[wm.Aquifer.fillna("").str.contains("Ogallala", case=False)]
wm["lat"], wm["lon"] = pd.to_numeric(wm.LatitudeDD, errors="coerce"), pd.to_numeric(wm.LongitudeDD, errors="coerce")
wm = wm[wm.lat.notna() & wm.lon.notna()]
lv = lv[lv.StateWellNumber.isin(wm.StateWellNumber)]
lv["depth_ft"] = pd.to_numeric(lv.DepthFromLSD, errors="coerce")
lv["date"] = pd.to_datetime(lv.MeasurementDate, errors="coerce")
log(f"TWDB: {len(wm)} Ogallala wells, {len(lv)} measurements; without a full date {int(lv.date.isna().sum())}, "
    f"without a depth {int(lv.depth_ft.isna().sum())}; status {lv.Status.value_counts().to_dict()}")
lv["year"] = pd.to_numeric(lv.MeasurementYear, errors="coerce")
lv = lv[lv.depth_ft.notna() & lv.year.notna()]
wm = wm[wm.StateWellNumber.isin(lv.StateWellNumber.unique())]
tx_w = pd.DataFrame(dict(
    well_id="TWDB:" + wm.StateWellNumber, source="TWDB", lat=wm.lat, lon=wm.lon,
    well_depth_ft=pd.to_numeric(wm.WellDepth, errors="coerce"), state="Texas", county=wm.County.fillna(""),
    name=wm.Owner.fillna(""), agency=wm.ReportingAgency.fillna(""), aquifer=wm.Aquifer.fillna(""),
    usgs_site_no=wm.USGSSiteNumber.fillna(""), key=wm.StateWellNumber))
tx_m = pd.DataFrame(dict(key=lv.StateWellNumber, date=lv.date, year=lv.year, depth_ft=lv.depth_ft,
                         status=(lv.Status.fillna("") + "; " + lv.Remarks.fillna("").str.strip()).str.strip("; "),
                         approval="", agency=lv.MeasuringAgency.fillna(""), method=lv.MethodOfMeasurement.fillna(""),
                         source="TWDB"))

# ───────────────────────── Nebraska (UNL CSD) ─────────────────────────
ne = RAW_STATE / "NE"
nw = pd.read_csv(ne / "Well_Info.csv", dtype=str)
nl = pd.read_csv(ne / "Water_Level_Data.csv", dtype=str)
nw["lat"], nw["lon"] = pd.to_numeric(nw.LatDD, errors="coerce"), pd.to_numeric(nw.LongDD, errors="coerce")
nw.loc[nw.lon > 0, "lon"] *= -1
nw = nw[nw.lat.notna() & nw.lon.notna()]
nl["depth_ft"] = pd.to_numeric(nl.WatLevel, errors="coerce")
nl["date"] = pd.to_datetime(nl.DateMsr, errors="coerce")
log(f"NECSD: {len(nw)} wells with coordinates, {len(nl)} measurements; without a date {int(nl.date.isna().sum())}, "
    f"without a depth {int(nl.depth_ft.isna().sum())}; seasons {nl.Season.value_counts().head(6).to_dict()}")
nl["year"] = pd.to_numeric(nl.YearMsr, errors="coerce")
nl.loc[nl.date.notna(), "year"] = nl.date.dt.year[nl.date.notna()]
nl = nl[nl.depth_ft.notna() & nl.year.notna() & nl.CSD_ID.isin(nw.CSD_ID)]
nw = nw[nw.CSD_ID.isin(nl.CSD_ID.unique())]
ne_w = pd.DataFrame(dict(
    well_id="NECSD:" + nw.CSD_ID, source="NECSD", lat=nw.lat, lon=nw.lon,
    well_depth_ft=pd.to_numeric(nw.Well_Dpth, errors="coerce"), state="Nebraska", county=nw.County.fillna(""),
    name=nw.Legal.fillna(""), agency=nw.Agency.fillna(""), aquifer="", usgs_site_no="", key=nw.CSD_ID))
ne_m = pd.DataFrame(dict(key=nl.CSD_ID, date=nl.date, year=nl.year, depth_ft=nl.depth_ft,
                         status=nl.Season.fillna("").str.strip().str.capitalize(),
                         approval="", agency="", method="", source="NECSD"))

# ───────────────────────── combine ─────────────────────────
wells = pd.concat([usgs_w, tx_w, ne_w], ignore_index=True)
wells["in_polygon"] = shapely.contains_xy(roi, wells.lon.values, wells.lat.values).astype("i1")
before = wells.groupby("source").size().to_dict()
wells = wells[(wells.in_polygon == 1) | (wells.source != "NECSD")].reset_index(drop=True)   # NE has no aquifer field
log(f"wells with measurements by source {before}; kept (Nebraska database limited to the polygon) "
    f"{wells.groupby('source').size().to_dict()}; inside the polygon {wells.groupby('source').in_polygon.sum().to_dict()}")
wells["station"] = np.arange(len(wells))
meas = pd.concat([usgs_m, tx_m, ne_m], ignore_index=True)
meas = meas.merge(wells[["source", "key", "station"]], on=["source", "key"], how="inner")
n0 = len(meas)
meas = meas.drop_duplicates(["station", "date", "year", "depth_ft", "status"])
log(f"measurements {n0}; exact repeats inside a source removed {n0 - len(meas)}")

# wells of different sources within MATCH_M metres -> one group
pts = gpd.GeoSeries(gpd.points_from_xy(wells.lon, wells.lat), crs=4326).to_crs(5070)
xy = np.column_stack([pts.x.values, pts.y.values])
tree = cKDTree(xy)
parent = np.arange(len(wells))


def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


src = wells.source.values
for i, j in tree.query_pairs(MATCH_M):
    if src[i] != src[j]:
        parent[find(j)] = find(i)
wells["well_group"] = pd.factorize(np.array([find(i) for i in range(len(wells))]))[0]
rank = {"USGS": 0, "TWDB": 1, "NECSD": 2}
meas["group"] = wells.well_group.values[meas.station.values]
meas["rank"] = meas.source.map(rank)
meas = meas.sort_values(["group", "year", "date", "rank"]).reset_index(drop=True)
same = (meas.group.values[1:] == meas.group.values[:-1]) & (meas.date.values[1:] == meas.date.values[:-1]) \
       & meas.date.notna().values[1:] \
       & (np.abs(meas.depth_ft.values[1:] - meas.depth_ft.values[:-1]) <= 0.05) \
       & (meas["rank"].values[1:] > meas["rank"].values[:-1])
meas["dup"] = np.concatenate([[False], same]).astype("i1")
shared = wells.groupby("well_group").source.nunique()
log(f"well groups {wells.well_group.nunique()} from {len(wells)} wells; groups seen by more than one source "
    f"{int((shared > 1).sum())}; records repeating an earlier source {int(meas.dup.sum())} "
    f"({meas[meas.dup == 1].source.value_counts().to_dict()})")
meas = meas.sort_values(["station", "year", "date"]).reset_index(drop=True)

# ───────────────────────── write ─────────────────────────
OUT.parent.mkdir(parents=True, exist_ok=True)
now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
with nc4.Dataset(OUT, "w", format="NETCDF4") as ds:
    ds.setncatts(dict(
        Conventions="CF-1.9", featureType="timeSeries",
        title="High Plains aquifer wells: depth-to-water measurements, whole record, from USGS, Texas and Nebraska",
        summary="Every depth-to-water field measurement obtained for wells of the High Plains aquifer from the USGS "
                "Water Data API (national aquifer code N100HGHPLN, parameter 72019), the Texas Water Development "
                "Board Groundwater Database (wells whose aquifer names Ogallala) and the University of Nebraska "
                "Conservation and Survey Division water-level database (wells inside hp_bound2010). Values are "
                "copied unchanged, in feet below land surface. No quality filter is applied: read status.",
        institution="USGS, TWDB, UNL Conservation and Survey Division (source data); packaged by the MSU "
                    "groundwater project",
        source="https://api.waterdata.usgs.gov/ogcapi/v0 ; https://www.twdb.texas.gov/groundwater/data/"
               "GWDBDownload.zip ; https://csdportal1.unl.edu/webdata/DB/WLDB.zip",
        partitions="The file holds three sources one after another on both axes. The source axis names each one and gives its first row and row count on the station axis and on the obs axis, so a source can be read as one slice.",
        comment="The same well can appear under two sources. well_group joins wells of different sources that lie "
                f"within {MATCH_M:.0f} m; duplicate_of_other_source marks records that repeat one from a source "
                "listed earlier. The Nebraska database has no aquifer field, so its wells are selected by position "
                "only. Kansas (WIZARD), Colorado, Oklahoma, New Mexico, Wyoming and South Dakota state databases "
                "are not included; their measurements are present only as far as the agencies share them with USGS.",
        history=f"{now} created by build_hpa_water_level_timeseries.py",
        date_downloaded="2026-10-07"))
    ds.createDimension("station", len(wells))
    ds.createDimension("obs", len(meas))

    def svar(name, values, long_name, dim, **a):
        v = ds.createVariable(name, str, (dim,))
        v[:] = np.array(["" if x is None or (isinstance(x, float) and np.isnan(x)) else str(x) for x in values],
                        dtype=object)
        v.setncatts(dict(long_name=long_name, **a))

    def nvar(name, values, dtype, long_name, dim, fill=None, **a):
        v = ds.createVariable(name, dtype, (dim,), fill_value=fill, **COMP)
        arr = np.asarray(values, dtype="f8" if np.dtype(dtype).kind == "f" else dtype)
        if fill is not None and np.dtype(dtype).kind == "f":
            arr = np.where(np.isnan(arr), fill, arr)
        v[:] = arr.astype(dtype)
        v.setncatts(dict(long_name=long_name, **a))

    svar("well_id", wells.well_id, "Well identifier: source, colon, the source's own well number", "station",
         cf_role="timeseries_id")
    # ---- one row per source: what it is, where it came from, and which rows of the two axes it occupies ----
    SRC = ["USGS", "TWDB", "NECSD"]
    info = dict(
        USGS=("U.S. Geological Survey Water Data API (National Water Information System)",
              "https://api.waterdata.usgs.gov/ogcapi/v0/collections/field-measurements",
              "wells with national aquifer code N100HGHPLN, parameter 72019 (depth to water below land surface); "
              "includes measurements that state agencies share with USGS",
              "U.S. Geological Survey, National Water Information System, https://doi.org/10.5066/F7P55KJN"),
        TWDB=("Texas Water Development Board Groundwater Database",
              "https://www.twdb.texas.gov/groundwater/data/GWDBDownload.zip",
              "wells whose aquifer name contains Ogallala; all four WaterLevels tables",
              "Texas Water Development Board, Groundwater Database (GWDB) Reports"),
        NECSD=("University of Nebraska-Lincoln Conservation and Survey Division, Nebraska Statewide Groundwater-Level "
               "Monitoring Program database", "https://csdportal1.unl.edu/webdata/DB/WLDB.zip",
               "wells inside the aquifer boundary hp_bound2010 (the database has no aquifer field); most records "
               "carry a year and a season but no date",
               "Conservation and Survey Division, School of Natural Resources, University of Nebraska-Lincoln"))
    ds.createDimension("source", len(SRC))
    st_src, ob_src = wells.source.values, meas.source.values
    for arr, what in ((st_src, "stations"), (ob_src, "records")):
        pos = [np.flatnonzero(arr == k) for k in SRC]
        if any(q.size and (q[-1] - q[0] + 1 != q.size) for q in pos):
            raise RuntimeError(f"{what} of a source are not contiguous")
    svar("source_name", SRC, "Short name of the source database", "source")
    svar("source_long_name", [info[k][0] for k in SRC], "Full name of the source database", "source")
    svar("source_url", [info[k][1] for k in SRC], "Address the data were downloaded from", "source")
    svar("source_selection", [info[k][2] for k in SRC], "How the wells of this source were selected", "source")
    svar("source_citation", [info[k][3] for k in SRC], "Credit line of the source", "source")
    nvar("source_first_station", [int(np.flatnonzero(st_src == k)[0]) for k in SRC], "i4",
         "First row of the station axis that belongs to this source", "source")
    nvar("source_station_count", [int((st_src == k).sum()) for k in SRC], "i4",
         "Number of consecutive station rows of this source", "source")
    nvar("source_first_obs", [int(np.flatnonzero(ob_src == k)[0]) for k in SRC], "i4",
         "First row of the obs axis that belongs to this source", "source")
    nvar("source_obs_count", [int((ob_src == k).sum()) for k in SRC], "i4",
         "Number of consecutive obs rows of this source", "source")
    nvar("source_first_year", [int(meas.year[ob_src == k].min()) for k in SRC], "i2",
         "Earliest measurement year of this source", "source")
    nvar("source_last_year", [int(meas.year[(ob_src == k) & (meas.year <= 2100)].max()) for k in SRC], "i2",
         "Latest measurement year of this source (years after 2100 in the source are typing errors, kept as given)",
         "source")
    nvar("station_source", pd.Series(st_src).map({k: i for i, k in enumerate(SRC)}), "i1",
         "Row of the source axis this well belongs to", "station", flag_values=np.arange(len(SRC), dtype="i1"),
         flag_meanings=" ".join(SRC))
    svar("source", wells.source, "Database the well record comes from (USGS, TWDB, NECSD)", "station")
    svar("agency", wells.agency, "Agency code or name given by the source", "station")
    svar("well_name", wells.name, "Well name, owner or legal description as given by the source", "station")
    svar("state", wells.state, "State", "station")
    svar("county", wells.county, "County", "station")
    svar("aquifer", wells.aquifer, "Aquifer code or name as given by the source (blank for Nebraska)", "station")
    svar("usgs_site_no", wells.usgs_site_no, "USGS site number where the source gives one", "station")
    nvar("lat", wells.lat, "f8", "Latitude of well", "station", standard_name="latitude", units="degrees_north")
    nvar("lon", wells.lon, "f8", "Longitude of well", "station", standard_name="longitude", units="degrees_east")
    nvar("well_depth", wells.well_depth_ft, "f8", "Well depth", "station", fill=-9999.0, units="ft",
         coordinates="lat lon")
    nvar("in_polygon", wells.in_polygon, "i1", "Well lies inside the aquifer boundary hp_bound2010", "station",
         flag_values=np.array([0, 1], "i1"), flag_meanings="no yes")
    nvar("well_group", wells.well_group, "i4", f"Group of wells of different sources within {MATCH_M:.0f} m", "station")
    nvar("station_index", meas.station, "i4", "Row of the station axis of this record", "obs",
         instance_dimension="station")
    nvar("time", (meas.date - EPOCH).dt.days, "f8", "Date of the measurement", "obs", fill=-9999.0,
         standard_name="time", units="days since 1900-01-01 00:00:00", calendar="standard",
         comment="Empty where the source gives only a year (see year and status); no date is invented.")
    nvar("year", meas.year, "i2", "Calendar year of the measurement as given by the source", "obs")

    def cvar(name, values, long_name, **a):
        """Text with few distinct values: stored as an index into a lookup table (keeps the file small)."""
        codes, table = pd.factorize(pd.Series(values).fillna("").astype(str), sort=True)
        ds.createDimension(name + "_n", len(table))
        t = ds.createVariable(name + "_text", str, (name + "_n",))
        t[:] = np.array(list(table), dtype=object)
        t.long_name = f"Lookup table of {name}"
        nvar(name, codes, "i4", long_name + f" (index into {name}_text)", "obs", **a)

    nvar("depth_to_water", meas.depth_ft, "f8", "Depth to water below land surface", "obs", fill=-9999.0, units="ft",
         positive="down", coordinates="time", comment="Negative values are water levels above land surface.")
    cvar("status", meas.status, "Status of the measurement as given by the source",
         comment="USGS: qualifier (for example Static, Pumping, Recently pumped, Dry, Obstructed). TWDB: "
                 "publication status and remark. NECSD: season of the measurement.")
    cvar("approval", meas.approval, "USGS approval status (blank for other sources)")
    cvar("measuring_agency", meas.agency, "Agency that made the measurement, as given by the source")
    cvar("method", meas.method, "Method of measurement, as given by the source")
    cvar("record_source", meas.source, "Database the record comes from")
    nvar("duplicate_of_other_source", meas.dup, "i1", "Record repeats one from a source listed earlier", "obs",
         flag_values=np.array([0, 1], "i1"), flag_meanings="no yes")

yr = meas.year
log(f"wrote {OUT.name} ({OUT.stat().st_size / 1e6:.1f} MB): {len(wells)} wells, {len(meas)} records, "
    f"{meas.date.min().date()} to {meas.date.max().date()}")
t = meas[meas.dup == 0].assign(year=yr).query("year >= 2000").groupby(["year", "source"]).size().unstack(fill_value=0)
print(t.to_string())
