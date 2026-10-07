"""
Rebuild the USGS High Plains aquifer products at their NATIVE resolution, with every column USGS publishes.

Sources (read-only), all under SRC:
  F01 hp_wlcpd19t.zip            raster, water-level change predevelopment (about 1950) to 2019, ft, 500 m, EPSG:5070
  F04 hp_wlc1719t.zip            raster, water-level change 2017 to 2019, ft, same grid
  F02/F03/F05 *_wells_A83.zip    well shapefiles (primary, supplemental, 2017-19 networks)
  ofr98-548.e00.gz               hydraulic-conductivity polygons (USGS Open-File Report 98-548)
  USGS data release doi:10.5066/P9WPP01S (McGuire and Strauch), companion to SIR 2023-5143.

Outputs (new files only) in OUT:
  USGS_HPA_water_level_change_500m.nc   both change rasters on USGS's own 500 m Albers grid, cell values
                                        copied bit for bit, cells whose centre is outside hp_bound2010 set to
                                        the raster's own no-data value; aquifer mask; conductivity class limits
                                        rasterised onto the same grid
  USGS_HPA_hydraulic_conductivity_polygons.gpkg   the conductivity polygons themselves, all attributes
  USGS_HPA_wells.nc                     every well, every water-level record and every attribute column of
                                        the three shapefiles

Nothing is resampled, averaged or converted to other units. Differences from the earlier derived_usgs files:
  * rasters stay on the 500 m grid (no 4 km averaging, so USGS's zero areas and contour plateaus are intact)
  * conductivity uses the RANGE text of each polygon; no-data polygons (RANGE -1) stay empty instead of 0
  * wells keep USGS's use codes, source codes, extra change columns, measurement times and partial dates;
    no date is invented (a record without a date has no time value)
"""
import gzip
import shutil
import tempfile
import warnings
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import netCDF4 as nc4
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize

warnings.filterwarnings("ignore")

SRC = Path(r"H:/USGS GW dataset")
OUT = Path(r"G:/MSU_GWB/datasets/merged_datasets/USGS_HPA_native")
BOUND = "zip://G:/MSU_GWB/datasets/high_plains_quifer.zip!high_plains_quifer/hp_bound2010.shp"
TMP = Path(tempfile.mkdtemp(prefix="usgs_hpa_"))
OUT.mkdir(parents=True, exist_ok=True)

RELEASE = ("McGuire, V.L., and Strauch, K.R., Data from maps of water-level changes in the High Plains aquifer, "
           "predevelopment (about 1950) to 2019 and 2017 to 2019: U.S. Geological Survey data release, "
           "https://doi.org/10.5066/P9WPP01S (companion to Scientific Investigations Report 2023-5143)")
K_REF = ("Cederstrand, J.R., and Becker, M.F., 1998, Digital map of hydraulic conductivity for the High Plains "
         "aquifer: U.S. Geological Survey Open-File Report 98-548 (digitised from Gutentag and others, 1984, "
         "Professional Paper 1400-B, 1:1,000,000)")
NOW = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
COMP = dict(zlib=True, complevel=4, shuffle=True)


def log(msg):
    print(f"[usgs-native] {msg}")


def one(pattern):
    hits = sorted(SRC.glob(pattern))
    if len(hits) != 1:
        raise FileNotFoundError(f"expected exactly one match for {pattern}, found {len(hits)}")
    return hits[0]


def unzip(pattern, suffix):
    z = one(pattern)
    dest = TMP / z.stem
    with zipfile.ZipFile(z) as f:
        f.extractall(dest)
    return next(dest.rglob(f"*{suffix}"))


# ════════════════════════════════════════════════════════════════════════════
# 1. Rasters, mask and conductivity on USGS's own 500 m grid
# ════════════════════════════════════════════════════════════════════════════
def build_rasters():
    tifs = {"dWL_predev_to_2019": unzip("F01_*/hp_wlcpd19t.zip", ".tif"),
            "dWL_2017_to_2019": unzip("F04_*/hp_wlc1719t.zip", ".tif")}
    texts = {"dWL_predev_to_2019": "Mapped water-level change, predevelopment (about 1950) to 2019",
             "dWL_2017_to_2019": "Mapped water-level change, 2017 to 2019"}
    data, ref = {}, None
    for name, path in tifs.items():
        with rasterio.open(path) as src:
            meta = dict(transform=src.transform, crs=src.crs, shape=src.shape, nodata=src.nodata, dtype=src.dtypes[0])
            if ref is None:
                ref = meta
            elif (meta["transform"], meta["shape"], meta["nodata"]) != (ref["transform"], ref["shape"], ref["nodata"]):
                raise ValueError("the two USGS rasters are not on the same grid")
            data[name] = src.read(1)
    if ref["crs"].to_epsg() != 5070:
        raise ValueError(f"unexpected raster CRS {ref['crs']}")
    H, W = ref["shape"]
    tr = ref["transform"]
    nodata = np.float32(ref["nodata"])
    x = tr.c + tr.a * (np.arange(W) + 0.5)
    y = tr.f + tr.e * (np.arange(H) + 0.5)

    bound = gpd.read_file(BOUND).to_crs(5070)
    mask = rasterize(((g, 1) for g in bound.geometry), out_shape=(H, W), transform=tr, fill=0,
                     dtype="uint8", all_touched=False)                 # 1 where the cell centre is in the polygon
    stats = {}
    for name, a in data.items():
        valid = a != nodata
        stats[name] = dict(raw_valid=int(valid.sum()), dropped=int((valid & (mask == 0)).sum()),
                           polygon_cells_without_value=int((~valid & (mask == 1)).sum()))
        data[name] = np.where(mask == 1, a, nodata).astype("float32")
        log(f"{name}: raw valid cells {stats[name]['raw_valid']}, outside the polygon and set to no-data "
            f"{stats[name]['dropped']}, polygon cells USGS left empty {stats[name]['polygon_cells_without_value']}")

    # ---- conductivity polygons ----
    e00 = TMP / "ofr98-548.e00"
    with gzip.open(one("Digital map of hydraulic conductivity*/ofr98-548.e00.gz"), "rb") as fi, open(e00, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    pal = gpd.read_file(e00, layer="PAL")
    pal["RANGE"] = pal["RANGE"].astype(str).str.strip()
    lo, hi = [], []
    for r in pal["RANGE"]:
        parts = r.split(" to ")
        if len(parts) == 2:
            lo.append(int(parts[0]))
            hi.append(int(parts[1]))
        else:                                            # "-1": USGS's code for no conductivity class
            lo.append(-1)
            hi.append(-1)
    pal["K_lower_ftday"], pal["K_upper_ftday"] = lo, hi
    pal["limits_conflict"] = ((pal["K_lower_ftday"] >= 0)
                              & ((pal["MAJOR1"] != pal["K_lower_ftday"]) | (pal["MINOR1"] != pal["K_upper_ftday"]))
                              ).astype("int8")
    pal5070 = pal.to_crs(5070)
    log(f"conductivity: {len(pal)} polygons, {(pal.K_lower_ftday < 0).sum()} without a class (RANGE -1), "
        f"{int(pal.limits_conflict.sum())} whose RANGE text disagrees with MAJOR1/MINOR1")
    gpkg = OUT / "USGS_HPA_hydraulic_conductivity_polygons.gpkg"
    if gpkg.exists():
        gpkg.unlink()
    pal5070.to_file(gpkg, layer="hydraulic_conductivity", driver="GPKG")

    def burn(col, fill, dtype):
        real = pal5070[pal5070.K_lower_ftday >= 0]
        out = rasterize(((g, int(v)) for g, v in zip(real.geometry, real[col])), out_shape=(H, W), transform=tr,
                        fill=fill, dtype=dtype, all_touched=False)
        return np.where(mask == 1, out, fill).astype(dtype)

    k_lo, k_hi = burn("K_lower_ftday", -1, "int16"), burn("K_upper_ftday", -1, "int16")
    k_conf = burn("limits_conflict", -1, "int8")

    path = OUT / "USGS_HPA_water_level_change_500m.nc"
    with nc4.Dataset(path, "w", format="NETCDF4") as ds:
        ds.setncatts(dict(
            Conventions="CF-1.9",
            title="USGS High Plains aquifer: mapped water-level change and hydraulic conductivity on the native "
                  "500 m grid",
            summary="USGS's two water-level change rasters copied cell for cell onto their own 500 m NAD83 / Conus "
                    "Albers grid (no resampling, no unit conversion). Cells whose centre is outside the aquifer "
                    "boundary hp_bound2010 hold the raster's own no-data value. The conductivity class limits of "
                    "Open-File Report 98-548 are rasterised onto the same grid by cell centre.",
            institution="U.S. Geological Survey (source data); packaged by the MSU groundwater project",
            source=f"{tifs['dWL_predev_to_2019'].name}, {tifs['dWL_2017_to_2019'].name}, ofr98-548.e00.gz",
            references=RELEASE + "; " + K_REF,
            comment="The change rasters are contour-based maps: USGS's own files are exactly 0 over 59 % "
                    "(predevelopment to 2019) and 74 % (2017 to 2019) of the aquifer and hold plateaus at contour "
                    "values. They are not continuous surfaces.",
            history=f"{NOW} created by build_usgs_hpa_native.py; sources read-only",
            license="USGS data are in the public domain"))
        ds.createDimension("y", H)
        ds.createDimension("x", W)
        vy = ds.createVariable("y", "f8", ("y",))
        vx = ds.createVariable("x", "f8", ("x",))
        vy[:], vx[:] = y, x
        vy.setncatts(dict(standard_name="projection_y_coordinate", long_name="northing of cell centre", units="m",
                          axis="Y"))
        vx.setncatts(dict(standard_name="projection_x_coordinate", long_name="easting of cell centre", units="m",
                          axis="X"))
        crs = ds.createVariable("crs", "i4")
        crs.setncatts(dict(grid_mapping_name="albers_conical_equal_area", standard_parallel=[29.5, 45.5],
                           latitude_of_projection_origin=23.0, longitude_of_central_meridian=-96.0,
                           false_easting=0.0, false_northing=0.0, semi_major_axis=6378137.0,
                           inverse_flattening=298.257222101, epsg_code="EPSG:5070",
                           crs_wkt=ref["crs"].to_wkt(), GeoTransform=" ".join(str(v) for v in tr.to_gdal())))
        ck = (min(512, H), min(512, W))
        for name, a in data.items():
            v = ds.createVariable(name, "f4", ("y", "x"), fill_value=nodata, chunksizes=ck, **COMP)
            v[:] = a
            v.setncatts(dict(long_name=texts[name], units="ft", grid_mapping="crs",
                             comment="Negative = water-level decline, positive = rise. Values copied unchanged "
                                     "from USGS's raster.",
                             source_file=tifs[name].name, **{k: np.int32(s) for k, s in stats[name].items()}))
        v = ds.createVariable("aquifer_mask", "i1", ("y", "x"), chunksizes=ck, **COMP)
        v[:] = mask.astype("i1")
        v.setncatts(dict(long_name="Cell centre lies inside the High Plains aquifer boundary (hp_bound2010)",
                         flag_values=np.array([0, 1], "i1"), flag_meanings="outside inside", grid_mapping="crs"))
        for name, a, text in (("K_lower", k_lo, "Lower limit of the hydraulic-conductivity class"),
                              ("K_upper", k_hi, "Upper limit of the hydraulic-conductivity class")):
            v = ds.createVariable(name, "i2", ("y", "x"), fill_value=np.int16(-1), chunksizes=ck, **COMP)
            v[:] = a
            v.setncatts(dict(long_name=text, units="ft day-1", grid_mapping="crs",
                             comment="From the RANGE text of the polygon holding the cell centre (USGS classes: "
                                     "0 to 25, 25 to 50, 50 to 100, 100 to 200, 200 to 300, 300 to 400, 400 to "
                                     "500 ft/day; one polygon reads '10 to 25'). Empty where USGS gives no class "
                                     "(RANGE -1) or no polygon. Source scale 1:1,000,000."))
        v = ds.createVariable("K_limits_conflict", "i1", ("y", "x"), fill_value=np.int8(-1), chunksizes=ck, **COMP)
        v[:] = k_conf
        v.setncatts(dict(long_name="RANGE text of the conductivity polygon disagrees with its MAJOR1/MINOR1 fields",
                         flag_values=np.array([0, 1], "i1"), flag_meanings="consistent conflict", grid_mapping="crs",
                         comment="USGS's table has 9 such polygons. K_lower / K_upper follow the RANGE text."))
    log(f"wrote {path.name} ({path.stat().st_size / 1e6:.1f} MB), grid {H} x {W}")
    return path, tifs, nodata


# ════════════════════════════════════════════════════════════════════════════
# 2. Wells: every record and every column
# ════════════════════════════════════════════════════════════════════════════
SENTINELS = (-9999.0, -999.0)
CAMPAIGNS = ["pd", "80", "2015", "2016", "2017", "2018", "2019"]          # USGS column suffixes
STATION_FLOAT = {"F02": ["deltapd_19"],
                 "F03": ["deltapd_15", "deltapd_16", "deltapd_17", "deltapd_18", "dpd_19est", "delta80_19",
                         "PD80_CNTR"],
                 "F05": ["delta17_19"]}
STATION_TEXT = {"F02": ["dpd_19src"], "F03": ["dpd_19src", "PD80_RANGE"], "F05": []}
STATION_LONG = {
    "deltapd_19": "Water-level change, predevelopment to 2019 (primary network; for 71 New Mexico wells to the "
                  "latest level of 2015-18, see dpd_19src)",
    "deltapd_15": "Water-level change, predevelopment to 2015 (supplemental network)",
    "deltapd_16": "Water-level change, predevelopment to 2016 (supplemental network)",
    "deltapd_17": "Water-level change, predevelopment to 2017 (supplemental network)",
    "deltapd_18": "Water-level change, predevelopment to 2018 (supplemental network)",
    "dpd_19est": "ESTIMATED water-level change, predevelopment to 2019 (supplemental network): either from a "
                 "post-predevelopment level or from the 1980 contour value plus the 1980-2019 change, see dpd_19src",
    "delta80_19": "Water-level change, 1980 to 2019 (supplemental network)",
    "PD80_CNTR": "Beginning contour value of the predevelopment-to-1980 change polygon holding the well",
    "delta17_19": "Water-level change, 2017 to 2019",
    "dpd_19src": "USGS code for how the predevelopment-to-2019 change of this well was obtained",
    "PD80_RANGE": "Contour range of the predevelopment-to-1980 change polygon holding the well",
}


def text(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    if isinstance(v, (float, np.floating)) and float(v).is_integer():
        v = int(v)
    return str(v).strip().strip("'\"")


def number(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return np.nan
    return np.nan if (np.isnan(f) or f in SENTINELS) else f


def parse_date(s):
    """USGS dates are YYYYMMDD, YYYYMM or YYYY. Returns (days since 1900-01-01 or NaN, precision flag):
    0 = day, 1 = month (first day of the month stored), 2 = year (1 January stored), 3 = no date."""
    try:
        if len(s) == 8 and s.isdigit():
            d, p = datetime(int(s[:4]), int(s[4:6]), int(s[6:])), 0
        elif len(s) == 6 and s.isdigit():
            d, p = datetime(int(s[:4]), int(s[4:6]), 1), 1
        elif len(s) == 4 and s.isdigit():
            d, p = datetime(int(s), 1, 1), 2
        else:
            return np.nan, 3
    except ValueError:
        return np.nan, 3
    return float((d - datetime(1900, 1, 1)).days), p


def build_wells():
    shp = {"F02": unzip("F02_*/hp_wlcpd19_wells_A83.zip", ".shp"),
           "F03": unzip("F03_*/hp_wlcpd19supwells_A83.zip", ".shp"),
           "F05": unzip("F05_*/hp_wlc1719_wells_A83.zip", ".shp")}
    tabs = {k: gpd.read_file(p) for k, p in shp.items()}
    for k, g in tabs.items():
        g["badge"] = g["sitebadge"].map(text)
        if g["badge"].duplicated().any() or (g["badge"] == "").any():
            raise ValueError(f"{k}: blank or repeated sitebadge")
        log(f"{k}: {len(g)} wells, {len(g.columns) - 2} attribute columns")

    # ---- stations ----
    static_cols = ["usgs_id", "station_na", "state", "county", "lat_nad83", "long_nad83", "well_depth"]
    allst = pd.concat([g[["badge"] + static_cols].assign(file=k) for k, g in tabs.items()], ignore_index=True)
    for c in ("usgs_id", "station_na", "state", "county"):
        allst[c] = allst[c].map(text)
    allst["well_depth"] = allst["well_depth"].map(number)
    conflicts = {c: int((allst.groupby("badge")[c].nunique(dropna=True) > 1).sum()) for c in static_cols}
    log(f"wells listed in more than one file with differing static values: {conflicts}")
    st = allst.drop_duplicates("badge").sort_values("badge").reset_index(drop=True)      # F02 > F03 > F05
    index = {b: i for i, b in enumerate(st["badge"])}
    n_st = len(st)

    # ---- observations: one row per file x well x campaign with a water level ----
    rows = []
    for rank, (k, g) in enumerate(tabs.items()):
        for c in CAMPAIGNS:
            if "wl" + c not in g.columns:
                continue
            wl = g["wl" + c].map(number)
            for i in np.flatnonzero(wl.notna().values):
                r = g.iloc[i]
                ds_ = text(r.get("date" + c)) if "date" + c in g.columns else ""
                if ds_ == "0":
                    ds_ = ""                                              # USGS: date not available
                t, prec = parse_date(ds_)
                use = r.get("use" + c) if "use" + c in g.columns else None
                rows.append(dict(station=index[r["badge"]], rank=rank, file=k,
                                 campaign="1980" if c == "80" else c, wl=float(wl.iloc[i]), date_string=ds_,
                                 time=t, precision=prec,
                                 time_of_day=text(r.get("time" + c)) if "time" + c in g.columns else "",
                                 use=int(use) if use is not None and not pd.isna(use) else -32768))
    ob = pd.DataFrame(rows)
    order = {"pd": 0, "1980": 1, "2015": 2, "2016": 3, "2017": 4, "2018": 5, "2019": 6}
    ob["corder"] = ob["campaign"].map(order)
    ob = ob.sort_values(["station", "corder", "rank"]).reset_index(drop=True)
    ob["repeat_file"] = ob.duplicated(["station", "campaign", "wl"], keep="first").astype("int8")
    dated = ob["date_string"].str.len() == 8
    ob["repeat_campaign"] = (dated & ob.duplicated(["station", "file", "date_string", "wl"], keep=False)).astype("int8")
    log(f"water-level records {len(ob)} ({int((ob.repeat_file == 0).sum())} after removing repeats between files); "
        f"date precision day/month/year/none = {ob.precision.value_counts().sort_index().to_dict()}")

    path = OUT / "USGS_HPA_wells.nc"
    with nc4.Dataset(path, "w", format="NETCDF4") as ds:
        ds.setncatts(dict(
            Conventions="CF-1.9", featureType="timeSeries",
            title="USGS High Plains aquifer wells: water levels and water-level changes, predevelopment to 2019",
            summary="Every well, every water-level record and every attribute column of the three USGS well "
                    "shapefiles: primary network (F02), supplemental network (F03) and the 2017-19 network (F05). "
                    "Values are copied unchanged, in feet. Records are not merged between files; "
                    "repeat_of_other_file marks the repeats.",
            institution="U.S. Geological Survey (source data); packaged by the MSU groundwater project",
            source=", ".join(p.name for p in shp.values()), references=RELEASE,
            comment="Read use_code before using an early level: in the supplemental network only use_code 1 means "
                    "a true predevelopment level; -2 is a later (post-predevelopment) level and -5 was not used by "
                    "USGS. USGS's missing-data codes -9999 and -999 are stored as _FillValue. No date is invented: "
                    "a record without a date has no time value (date_precision 3).",
            history=f"{NOW} created by build_usgs_hpa_native.py; sources read-only",
            license="USGS data are in the public domain"))
        ds.createDimension("station", n_st)
        ds.createDimension("obs", len(ob))

        def svar(name, values, long_name, dim="station", **attrs):
            v = ds.createVariable(name, str, (dim,))
            v[:] = np.array(list(values), dtype=object)
            v.setncatts(dict(long_name=long_name, **{k: w for k, w in attrs.items() if w is not None}))
            return v

        def nvar(name, values, dtype, long_name, dim="station", fill=None, **attrs):
            v = ds.createVariable(name, dtype, (dim,), fill_value=fill, **COMP)
            a = np.asarray(values, dtype="f8" if np.dtype(dtype).kind == "f" else dtype)
            if fill is not None and np.dtype(dtype).kind == "f":
                a = np.where(np.isnan(a), fill, a)
            v[:] = a.astype(dtype)
            v.setncatts(dict(long_name=long_name, **{k: w for k, w in attrs.items() if w is not None}))
            return v

        svar("sitebadge", st["badge"], "USGS well identifier (agency code : well number)", cf_role="timeseries_id")
        svar("usgs_id", st["usgs_id"], "USGS site identifier (blank where USGS gives none)")
        svar("station_na", st["station_na"], "Station name")
        svar("state", st["state"], "State FIPS code")
        svar("county", st["county"], "County FIPS code")
        nvar("lat", st["lat_nad83"], "f8", "Latitude of well (NAD83)", standard_name="latitude", units="degrees_north")
        nvar("lon", st["long_nad83"], "f8", "Longitude of well (NAD83)", standard_name="longitude",
             units="degrees_east")
        nvar("well_depth", st["well_depth"], "f8", "Well depth", fill=-9999.0, units="ft", coordinates="lat lon")
        for k, g in tabs.items():
            member = np.zeros(n_st, "i1")
            pos = g["badge"].map(index).values
            member[pos] = 1
            nvar(f"in_{k.lower()}", member, "i1", f"Well is listed in USGS file {k} ({shp[k].name})",
                 flag_values=np.array([0, 1], "i1"), flag_meanings="no yes")
            for c in STATION_FLOAT[k]:
                a = np.full(n_st, np.nan)
                a[pos] = g[c].map(number).values
                nvar(f"{k.lower()}_{c}", a, "f8", STATION_LONG[c], fill=-9999.0, units="ft", coordinates="lat lon",
                     source_column=f"{k} {c}", comment="Negative = decline, positive = rise.")
            for c in STATION_TEXT[k]:
                a = np.full(n_st, "", dtype=object)
                a[pos] = g[c].map(text).values
                svar(f"{k.lower()}_{c}", a, STATION_LONG[c], source_column=f"{k} {c}",
                     comment="pd19: predevelopment to 2019. pd15..pd18: predevelopment to that year, used as the "
                             "2019 estimate. pd19est: from a post-predevelopment level (dated 15 June 1978 or "
                             "earlier) to 2019. RASA19: predevelopment-to-1980 contour value plus the 1980-2019 "
                             "change." if c == "dpd_19src" else None)

        nvar("station_index", ob["station"], "i4", "Row of the station axis of this record", dim="obs",
             instance_dimension="station")
        svar("source_file", ob["file"], "USGS file the record comes from (F02 primary, F03 supplemental, F05 2017-19)",
             dim="obs")
        svar("campaign", ob["campaign"], "USGS water-level column of the record", dim="obs",
             comment="pd = column wlpd (the earliest or one of the earliest levels of the well; see use_code for "
                     "whether USGS classes it as predevelopment); 1980 = wl80; otherwise the water year of the "
                     "pre-irrigation-season measurement. For New Mexico wells a level of 2015-18 is also the 2019 "
                     "estimate.")
        t = nvar("time", ob["time"], "f8", "Date of the water-level measurement", dim="obs", fill=-9999.0,
                 standard_name="time", units="days since 1900-01-01 00:00:00", calendar="standard",
                 comment="From USGS's date field. Where USGS gives only a year-month or a year, the first day of "
                         "that month or year is stored and date_precision says so. Empty where USGS gives no date "
                         "(all 1980 levels).")
        nvar("date_precision", ob["precision"], "i1", "Precision of USGS's measurement date", dim="obs",
             flag_values=np.array([0, 1, 2, 3], "i1"), flag_meanings="day month year no_date")
        svar("date_string", ob["date_string"], "USGS's measurement date as written (YYYYMMDD, YYYYMM or YYYY)",
             dim="obs")
        svar("time_of_day", ob["time_of_day"], "USGS's measurement time as written (HHMM), blank where not given",
             dim="obs")
        nvar("water_level", ob["wl"], "f8", "Depth to water below land surface", dim="obs", fill=-9999.0, units="ft",
             positive="down", coordinates="time",
             comment="Negative values are water levels above land surface.")
        nvar("use_code", ob["use"], "i2", "USGS use code of this water level (columns usepd, use80, use2015 ... "
                                           "use2019)", dim="obs", fill=np.int16(-32768),
             comment="Greater than 0: used by USGS in the change calculation. Early levels (campaign pd): 1 = "
                     "predevelopment level; -2 = post-predevelopment level, used for an estimate; -5 = not used. "
                     "-999 = not applicable. In F02, use2019 values 3 to 6 mean the 2019 value is the level of "
                     "2015 to 2018 (New Mexico). See USGS's metadata for the full definitions.")
        nvar("repeat_of_other_file", ob["repeat_file"], "i1",
             "Same well, campaign and level as a record from an earlier-listed file", dim="obs",
             flag_values=np.array([0, 1], "i1"), flag_meanings="no yes",
             comment="Keep the records with 0 to get each measurement once.")
        nvar("same_measurement_as_other_campaign", ob["repeat_campaign"], "i1",
             "USGS lists this dated measurement under more than one campaign column of the same file", dim="obs",
             flag_values=np.array([0, 1], "i1"), flag_meanings="no yes")
        del t
    log(f"wrote {path.name} ({path.stat().st_size / 1e6:.2f} MB): {n_st} wells, {len(ob)} records")
    return path, tabs, ob


# ════════════════════════════════════════════════════════════════════════════
# 3. Verification against the sources
# ════════════════════════════════════════════════════════════════════════════
def verify(raster_path, tifs, nodata, wells_path, tabs):
    with nc4.Dataset(raster_path) as ds:
        ds.set_auto_maskandscale(False)
        mask = ds["aquifer_mask"][:] == 1
        for name, tif in tifs.items():
            with rasterio.open(tif) as src:
                a = src.read(1)
            b = ds[name][:]
            same_in = a[mask].tobytes() == b[mask].tobytes()
            out_fill = bool(np.all(b[~mask] == nodata))
            log(f"VERIFY {name}: cells inside the polygon identical to USGS's raster bit for bit = {same_in}; "
                f"cells outside all no-data = {out_fill}")
            if not (same_in and out_fill):
                raise AssertionError(name)
        lo, hi = ds["K_lower"][:], ds["K_upper"][:]
        u = sorted({(int(p), int(q)) for p, q in zip(lo[lo >= 0], hi[lo >= 0])})
        log(f"VERIFY conductivity classes on the grid (ft/day): {u}; cells with a class {int((lo >= 0).sum())}, "
            f"polygon cells without one {int(((lo < 0) & mask).sum())}, zero-valued cells {int((hi == 0).sum())}")
    with nc4.Dataset(wells_path) as ds:
        ds.set_auto_maskandscale(False)
        badge = np.array([str(v) for v in ds["sitebadge"][:]])
        pos = {b: i for i, b in enumerate(badge)}
        si = ds["station_index"][:]
        files = np.array([str(v) for v in ds["source_file"][:]])
        camp = np.array([str(v) for v in ds["campaign"][:]])
        wl = ds["water_level"][:]
        checked = bad = 0
        for k, g in tabs.items():
            idx = g["badge"].map(pos).values
            for c in STATION_FLOAT[k]:
                want = g[c].map(number).values
                got = ds[f"{k.lower()}_{c}"][:][idx]
                got = np.where(got == -9999.0, np.nan, got)
                bad += int((~((want == got) | (np.isnan(want) & np.isnan(got)))).sum())
                checked += len(want)
            for c in CAMPAIGNS:
                if "wl" + c not in g.columns:
                    continue
                want = g["wl" + c].map(number)
                sel = (files == k) & (camp == ("1980" if c == "80" else c))
                got = pd.Series(wl[sel], index=badge[si[sel]])
                w = want[want.notna()]
                w.index = g["badge"][want.notna()].values
                bad += int(len(w) != len(got)) + int((got.reindex(w.index).values != w.values).sum())
                checked += len(w)
        log(f"VERIFY wells: {checked} values compared with the shapefiles, {bad} differ")
        if bad:
            raise AssertionError("wells")


if __name__ == "__main__":
    try:
        rp, tifs, nodata = build_rasters()
        wp, tabs, ob = build_wells()
        verify(rp, tifs, nodata, wp, tabs)
        log("DONE")
    finally:
        shutil.rmtree(TMP, ignore_errors=True)
