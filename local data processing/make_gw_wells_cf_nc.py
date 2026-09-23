"""
Convert raw USGS High Plains groundwater WELL data (F02/F03/F05 shapefiles)
into ONE CF-1.8 compliant NetCDF (discrete sampling geometry, featureType="point").

Sources (READ-ONLY):
  F02 hp_wlcpd19_wells_A83.shp      : primary wells, predev(~1950)->2019 campaign
  F03 hp_wlcpd19supwells_A83.shp    : supplemental wells (adds 1980 epoch, est. deltas)
  F05 hp_wlc1719_wells_A83.shp      : dense network, 2017->2019 campaign

Output (NEW file only):
  derived_usgs/GW_wells_HPA_CF.nc

CF design
---------
featureType = "point"
dim station : unique wells (deduped by site badge across the three files)
dim obs     : one record per dated water-level measurement
coords  : lat(station), lon(station)            [degrees_north/east, NAD83 datum noted]
          time(obs)                             [days since 1900-01-01]
link    : station_index(obs) -> station
data    : water_level_ft(obs)                   depth-to-water below land surface, ft
static  : well_depth_ft(station),
          dWL_pd19_primary_ft(station)   (F02 deltapd_19)
          dWL_pd19_supplemental_ft(station) (F03 dpd_19est)
          dWL_1980_2019_ft(station)      (F03 delta80_19)
          dWL_2017_2019_ft(station)      (F05 delta17_19)
ids     : site_badge(station), usgs_id(station), station_name(station),
          state(station), county(station), source_dataset(station)
Sentinels -9999 / -999 -> NaN (_FillValue).
"""
import warnings
warnings.filterwarnings("ignore")
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr

SRC = Path(r"G:/USGS GW dataset")
TMP = Path(r"C:/Users/AlienX/AppData/Local/Temp/opencode/usgs_raw")
OUT = Path(r"G:/MSU_GWB/datasets/derived_usgs/GW_wells_HPA_CF.nc")
OUT.parent.mkdir(parents=True, exist_ok=True)

def log(m): print(f"[wells] {m}")

# ------------------------------------------------------------------
# extract zips (idempotent)
# ------------------------------------------------------------------
def ensure(shp_zip_glob, dest, member_dir):
    shp = next(dest.rglob("*.shp")) if dest.exists() and any(dest.rglob("*.shp")) else None
    if shp is None:
        zp = next(SRC.glob(shp_zip_glob))
        dest.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zp) as z:
            z.extractall(dest)
        shp = next(dest.rglob("*.shp"))
    return shp

f02 = ensure("F02_*/*.zip", TMP/"f02", None)
f03 = ensure("F03_*/*.zip", TMP/"f03", None)
f05 = ensure("F05_*/*.zip", TMP/"f05", None)
log(f"F02 {f02.name} | F03 {f03.name} | F05 {f05.name}")

def clean_str(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    return str(v).strip().strip("'\"")

def num(v, sentinels=(-9999.0, -999.0, -9999)):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return np.nan
    return np.nan if f in sentinels else f

def read_wells(path, kind):
    g = gpd.read_file(path)
    recs = []
    for _, r in g.iterrows():
        base = dict(
            site_badge=clean_str(r.get("sitebadge")),
            usgs_id=clean_str(r.get("usgs_id")),
            station_name=clean_str(r.get("station_na")),
            state=clean_str(r.get("state")),
            county=clean_str(r.get("county")),
            lat=float(r["lat_nad83"]),
            lon=float(r["long_nad83"]),
            well_depth_ft=num(r.get("well_depth")),
            source=kind,
        )
        # per-file water level epochs: (column_wl, column_date, campaign_label)
        if kind == "F02":
            epochs = [(f"wl{y}", f"date{y}", y) for y in
                      ["pd", "2015", "2016", "2017", "2018", "2019"]]
        elif kind == "F03":
            epochs = [("wl80", None, "1980")] + \
                     [(f"wl{y}", f"date{y}", y) for y in
                      ["pd", "2015", "2016", "2017", "2018", "2019"]]
        else:  # F05
            epochs = [("wl2017", "date2017", "2017"), ("wl2019", "date2019", "2019")]
        for wl_col, dt_col, camp in epochs:
            wl = num(r.get(wl_col))
            if not np.isfinite(wl):
                continue
            d = None
            if dt_col is not None:
                raw = r.get(dt_col)
                s = "" if raw is None else str(raw).strip()
                if len(s) == 8 and s.isdigit():
                    d = pd.to_datetime(s, format="%Y%m%d", errors="coerce")
            if d is None or pd.isna(d):
                # fallback mid-epoch nominal date (documented approximation for
                # records where the measurement date field is blank)
                nom = {"1980": "1980-04-01", "pd": "1950-01-01",
                       "2015": "2015-04-01", "2016": "2016-04-01",
                       "2017": "2017-04-01", "2018": "2018-04-01",
                       "2019": "2019-04-01"}[str(camp)]
                d = pd.to_datetime(nom)
            recs.append({**base, "campaign": str(camp), "time": d,
                         "water_level_ft": wl})
        # stash static deltas on first record holder
        base.update(dict(
            dWL_pd19_primary_ft=num(r.get("deltapd_19")),
            dWL_pd19_supplemental_ft=num(r.get("dpd_19est")),
            dWL_1980_2019_ft=num(r.get("delta80_19")),
            dWL_2017_2019_ft=num(r.get("delta17_19")),
            geometry_valid=g.geometry is not None,
        ))
        static_recs.append(base)

static_recs = []
df02 = gpd.read_file(f02); df03 = gpd.read_file(f03); df05 = gpd.read_file(f05)
log(f"rows: F02 {len(df02)}, F03 {len(df03)}, F05 {len(df05)}")

def collect(gdf, kind):
    out_static, out_obs = [], []
    for _, r in gdf.iterrows():
        b = dict(
            site_badge=clean_str(r.get("sitebadge")),
            usgs_id=clean_str(r.get("usgs_id")),
            station_name=clean_str(r.get("station_na")),
            state=clean_str(r.get("state")),
            county=clean_str(r.get("county")),
            lat=float(r["lat_nad83"]), lon=float(r["long_nad83"]),
            well_depth_ft=num(r.get("well_depth")), source=kind)
        if kind == "F02":
            epochs = [("wlpd", "datepd", "pd")] + \
                     [(f"wl{y}", f"date{y}", y) for y in ["2015","2016","2017","2018","2019"]]
            b["dWL_pd19_primary_ft"] = num(r.get("deltapd_19"))
        elif kind == "F03":
            epochs = [("wl80", None, "1980"), ("wlpd", "datepd", "pd")] + \
                     [(f"wl{y}", f"date{y}", y) for y in ["2015","2016","2017","2018","2019"]]
            b["dWL_pd19_supplemental_ft"] = num(r.get("dpd_19est"))
            b["dWL_1980_2019_ft"] = num(r.get("delta80_19"))
        else:
            epochs = [("wl2017", "date2017", "2017"), ("wl2019", "date2019", "2019")]
            b["dWL_2017_2019_ft"] = num(r.get("delta17_19"))
        for wl_col, dt_col, camp in epochs:
            wl = num(r.get(wl_col))
            if not np.isfinite(wl):
                continue
            d = pd.NaT
            if dt_col is not None:
                s = "" if r.get(dt_col) is None else str(r.get(dt_col)).strip().strip("'\"")
                if len(s) == 8 and s.isdigit():
                    d = pd.to_datetime(s, format="%Y%m%d", errors="coerce")
            if pd.isna(d):
                nom = {"1980":"1980-04-01","pd":"1950-01-01","2015":"2015-04-01",
                       "2016":"2016-04-01","2017":"2017-04-01","2018":"2018-04-01",
                       "2019":"2019-04-01"}[str(camp)]
                d = pd.to_datetime(nom)
            out_obs.append((b["site_badge"], str(camp), d, wl))
        out_static.append(b)
    return out_static, out_obs

S, O = [], []
for gdf, kind in [(df02,"F02"),(df03,"F03"),(df05,"F05")]:
    s_, o_ = collect(gdf, kind)
    S += s_; O += o_

stat = pd.DataFrame(S)
obs  = pd.DataFrame(O, columns=["site_badge","campaign","time","water_level_ft"])
log(f"stations {len(stat)} (raw), obs {len(obs)}")

# dedup stations by site_badge: prefer F02 > F03 > F05 for static attrs;
# combine delta columns from whichever file provides them.
prio = {"F02":0, "F03":1, "F05":2}
stat["_p"] = stat["source"].map(prio)
stat = stat.sort_values("_p")
delta_cols = ["dWL_pd19_primary_ft","dWL_pd19_supplemental_ft",
              "dWL_1980_2019_ft","dWL_2017_2019_ft"]
agg = stat.groupby("site_badge", sort=True).agg(
    usgs_id=("usgs_id","first"),
    station_name=("station_name","first"),
    state=("state","first"),
    county=("county","first"),
    lat=("lat","first"), lon=("lon","first"),
    well_depth_ft=("well_depth_ft","first"),
    source=("source","first"),
    **{c:(c,"first") for c in delta_cols},
).reset_index()
# fill deltas across sources (first-nonnull within group already handled by sort order)
log(f"unique stations after dedup {len(agg)}")

# map station -> index; drop obs whose station vanished (shouldn't happen)
sid = {sb:i for i,sb in enumerate(agg["site_badge"].tolist())}
obs = obs[obs.site_badge.isin(sid)]
# dedupe identical measurements reported by multiple source files
n_before = len(obs)
obs = obs.drop_duplicates(subset=["site_badge", "time", "water_level_ft"])
log(f"dropped {n_before - len(obs)} cross-file duplicate measurements")
obs["station_index"] = obs.site_badge.map(sid).astype("int32")
obs = obs.sort_values(["station_index","time"]).reset_index(drop=True)
log(f"final obs {len(obs)}; campaigns {obs.campaign.value_counts().to_dict()}")

# ------------------------------------------------------------------
# build xarray Dataset (CF-1.8 point DSG)
# ------------------------------------------------------------------
def s_arr(series, itemsize):
    """fixed-width char array ('S') for max CF/netCDF3-style compat"""
    return series.fillna("").astype(str).str.encode("ascii","replace").to_numpy().astype(f"S{itemsize}")

st_ids   = s_arr(agg.site_badge, 24)
st_usgs  = s_arr(agg.usgs_id, 32)
st_names = s_arr(agg.station_name, 32)
st_state = s_arr(agg.state, 4)
st_county= s_arr(agg.county, 6)
st_src   = s_arr(agg.source, 8)
ob_camp  = s_arr(obs.campaign, 8)

ds = xr.Dataset(
    data_vars=dict(
        site_badge=(("station",), st_ids),
        usgs_id=(("station",), st_usgs),
        station_name=(("station",), st_names),
        state=(("station",), st_state),
        county=(("station",), st_county),
        source_dataset=(("station",), st_src),
        well_depth_ft=(("station",), agg.well_depth_ft.astype("float32").values,
                       {"long_name":"total well depth below land surface",
                        "units":"ft", "coordinates":"lat lon"}),
        dWL_pd19_primary_ft=(("station",), agg.dWL_pd19_primary_ft.astype("float32").values,
                       {"long_name":"mapped water-level change, predevelopment to 2019 (primary network)",
                        "units":"ft", "coordinates":"lat lon",
                        "comment":"negative = decline; source F02 attribute deltapd_19"}),
        dWL_pd19_supplemental_ft=(("station",), agg.dWL_pd19_supplemental_ft.astype("float32").values,
                       {"long_name":"estimated water-level change, predevelopment to 2019 (supplemental network)",
                        "units":"ft", "coordinates":"lat lon", "comment":"source F03 dpd_19est"}),
        dWL_1980_2019_ft=(("station",), agg.dWL_1980_2019_ft.astype("float32").values,
                       {"long_name":"water-level change, 1980 to 2019",
                        "units":"ft", "coordinates":"lat lon", "comment":"source F03 delta80_19"}),
        dWL_2017_2019_ft=(("station",), agg.dWL_2017_2019_ft.astype("float32").values,
                       {"long_name":"water-level change, 2017 to 2019",
                        "units":"ft", "coordinates":"lat lon", "comment":"source F05 delta17_19"}),

        station_index=(("obs",), obs.station_index.values.astype("int32"),
                       {"long_name":"station index for this measurement",
                        "instance_dimension":"station"}),
        campaign=(("obs",), ob_camp,
                  {"long_name":"measurement campaign label (epoch of the water-level program)",
                   "comment":"pd=predevelopment(~1950); else calendar year; 1980 epoch only in F03"}),
        water_level_ft=(("obs",), obs.water_level_ft.astype("float32").values,
                        {"long_name":"depth to water below land surface (groundwater level)",
                         "standard_name":"depth_below_geoid",  # closest CF std name; see comment
                         "units":"ft", "positive":"down",
                         "comment":"NOT a CF-recognized standard usage for gw level; kept for searchability. True meaning: depth from land surface down to water table, per USGS field forms."}),
    ),
    coords=dict(
        station=np.arange(len(agg), dtype="int32"),
        obs=np.arange(len(obs), dtype="int32"),
        lat=(("station",), agg.lat.astype("float64").values,
             {"standard_name":"latitude","long_name":"latitude of well",
              "units":"degrees_north","axis":"Y",
              "comment":"datum NAD83 (read from lat_nad83/long_nad83 attributes); ~1 m from WGS84"}),
        lon=(("station",), agg.lon.astype("float64").values,
             {"standard_name":"longitude","long_name":"longitude of well",
              "units":"degrees_east","axis":"X",
              "comment":"datum NAD83"}),
        time=(("obs",), obs.time.values,
              {"standard_name":"time","long_name":"date of water-level measurement"}),
    ),
    attrs=dict(
        Conventions="CF-1.8",
        featureType="point",
        title="High Plains Aquifer groundwater well water levels (USGS McGuire 2017/2019 campaigns)",
        institution="USGS (source data); packaged by MSU_GWB pipeline",
        summary=("Point observations of depth-to-water for High Plains aquifer wells, "
                 "compiled from USGS water-level change releases: F02 primary predev->2019 "
                 "network, F03 supplemental network (incl. 1980 epoch), F05 2017->2019 network."),
        source=("G:\\USGS GW dataset\\F02_hpwicpd19pt_Water-level change data ... ; "
                "F03_hpw1cpd19sp_Supplemental ... ; F05_hpwlcp1719pt_... (shapefiles, EPSG:5070; "
                "native lat_nad83/long_nad83 attribute columns used for coordinates)"),
        history=("Created by preprocess script from USGS shapefiles; sentinels -9999/-999 -> "
                 "NaN; measurement-date blanks in some epochs filled with mid-season nominal date "
                 "(documented approximation)."),
        geospatial_lat_min=f"{agg.lat.min():.4f}", geospatial_lat_max=f"{agg.lat.max():.4f}",
        geospatial_lon_min=f"{agg.lon.min():.4f}", geospatial_lon_max=f"{agg.lon.max():.4f}",
        time_coverage_start=str(obs.time.min().date()),
        time_coverage_end=str(obs.time.max().date()),
        license="USGS data are public domain",
        note="Derived file; original USGS zips/shapefiles unmodified.",
    ),
)

# CF time encoding
enc = {
    "time": {"units": "days since 1900-01-01 00:00:00", "calendar": "gregorian",
             "dtype": "float64", "_FillValue": None},
    "lat": {"dtype":"float64"}, "lon": {"dtype":"float64"},
    "water_level_ft": {"dtype":"float32","_FillValue": np.float32(9.96921e36)},
}
ds.to_netcdf(OUT, engine="netcdf4", format="NETCDF4", encoding=enc)
log(f"wrote {OUT} ({OUT.stat().st_size/1e6:.2f} MB)")

# ------------------------------------------------------------------
# verification pass: reopen and sanity check
# ------------------------------------------------------------------
chk = xr.open_dataset(OUT, decode_cf=True)
assert chk.attrs.get("Conventions") == "CF-1.8"
assert chk.attrs.get("featureType") == "point"
assert "days since 1900" in str(chk.time.encoding.get("units",""))
n_st = chk.sizes["station"]; n_ob = chk.sizes["obs"]
tmin = str(pd.to_datetime(chk.time.values.min()).date())
tmax = str(pd.to_datetime(chk.time.values.max()).date())
si_ok = bool(((chk.station_index.values >= 0) & (chk.station_index.values < n_st)).all())
wl_nan = float(np.isnan(chk.water_level_ft.values).mean())
print(chk)
print(f"sanity: stations={n_st}, obs={n_ob}, time {tmin}..{tmax}, "
      f"station_index_in_range={si_ok}, water_level NaN frac={wl_nan:.3%}")
print("sample:", chk.isel(obs=slice(0,3)).to_dataframe()[["station_index","campaign","water_level_ft","time"]])
chk.close()
log("VERIFY OK")
