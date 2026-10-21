"""
Preprocess raw USGS GW datasets -> the gridMET lab grid (1/24 degree, EPSG:4326).

Reprocessed October 2026 (data audit of 2026-10-20). Changes from the first version:
  * Target grid = GRIDMET_Ogallala_lab_latest.nc (287 x 233 cells of 1/24 degree, the Climatology
    Lab's own grid). The first version used the old Earth Engine gridMET grid (287 x 181, x step
    0.05358 degree), which no longer lines up with the gridMET/drought files.
  * Output dims are lat/lon with the same values as the gridMET lab files, so cells match 1:1.
  * File names carry the generation date (..._reprocessed_YYYYMMDD.nc); date_created and
    generation_date are global attributes. The old files are left untouched.
  * Metadata: CF-1.8 Conventions, a latitude_longitude grid mapping, numeric flag_values for
    the mask, no _FillValue on coordinates, `source` points to the USGS release files.
  * The aquifer mask is applied to K too (the old K file had 427 cells outside the mask).
  * Sources are the in-repo copies in merged_datasets/USGS_HPA_native/ (built cell for cell from
    USGS's hp_wlcpd19t.tif, hp_wlc1719t.tif and ofr98-548.e00.gz by build_usgs_hpa_native.py);
    the raw release is no longer on G:. Those rasters are already aquifer-masked at 500 m, so a
    4 km cell on the boundary now averages only its aquifer part.
  * K uses each zone's class limits from its RANGE text (as in USGS_HPA_native) and leaves zones
    without a class (RANGE -1, 37 polygons) empty; the first version painted them as K = 0.
Methods otherwise unchanged: water-level change = area average (rasterio 'average') of the 500 m
rasters; K = class midpoint painted by cell centre (largest zones first, smaller zones
overwrite); mask = cell centre inside hp_bound2010.

READ-ONLY on sources. Outputs (new files only) in
G:/MSU_GWB/datasets/merged_datasets/USGS data/derived_usgs/:
   dWL_predev_to_2019_ft_4km_reprocessed_<date>.nc   (lat,lon) float32, feet
   dWL_2017_to_2019_ft_4km_reprocessed_<date>.nc     (lat,lon) float32, feet
   K_hydraulic_conductivity_mday_4km_reprocessed_<date>.nc  (lat,lon) float32, m/day
   aquifer_mask_4km_reprocessed_<date>.nc            (lat,lon) uint8 1=inside hp_bound2010
   preprocess_log_reprocessed_<date>.txt

Run with the MSUGW conda env (rasterio, geopandas, netCDF4).
"""
import warnings
warnings.filterwarnings("ignore")

from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import netCDF4
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import Affine
import geopandas as gpd
from rasterio.features import rasterize

MERGED = Path(r"G:/MSU_GWB/datasets/merged_datasets")
NATIVE_NC = MERGED / "USGS_HPA_native" / "USGS_HPA_water_level_change_500m.nc"
K_GPKG = MERGED / "USGS_HPA_native" / "USGS_HPA_hydraulic_conductivity_polygons.gpkg"
OUT = MERGED / "USGS data" / "derived_usgs"
OUT.mkdir(parents=True, exist_ok=True)
GRIDMET = MERGED / "GRIDMET_Ogallala_lab_latest.nc"
BOUND_SHP = Path(r"G:/MSU_GWB/datasets/HPA_polygon/hp_bound2010.shp")

NOW = datetime.now(timezone.utc)
GEN = NOW.strftime("%Y%m%d")
STAMP = NOW.strftime("%Y-%m-%dT%H:%M:%SZ")
PROCESSING_VERSION = "reprocessed-2026-10 (audit 2026-10-20: gridMET lab grid)"

log_lines = []
def log(msg):
    print(f"[prep] {msg}")
    log_lines.append(msg)

# ------------------------------------------------------------------
# 0. Target grid = the gridMET lab grid (cell centres and GeoTransform from the file)
# ------------------------------------------------------------------
with netCDF4.Dataset(GRIDMET) as g:
    LAT = g["lat"][:].astype("f8"); LON = g["lon"][:].astype("f8")
    LAT_BNDS = g["lat_bnds"][:].astype("f8"); LON_BNDS = g["lon_bnds"][:].astype("f8")
    CRS_ATTRS = {k: g["crs"].getncattr(k) for k in g["crs"].ncattrs()}
    gt = [float(v) for v in str(g["crs"].GeoTransform).split()]
dst_transform = Affine(gt[1], gt[2], gt[0], gt[4], gt[5], gt[3])
H, W = len(LAT), len(LON)
assert np.allclose(LON, gt[0] + (np.arange(W) + 0.5) * gt[1]) and np.allclose(LAT, gt[3] + (np.arange(H) + 0.5) * gt[5])
log(f"Target grid = {GRIDMET.name}: {W} x {H} cells, step {gt[1]:.9f} x {gt[5]:.9f} deg, origin ({gt[0]}, {gt[3]})")


def write_grid(path, var, data, attrs, dtype, fill, title, extra_global=None):
    """Write one 2-D layer on the gridMET lab grid as a CF-1.8 NetCDF."""
    with netCDF4.Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("lat", H); ds.createDimension("lon", W); ds.createDimension("nv", 2)
        la = ds.createVariable("lat", "f8", ("lat",)); la[:] = LAT
        la.units = "degrees_north"; la.standard_name = "latitude"; la.axis = "Y"; la.bounds = "lat_bnds"
        lo = ds.createVariable("lon", "f8", ("lon",)); lo[:] = LON
        lo.units = "degrees_east"; lo.standard_name = "longitude"; lo.axis = "X"; lo.bounds = "lon_bnds"
        ds.createVariable("lat_bnds", "f8", ("lat", "nv"))[:] = LAT_BNDS
        ds.createVariable("lon_bnds", "f8", ("lon", "nv"))[:] = LON_BNDS
        crs = ds.createVariable("crs", "i4")
        for k, v in CRS_ATTRS.items():
            setattr(crs, k, v)
        v = ds.createVariable(var, dtype, ("lat", "lon"), zlib=True, complevel=4,
                              fill_value=None if fill is None else np.dtype(dtype).type(fill))
        v[:] = data
        for k, val in attrs.items():
            setattr(v, k, val)
        v.grid_mapping = "crs"
        ds.Conventions = "CF-1.8"
        ds.title = title
        ds.institution = "U.S. Geological Survey (source data); packaged by the MSU groundwater-buffering project"
        ds.date_created = STAMP
        ds.generation_date = GEN
        ds.processing_version = PROCESSING_VERSION
        ds.history = f"{STAMP} created by local data processing/preprocess_usgs_gw.py; sources read-only"
        ds.grid = f"gridMET lab grid of {GRIDMET.name} (Climatology Lab 1/24 degree, cells match 1:1)"
        for k, val in (extra_global or {}).items():
            setattr(ds, k, val)
    log(f"wrote {path.name}")


def warp_native(var):
    """Area-average one 500 m layer of the native USGS file onto the target grid (NaN where no data)."""
    with netCDF4.Dataset(NATIVE_NC) as n:
        n.set_auto_maskandscale(False)
        a = n[var][:].astype("float32")
        fill = np.float32(n[var]._FillValue)
        gt5070 = [float(v) for v in str(n["crs"].GeoTransform).split()]
        src_crs = n["crs"].crs_wkt
        src_file = n[var].source_file
    a = np.where(a == fill, np.nan, a)
    src_tr = Affine(gt5070[1], gt5070[2], gt5070[0], gt5070[4], gt5070[5], gt5070[3])
    dst = np.full((H, W), np.nan, dtype="float32")
    reproject(source=a, destination=dst, src_transform=src_tr, src_crs=src_crs,
              dst_transform=dst_transform, dst_crs="EPSG:4326",
              resampling=Resampling.average, src_nodata=np.nan, dst_nodata=np.nan)
    return dst, src_file


# ------------------------------------------------------------------
# 1. Aquifer boundary mask (cell centre inside hp_bound2010)
# ------------------------------------------------------------------
bnd = gpd.read_file(BOUND_SHP).to_crs("EPSG:4326")
mask = rasterize(((g, 1) for g in bnd.geometry), out_shape=(H, W),
                 transform=dst_transform, fill=0, dtype="uint8", all_touched=False)
inside = mask == 1
log(f"aquifer mask: {int(mask.sum())} cells inside ({100 * mask.mean():.1f}% of grid)")
write_grid(OUT / f"aquifer_mask_4km_reprocessed_{GEN}.nc", "aquifer_mask", mask,
           {"long_name": "cell centre inside the High Plains aquifer boundary (hp_bound2010)",
            "flag_values": np.array([0, 1], dtype="uint8"), "flag_meanings": "outside inside",
            "units": "1", "source": str(BOUND_SHP)},
           "uint8", None, "High Plains aquifer mask on the gridMET lab grid")

# ------------------------------------------------------------------
# 2. Water-level change (area average of the 500 m rasters), masked to the aquifer
# ------------------------------------------------------------------
for var, long_name, base in (
        ("dWL_predev_to_2019", "Mapped water-level change, predevelopment (~1950) to 2019", "dWL_predev_to_2019_ft_4km"),
        ("dWL_2017_to_2019", "Mapped water-level change, 2017 to 2019", "dWL_2017_to_2019_ft_4km")):
    arr, src_file = warp_native(var)
    arr = np.where(inside, arr, np.nan).astype("float32")
    write_grid(OUT / f"{base}_reprocessed_{GEN}.nc", var, arr,
               {"long_name": long_name, "units": "ft",
                "source": f"{NATIVE_NC.relative_to(MERGED).as_posix()} variable {var} (USGS {src_file}, 500 m, EPSG:5070)",
                "cell_methods": "area: mean",
                "regrid": "rasterio.warp.reproject Resampling.average onto the gridMET lab grid, then aquifer mask",
                "comment": ("Negative = decline. USGS's contour-based maps are exactly 0 over much of the aquifer; "
                            "zeros are data. The 500 m source is aquifer-masked, so boundary cells average their "
                            "aquifer part only.")},
               "float32", np.nan, f"USGS {long_name} on the gridMET lab grid")
    v = arr[np.isfinite(arr)]
    log(f"{base}: valid {v.size} cells, min {v.min():.2f} max {v.max():.2f} mean {v.mean():.2f} ft")

# ------------------------------------------------------------------
# 3. Hydraulic conductivity: class midpoint of each OFR 98-548 zone
# ------------------------------------------------------------------
pal = gpd.read_file(K_GPKG, layer="hydraulic_conductivity")
n_noclass = int((pal["K_lower_ftday"] < 0).sum())
pal = pal[pal["K_lower_ftday"] >= 0].copy()          # RANGE -1: USGS gives no class -> left empty
pal["K_mid_md"] = (pal["K_lower_ftday"] + pal["K_upper_ftday"]) / 2.0 * 0.3048      # ft/day -> m/day
pal_wgs = pal.to_crs("EPSG:4326").sort_values("AREA", ascending=False)
log(f"K: {len(pal)} zones with a class, {n_noclass} without (left empty), "
    f"{int(pal['limits_conflict'].sum())} whose RANGE text disagrees with MAJOR1/MINOR1 (RANGE used)")
K = rasterize(((geom, float(val)) for geom, val in zip(pal_wgs.geometry, pal_wgs["K_mid_md"])),
              out_shape=(H, W), transform=dst_transform, fill=np.nan, dtype="float32", all_touched=False)
log(f"K before mask: {int(np.isfinite(K).sum())} cells, of which {int((np.isfinite(K) & ~inside).sum())} outside the aquifer")
K = np.where(inside, K, np.nan).astype("float32")
write_grid(OUT / f"K_hydraulic_conductivity_mday_4km_reprocessed_{GEN}.nc", "K_mday", K,
           {"long_name": "Horizontal hydraulic conductivity (class midpoint)", "units": "m day-1",
            "source": f"{K_GPKG.relative_to(MERGED).as_posix()} (USGS OFR 98-548, Cederstrand & Becker 1998, 1:1,000,000)",
            "comment": ("Class limits from each zone's RANGE text (ft/day), midpoint x0.3048 to m/day; painted by cell "
                        "centre, largest zones first so smaller zones overwrite; zones without a class (RANGE -1) and "
                        "cells outside the aquifer are empty.")},
           "float32", np.nan, "USGS hydraulic conductivity zones on the gridMET lab grid")
v = K[np.isfinite(K)]
log(f"K: valid {v.size} cells, min {v.min():.2f} mean {v.mean():.2f} p90 {np.percentile(v, 90):.2f} m/d")

with open(OUT / f"preprocess_log_reprocessed_{GEN}.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log_lines))
print("PREP DONE")
