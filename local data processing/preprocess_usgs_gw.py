"""
Preprocess raw USGS GW datasets -> common GRIDMET ~4km EPSG:4326 grid.

READ-ONLY on sources in "G:\\USGS GW dataset". All derived outputs written to
G:/MSU_GWB/datasets/derived_usgs/ (new files only).

Inputs processed
----------------
1. F01 hp_wlcpd19t.tif : Water-level change predevelopment(~1950)->2019, ft
   EPSG:5070 (NAD83 Conus Albers), 500 m, nodata=-3.4028e38
   Source: McGuire 2017-style USGS SIR raster (zip inside folder F01)
2. F04 hp_wlc1719t.tif : Water-level change 2017->2019, ft (same grid)
3. OFR98-548 .e00      : Digital map of hydraulic conductivity (K), High Plains
   AVCE00 vector, Albers, zones with RANGE strings "25 to 50" (ft/day).
   Layers ARC/CNT/LAB/PAL. We rasterize PAL polygons using midpoint of
   [MINOR1, MAJOR1] as K (ft/day -> m/day x0.3048).

Outputs (derived_usgs/)
-----------------------
   dWL_predev_to_2019_ft_4km.nc   (y,x) float32, feet
   dWL_2017_to_2019_ft_4km.nc     (y,x) float32, feet
   K_hydraulic_conductivity_mday_4km.nc  (y,x) m/day
   aquifer_mask_4km.nc            (y,x) uint8 1=inside hp_bound2010
   preprocess_log.txt

Common target grid = GRIDMET_Merged_Ogallala.nc coords (287x181, ~4 km).
"""
import warnings
warnings.filterwarnings("ignore")

import os, zipfile, gzip, shutil
from pathlib import Path
import numpy as np
import xarray as xr
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import array_bounds
import geopandas as gpd
from rasterio.features import rasterize

SRC_DIR = Path(r"G:/USGS GW dataset")
TMP = Path(r"C:/Users/AlienX/AppData/Local/Temp/opencode/usgs_raw")
OUT = Path(r"G:/MSU_GWB/datasets/derived_usgs")
OUT.mkdir(parents=True, exist_ok=True)

MERGED = Path(r"G:/MSU_GWB/datasets/merged_datasets")
BOUND_SHP = Path(r"G:/MSU_GWB/datasets/high_plains_quifer/hp_bound2010.shp")

log_lines = []
def log(msg):
    print(f"[prep] {msg}")
    log_lines.append(msg)

# ------------------------------------------------------------------
# 0. Common target grid from GRIDMET
# ------------------------------------------------------------------
ds_gridmet = xr.open_dataset(MERGED / "GRIDMET_Merged_Ogallala.nc",
                             chunks={"time": 1}, mask_and_scale=False)
ty, tx = ds_gridmet["y"].values, ds_gridmet["x"].values
res_y = abs(float(ty[1] - ty[0])); res_x = abs(float(tx[1] - tx[0]))
dst_transform = rasterio.transform.from_origin(
    float(tx[0]) - res_x/2, float(ty[0]) + res_y/2,
    res_x, res_y)
# note: ty decreasing; from_origin takes top-left corner. top edge center y[0]+res/2
H, W = len(ty), len(tx)
dst_bounds_4326 = array_bounds(H, W, dst_transform)
log(f"Target grid {W}x{H} res ({res_x:.5f},{res_y:.5f}) deg bounds {np.round(dst_bounds_4326,3)}")

def warp_to_gridmet(src_path, out_name, long_name, units, resamp=Resampling.average):
    """Reproject any raster to GRIDMET grid, return xr.DataArray."""
    with rasterio.open(src_path) as src:
        src_data = src.read(1)
        src_nodata = src.nodata if src.nodata is not None else -3.4028230607370965e+38
        src_data = np.where(src_data == src_nodata, np.nan, src_data)
        dst = np.full((H, W), np.nan, dtype="float32")
        reproject(
            source=src_data,
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=dst_transform,
            dst_crs="EPSG:4326",
            resampling=resamp,
            src_nodata=np.nan,
            dst_nodata=np.nan,
        )
    da = xr.DataArray(dst, coords={"y": ty, "x": tx}, dims=("y", "x"),
                      name=out_name,
                      attrs={"long_name": long_name, "units": units,
                             "source": str(src_path),
                             "source_crs": str(src.crs), "source_res_m": str(src.res),
                             "regrid": "rasterio.warp.reproject average to GRIDMET 4km EPSG:4326"})
    da.attrs["note"] = "Derived file; original USGS zip/tif unmodified."
    return da

# ------------------------------------------------------------------
# 1. Extract zips (to temp) if not already
# ------------------------------------------------------------------
def ensure_extract(zip_path, dest_dir):
    marker = dest_dir / ".done"
    if marker.exists():
        return next(dest_dir.rglob("*.tif"))
    dest_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as z:
        z.extractall(dest_dir)
    marker.write_text("ok")
    return next(dest_dir.rglob("*.tif"))

f01_zip = next(SRC_DIR.glob("F01_*/hp_wlcpd19t.zip"))
f04_zip = next(SRC_DIR.glob("F04_*/hp_wlc1719t.zip"))
f01_tif = ensure_extract(f01_zip, TMP / "f01")
f04_tif = ensure_extract(f04_zip, TMP / "f04")
log(f"F01 tif {f01_tif}")
log(f"F04 tif {f04_tif}")

dwl_pre = warp_to_gridmet(f01_tif, "dWL_predev_to_2019",
                          "Mapped water-level change, predevelopment (~1950) to 2019",
                          "feet", Resampling.average)
dwl_1719 = warp_to_gridmet(f04_tif, "dWL_2017_to_2019",
                           "Mapped water-level change, 2017 to 2019",
                           "feet", Resampling.average)
for da, nm in [(dwl_pre, "dWL_predev_to_2019_ft_4km"), (dwl_1719, "dWL_2017_to_2019_ft_4km")]:
    ds = da.to_dataset()
    ds.attrs.update(da.attrs)
    ds.to_netcdf(OUT / nm, engine="netcdf4")
    log(f"wrote {nm}: min {float(np.nanmin(da.values)):.2f} max {float(np.nanmax(da.values)):.2f} "
        f"mean {float(np.nanmean(da.values)):.2f} ft, NaN% {100*np.isnan(da.values).mean():.1f}")

# ------------------------------------------------------------------
# 2. Hydraulic conductivity from E00 (AVCE00 vector)
# ------------------------------------------------------------------
e00_src = next(SRC_DIR.glob("Digital map of hydraulic conductivity*/ofr98-548.e00.gz"))
e00 = TMP / "ofr98-548.e00"
if not e00.exists():
    e00.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(e00_src, "rb") as fin, open(e00, "wb") as fout:
        shutil.copyfileobj(fin, fout)
log(f"E00 {e00}")

pal = gpd.read_file(e00, layer="PAL")           # polygons w/ RANGE MAJOR1 MINOR1
lab = gpd.read_file(e00, layer="LAB")           # label points (sanity)
log(f"PAL polygons {len(pal)}, LAB labels {len(lab)}")

pal = pal[pal["MINOR1"] >= 0].copy()            # drop -1 background
pal["K_mid_ft"] = (pal["MAJOR1"] + pal["MINOR1"]) / 2.0
pal["K_mid_md"] = pal["K_mid_ft"] * 0.3048      # ft/day -> m/day
log(f"K classes ft/d: {sorted(pal['RANGE'].unique())}")
log(f"K mid m/d range {pal['K_mid_md'].min():.2f}-{pal['K_mid_md'].max():.2f}")

# Reproject polygons to 4326 and clip to GRIDMET bbox
pal_wgs = pal.to_crs("EPSG:4326")
xmin, ymin, xmax, ymax = dst_bounds_4326
pal_clip = pal_wgs.cx[xmin:xmax, ymin:ymax]
log(f"polygons after ROI clip: {len(pal_clip)}")
# Paint order: largest area first so small zones overwrite (mode-like approx)
pal_clip = pal_clip.sort_values("AREA", ascending=False)
shapes = ((geom, float(val)) for geom, val in zip(pal_clip.geometry, pal_clip["K_mid_md"]))
K_grid = rasterize(shapes, out_shape=(H, W), transform=dst_transform,
                   fill=np.nan, dtype="float32", all_touched=False)

K_da = xr.DataArray(K_grid, coords={"y": ty, "x": tx}, dims=("y", "x"),
                    name="K_mday",
                    attrs={"long_name": "Horizontal hydraulic conductivity (zone midpoint)",
                           "units": "m/day",
                           "source": "USGS OFR 98-548 digital map of hydraulic conductivity, High Plains aquifer (Cederstrand & Becker 1998), E00 PAL zones",
                           "note": ("Zone RANGE strings (ft/day) converted using interval midpoint, "
                                    "x0.3048 to m/day; rasterized onto GRIDMET 4km by area-descending paint "
                                    "(smaller zones overwrite larger). Original .e00.gz unmodified.")})
ds = K_da.to_dataset(); ds.attrs.update(K_da.attrs)
ds.to_netcdf(OUT / "K_hydraulic_conductivity_mday_4km.nc", engine="netcdf4")
log(f"K grid: valid px {int(np.isfinite(K_grid).sum())}/{H*W}, "
    f"mean {float(np.nanmean(K_grid)):.2f} m/d, p90 {float(np.nanpercentile(K_grid,90)):.2f}")

# ------------------------------------------------------------------
# 3. Aquifer boundary mask (hp_bound2010.shp) on GRIDMET grid
# ------------------------------------------------------------------
bnd = gpd.read_file(BOUND_SHP)
if bnd.crs is None or bnd.crs.to_epsg() != 4326:
    bnd = bnd.to_crs("EPSG:4326")
mask = rasterize(((g, 1) for g in bnd.geometry), out_shape=(H, W),
                 transform=dst_transform, fill=0, dtype="uint8")
mask_da = xr.DataArray(mask, coords={"y": ty, "x": tx}, dims=("y", "x"),
                       name="aquifer_mask",
                       attrs={"long_name": "High Plains aquifer boundary (hp_bound2010)",
                              "flag_values": "0 outside, 1 inside",
                              "units": "1"})
ds = mask_da.to_dataset(); ds.attrs.update(mask_da.attrs)
ds.to_netcdf(OUT / "aquifer_mask_4km.nc", engine="netcdf4")
log(f"aquifer mask: {int(mask.sum())} px inside ({100*mask.mean():.1f}% of grid)")

# Apply mask to derived layers (spatially consistent: aquifer-only)
for da, nm in [(dwl_pre, "dWL_predev_to_2019_ft_4km"), (dwl_1719, "dWL_2017_to_2019_ft_4km"),
               (K_da, "K_hydraulic_conductivity_mday_4km")]:
    masked = da.where(mask_da == 1)
    ds = masked.to_dataset(); ds.attrs.update(da.attrs)
    ds.to_netcdf(OUT / nm, mode="w", engine="netcdf4")
    log(f"masked {nm}: valid {int(np.isfinite(masked.values).sum())} px")

with open(OUT / "preprocess_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log_lines))

ds_gridmet.close()
print("PREP DONE")
