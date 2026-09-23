"""
Process gSSURGO CONUS gdb -> REAL Plant-Available-Water (PAW) on the common
GRIDMET ~4 km grid, saved as a CF-1.8 compliant NetCDF.

Chain (sources strictly READ-ONLY):
  1. Embedded raster `MURASTER_30m` (30 m mukey grid, EPSG:5070, int32,
     overviews [2..513]) opened via GDAL OpenFileGDB driver (rasterio).
  2. Window read over the Ogallala ROI at overview factor 4 (~120 m effective)
     -> mukey samples.
  3. Join mukey -> Valu1.rootznaws  ("Root Zone Available Water Storage",
     USDA-NRCS gSSURGO Valu1 table; units cm) via vectorized lookup.
  4. Block-mean aggregate of rootznaws onto the GRIDMET 287x181 cells
     (area-weighted mean of all valid 120 m samples inside each cell).
  5. Mask to High Plains aquifer boundary (derived_usgs/aquifer_mask_4km.nc).
  6. Write CF-1.8 NetCDF with full provenance.

Output: derived_usgs/PAW_rootznaws_4km.nc   (+ .png quicklook in tbi_maps_v3_usgs/)
"""
import warnings
warnings.filterwarnings("ignore")
import json, time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import pyogrio
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import from_bounds as win_from_bounds
from rasterio.warp import transform_bounds
from rasterio.features import rasterize
import geopandas as gpd

GDB = Path(r"G:/USGS GW dataset/gSSURGO_CONUS.gdb")
MERGED = Path(r"G:/MSU_GWB/datasets/merged_datasets")
DERIVED = Path(r"G:/MSU_GWB/datasets/derived_usgs")
V3DIR = Path(r"G:/MSU_GWB/datasets/tbi_maps_v3_usgs")
V3DIR.mkdir(parents=True, exist_ok=True)
OUT_NC = DERIVED / "PAW_rootznaws_4km.nc"

OVERVIEW = 4            # 30 m base * 4 = ~120 m effective sampling
def log(m): print(f"[paw] {m}", flush=True)

t_start = time.time()

# ------------------------------------------------------------------
# 1. Common grid + aquifer mask
# ------------------------------------------------------------------
ds_gridmet = xr.open_dataset(MERGED / "GRIDMET_Merged_Ogallala.nc",
                             chunks={"time": 1}, mask_and_scale=False)
ty = ds_gridmet["y"].values; tx = ds_gridmet["x"].values
H, W = len(ty), len(tx)
res_y = abs(float(ty[1]-ty[0])); res_x = abs(float(tx[1]-tx[0]))
dst_transform = rasterio.transform.from_origin(
    float(tx[0])-res_x/2, float(ty[0])+res_y/2, res_x, res_y)

mask_da = xr.open_dataset(DERIVED/"aquifer_mask_4km.nc", chunks="auto")["aquifer_mask"]
mask_np = mask_da.values.astype("uint8")
log(f"grid {W}x{H}; aquifer px {int(mask_np.sum())}")

roi_wgs = (float(tx[0])-res_x/2, float(ty[-1])-res_y/2,
           float(tx[-1])+res_x/2, float(ty[0])+res_y/2)   # (W, S, E, N)
roi_5070 = transform_bounds("EPSG:4326", "EPSG:5070", *roi_wgs)
log(f"ROI wgs {tuple(round(v,3) for v in roi_wgs)}")
log(f"ROI 5070 {tuple(round(v) for v in roi_5070)}")

# ------------------------------------------------------------------
# 2. Open embedded raster & windowed decimated read
# ------------------------------------------------------------------
src = rasterio.open(GDB, driver="OpenFileGDB")
log(f"raster {src.width}x{src.height} crs {src.crs} dt {src.dtypes} "
    f"nodata {src.nodata} tags {src.tags()}")
ovs = src.overviews(1)
assert OVERVIEW in ovs, f"overview {OVERVIEW} not in {ovs}"

win = win_from_bounds(*roi_5070, transform=src.transform).round_offsets().round_lengths()
full = rasterio.windows.Window(0, 0, src.width, src.height)
win = win.intersection(full)   # clamp to raster extent
log(f"window(full-res px): cols {win.col_off:.0f}+{win.width:.0f}, rows {win.row_off:.0f}+{win.height:.0f}")

out_h = max(1, int(win.height // OVERVIEW))
out_w = max(1, int(win.width // OVERVIEW))
log(f"decimated read at overview {OVERVIEW} -> {out_w}x{out_h} ({out_w*out_h/1e6:.0f} M px)")
t0=time.time()
mukey = src.read(1, window=win, out_shape=(out_h, out_w),
                 resampling=Resampling.nearest)      # ids: keep exact values
# transform OF THE DECIMATED ARRAY:
dec_transform = rasterio.transform.from_origin(
    win.col_off*src.transform.a + src.transform.c,
    win.row_off*src.transform.e + src.transform.f,
    abs(src.transform.a)*OVERVIEW, abs(src.transform.e)*OVERVIEW)
log(f"read done in {time.time()-t0:.1f}s; nodata={src.nodata}")
if src.nodata is not None:
    mukey = np.where(mukey == np.int32(src.nodata), np.int32(-1), mukey)

# ------------------------------------------------------------------
# 3. mukey -> rootznaws lookup (Valu1 table)
# ------------------------------------------------------------------
t0=time.time()
val = pyogrio.read_dataframe(GDB, layer="Valu1",
                             columns=["mukey", "rootznaws"], read_geometry=False)
val["mukey"] = pd.to_numeric(val["mukey"], errors="coerce")
val["rootznaws"] = pd.to_numeric(val["rootznaws"], errors="coerce")
val = val.dropna()
val = val[val.mukey >= 0]
val["mukey"] = val.mukey.astype(np.int64)
max_mk = int(np.maximum(val.mukey.max(), mukey.max()))
lut = np.full(max_mk+1, np.nan, dtype="float32")
lut[val.mukey.to_numpy(dtype=np.int64)] = val.rootznaws.to_numpy(dtype="float32")
log(f"Valu1 rows {len(val)}; lut built in {time.time()-t0:.1f}s "
    f"(rootznaws cm: min {np.nanmin(lut[lut>=0]):.2f} med {np.nanmedian(lut[lut>=0]):.2f} "
    f"max {np.nanmax(lut[lut>=0]):.2f})")

paw_src = lut[mukey.astype(np.int64)]          # NaN where unmapped/nodata
valid = np.isfinite(paw_src)
log(f"mapped source samples: {valid.sum()/valid.size:.1%}")

# ------------------------------------------------------------------
# 4. aggregate to 4 km grid (mean of valid samples per cell)
# ------------------------------------------------------------------
xs_c = dec_transform.c + (np.arange(out_w)+0.5)*dec_transform.a
ys_c = dec_transform.f + (np.arange(out_h)+0.5)*dec_transform.e
X, Y = np.meshgrid(xs_c, ys_c)
from rasterio.warp import transform as warp_xy
XL, YL = warp_xy("EPSG:5070", "EPSG:4326", X.ravel(), Y.ravel())
XL = np.asarray(XL).reshape(X.shape); YL = np.asarray(YL).reshape(Y.shape)
col_t = np.floor((XL - (tx[0]-res_x/2)) / res_x).astype(np.int64)
row_t = np.floor(((ty[0]+res_y/2) - YL) / res_y).astype(np.int64)
inb = (row_t >= 0) & (row_t < H) & (col_t >= 0) & (col_t < W) & valid

flat = row_t[inb]*W + col_t[inb]
vals = paw_src[inb]
sums = np.bincount(flat, weights=vals.astype("float64"), minlength=H*W)
cnts = np.bincount(flat, minlength=H*W).astype("float64")
paw_grid = np.full(H*W, np.nan, dtype="float32")
ok = cnts > 0
paw_grid[ok] = (sums[ok]/cnts[ok]).astype("float32")
paw_grid = paw_grid.reshape(H, W)
n_per = np.full(H*W, 0.0); n_per[ok]=cnts[ok]; n_per=n_per.reshape(H,W)
log(f"cells with PAW: {(np.isfinite(paw_grid)).sum()}; "
    f"samples/cell median {np.nanmedian(n_per):.0f}")

# mask to aquifer
paw_grid = np.where(mask_np == 1, paw_grid, np.nan)
v = paw_grid[np.isfinite(paw_grid)]
log(f"aquifer PAW stats (cm): min {v.min():.2f} p25 {np.percentile(v,25):.2f} "
    f"med {np.median(v):.2f} p75 {np.percentile(v,75):.2f} max {v.max():.2f}")

# ------------------------------------------------------------------
# 5. write CF-compliant NetCDF
# ------------------------------------------------------------------
meta = {}
mp = DERIVED/"gssurgo_iso19115_summary.json"
if mp.exists():
    meta = json.loads(mp.read_text(encoding="utf-8"))

crs_var = xr.DataArray(
    np.int32(4326),
    attrs=dict(grid_mapping_name="latitude_longitude",
               longitude_of_prime_meridian=0.0, semi_major_axis=6378137.0,
               inverse_flattening=298.257223563, epsg_code="EPSG:4326"))

da = xr.DataArray(
    paw_grid.reshape(1, H, W),
    coords={"time": np.array([np.datetime64("2019-07-01")]),
            "y": ty, "x": tx},
    dims=("time","y","x"), name="root_zone_available_water_storage")

ds = da.to_dataset()
ds["crs"] = crs_var
ds["root_zone_available_water_storage"] = ds["root_zone_available_water_storage"].assign_attrs(
    long_name="Plant available water stored in the soil root zone (PAW stock)",
    standard_name="soil_moisture_content_at_field_capacity",  # nearest CF token; see comment
    units="cm",
    grid_mapping="crs",
    coordinates="time y x",
    cell_methods="area: mean",
    comment=("SSURGO 'root zone available water storage' (Valu1.rootznaws) joined by mapunit key; "
             "depth-integrated to the soil-mapunit root-zone limiting depth. Units cm of water. "
             "standard_name token approximates CF vocabulary; semantic definition is SSURGO AWS."),
    source_raster="gSSURGO_CONUS.gdb embedded raster MURASTER_30m (30 m mukey grid, EPSG:5070)",
    attribute_table="gSSURGO_CONUS.gdb Valu1 table (mukey -> rootznaws)",
    aggregation=(f"window read at overview {OVERVIEW} (~{30*OVERVIEW} m effective), nearest resampling "
                 "for keys, then area-mean of mapped samples into each ~4 km cell; aquifer-masked"),
    valid_px=int(np.isfinite(paw_grid).sum()),
    samples_per_cell_median=float(np.nanmedian(n_per)),
)
ds.attrs.update(dict(
    Conventions="CF-1.8",
    title="gSSURGO Root-Zone Available Water Storage (PAW) aggregated to Ogallala 4 km analysis grid",
    institution="USDA-NRCS (source data); packaged by MSU_GWB pipeline",
    summary=("Real PAW for VVIP Soil Buffer Index: gSSURGO MURASTER_30m mukey grid joined to "
             "Valu1.rootznaws, averaged onto GRIDMET-matched 287x181 EPSG:4326 grid and masked "
             "to the High Plains aquifer boundary."),
    iso_metadata_titles=" | ".join(meta.get("titles", [])[:3]),
    iso_metadata_abstract_snippet=meta.get("abstract_snippet","")[:400],
    iso_bounding_box=json.dumps(meta.get("bounding_boxes",[{}])[0]),
    history="Created by process_gssurgo_paw.py; sources read-only.",
    note="Derived file; original gdb unmodified.",
))
enc = {"root_zone_available_water_storage":
       {"dtype":"float32", "zlib":True, "complevel":3,
        "_FillValue": np.float32(9.96921e36)}}
ds.to_netcdf(OUT_NC, engine="netcdf4", encoding=enc)
log(f"wrote {OUT_NC} ({OUT_NC.stat().st_size/1e6:.2f} MB)")

# quicklook png
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6.5,7))
extent=[float(tx[0]),float(tx[-1]),float(ty[-1]),float(ty[0])]
im=ax.imshow(np.ma.masked_invalid(paw_grid), extent=extent, origin="upper",
             cmap="YlGnBu", interpolation="nearest")
ax.set_title("REAL PAW from gSSURGO rootznaws (cm)\n30 m raster -> 120 m samples -> 4 km mean, aquifer-masked", fontsize=9)
ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
fig.colorbar(im, ax=ax, shrink=.85, pad=.02, label="PAW (cm)")
fig.tight_layout(); fig.savefig(V3DIR/"PAW_static.png", dpi=180); plt.close(fig)

src.close(); ds_gridmet.close()
print(f"PAW PREP DONE in {time.time()-t_start:.1f}s")
