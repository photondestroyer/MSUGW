"""
Process gSSURGO CONUS gdb -> REAL Plant-Available-Water (PAW) on the gridMET lab grid
(1/24 degree), saved as a CF-1.8 compliant NetCDF.

Reprocessed October 2026 (data audit of 2026-10-20): the target grid is now the grid of
GRIDMET_Ogallala_lab_latest.nc (287 x 233, lat/lon dims, cells match the gridMET/drought files
1:1) instead of the old Earth Engine gridMET grid (287 x 181); the mask is the reprocessed
aquifer_mask_4km_reprocessed_<date>.nc; the output name carries the generation date. The method
(120 m mukey samples, Valu1.rootznaws join, mean of samples per cell) is unchanged.

Chain (sources strictly READ-ONLY):
  1. gSSURGO CONUS raster `MURASTER_30m` (30 m mukey grid, EPSG:5070), read at overview 4 (~120 m,
     nearest) and joined to Valu1.rootznaws (USDA-NRCS "Root Zone Available Water Storage", cm).
     That join was stored on 2026-08-26 as `rootznaws` in
     merged_datasets/USGS data/derived_usgs/GSSURGO_Ogallala_all_attrs_120m.nc. The geodatabase on
     G: has since lost its catalog files and cannot be opened, so this script reads the stack.
  2. (The 120 m samples are the stack's `rootznaws` cells; it was clipped to the aquifer polygon.)
  4. Block-mean aggregate of rootznaws onto the gridMET lab 287x233 cells
     (area-weighted mean of all valid 120 m samples inside each cell).
  5. Mask to High Plains aquifer boundary (aquifer_mask_4km_reprocessed_<date>.nc).
  6. Write CF-1.8 NetCDF with full provenance.

Output: merged_datasets/USGS data/derived_usgs/PAW_rootznaws_4km_reprocessed_<date>.nc
        (+ .png quicklook in tbi_maps_v3_usgs/)
"""
import warnings
warnings.filterwarnings("ignore")
import json, time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import netCDF4
import rasterio

MERGED = Path(r"G:/MSU_GWB/datasets/merged_datasets")
DERIVED = MERGED / "USGS data" / "derived_usgs"
V3DIR = Path(r"G:/MSU_GWB/datasets/tbi_maps_v3_usgs")
V3DIR.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc)
GEN = NOW.strftime("%Y%m%d")
OUT_NC = DERIVED / f"PAW_rootznaws_4km_reprocessed_{GEN}.nc"

STACK = DERIVED / "GSSURGO_Ogallala_all_attrs_120m.nc"
OVERVIEW = 4            # 30 m base * 4 = ~120 m effective sampling (of the stack)
def log(m): print(f"[paw] {m}", flush=True)

t_start = time.time()

# ------------------------------------------------------------------
# 1. Common grid + aquifer mask
# ------------------------------------------------------------------
ds_gridmet = xr.open_dataset(MERGED / "GRIDMET_Ogallala_lab_latest.nc",
                             chunks={"time": 1}, mask_and_scale=False)
ty = ds_gridmet["lat"].values; tx = ds_gridmet["lon"].values
H, W = len(ty), len(tx)
res_y = abs(float(ty[1]-ty[0])); res_x = abs(float(tx[1]-tx[0]))
dst_transform = rasterio.transform.from_origin(
    float(tx[0])-res_x/2, float(ty[0])+res_y/2, res_x, res_y)

MASK_NC = sorted(DERIVED.glob("aquifer_mask_4km_reprocessed_*.nc"))[-1]   # newest reprocessed mask
mask_da = xr.open_dataset(MASK_NC)["aquifer_mask"]
assert np.array_equal(mask_da["lat"].values, ty) and np.array_equal(mask_da["lon"].values, tx), "mask grid differs"
mask_np = mask_da.values.astype("uint8")
log(f"grid {W}x{H}; aquifer px {int(mask_np.sum())} (mask {MASK_NC.name})")

# ------------------------------------------------------------------
# 2. 120 m rootznaws samples from the gSSURGO stack
# ------------------------------------------------------------------
# The gSSURGO_CONUS.gdb on G: is incomplete since Sep 2026 (6 table files, no catalog: GDAL cannot
# open it), so the MURASTER_30m -> Valu1.rootznaws join cannot be repeated. The same join was
# stored on 2026-08-26 in GSSURGO_Ogallala_all_attrs_120m.nc (mukey read at overview 4 = 120 m,
# nearest, joined to Valu1.rootznaws), which is exactly the sample set this script used before.
# Difference: that stack was clipped to the aquifer polygon, so a boundary 4 km cell now averages
# only the samples inside the polygon.
t0 = time.time()
with netCDF4.Dataset(STACK) as st:
    st.set_auto_maskandscale(False)
    v = st["rootznaws"]
    fillv = np.float32(v._FillValue)
    paw_src = v[:].astype("float32")
    gt = [float(x) for x in str(st["crs"].comment).split("(")[1].split(")")[0].split(",")]
    stack_history = str(st.history)
paw_src = np.where(paw_src == fillv, np.nan, paw_src)
dec_transform = rasterio.Affine(gt[0], gt[1], gt[2], gt[3], gt[4], gt[5])
out_h, out_w = paw_src.shape
valid = np.isfinite(paw_src)
log(f"read {STACK.name}: {out_w}x{out_h} samples of {abs(gt[0]):.0f} m in {time.time()-t0:.1f}s; "
    f"valid samples {valid.sum()/valid.size:.1%}; rootznaws cm min {np.nanmin(paw_src):.2f} "
    f"med {np.nanmedian(paw_src):.2f} max {np.nanmax(paw_src):.2f}")

# ------------------------------------------------------------------
# 4. aggregate to 4 km grid (mean of valid samples per cell)
# ------------------------------------------------------------------
xs_c = dec_transform.c + (np.arange(out_w)+0.5)*dec_transform.a
ys_c = dec_transform.f + (np.arange(out_h)+0.5)*dec_transform.e
from rasterio.warp import transform as warp_xy
# sums/counts per target cell, built in row blocks so memory stays at one block of samples
sums = np.zeros(H*W, dtype="float64"); cnts = np.zeros(H*W, dtype="float64")
ROWS = 500
for r0 in range(0, out_h, ROWS):
    r1 = min(out_h, r0 + ROWS)
    X, Y = np.meshgrid(xs_c, ys_c[r0:r1])
    XL, YL = warp_xy("EPSG:5070", "EPSG:4326", X.ravel(), Y.ravel())
    XL = np.asarray(XL).reshape(X.shape); YL = np.asarray(YL).reshape(Y.shape)
    col_t = np.floor((XL - (tx[0]-res_x/2)) / res_x).astype(np.int64)
    row_t = np.floor(((ty[0]+res_y/2) - YL) / res_y).astype(np.int64)
    inb = (row_t >= 0) & (row_t < H) & (col_t >= 0) & (col_t < W) & valid[r0:r1]
    flat = row_t[inb]*W + col_t[inb]
    sums += np.bincount(flat, weights=paw_src[r0:r1][inb].astype("float64"), minlength=H*W)
    cnts += np.bincount(flat, minlength=H*W)
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

crs_var = xr.DataArray(np.int32(0), attrs=dict(ds_gridmet["crs"].attrs))   # the gridMET lab grid mapping

da = xr.DataArray(
    paw_grid.reshape(1, H, W),
    coords={"time": np.array([np.datetime64("2019-07-01")]),
            "lat": ty, "lon": tx},
    dims=("time","lat","lon"), name="root_zone_available_water_storage")

ds = da.to_dataset()
ds["crs"] = crs_var
ds["root_zone_available_water_storage"] = ds["root_zone_available_water_storage"].assign_attrs(
    long_name="Plant available water stored in the soil root zone (PAW stock)",
    standard_name="soil_moisture_content_at_field_capacity",  # nearest CF token; see comment
    units="cm",
    grid_mapping="crs",
    coordinates="time lat lon",
    cell_methods="area: mean",
    comment=("SSURGO 'root zone available water storage' (Valu1.rootznaws) joined by mapunit key; "
             "depth-integrated to the soil-mapunit root-zone limiting depth. Units cm of water. "
             "standard_name token approximates CF vocabulary; semantic definition is SSURGO AWS."),
    source_raster=("gSSURGO_CONUS.gdb raster MURASTER_30m (30 m mukey grid, EPSG:5070) at overview 4 (~120 m), "
                   "joined to Valu1.rootznaws on 2026-08-26"),
    source_file=f"{STACK.relative_to(MERGED).as_posix()} variable rootznaws",
    attribute_table="gSSURGO_CONUS.gdb Valu1 table (mukey -> rootznaws)",
    aggregation=("area-mean of the mapped 120 m rootznaws samples whose centre lies in each cell, "
                 "then aquifer mask; the 120 m stack was already clipped to the aquifer polygon"),
    valid_px=int(np.isfinite(paw_grid).sum()),
    samples_per_cell_median=float(np.nanmedian(n_per)),
)
ds.attrs.update(dict(
    Conventions="CF-1.8",
    title="gSSURGO Root-Zone Available Water Storage (PAW) aggregated to Ogallala 4 km analysis grid",
    institution="USDA-NRCS (source data); packaged by MSU_GWB pipeline",
    summary=("Real PAW for VVIP Soil Buffer Index: gSSURGO MURASTER_30m mukey grid joined to "
             "Valu1.rootznaws, averaged onto the gridMET lab 287x233 1/24-degree grid and masked "
             "to the High Plains aquifer boundary."),
    iso_metadata_titles=" | ".join(meta.get("titles", [])[:3]),
    iso_metadata_abstract_snippet=meta.get("abstract_snippet","")[:400],
    iso_bounding_box=json.dumps(meta.get("bounding_boxes",[{}])[0]),
    source_stack_history=stack_history,
    history=f"{NOW:%Y-%m-%dT%H:%M:%SZ} created by local data processing/process_gssurgo_paw.py; sources read-only.",
    date_created=NOW.strftime("%Y-%m-%dT%H:%M:%SZ"),
    generation_date=GEN,
    processing_version="reprocessed-2026-10 (audit 2026-10-20: gridMET lab grid)",
    aquifer_mask=MASK_NC.name,
    note="Derived file; original gdb unmodified.",
))
enc = {"root_zone_available_water_storage":
       {"dtype":"float32", "zlib":True, "complevel":3,
        "_FillValue": np.float32(9.96921e36)},
       "lat": {"_FillValue": None}, "lon": {"_FillValue": None},
       "time": {"units": "days since 1970-01-01", "dtype": "int32"}}
ds["lat"].attrs.update(units="degrees_north", standard_name="latitude", axis="Y")
ds["lon"].attrs.update(units="degrees_east", standard_name="longitude", axis="X")
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

ds_gridmet.close()
print(f"PAW PREP DONE in {time.time()-t_start:.1f}s")
