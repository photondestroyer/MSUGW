"""
Export ALL gSSURGO attribute data clipped to the Ogallala aquifer as ONE
CF-1.8 compliant NetCDF.

Sources (READ-ONLY):
  G:/USGS GW dataset/gSSURGO_CONUS.gdb
      - embedded raster MURASTER_30m  (30 m mukey grid, EPSG:5070, pyramids)
      - Valu1 attribute table         (56 thematic columns + mukey, 304,834 rows)
  G:/MSU_GWB/datasets/high_plains_quifer/hp_bound2010.shp   (aquifer clip polygon)

Method:
  1. Window-read MURASTER_30m over the aquifer bounding box at OVERVIEW=4
     (~120 m effective; nearest keeps exact mukeys).
  2. Rasterize hp_bound2010 polygon onto the same 120 m grid -> clip mask.
  3. For every Valu1 column build a mukey->value LUT, map onto the grid,
     apply clip mask, stream-write one compressed variable per column.
  4. CF projected-grid encoding: 1-D x/y in metres (EPSG:5070),
     `crs` grid_mapping variable = albers_conical_equal_area (NAD83).

Output: derived_usgs/GSSURGO_Ogallala_all_attrs_120m.nc   (NEW file only)
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
from rasterio.warp import transform_bounds, transform as warp_xy
from rasterio.features import rasterize
import geopandas as gpd
import netCDF4

GDB       = Path(r"G:/USGS GW dataset/gSSURGO_CONUS.gdb")
BOUND_SHP = Path(r"G:/MSU_GWB/datasets/high_plains_quifer/hp_bound2010.shp")
OUT_NC    = Path(r"G:/MSU_GWB/datasets/derived_usgs/GSSURGO_Ogallala_all_attrs_120m.nc")
OUT_NC.parent.mkdir(parents=True, exist_ok=True)

OVERVIEW = 4                      # 30 m base -> ~120 m effective
FILL_F   = np.float32(9.96921e36) # CF default-style float fill
FILL_I   = np.int32(-9999)

def log(m): print(f"[gssurgo-all] {m}", flush=True)

t_start = time.time()

# ------------------------------------------------------------------
# 1. Clip polygon -> ROI bounds (EPSG:5070) and 120 m clip mask later
# ------------------------------------------------------------------
bnd = gpd.read_file(BOUND_SHP)
if bnd.crs is None or bnd.crs.to_epsg() != 5070:
    bnd = bnd.to_crs("EPSG:5070")
bminx, bminy, bmaxx, bmaxy = bnd.total_bounds
log(f"aquifer bounds 5070: {bminx:.0f} {bminy:.0f} {bmaxx:.0f} {bmaxy:.0f}")

# ------------------------------------------------------------------
# 2. Window read of mukey raster
# ------------------------------------------------------------------
src = rasterio.open(GDB, driver="OpenFileGDB")
pad = 2000.0                       # 2 km padding around bbox
roi = (bminx-pad, bminy-pad, bmaxx+pad, bmaxy+pad)
win = win_from_bounds(*roi, transform=src.transform).round_offsets().round_lengths()
win = win.intersection(rasterio.windows.Window(0, 0, src.width, src.height))
out_w = max(1, int(win.width  // OVERVIEW))
out_h = max(1, int(win.height // OVERVIEW))
log(f"window {int(win.width)}x{int(win.height)} full-res px -> "
    f"{out_w}x{out_h} @~{30*OVERVIEW} m")

t0=time.time()
mukey = src.read(1, window=win, out_shape=(out_h,out_w),
                 resampling=Resampling.nearest).astype(np.int32)
log(f"read ok in {time.time()-t0:.1f}s")

dec_transform = rasterio.transform.from_origin(
    win.col_off*src.transform.a + src.transform.c,
    win.row_off*src.transform.e + src.transform.f,
    abs(src.transform.a)*OVERVIEW, abs(src.transform.e)*OVERVIEW)

# ------------------------------------------------------------------
# 3. Clip mask from shapefile on this exact grid
# ------------------------------------------------------------------
clip_mask = rasterize(((g,1) for g in bnd.geometry), out_shape=(out_h,out_w),
                      transform=dec_transform, fill=0, dtype="uint8").astype(bool)
n_inside = int(clip_mask.sum())
log(f"clip mask px inside aquifer: {n_inside} ({n_inside/(out_h*out_w):.1%})")

# ------------------------------------------------------------------
# 4. Valu1 table -> LUTs
# ------------------------------------------------------------------
val = pyogrio.read_dataframe(GDB, layer="Valu1", read_geometry=False)
val["mukey"] = pd.to_numeric(val["mukey"], errors="coerce")
val = val.dropna(subset=["mukey"])
val = val[val.mukey >= 0]
val["mukey"] = val.mukey.astype(np.int64)
attr_cols = [c for c in val.columns if c != "mukey"]
log(f"Valu1 rows {len(val)}, attribute columns {len(attr_cols)}")

max_mk = int(max(val.mukey.max(), int(mukey.max())))
lut_idx = np.zeros(max_mk+1, dtype=bool)
lut_vals = {}
mk = val.mukey.to_numpy()
for c in attr_cols:
    v = pd.to_numeric(val[c], errors="coerce").to_numpy(dtype="float32")
    lut_vals[c] = v

# ------------------------------------------------------------------
# 5. Create CF NetCDF skeleton (netCDF4 direct for streaming writes)
# ------------------------------------------------------------------
if OUT_NC.exists():
    OUT_NC.unlink()
nc = netCDF4.Dataset(OUT_NC, "w", format="NETCDF4")

# dims
nc.createDimension("y", out_h)
nc.createDimension("x", out_w)

# coordinates
xv = nc.createVariable("x", "f8", ("x",))
yv = nc.createVariable("y", "f8", ("y",))
xs = dec_transform.c + (np.arange(out_w)+0.5)*dec_transform.a
ys = dec_transform.f + (np.arange(out_h)+0.5)*dec_transform.e
xv[:] = xs; yv[:] = ys
xv.standard_name = "projection_x_coordinate"
yv.standard_name = "projection_y_coordinate"
xv.units = yv.units = "m"
xv.axis = "X"; yv.axis = "Y"
xv.long_name = "easting (NAD83 Conus Albers)"; yv.long_name = "northing"

crsv = nc.createVariable("crs", "i4")
crsv.grid_mapping_name = "albers_conical_equal_area"
crsv.latitude_of_projection_origin = 23.0
crsv.longitude_of_central_meridian = -96.0
crsv.standard_parallel = np.array([29.5, 45.5], dtype="f8")
crsv.false_easting = 0.0
crsv.false_northing = 0.0
crsv.semi_major_axis = 6378137.0
crsv.inverse_flattening = 298.257222101
crsv.datum = "NAD83"
crsv.epsg_code = "EPSG:5070"
crsv.comment = "GeoTransform of this grid: " + str(tuple(dec_transform)[:6])

# mukey band first
mv = nc.createVariable("mukey", "i4", ("y","x"), zlib=True, complevel=4,
                       chunksizes=(512,512), fill_value=FILL_I)
mukey_masked = np.where(clip_mask, mukey, FILL_I)
mv[:] = mukey_masked
mv.long_name = "USDA-NRCS gSSURGO map unit key (mukey)"
mv.grid_mapping = "crs"; mv.coordinates = "y x"
mv.comment = ("Read from embedded raster MURASTER_30m at overview "
              f"{OVERVIEW} (~{30*OVERVIEW} m effective); pixels outside the "
              "High Plains aquifer polygon are fill.")
del mukey_masked

# global attrs (fill more after writing)
nc.Conventions = "CF-1.8"
nc.title = ("gSSURGO CONUS thematic attribute stack clipped to the High Plains "
            "(Ogallala) aquifer")
nc.institution = "USDA-NRCS (source data); packaged by MSU_GWB pipeline"
nc.summary = ("All numeric Valu1 attribute-table columns of gSSURGO CONUS joined "
              "through the mapunit key raster (MURASTER_30m) and clipped to "
              "hp_bound2010 aquifer polygon, on a uniform ~120 m Albers grid.")
nc.source = (f"{GDB} :: raster 'MURASTER_30m' + table 'Valu1' ; clip polygon "
             f"{BOUND_SHP}")
nc.history = (f"Created {time.strftime('%Y-%m-%d %H:%M:%S')} by "
              "make_gssurgo_ogallala_stack.py; sources strictly read-only.")
nc.geospatial_bounds_crs = "EPSG:5070"
wgs = transform_bounds("EPSG:5070", "EPSG:4326", *roi)
nc.geospatial_lon_min, nc.geospatial_lat_min, \
    nc.geospatial_lon_max, nc.geospatial_lat_max = [float(v) for v in wgs]
nc.overview_factor = OVERVIEW
nc.effective_resolution_m = 30*OVERVIEW
nc.pixels_inside_aquifer = np.int64(n_inside)

# ------------------------------------------------------------------
# 6. Stream every attribute column
# ------------------------------------------------------------------
report = []
for i, c in enumerate(attr_cols, start=1):
    t0=time.time()
    lut = np.full(max_mk+1, np.nan, dtype="float32")
    lut[mk] = lut_vals[c]
    grid = lut[mukey]                       # float32 map
    grid = np.where(clip_mask, grid, np.nan)
    var = nc.createVariable(c, "f4", ("y","x"), zlib=True, complevel=4,
                            chunksizes=(512,512), fill_value=FILL_F)
    var[:] = np.where(np.isfinite(grid), grid, FILL_F).astype("float32")
    var.long_name = f"gSSURGO Valu1.{c}"
    var.units = "cm" if c.startswith(("aws","rootzn")) else \
                ("g/m2" if c.startswith(("soc","tk")) else "1")
    var.units_comment = ("units follow SSURGO conventions: AWS/root-zone storage in cm "
                         "of water; SOC stocks g/m2 (per depth interval prefix); "
                         "tk thickness cm; index scores unitless. See SSURGO docs.")
    var.grid_mapping = "crs"; var.coordinates = "y x"
    vv = grid[np.isfinite(grid)]
    report.append((c, int(vv.size),
                   float(vv.min()) if vv.size else np.nan,
                   float(np.median(vv)) if vv.size else np.nan,
                   float(vv.max()) if vv.size else np.nan))
    del lut, grid, vv
    if i % 10 == 0 or i == len(attr_cols):
        log(f"  [{i}/{len(attr_cols)}] wrote '{c}' ({time.time()-t0:.1f}s)")

nc.close()
src.close()

size_mb = OUT_NC.stat().st_size/1e6
log(f"wrote {OUT_NC} ({size_mb:.1f} MB) in {time.time()-t_start:.0f}s")

# save write-report next to file
rep = pd.DataFrame(report, columns=["column","n_valid","min","median","max"])
rep.to_csv(OUT_NC.with_suffix(".report.csv"), index=False)
print(rep.to_string(index=False))

# ------------------------------------------------------------------
# 7. verification pass
# ------------------------------------------------------------------
chk = xr.open_dataset(OUT_NC, decode_cf=True)
assert chk.attrs.get("Conventions") == "CF-1.8"
need = ["mukey","x","y","crs"] + attr_cols
missing = [v for v in need if v not in chk]
rz = chk["rootznaws"].values
rz_valid = rz[np.isfinite(rz)]
print(f"verify: missing_vars={missing}")
print(f"        dims={dict(chk.sizes)}  vars={len(chk.data_vars)}")
print(f"        rootznaws valid px={rz_valid.size} ({100*rz_valid.size/rz.size:.1f}%), "
      f"min {rz_valid.min():.1f} med {np.median(rz_valid):.1f} max {rz_valid.max():.1f} cm "
      f"(4 km product median was 203.5)")
print(f"        mukey inside-aquifer px={int((chk['mukey'].values != FILL_I).sum())}")
print(f"        crs: {chk['crs'].grid_mapping_name}, epsg {chk['crs'].epsg_code}")
assert not missing and rz_valid.size > 0
chk.close()
log("VERIFY OK")
