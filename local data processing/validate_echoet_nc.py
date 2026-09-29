"""Sanity + corruption audit for ECHOET_Ogallala.nc (read-only, strided sampling)."""
import os, sys
import numpy as np
import netCDF4 as nc4
from datetime import datetime, timedelta

PATH = sys.argv[1] if len(sys.argv) > 1 else r"D:\Downloads\ECHOET_Ogallala.nc"

print("=" * 70)
print("FILE:", PATH)
print("SIZE: %.2f MiB" % (os.path.getsize(PATH) / 1024 ** 2))
ds = nc4.Dataset(PATH, "r")
print("FORMAT:", ds.file_format)
print("GLOBAL Conventions:", getattr(ds, "Conventions", "<missing>"))
for a in ("title", "source", "references", "history", "comment", "institution"):
    v = getattr(ds, a, "<missing>")
    print("  %s: %s" % (a, str(v)[:130]))
print("DIMENSIONS:", {k: len(v) for k, v in ds.dimensions.items()})
print("UNLIMITED:", [k for k, v in ds.dimensions.items() if v.isunlimited()])

problems, notes, improvements = [], [], []

# ---- variable inventory ----
print("\n[VARIABLES]")
for name, v in ds.variables.items():
    fill = getattr(v, "_FillValue", "<none>")
    print("  %-12s dims=%s dtype=%s fill=%s" % (name, v.dimensions, v.dtype, fill))
    for attr in ("units", "standard_name", "long_name", "calendar", "axis",
                 "coordinates", "grid_mapping", "cell_methods"):
        if hasattr(v, attr):
            print("      %s = %s" % (attr, getattr(v, attr)))
    try:
        print("      filters =", v.filters(), "| chunksizes =", v.chunking())
    except Exception as e:
        print("      filters/chunking unreadable:", e)

nt = len(ds.variables["time"])
ny = len(ds.variables["y"])
nx = len(ds.variables["x"])
print("\nGRID: nt=%d ny=%d nx=%d" % (nt, ny, nx))

# ---- time: hourly regularity ----
t = ds.variables["time"]
tvals = np.array(t[:], dtype=np.float64)
print("\n[TIME] n=%d units=%s calendar=%s" % (nt, t.units, getattr(t, "calendar", "?")))
dts = nc4.num2date(tvals, t.units, getattr(t, "calendar", "proleptic_gregorian"))
dts = [datetime(d.year, d.month, d.day, d.hour, d.minute, d.second) for d in dts]
print("  range: %s -> %s" % (dts[0], dts[-1]))
if np.any(np.diff(tvals) <= 0):
    problems.append("time axis not strictly increasing")
else:
    print("  monotonic: OK")
if len(set(tvals.tolist())) != nt:
    problems.append("duplicate timesteps present")
else:
    print("  duplicates: none")
diffs = np.diff(tvals)
step_unit = {"seconds": 1, "minutes": 60, "hours": 3600, "days": 86400}
unit = str(t.units).split()[0].lower().rstrip("s") + "s"
sec = diffs * step_unit.get(unit, float("nan"))
n_hourly = int((sec == 3600).sum())
print("  steps exactly 1h: %d/%d (%.1f%%)" % (n_hourly, len(diffs), 100.0 * n_hourly / max(1, len(diffs))))
if n_hourly != len(diffs):
    bad = np.where(sec != 3600)[0][:15]
    for i in bad:
        print("    gap %s -> %s = %.1f h" % (dts[i], dts[i + 1], sec[i] / 3600.0))
    ngaps = int((sec != 3600).sum())
    problems.append("%d/%d steps are not exactly 1-hourly (gaps/duplicates)" % (ngaps, len(diffs)))
    # missing-hours estimate
    span_h = (dts[-1] - dts[0]).total_seconds() / 3600.0
    print("  span %.1f h vs %d steps -> ~%d hourly slots missing" % (span_h, nt, int(round(span_h)) + 1 - nt))
else:
    print("  hourly continuity: PERFECT (1000 consecutive hours)")

# ---- coordinates ----
for c in ("y", "x"):
    a = ds.variables[c][:]
    print("\n[COORD %s] n=%d finite %d/%d range %.4f..%.4f monotonic_inc=%s" % (
        c, a.size, np.isfinite(a).sum(), a.size, np.nanmin(a), np.nanmax(a),
        bool(np.all(np.diff(a) > 0))))
    if not np.all(np.isfinite(a)):
        problems.append("non-finite values in coordinate %s" % c)
lat = ds.variables["lat"][:]
lon = ds.variables["lon"][:]
print("\n[LAT] shape=%s finite=%d/%d range %.3f..%.3f" % (
    lat.shape, np.isfinite(lat).sum(), lat.size, np.nanmin(lat), np.nanmax(lat)))
print("[LON] shape=%s finite=%d/%d range %.3f..%.3f" % (
    lon.shape, np.isfinite(lon).sum(), lon.size, np.nanmin(lon), np.nanmax(lon)))
if lat.shape != (ny, nx) or lon.shape != (ny, nx):
    problems.append("lat/lon shape %s/%s != grid (y,x)=( %d,%d)" % (lat.shape, lon.shape, ny, nx))
if not (np.all(np.isfinite(lat)) and np.all(np.isfinite(lon))):
    problems.append("non-finite lat/lon cells")
# 2D lat/lon should be consistent with 1D y/x for a regular grid
if np.all(np.isfinite(lat)):
    rowsp = np.abs(np.diff(lat, axis=0)).max()
    colsp = np.abs(np.diff(lat, axis=1)).max()
    print("  lat varies along y (max row-diff %.4f) vs along x (max col-diff %.6f)" % (rowsp, colsp))
if lon.min() < -115 or lon.max() > -85 or lat.min() < 25 or lat.max() > 50:
    problems.append("grid bbox outside expected Ogallala region")
else:
    print("  bbox inside Ogallala region: OK")

# ---- roi_mask ----
mask = ds.variables["roi_mask"][:]
uniq = np.unique(mask)
inside = float((mask == 1).sum()) / mask.size
print("\n[MASK] unique=%s inside=%.1f%% (%d/%d)" % (uniq, 100 * inside, int((mask == 1).sum()), mask.size))
if not set(uniq.tolist()) <= {0, 1}:
    problems.append("roi_mask values outside {0,1}: %s" % uniq)

# ---- crs ----
if "crs" not in ds.variables:
    problems.append("missing crs variable")
else:
    c = ds.variables["crs"]
    print("\n[CRS] attrs:", {a: str(getattr(c, a))[:90] for a in c.ncattrs()})
    if not hasattr(c, "crs_wkt") and not hasattr(c, "spatial_ref"):
        notes.append("crs has no WKT (crs_wkt/spatial_ref) — weak for GIS readers")

# ---- ET variable sanity (strided + coarse full scan) ----
ET = ds.variables["ET"]
fill = getattr(ET, "_FillValue", None)
print("\n[ET] dtype=%s fill=%s units=%s" % (ET.dtype, fill, getattr(ET, "units", "?")))
sy, sx = max(1, ny // 200), max(1, nx // 200)
idx = sorted(set(np.linspace(0, nt - 1, min(12, nt)).astype(int)))
fin_all, neg_all, samp_all = 0, 0, 0
gmin, gmax, gsum, gsum2, gn = np.inf, -np.inf, 0.0, 0.0, 0
t_empty_samp, t_const = 0, 0
for ti in idx:
    a = ET[ti, ::sy, ::sx].astype(np.float64)
    if fill is not None:
        a[a == fill] = np.nan
    a[~np.isfinite(a)] = np.nan
    fin = a[np.isfinite(a)]
    samp_all += a.size
    if fin.size == 0:
        t_empty_samp += 1
        continue
    fin_all += fin.size
    neg_all += int((fin < 0).sum())
    gmin = min(gmin, float(fin.min())); gmax = max(gmax, float(fin.max()))
    gsum += float(fin.sum()); gsum2 += float((fin ** 2).sum()); gn += fin.size
    if fin.size >= 20:
        vals, counts = np.unique(fin, return_counts=True)
        if counts.max() / fin.size > 0.995:
            t_const += 1
mean = gsum / gn if gn else float("nan")
std = (gsum2 / gn - mean ** 2) ** 0.5 if gn else float("nan")
print("  sample: finite %.1f%% | neg %.3f%% | min %.4g max %.4g mean %.4g std %.4g"
      % (100 * fin_all / max(1, samp_all), 100 * neg_all / max(1, fin_all),
         gmin, gmax, mean, std))
print("  sampled empty steps: %d/%d | const-plane steps: %d/%d" % (t_empty_samp, len(idx), t_const, len(idx)))
if neg_all:
    problems.append("negative ET values present (%.3f%% of sample; hourly ET should be >= 0)" % (100 * neg_all / max(1, fin_all)))
if t_const:
    problems.append("constant-plane artifact in %d sampled steps" % t_const)
# units inference helper
print("  units hint: max %.3g, mean %.3g -> %s" % (
    gmax, mean,
    "Kelvin-like (LST, NOT mm ET)" if gmax > 100 else
    "mm-like (plausible hourly ET)" if gmax < 30 else "AMBIGUOUS — confirm units"))
# coarse full-cube scan: empty steps + outside-mask leakage on 3 steps
my, mx = max(1, ny // 100), max(1, nx // 100)
mm = mask[::my, ::mx]
full_empty, leak_worst = 0, 0.0
for ti in range(nt):
    a = ET[ti, ::my, ::mx].astype(np.float64)
    if fill is not None:
        a[a == fill] = np.nan
    if np.isfinite(a).sum() == 0:
        full_empty += 1
print("  full-scan empty timesteps: %d/%d" % (full_empty, nt))
if full_empty == nt:
    problems.append("ET is entirely empty (all %d steps)" % nt)
elif full_empty:
    notes.append("ET: %d/%d timesteps empty" % (full_empty, nt))
for ti in (0, nt // 2, nt - 1):
    a = ET[ti, ::my, ::mx].astype(np.float64)
    if fill is not None:
        a[a == fill] = np.nan
    fin = np.isfinite(a)
    leak = float(((fin) & (mm == 0)).sum()) / max(1, fin.sum())
    leak_worst = max(leak_worst, leak)
print("  worst outside-mask finite (steps 0/mid/last): %.2f%%" % (100 * leak_worst))
if leak_worst > 0.02:
    problems.append("%.1f%% of finite pixels fall outside roi_mask" % (100 * leak_worst))

ds.close()
print("\n" + "=" * 70)
if problems:
    print("VERDICT: CORRUPT / SUSPECT — %d problem(s):" % len(problems))
    for p in problems:
        print("  - " + p)
else:
    print("VERDICT: CORRECT (no corruption detected in audit)")
if notes:
    print("Notes:")
    for n in notes:
        print("  * " + n)
print("=" * 70)
