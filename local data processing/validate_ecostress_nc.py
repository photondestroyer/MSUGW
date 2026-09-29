"""Corruption audit for ECO_L3T_JET_Ogallala_ET.nc (read-only, strided sampling).

Usage (local):  python validate_ecostress_nc.py ["path/to/file.nc"]
Usage (Colab):  edit DEFAULT_PATH below, then Run.
"""
import os
import sys
import numpy as np
import netCDF4 as nc4
from datetime import datetime

DEFAULT_PATH = r"D:\Downloads\ECO_L3T_JET_Ogallala_ET.nc"
PATH = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_PATH
VALID = {
    "ETdaily": (-1.0, 30.0),
    "PTJPLSMinst": (-50.0, 1200.0),
    "STICinst": (-50.0, 1200.0),
    "MOD16inst": (-50.0, 1200.0),
    "BESSinst": (-50.0, 1200.0),
    "ETinstUncertainty": (0.0, 1000.0),
}
MAX_TIMESTEPS_SAMPLE = 8
CELL_CAP = 250_000

print("=" * 70)
print("FILE:", PATH)
print("SIZE: %.2f MiB" % (os.path.getsize(PATH) / 1024 ** 2))
ds = nc4.Dataset(PATH, "r")
print("GLOBAL Conventions:", getattr(ds, "Conventions", "<missing>"))
for a in ("title", "source", "history", "references", "comment", "roi"):
    print("  %s: %s" % (a, str(getattr(ds, a, "<missing>"))[:110]))
print("DIMENSIONS:", {k: len(v) for k, v in ds.dimensions.items()})
print("VARIABLES:", list(ds.variables.keys()))

problems = []
notes = []

# --- time ---
t = ds.variables["time"]
tvals = t[:]
print("\n[TIME] n=%d units=%s calendar=%s" % (len(tvals), t.units, getattr(t, "calendar", "?")))
dts = nc4.num2date(tvals, t.units, getattr(t, "calendar", "proleptic_gregorian"))
dts = [datetime(d.year, d.month, d.day, d.hour, d.minute, d.second) for d in dts]
print("  range: %s -> %s" % (dts[0], dts[-1]))
if np.any(np.diff(tvals) <= 0):
    problems.append("time axis not strictly increasing")
    print("  FAIL: non-monotonic time")
else:
    print("  monotonic: OK")
if len(set(tvals)) != len(tvals):
    problems.append("duplicate timesteps")
    print("  FAIL: duplicates")
else:
    print("  duplicates: none")
gaps = np.diff(sorted(tvals))
if len(gaps):
    # gaps in days
    gap_days = gaps / 86400.0
    print("  median gap: %.2f days | max gap: %.1f days" % (float(np.median(gap_days)), float(gap_days.max())))
    if gap_days.max() > 90:
        notes.append("max overpass gap %.0f days (ISS irregularity/MSU outages expected)" % gap_days.max())

# --- coverage: timesteps per month + missing months + big gaps ---
print("\n[COVERAGE] timesteps per month (YYYY-MM: count; '--' = zero, i.e. missing month)")
from collections import Counter
ym_counts = Counter((d.year, d.month) for d in dts)
y0, m0 = dts[0].year, dts[0].month
y1, m1 = dts[-1].year, dts[-1].month
missing_months = []
yy, mm = y0, m0
while (yy, mm) <= (y1, m1):
    c = ym_counts.get((yy, mm), 0)
    if c == 0:
        missing_months.append("%04d-%02d" % (yy, mm))
    mm += 1
    if mm > 12:
        mm, yy = 1, yy + 1
# compact year rows
yy, mm = y0, m0
while (yy, mm) <= (y1, m1):
    row_year = yy
    row = []
    while mm <= 12 and (yy, mm) <= (y1, m1):
        c = ym_counts.get((yy, mm), 0)
        row.append("%02d:%s" % (mm, c if c else "--"))
        mm += 1
    print("  %d  %s" % (row_year, " ".join(row)))
    if mm > 12:
        mm, yy = 1, yy + 1
print("  months in span: %d | months with data: %d | missing months: %d"
      % (len(ym_counts) + len(missing_months), len(ym_counts), len(missing_months)))
if missing_months:
    print("  MISSING MONTHS (%d): %s" % (len(missing_months),
          ", ".join(missing_months[:40]) + (" ..." if len(missing_months) > 40 else "")))
    notes.append("%d of %d months have zero overpasses" % (len(missing_months), len(ym_counts) + len(missing_months)))
# big gaps (>14 d) with bounding overpass dates
BIG = 14.0
big_gaps = [(dts[i], dts[i + 1], float(gap_days[i]))
            for i in range(len(gap_days)) if gap_days[i] > BIG]
print("  gaps > %.0f d: %d" % (BIG, len(big_gaps)))
for a, b, g in big_gaps[:20]:
    print("    %.1f d missing: %s -> %s" % (g, a.strftime("%Y-%m-%d %H:%M"), b.strftime("%Y-%m-%d %H:%M")))
if len(big_gaps) > 20:
    print("    ... and %d more" % (len(big_gaps) - 20))

# --- coords ---
lat = ds.variables["lat"][:]
lon = ds.variables["lon"][:]
print("\n[COORDS] lat shape %s finite %d/%d range %.3f..%.3f | lon shape %s finite %d/%d range %.3f..%.3f" % (
    lat.shape, np.isfinite(lat).sum(), lat.size, np.nanmin(lat), np.nanmax(lat),
    lon.shape, np.isfinite(lon).sum(), lon.size, np.nanmin(lon), np.nanmax(lon)))
if not np.all(np.isfinite(lat)) or not np.all(np.isfinite(lon)):
    problems.append("non-finite lat/lon coordinates")
if not (np.all(np.diff(lat) < 0)):
    notes.append("lat not strictly descending (check grid orientation)")
if not (np.all(np.diff(lon) > 0)):
    problems.append("lon not strictly ascending")
# Ogallala bbox sanity: lon ~ -106..-96, lat ~ 31..44
if not (lon.min() > -115 and lon.max() < -85 and lat.min() > 25 and lat.max() < 50):
    problems.append("grid bbox outside expected Ogallala region")

# --- roi_mask ---
mask = ds.variables["roi_mask"][:]
uniq = np.unique(mask)
inside = float((mask == 1).sum()) / mask.size
print("\n[MASK] unique=%s inside=%.1f%% (%d/%d)" % (uniq, 100 * inside, int((mask == 1).sum()), mask.size))
if not set(uniq.tolist()) <= {0, 1}:
    problems.append("roi_mask has values outside {0,1}: %s" % uniq)
if not (0.05 < inside < 0.95):
    notes.append("inside fraction %.1f%% unusual (SMAP ROI was ~40%%)" % (100 * inside))

# --- crs ---
if "crs" not in ds.variables:
    problems.append("missing crs grid-mapping variable")
else:
    print("\n[CRS] grid_mapping_name=%s" % getattr(ds.variables["crs"], "grid_mapping_name", "?"))
for b in VALID:
    if b in ds.variables:
        gm = getattr(ds.variables[b], "grid_mapping", None)
        co = getattr(ds.variables[b], "coordinates", None)
        if gm != "crs":
            problems.append("%s grid_mapping=%r (expected 'crs')" % (b, gm))

# --- per-band strided audit ---
nt, ny, nx = (len(tvals), len(lat), len(lon))
print("\nGRID: nt=%d ny=%d nx=%d (%.2fM cells/band full cube)" % (nt, ny, nx, nt * ny * nx / 1e6))
idx = sorted(set(np.linspace(0, max(0, nt - 1), min(MAX_TIMESTEPS_SAMPLE, nt)).astype(int)))
sy = max(1, ny // int(np.sqrt(CELL_CAP * max(1, len(idx)) / max(1, nx)) + 1)) if nx else 1
sx = max(1, nx // int(np.sqrt(CELL_CAP * max(1, len(idx)) / max(1, ny)) + 1)) if ny else 1
# simpler: stride so idx*subsample <= ~250k per band total
sy = max(1, int(np.ceil(ny / 250)))
sx = max(1, int(np.ceil(nx / 250)))
print("SAMPLE: %d timesteps %s, spatial stride y=%d x=%d" % (len(idx), idx, sy, sx))

print("\n[BANDS]")
all_empty_bands = []
for band, (vmin, vmax) in VALID.items():
    if band not in ds.variables:
        print("  %-18s MISSING from file" % band)
        notes.append("band %s absent" % band)
        continue
    v = ds.variables[band]
    fill = getattr(v, "_FillValue", None)
    print("  %-18s dtype=%s fill=%s units=%s" % (band, v.dtype, fill, getattr(v, "units", "?")))
    finite_total, range_ok_total, samp_total = 0, 0, 0
    t_empty, t_const = 0, 0
    bmin, bmax, bsum, bsum2, bn = np.inf, -np.inf, 0.0, 0.0, 0
    for ti in idx:
        a = v[ti, ::sy, ::sx].astype(np.float64)
        a[~np.isfinite(a)] = np.nan
        fin = a[np.isfinite(a)]
        samp_total += a.size
        if fin.size == 0:
            t_empty += 1
            continue
        finite_total += fin.size
        inr = fin[(fin >= vmin) & (fin <= vmax)]
        range_ok_total += inr.size
        bmin = min(bmin, float(fin.min())); bmax = max(bmax, float(fin.max()))
        bsum += float(fin.sum()); bsum2 += float((fin ** 2).sum()); bn += fin.size
        if inr.size >= 20:
            vals, counts = np.unique(inr, return_counts=True)
            if counts.max() / inr.size > 0.995:
                t_const += 1
    # full-cube empty-timestep scan (cheap: per-timestep finite count via stride-4 scan)
    full_empty = 0
    for ti in range(nt):
        a = v[ti, ::max(1, ny // 60), ::max(1, nx // 60)]
        if np.isfinite(a).sum() == 0:
            full_empty += 1
    frac = finite_total / max(1, samp_total)
    rok = range_ok_total / max(1, finite_total)
    mean = bsum / bn if bn else float("nan")
    std = (bsum2 / bn - mean ** 2) ** 0.5 if bn else float("nan")
    flag = []
    if full_empty:
        flag.append("EMPTY_TIMESTEPS=%d/%d" % (full_empty, nt))
    if frac < 0.01:
        flag.append("NEAR-EMPTY sample finite=%.2f%%" % (100 * frac))
    if rok < 0.99:
        flag.append("OUT-OF-RANGE %.2f%%" % (100 * (1 - rok)))
    if t_const:
        flag.append("CONST-PLANE in %d sampled steps" % t_const)
    status = "OK" if not flag else "FAIL: " + "; ".join(flag)
    print("    sample finite %.1f%% | in-range %.2f%% | min %.3g max %.3g mean %.3g std %.3g | empty-steps %d/%d | %s"
          % (100 * frac, 100 * rok, bmin, bmax, mean, std, full_empty, nt, status))
    if full_empty == nt:
        problems.append("%s: all %d timesteps empty" % (band, nt))
        all_empty_bands.append(band)
    elif full_empty:
        notes.append("%s: %d/%d timesteps empty" % (band, full_empty, nt))
    if rok < 0.99:
        problems.append("%s: %.2f%% of finite sample outside [%.3g, %.3g]" % (band, 100 * (1 - rok), vmin, vmax))
    if t_const:
        problems.append("%s: constant-plane artifact in %d sampled steps" % (band, t_const))

# --- ROI masking check: finite pixels should sit inside mask (allow small edge tolerance) ---
print("\n[ROI-MASK CONSISTENCY] (one middle timestep, coarse scan)")
ti = nt // 2
my, mx = max(1, ny // 120), max(1, nx // 120)
mm = mask[::my, ::mx]
for band in VALID:
    if band not in ds.variables:
        continue
    a = ds.variables[band][ti, ::my, ::mx]
    fin = np.isfinite(a)
    outside = fin & (mm == 0)
    if fin.sum() and outside.sum() / fin.sum() > 0.02:
        problems.append("%s: %.1f%% of finite pixels fall OUTSIDE roi_mask" % (band, 100 * outside.sum() / fin.sum()))
        print("  %-18s FAIL: %.1f%% finite outside mask" % (band, 100 * outside.sum() / fin.sum()))
    else:
        print("  %-18s OK (outside-mask finite: %.2f%%)" % (band, 100 * outside.sum() / max(1, fin.sum())))

ds.close()
print("\n" + "=" * 70)
if problems:
    print("VERDICT: CORRUPT / SUSPECT — %d problem(s):" % len(problems))
    for p in problems:
        print("  - " + p)
else:
    print("VERDICT: CORRECT (no corruption detected in strided audit)")
if notes:
    print("Notes (not failures):")
    for n in notes:
        print("  * " + n)
print("=" * 70)
