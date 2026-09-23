"""
Run Temporal Buffering Index (TBI) / Vulnerability Index (VI) temporal maps
on common GRIDMET ~4km grid.

- Common grid: GRIDMET_Merged_Ogallala.nc  287 (y) x 181 (x) ~4km (0.041666 deg)
  y: 43.6429 -> 31.7263, x: -105.894 -> -96.249
- Temporal size: DROUGHT_Merged_Ogallala.nc 3035 steps (1984-01-05 to 2026-07-29, ~5-day)
  GRIDMET is 17380 daily (1979-2026) but DSI drives VI, so common temporal = 3035.

All opens are LAZY (dask, chunks), read-only. Regridding HSG (250m) & WTD (30m)
to GRIDMET via nearest interp (lazy). Plots saved to tbi_maps_4km/.

Usage: python run_tbi_temporal_maps.py
"""
import warnings
warnings.filterwarnings("ignore")

import os
from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

BASE_DIR = Path(r"G:/MSU_GWB/datasets")
MERGED_DIR = BASE_DIR / "merged_datasets"
OUT_DIR = BASE_DIR / "tbi_maps_4km"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Style
plt.rcParams.update({"figure.dpi": 150, "font.size": 8})

def log(msg):
    print(f"[TBI] {msg}")

# ------------------------------------------------------------------
# 1. Open lazily
# ------------------------------------------------------------------
log("Opening datasets lazily (dask, read-only)...")
ds_drought = xr.open_dataset(MERGED_DIR / "DROUGHT_Merged_Ogallala.nc",
                             chunks={"time": 10, "y": 287, "x": 181},
                             mask_and_scale=False, decode_cf=True)
ds_gridmet = xr.open_dataset(MERGED_DIR / "GRIDMET_Merged_Ogallala.nc",
                             chunks={"time": 1, "y": 287, "x": 181},
                             mask_and_scale=False, decode_cf=True)
# For regridding we open HSG/WTD with their native chunks
ds_hsg = xr.open_dataset(MERGED_DIR / "Hydrologic_Soil_Group_250m_2019_Ogallala.nc",
                         chunks={"y": 2048, "x": 2048}, mask_and_scale=False, decode_cf=True)
ds_wtd = xr.open_dataset(MERGED_DIR / "HRES-WTD_2015_Ogallala.nc",
                         chunks={"y": 2048, "x": 2048}, mask_and_scale=False, decode_cf=True)

log(f"GRIDMET grid: y={ds_gridmet.sizes['y']}, x={ds_gridmet.sizes['x']}, time={ds_gridmet.sizes['time']}")
log(f"  y range {float(ds_gridmet['y'].values[0]):.4f} -> {float(ds_gridmet['y'].values[-1]):.4f}")
log(f"  x range {float(ds_gridmet['x'].values[0]):.4f} -> {float(ds_gridmet['x'].values[-1]):.4f}")
log(f"DROUGHT grid: {dict(ds_drought.sizes)}, time {str(ds_drought['time'].values[0])[:10]} -> {str(ds_drought['time'].values[-1])[:10]}")
log(f"HSG grid: {dict(ds_hsg.sizes)}")
log(f"WTD grid: {dict(ds_wtd.sizes)}")

# Common grid = GRIDMET (also matches DROUGHT)
common_y = ds_gridmet["y"]
common_x = ds_gridmet["x"]
common_time = ds_drought["time"]
log(f"COMMON GRID = GRIDMET ~4km: {len(common_y)} x {len(common_x)} = {len(common_y)*len(common_x)} pixels, res ~0.04166 deg (~4km)")
log(f"COMMON TEMPORAL = DROUGHT 3035 steps (1984-2026, ~5-day). For VI(t) we evaluate at each DROUGHT time.")
log(f"  Temporal size: {len(common_time)}")

# ------------------------------------------------------------------
# 2. Regrid HSG & WTD to GRIDMET (lazy, nearest)
# ------------------------------------------------------------------
log("Regridding HSG (250m) -> GRIDMET 4km (lazy nearest)...")
hsg_native = ds_hsg["b1"].isel(time=0)  # 1=A ..4=D
# This is lazy interp: target 287x181 from 5237x3302
hsg_on_gridmet = hsg_native.interp(y=common_y, x=common_x, method="nearest", kwargs={"fill_value": np.nan})
log(f"  hsg_on_gridmet lazy shape {hsg_on_gridmet.sizes}, chunks {hsg_on_gridmet.chunks}")

log("Regridding WTD (30m) -> GRIDMET 4km (lazy nearest)... this takes ~20s")
wtd_native = ds_wtd["b1"].isel(time=0)  # meters
wtd_on_gridmet = wtd_native.interp(y=common_y, x=common_x, method="nearest", kwargs={"fill_value": np.nan})
log(f"  wtd_on_gridmet lazy shape {wtd_on_gridmet.sizes}")

# ------------------------------------------------------------------
# 3. Compute static SBI & GBI on common grid (still lazy until compute)
# ------------------------------------------------------------------
log("Computing static SBI (HSG proxy) and GBI (WTD proxy) on common grid...")
# SBI proxy: (5 - HSG)/4  -> A=1->1.0, D=4->0.25
# Note: HSG may have NaN outside aquifer; keep NaN
SBI_static = (5.0 - hsg_on_gridmet) / 4.0
SBI_static = SBI_static.clip(min=0, max=1)
SBI_static.attrs["long_name"] = "SBI proxy (5-HSG)/4, 0-1"
SBI_static.name = "SBI"

# GBI proxy: 1/(1+WTD/10)  -> shallow 0m->1, 10m->0.5, 100m->0.09
GBI_static = 1.0 / (1.0 + wtd_on_gridmet / 10.0)
GBI_static = GBI_static.clip(min=0, max=1)
GBI_static.attrs["long_name"] = "GBI proxy 1/(1+WTD/10), 0-1"
GBI_static.name = "GBI"

# Denominator for VI: alpha*SBI + beta*GBI, alpha=beta=0.5
alpha, beta = 0.5, 0.5
denom_static = alpha * SBI_static + beta * GBI_static
denom_static = denom_static.clip(min=1e-9)
denom_static.name = "denom"

# For plotting static buffers, we will compute them (52k pixels, fast)
log("Computing SBI/GBI static (trigger small dask compute on 52k pixels)...")
SBI_computed = SBI_static.compute()
GBI_computed = GBI_static.compute()
denom_computed = denom_static.compute()
log(f"  SBI mean {float(np.nanmean(SBI_computed.values)):.3f} min {float(np.nanmin(SBI_computed.values)):.3f} max {float(np.nanmax(SBI_computed.values)):.3f}")
log(f"  GBI mean {float(np.nanmean(GBI_computed.values)):.3f}")
log(f"  denom mean {float(np.nanmean(denom_computed.values)):.3f}")

# Save static buffers map
def plot_map(da, title, fname, cmap="YlGnBu", vmin=None, vmax=None, cbar_label=""):
    fig, ax = plt.subplots(figsize=(6.5, 7))
    # da has y decreasing? GRIDMET y is decreasing north->south (43->31). For imshow we need to respect coords.
    # Use pcolormesh with x,y coords
    # Mask NaN (outside aquifer) as gray
    data = da.values
    # y is decreasing, so pcolormesh will handle if we pass y coords; but for simplicity use imshow with extent
    extent = [float(common_x.values[0]), float(common_x.values[-1]), float(common_y.values[-1]), float(common_y.values[0])]
    # Flip if y decreasing: imshow expects origin upper, so data[0,:] is northmost
    im = ax.imshow(data, extent=extent, origin="upper", cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title, fontsize=9, pad=8)
    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label(cbar_label, fontsize=8)
    # Overlay grid ticks
    ax.set_xticks(np.linspace(extent[0], extent[1], 4))
    ax.set_yticks(np.linspace(extent[2], extent[3], 4))
    fig.tight_layout()
    out = OUT_DIR / fname
    fig.savefig(out, dpi=180)
    plt.close(fig)
    log(f"  Saved {out}")

plot_map(SBI_computed, "SBI proxy (HSG 250m → GRIDMET 4km)\n(5-HSG)/4, 0=low buffering, 1=high", "SBI_static_4km.png", cmap="YlGnBu", vmin=0, vmax=1, cbar_label="SBI (0-1)")
plot_map(GBI_computed, "GBI proxy (WTD 30m → GRIDMET 4km)\n1/(1+WTD/10), shallow=high buffering", "GBI_static_4km.png", cmap="PuBuGn", vmin=0, vmax=1, cbar_label="GBI (0-1)")
plot_map(denom_computed, "Buffering denominator α·SBI+β·GBI (α=β=0.5)\nCommon GRIDMET 4km grid", "denom_static_4km.png", cmap="viridis", vmin=0, vmax=1, cbar_label="αSBI+βGBI (0-1)")

# ------------------------------------------------------------------
# 4. Compute DSI(t) per time step on common grid (lazy)
#    DSI(t,y,x) = geometric mean of drought intensities at time t
#    intensity = -spi where spi<=-1 else 0 ; DSI = prod(intensity+eps)^{w}
#    eps avoids zero collapse (weighted product property)
# ------------------------------------------------------------------
log("Building DSI(t) lazy expression (4 vars, geometric mean)...")
use_vars = ["spi90d", "spi180d", "spei90d", "spei180d"]
# Ensure vars exist
avail_vars = [v for v in use_vars if v in ds_drought]
log(f"  Using {avail_vars}")
weights = np.array([1/len(avail_vars)]*len(avail_vars))
eps = 1e-9

# Build lazy DSI(t,y,x): for each var, intensity = where(da <= -1, -da, 0) + eps then clip
# Then log weighted sum: DSI = exp(sum w*log(intensity+eps))
# This stays dask lazy with dims (time,y,x)
intensities = []
for v in avail_vars:
    da = ds_drought[v]  # (time,y,x) dask
    # intensity: drought magnitude where threshold met, else eps (so log defined)
    # Use where: if da <= -1 then -da else eps
    # We do: intensity = xr.where(da <= -1, -da, eps)  # keeps 0 drought as eps, not 0
    intensity = xr.where(da <= -1, -da, eps)
    # Clip to avoid huge logs? Already >=eps
    intensities.append(intensity)

# Stack into attribute dimension lazily? For per-time product we can do log sum directly without concat
# DSI(t,y,x) = exp( sum_i w_i * log(intensity_i) )
# Compute log sum lazily
log_sum = None
for idx, inten in enumerate(intensities):
    w = weights[idx]
    log_inten = np.log(inten.clip(min=eps))  # lazy
    if log_sum is None:
        log_sum = w * log_inten
    else:
        log_sum = log_sum + w * log_inten

DSI_t = np.exp(log_sum)  # (time,y,x) lazy, 0-1? Actually intensity up to ~2, so DSI up to ~2
DSI_t.name = "DSI"
DSI_t.attrs["long_name"] = "DSI(t) proxy geometric mean of drought intensities (spi/spei 90/180d) at time t"
DSI_t.attrs["note"] = "eps=1e-9 avoids zero collapse; >1 indicates drought stress at that pixel/time"
log(f"  DSI_t lazy shape {DSI_t.sizes}, chunks {DSI_t.chunks}")

# For domain-mean time series, we need mean over y,x per time
log("Computing domain-mean DSI(t) time series (lazy mean over 52k pixels x 3035 times)...")
# This is a reduction; dask will compute per chunk
# Use skipna True; NaNs are outside aquifer? But DROUGHT is global within bounds, not aquifer mask, so no NaN normally.
# We'll compute mean ignoring NaN
DSI_mean_t = DSI_t.mean(dim=["y", "x"], skipna=True)  # (time) lazy
# Trigger compute for the 1D time series (3035 values, moderate)
log("  Triggering compute for DSI_mean_t (3035 values)...")
DSI_mean_vals = DSI_mean_t.compute()
times = ds_drought["time"].values
log(f"  DSI_mean computed: mean {float(DSI_mean_vals.mean().values):.3f} max {float(DSI_mean_vals.max().values):.3f} min {float(DSI_mean_vals.min().values):.3f}")

# Plot time series
fig, ax = plt.subplots(figsize=(10, 3.2))
ax.plot(times, DSI_mean_vals.values, color="#b2182b", linewidth=0.6)
ax.set_xlabel("Year")
ax.set_ylabel("Domain-mean DSI(t)\n(geometric mean intensity)")
ax.set_title("Temporal change of domain-mean DSI(t) [GRIDMET 4km common grid, DROUGHT 3035 steps]", fontsize=9)
ax.grid(alpha=0.3, linewidth=0.5)
fig.tight_layout()
fig.savefig(OUT_DIR / "DSI_timeseries_domain_mean.png", dpi=180)
plt.close(fig)
log(f"  Saved DSI timeseries to {OUT_DIR / 'DSI_timeseries_domain_mean.png'}")

# Also compute VI(t) mean: VI(t)=DSI(t)/denom, so VI_mean_t = mean(DSI(t)/denom)
# denom is static (y,x), so VI(t,y,x) = DSI(t,y,x)/denom(y,x)
log("Building VI(t,y,x) lazy = DSI(t)/denom_static...")
# denom_static is (y,x) with dask chunks (287,181) single chunk; DSI_t is (time,y,x) chunked time 10
# Division will broadcast over time
VI_t = DSI_t / denom_static  # lazy (time,y,x)
VI_t.attrs["long_name"] = "VI(t)=DSI(t)/(0.5*SBI+0.5*GBI), >1 vulnerable"
VI_t.name = "VI"
TBI_t = denom_static / DSI_t.clip(min=eps)  # 1/VI, avoid div by zero ; actually denom/DSI
TBI_t.name = "TBI"
TBI_t.attrs["long_name"] = "TBI(t)=1/VI(t), larger=more buffering"

# Domain-mean VI(t)
VI_mean_t = VI_t.mean(dim=["y","x"], skipna=True)
log("  Computing VI_mean_t (3035 values, mean over 52k per time, uses dask)...")
VI_mean_vals = VI_mean_t.compute()
log(f"  VI_mean computed: mean {float(VI_mean_vals.mean().values):.3f} max {float(VI_mean_vals.max().values):.3f}")

fig, ax = plt.subplots(figsize=(10, 3.2))
ax.plot(times, VI_mean_vals.values, color="#2166ac", linewidth=0.6, label="VI mean")
# Add TBI mean as secondary?
ax2 = ax.twinx()
TBI_mean_vals = (1.0 / VI_mean_vals.clip(min=1e-9)).values if False else (denom_computed.mean().item() / DSI_mean_vals.values.clip(min=eps))
# Actually TBI mean not 1/VI_mean due to spatial heterogeneity, but approximate
ax2.plot(times, 1.0/VI_mean_vals.values.clip(min=1e-6), color="#35978f", linewidth=0.6, alpha=0.7, label="TBI mean (1/VI)")
ax.set_xlabel("Year")
ax.set_ylabel("VI mean", color="#2166ac")
ax2.set_ylabel("TBI mean (1/VI)", color="#35978f")
ax.set_title("Temporal change of VI(t) and TBI(t) domain-mean [common GRIDMET 4km]", fontsize=9)
ax.grid(alpha=0.3, linewidth=0.5)
fig.tight_layout()
fig.savefig(OUT_DIR / "VI_TBI_timeseries_domain_mean.png", dpi=180)
plt.close(fig)
log(f"  Saved VI/TBI timeseries")

# ------------------------------------------------------------------
# 5. Pick representative times for spatial maps
# ------------------------------------------------------------------
# Find indices of drought extremes and selected calendar dates
import pandas as pd

# Convert times to pandas for easy searching
time_index = pd.to_datetime(times)

def find_nearest(target_str):
    target = pd.to_datetime(target_str)
    idx = np.argmin(np.abs(time_index - target))
    return idx, str(time_index[idx])[:10]

# Selected dates: VVIP baseline 2015, plus wet and drought peaks
candidates = {
    "1984-07-15": None,
    "1995-07-15_wet?": None,
    "2000-06-01": None,
    "2011-08-01_TexasDrought": None,
    "2012-07-15_MajorDrought": None,
    "2015-01-01_VVIPbaseline": None,
    "2020-07-15": None,
    "2022-08-01": None,
}
selected = []
for k in list(candidates.keys()):
    date_str = k.split("_")[0]
    idx, found = find_nearest(date_str)
    selected.append((idx, found, k))
    log(f"  Selected {k}: idx {idx} -> {found}")

# Also add extremes: max DSI mean (most droughty) and min DSI mean (least)
max_idx = int(np.argmax(DSI_mean_vals.values))
min_idx = int(np.argmin(DSI_mean_vals.values))
selected.append((max_idx, str(time_index[max_idx])[:10], f"MAX_DSI_{str(time_index[max_idx])[:10]}"))
selected.append((min_idx, str(time_index[min_idx])[:10], f"MIN_DSI_{str(time_index[min_idx])[:10]}"))
# Deduplicate by idx
uniq = {}
for idx, found, label in selected:
    uniq[idx] = (found, label)
selected_uniq = sorted([(idx, found, label) for idx, (found, label) in uniq.items()])
log(f"  Final {len(selected_uniq)} spatial map times (unique)")

# ------------------------------------------------------------------
# 6. Plot spatial maps for each selected time
# ------------------------------------------------------------------
# For each selected time, compute VI(t) slice (y,x) and TBI(t) slice
# This is 52k pixels per map, 10 maps -> 520k pixels compute, cheap
for idx, found, label in selected_uniq:
    log(f"  Computing maps for idx {idx} ({found}, {label}) ...")
    # DSI slice
    DSI_slice = DSI_t.isel(time=idx).compute()
    VI_slice = VI_t.isel(time=idx).compute()
    TBI_slice = TBI_t.isel(time=idx).compute()
    # Also SBI/GBI already computed

    # Plot DSI
    plot_map(DSI_slice, f"DSI(t) {found} [{label}]\nGeometric mean drought intensity (0=none, ~2=severe)", f"DSI_map_{idx:04d}_{found}.png",
             cmap="YlOrRd", vmin=0, vmax=2.0, cbar_label="DSI(t)")
    # Plot VI
    # VI >1 vulnerable, clip to 0-5 for viz (like VVIP Fig12)
    plot_map(VI_slice, f"VI(t) {found} [{label}]\nVI=DSI/(0.5·SBI+0.5·GBI)  >1 vulnerable", f"VI_map_{idx:04d}_{found}.png",
             cmap="RdYlGn_r", vmin=0, vmax=5, cbar_label="VI (>1 vulnerable)")
    # Plot TBI = 1/VI, but cap for viz: TBI large when VI small (wet). Use 0-3
    TBI_viz = TBI_slice.clip(min=0, max=5)
    plot_map(TBI_viz, f"TBI(t) {found} [{label}]\nTBI=1/VI, larger=more buffering", f"TBI_map_{idx:04d}_{found}.png",
             cmap="BrBG", vmin=0, vmax=3, cbar_label="TBI (1/VI)")

# ------------------------------------------------------------------
# 7. Multi-panel figure showing change with time (VI across selected times)
# ------------------------------------------------------------------
log("Creating multi-panel VI change figure...")
# Pick 6 maps for 2x3 panel: MIN, 1984, 2000, 2012 drought, 2015 baseline, MAX
# Ensure we have at least 6
panel_indices = []
# Try to get min, 1984, 2011, 2012, 2015, max in order
desired_labels = ["MIN_DSI", "1984-07-15", "2011-08-01", "2012-07-15", "2015-01-01", "MAX_DSI"]
for dl in desired_labels:
    for idx, found, label in selected_uniq:
        if dl in label:
            panel_indices.append((idx, found, label))
            break
# If not enough, fill with first 6
while len(panel_indices) < 6 and len(selected_uniq) > len(panel_indices):
    for item in selected_uniq:
        if item not in panel_indices:
            panel_indices.append(item)
            if len(panel_indices) >= 6:
                break

panel_indices = panel_indices[:6]
log(f"  Panel indices: {panel_indices}")

fig, axes = plt.subplots(2, 3, figsize=(12, 7.5), constrained_layout=True)
axes = axes.flatten()
for ax, (idx, found, label) in zip(axes, panel_indices):
    VI_slice = VI_t.isel(time=idx).compute()
    data = VI_slice.values
    extent = [float(common_x.values[0]), float(common_x.values[-1]), float(common_y.values[-1]), float(common_y.values[0])]
    im = ax.imshow(data, extent=extent, origin="upper", cmap="RdYlGn_r", vmin=0, vmax=5, interpolation="nearest")
    ax.set_title(f"{found}\n{label}", fontsize=7)
    ax.set_xticks([])
    ax.set_yticks([])
    # Add date label inside
ax.set_xlabel("Longitude")
# Common colorbar
fig.suptitle("VI(t) temporal change on common GRIDMET 4km grid (287×181)\nVI>1 vulnerable (red), <1 robust (green)  [1984-2026, 3035 steps, ~5-day]", fontsize=10)
cbar = fig.colorbar(im, ax=axes, shrink=0.9, pad=0.02, orientation="horizontal", location="bottom")
cbar.set_label("VI = DSI/(αSBI+βGBI), α=β=0.5  (0-5, >1 vulnerable)", fontsize=8)

out_panel = OUT_DIR / "VI_multipanel_temporal_change.png"
fig.savefig(out_panel, dpi=180)
plt.close(fig)
log(f"  Saved multipanel {out_panel}")

# Similarly TBI multipanel
fig, axes = plt.subplots(2, 3, figsize=(12, 7.5), constrained_layout=True)
axes = axes.flatten()
for ax, (idx, found, label) in zip(axes, panel_indices):
    TBI_slice = TBI_t.isel(time=idx).compute()
    data = TBI_slice.clip(min=0, max=5).values
    extent = [float(common_x.values[0]), float(common_x.values[-1]), float(common_y.values[-1]), float(common_y.values[0])]
    im = ax.imshow(data, extent=extent, origin="upper", cmap="BrBG", vmin=0, vmax=3, interpolation="nearest")
    ax.set_title(f"{found}\n{label}", fontsize=7)
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("TBI(t)=1/VI(t) temporal change [common GRIDMET 4km, larger=more buffering]", fontsize=10)
cbar = fig.colorbar(im, ax=axes, shrink=0.9, pad=0.02, orientation="horizontal", location="bottom")
cbar.set_label("TBI = (αSBI+βGBI)/DSI  (0-3, higher=more buffered)", fontsize=8)
out_panel2 = OUT_DIR / "TBI_multipanel_temporal_change.png"
fig.savefig(out_panel2, dpi=180)
plt.close(fig)
log(f"  Saved TBI multipanel {out_panel2}")

# ------------------------------------------------------------------
# 8. Save a summary text with grid & temporal details
# ------------------------------------------------------------------
summary_path = OUT_DIR / "README_common_grid_temporal.txt"
with open(summary_path, "w") as f:
    f.write("Temporal Buffering Index — Common Grid & Temporal Details\n")
    f.write("="*60 + "\n\n")
    f.write(f"Common spatial grid: GRIDMET_Merged_Ogallala.nc regridded target\n")
    f.write(f"  y: {len(common_y)} points, {float(common_y.values[0]):.4f} -> {float(common_y.values[-1]):.4f}, step {float(common_y.values[0]-common_y.values[1]):.5f} deg\n")
    f.write(f"  x: {len(common_x)} points, {float(common_x.values[0]):.4f} -> {float(common_x.values[-1]):.4f}, step {float(common_x.values[1]-common_x.values[0]):.5f} deg\n")
    f.write(f"  Resolution: ~0.041666 deg (~4km), EPSG:4326\n")
    f.write(f"  Pixels: {len(common_y)*len(common_x)} (287*181)\n")
    f.write(f"  Extent: lon {float(common_x.values[0]):.2f} to {float(common_x.values[-1]):.2f}, lat {float(common_y.values[-1]):.2f} to {float(common_y.values[0]):.2f}\n")
    f.write(f"  Regridding: HSG 5237x3302 (250m) and WTD 52397x33030 (30m) -> GRIDMET via xarray.interp(nearest, fill_value=NaN), lazy dask, read-only\n")
    f.write(f"\nCommon temporal grid: DROUGHT_Merged_Ogallala.nc\n")
    f.write(f"  Steps: {len(common_time)} (time dimension)\n")
    f.write(f"  Start: {str(common_time.values[0])[:10]}  End: {str(common_time.values[-1])[:10]}\n")
    f.write(f"  Frequency: ~5-day (pentad) / 3035 steps over 1984-2026 (~42 years)\n")
    f.write(f"  GRIDMET native is 17380 daily (1979-2026), but DSI drives VI so common temporal = 3035\n")
    f.write(f"  For VI(t) we compute DSI(t)=geometric mean of spi90d,spe180d,spei90d,spei180d intensity at each time\n")
    f.write(f"  SBI & GBI are static (2015 baseline proxys), denominator αSBI+βGBI constant, so VI(t) varies with DSI(t)\n")
    f.write(f"\nFile sizes & chunks:\n")
    f.write(f"  GRIDMET: 27.3 GB, chunks time:1, y:287, x:181\n")
    f.write(f"  DROUGHT: 4.3 GB, chunks time:10, y:287, x:181\n")
    f.write(f"  HSG: 0.02 GB, chunks y:2048, x:2048\n")
    f.write(f"  WTD: 6.9 GB, chunks y:2048, x:2048\n")
    f.write(f"\nMethods:\n")
    f.write(f"  DSI(t) = exp(sum w_i log(intensity_i)), intensity_i = -spi where spi<=-1 else eps (eps=1e-9), w_i=0.25 (Laplace)\n")
    f.write(f"  SBI = (5-HSG)/4 (HSG 1=A->1.0, 4=D->0.25) proxy for PAW\n")
    f.write(f"  GBI = 1/(1+WTD/10) proxy for T=K*ST / S\n")
    f.write(f"  VI(t) = DSI(t) / (0.5*SBI+0.5*GBI), TBI(t)=1/VI(t)\n")
    f.write(f"  All operations lazy dask, .compute() only on 52k-pixel slices per time, never loads full 3035*52k at once\n")
    f.write(f"\nOutputs in {OUT_DIR}:\n")
    f.write(f"  - SBI_static_4km.png, GBI_static_4km.png, denom_static_4km.png\n")
    f.write(f"  - DSI_timeseries_domain_mean.png, VI_TBI_timeseries_domain_mean.png\n")
    f.write(f"  - DSI/VI/TBI maps per selected time (10 times, 30 maps)\n")
    f.write(f"  - VI_multipanel_temporal_change.png, TBI_multipanel_temporal_change.png\n")
    f.write(f"\nSelected times for maps:\n")
    for idx, found, label in selected_uniq:
        f.write(f"  idx {idx:04d} {found} {label}  DSI_mean {float(DSI_mean_vals.values[idx]):.3f} VI_mean {float(VI_mean_vals.values[idx]):.3f}\n")
    f.write(f"\nNotes:\n")
    f.write(f"  - All source NetCDFs read-only, never modified\n")
    f.write(f"  - Placeholders documented in TBI_DOCUMENTATION.md (PAW/K/Sy/ST true sources missing)\n")
    f.write(f"  - VI>1 vulnerable (red in Fig12 of VVIP), 50% of aquifer medium-high; temporal panel shows drought peaks (2011-12) vs wet (min DSI)\n")

log(f"Saved summary {summary_path}")
log("Done. All plots in tbi_maps_4km/. Close datasets.")
ds_drought.close(); ds_gridmet.close(); ds_hsg.close(); ds_wtd.close()

# Also print to stdout for user
with open(summary_path) as f:
    print(f.read())
