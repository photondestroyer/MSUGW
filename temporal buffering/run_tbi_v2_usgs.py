"""
TBI v2 — Temporal Buffering Index with preprocessed USGS groundwater data.

Changes vs v1 (run_tbi_temporal_maps.py):
  * GBI now uses REAL hydraulic conductivity K (USGS OFR 98-548 E00 zones -> 4 km)
    with Specific Yield / storage coefficient  S = 0.15  (user-specified constant):
        GBI_raw = K_norm^{wT} * S^{-wS},  wT=wS=1  =>  GBI_raw = K_norm / 0.15
    then min-max normalized to [0,1] over the aquifer (VVIP regional comparison mode).
    ST (saturated thickness) remains a documented PLACEHOLDER -> T = K*ST not yet
    computable; ranking currently K-driven (see TBI_DOCUMENTATION.md section 3b).
  * SBI unchanged proxy (5-HSG)/4 -- true gSSURGO PAW pending polygon extraction
    (CONUS gdb lacks spatial index; see docs). DAC/ADP neutral placeholders.
  * New diagnostics: USGS water-level change rasters (predev->2019, 2017->2019)
    mapped on the same grid for sustainability screening (Research Plan sec 6.3).
  * Everything masked to hp_bound2010 aquifer boundary (spatially consistent).

Common grid : GRIDMET/DROUGHT native 287x181 (~4 km, EPSG:4326)
Common time : DROUGHT 3035 steps, 1984-01-05 .. 2026-07-29 (~5-day)

Outputs -> G:/MSU_GWB/datasets/tbi_maps_v2_usgs/
Sources read-only (merged_datasets/*.nc lazy dask; derived_usgs/*.nc small).
"""
import warnings
warnings.filterwarnings("ignore")

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(r"G:/MSU_GWB/datasets")
MERGED = BASE / "merged_datasets"
DERIVED = BASE / "derived_usgs"
OUT = BASE / "tbi_maps_v2_usgs"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"figure.dpi": 150, "font.size": 8})

S_SY = 0.15          # <-- USER-SPECIFIED Specific Yield / storage coefficient
ALPHA = 0.5          # soil weight   (alpha+beta=1)
BETA = 0.5           # groundwater weight
EPS = 1e-9

def log(m): print(f"[tbi-v2] {m}")

# ------------------------------------------------------------------
# 1. Lazy opens
# ------------------------------------------------------------------
log("opening datasets lazily...")
ds_drought = xr.open_dataset(MERGED / "DROUGHT_Merged_Ogallala.nc",
                             chunks={"time": 10, "y": 287, "x": 181},
                             mask_and_scale=False)
ds_hsg = xr.open_dataset(MERGED / "Hydrologic_Soil_Group_250m_2019_Ogallala.nc",
                         chunks={"y": 2048, "x": 2048}, mask_and_scale=False)
ds_k = xr.open_dataset(DERIVED / "K_hydraulic_conductivity_mday_4km.nc", chunks="auto")
ds_mask = xr.open_dataset(DERIVED / "aquifer_mask_4km.nc", chunks="auto")
ds_dwl_pre = xr.open_dataset(DERIVED / "dWL_predev_to_2019_ft_4km.nc", chunks="auto")
ds_dwl_1719 = xr.open_dataset(DERIVED / "dWL_2017_to_2019_ft_4km.nc", chunks="auto")

ty = ds_drought["y"].values; tx = ds_drought["x"].values
times = ds_drought["time"].values
H, W = len(ty), len(tx)
mask_da = ds_mask["aquifer_mask"]
log(f"grid {W}x{H}; aquifer px {int(mask_da.sum().compute())}; time steps {len(times)}")

extent = [float(tx[0]), float(tx[-1]), float(ty[-1]), float(ty[0])]

def plot_map(da, title, fname, cmap, vmin, vmax, cbar_label):
    fig, ax = plt.subplots(figsize=(6.5, 7))
    data = np.ma.masked_invalid(np.asarray(da.values, dtype="float64"))
    im = ax.imshow(data, extent=extent, origin="upper", cmap=cmap,
                   vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    ax.set_title(title, fontsize=9, pad=8)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02); cb.set_label(cbar_label, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=180); plt.close(fig)
    log(f"  saved {fname}")

# ------------------------------------------------------------------
# 2. Static components (small computes)
# ------------------------------------------------------------------
log("computing SBI (HSG proxy) on 4km...")
hsg_native = ds_hsg["b1"].isel(time=0)
hsg_4km = hsg_native.interp(y=ds_drought["y"], x=ds_drought["x"],
                            method="nearest", kwargs={"fill_value": np.nan})
SBI = ((5.0 - hsg_4km) / 4.0).clip(0, 1).where(mask_da == 1)
SBI = SBI.compute()
log(f"  SBI mean {float(np.nanmean(SBI.values)):.3f} valid {int(np.isfinite(SBI.values).sum())}")

log(f"computing GBI from real K with S={S_SY} (Eq7, ST placeholder)...")
K = ds_k["K_mday"].where(mask_da == 1).compute()
# normalize K 0-1 across aquifer (min-max), then GBI_raw = K_norm / S ; renormalize for VI
K_valid = K.values[np.isfinite(K.values)]
kmin, kmax = float(np.nanmin(K_valid)), float(np.nanmax(K_valid))
K_norm = ((K - kmin) / (kmax - kmin)).clip(0, 1)
GBI_raw = K_norm / S_SY                       # wT=wS=1
gmin, gmax = float(np.nanmin(GBI_raw.values)), float(np.nanmax(GBI_raw.values))
GBI = ((GBI_raw - gmin) / (gmax - gmin)).clip(0, 1)
GBI.attrs.update({
    "formula": f"GBI = norm( norm(K)^{1} / {S_SY} )  [Eq7 with ST placeholder]",
    "K_source_range_md": f"{kmin:.2f}-{kmax:.2f}",
    "note": ("ST missing -> T=K*ST replaced by K alone; ranking K-driven. "
             "S=0.15 user-specified scalar (rescales linearly, removed by renorm).")})
log(f"  GBI mean {float(np.nanmean(GBI.values)):.3f}")

DENOM = (ALPHA * SBI + BETA * GBI).clip(min=EPS)
log(f"  DENOM mean {float(np.nanmean(DENOM.values)):.3f}")

plot_map(SBI, f"SBI proxy (HSG->4km)  alpha={ALPHA}\n(5-HSG)/4, aquifer mask applied",
         "SBI_static.png", "YlGnBu", 0, 1, "SBI (0-1)")
plot_map(GBI, f"GBI from USGS K (OFR98-548) & S={S_SY}\nEq7 with ST placeholder",
         "GBI_static.png", "PuBuGn", 0, 1, "GBI (0-1)")
plot_map(DENOM, "Buffering denominator\n0.5*SBI + 0.5*GBI",
         "denom_static.png", "viridis", 0, 1, "aSBI+bGBI")

# Diagnostics: water-level change maps
dwl_pre = ds_dwl_pre["dWL_predev_to_2019"].compute()
dwl_1719 = ds_dwl_1719["dWL_2017_to_2019"].compute()
plot_map(dwl_pre, "USGS water-level change predev(~1950)->2019\n(regridded 500m EPSG:5070 -> 4km)",
         "diag_dWL_predev_to_2019.png", "RdBu", -100, 100, "feet")
plot_map(dwl_1719, "USGS water-level change 2017->2019",
         "diag_dWL_2017_to_2019.png", "RdBu", -10, 10, "feet")

# ------------------------------------------------------------------
# 3. DSI(t) domain-mean time series (lazy, same as v1)
# ------------------------------------------------------------------
log("building lazy DSI(t,y,x) from GRIDMET-DROUGHT indices...")
use_vars = ["spi90d", "spi180d", "spei90d", "spei180d"]
w_each = 1.0 / len(use_vars)
log_sum = None
for v in use_vars:
    inten = xr.where(ds_drought[v] <= -1.0, -ds_drought[v], EPS)
    term = w_each * np.log(inten.clip(min=EPS))
    log_sum = term if log_sum is None else log_sum + term
DSI_t = np.exp(log_sum)

log("computing domain-mean DSI(t) (aquifer-masked)...")
DSI_masked_t = DSI_t.where(mask_da == 1)
DSI_mean = DSI_masked_t.mean(dim=["y", "x"], skipna=True).compute()

VI_mean = (DSI_masked_t / DENOM).mean(dim=["y", "x"], skipna=True).compute()
TBI_mean = 1.0 / VI_mean.copy(deep=True)

fig, ax = plt.subplots(figsize=(10, 3.4))
ax.plot(times, DSI_mean.values, color="#b2182b", lw=0.7)
ax.set_ylabel("Domain-mean DSI(t)")
ax.set_title("DSI(t) aquifer-mean | GRIDMET/DROUGHT spi/spei 90/180d | 3035 steps ~5-day")
ax.grid(alpha=.3, lw=.4)
fig.tight_layout(); fig.savefig(OUT / "DSI_timeseries.png", dpi=180); plt.close(fig)

fig, ax = plt.subplots(figsize=(10, 3.4))
ax.plot(times, VI_mean.values, color="#2166ac", lw=0.7, label="VI(t) mean")
ax.axhline(1.0, color="k", lw=.6, ls="--")
ax2 = ax.twinx()
ax2.plot(times, TBI_mean.values, color="#35978f", lw=0.7, alpha=.8)
ax.set_ylabel("VI mean"); ax2.set_ylabel("TBI mean (1/VI)", color="#35978f")
ax.set_title("VI(t) & TBI(t)=1/VI(t) aquifer-mean | SBI=HSG proxy, GBI=f(K,S=0.15)")
ax.grid(alpha=.3, lw=.4)
fig.tight_layout(); fig.savefig(OUT / "VI_TBI_timeseries.png", dpi=180); plt.close(fig)
log(f"  VI mean range {float(VI_mean.min()):.3f}..{float(VI_mean.max()):.3f}")

# ------------------------------------------------------------------
# 4. Selected-time spatial maps (fast per-index path)
# ------------------------------------------------------------------
time_index = pd.to_datetime(times)
def nearest_idx(date_str):
    t = pd.to_datetime(date_str)
    return int(np.argmin(np.abs((time_index - t).values)))

candidates = ["1984-07-13", "2000-06-03", "2011-08-03", "2012-07-13",
              "2014-12-31", "2022-08-03"]
imax = int(np.argmax(DSI_mean.values)); imin = int(np.argmin(DSI_mean.values))
sel = sorted(set([nearest_idx(c) for c in candidates] + [imin, imax]))
log(f"selected time idx {sel}")

for i in sel:
    d = str(time_index[i])[:10]
    inten_stack = []
    for v in use_vars:
        da = ds_drought[v].isel(time=i)
        inten_stack.append(xr.where(da <= -1.0, -da, EPS))
    logs = sum(w_each * np.log(s.clip(min=EPS)) for s in inten_stack)
    DSI_i = np.exp(logs).where(mask_da == 1).compute()
    VI_i = (DSI_i / DENOM).compute()
    TBI_i = (1.0 / VI_i.clip(min=EPS)).clip(max=5).compute()
    plot_map(DSI_i, f"DSI(t) {d}", f"DSI_map_{i:04d}_{d}.png", "YlOrRd", 0, 2, "DSI")
    plot_map(VI_i, f"VI(t) {d}  >1 vulnerable\nVI = DSI/(0.5 SBI + 0.5 GBI)",
             f"VI_map_{i:04d}_{d}.png", "RdYlGn_r", 0, 5, "VI (0-5)")
    plot_map(TBI_i, f"TBI(t) {d} = 1/VI  higher=buffered",
             f"TBI_map_{i:04d}_{d}.png", "BrBG", 0, 3, "TBI (0-3)")

# ------------------------------------------------------------------
# 5. Multipanels
# ------------------------------------------------------------------
panel = sel[:6]
fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
for ax, i in zip(axes.flat, panel):
    d = str(time_index[i])[:10]
    inten = sum(w_each*np.log(xr.where(ds_drought[v].isel(time=i) <= -1.0, -ds_drought[v].isel(time=i), EPS).clip(min=EPS)) for v in use_vars)
    VI_i = (np.exp(inten) / DENOM).compute()
    im = ax.imshow(np.ma.masked_invalid(VI_i.values), extent=extent, origin="upper",
                   cmap="RdYlGn_r", vmin=0, vmax=5, interpolation="nearest")
    ax.set_title(d, fontsize=8); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("VI(t) temporal change | GRIDMET 4km common grid, aquifer-masked\n"
             "SBI=(5-HSG)/4 proxy | GBI=norm(K/S=0.15), OFR98-548 K", fontsize=10)
cb = fig.colorbar(im, ax=list(axes.flat), shrink=.9, pad=.02, orientation="horizontal", location="bottom")
cb.set_label("VI (>1 vulnerable)", fontsize=8)
fig.savefig(OUT / "VI_multipanel.png", dpi=180); plt.close(fig)

# ------------------------------------------------------------------
# 6. Sustainability screening: VI vs dWL decline (Research Plan 6.3)
# ------------------------------------------------------------------
log("screening: binned mean dWL(predev->2019) vs long-term mean VI...")
VI_ltm = (DSI_masked_t / DENOM).mean(dim="time", skipna=True).compute()
df = pd.DataFrame({"vi": VI_ltm.values.ravel(),
                   "dwl": dwl_pre.values.ravel(),
                   "gb": GBI.values.ravel(), "sb": SBI.values.ravel()})
df = df.dropna()
bins = pd.qcut(df["vi"], 5, labels=["VLow","Low","Mid","High","VHigh"])
prof = df.groupby(bins, observed=True).agg(vi=("vi","mean"), dwl=("dwl","mean"),
                                           gbi=("gb","mean"), sbi=("sb","mean"))
print(prof)
fig, ax = plt.subplots(figsize=(6, 3.6))
colors = ["#1a9850","#91cf60","#fee08b","#fc8d59","#d73027"]
ax.bar(range(5), prof["dwl"].values, color=colors)
ax.set_xticks(range(5)); ax.set_xticklabels(prof.index)
ax.axhline(0, color="k", lw=.6)
ax.set_xlabel("long-term VI class (1984-2026)")
ax.set_ylabel("mean dWL predev->2019 (ft)")
ax.set_title("Sustainability screen: more vulnerable classes should show deeper declines\n(if mined-buffering hypothesis holds; Research Plan 6.3)", fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "screen_VI_vs_dWL_decline.png", dpi=180); plt.close(fig)
prof.to_csv(OUT / "screen_VI_vs_dWL_profile.csv")

with open(OUT / "README_v2.txt", "w", encoding="utf-8") as f:
    f.write(f"""TBI v2 outputs (this folder)
Common grid: GRIDMET 287x181 ~4km EPSG:4326, aquifer-masked (hp_bound2010)
Common time: DROUGHT 3035 steps 1984-01-05..2026-07-29 (~5-day)
SBI: (5-HSG)/4 proxy, DAC/ADP neutral placeholders
GBI: norm( norm(K)^1 / S ), K = USGS OFR98-548 zone midpoints (m/d, {kmin:.1f}-{kmax:.1f}),
     Specific yield S = {S_SY} (user-specified constant; rescale only, removed by renormalization)
     ST placeholder: true Eq7 needs T=K*ST (McGuire saturated-thickness raster still missing)
VI(t)=DSI(t)/(0.5*SBI+0.5*GBI); TBI=1/VI; VI>1 vulnerable
Diagnostics: diag_dWL_predev_to_2019.png ({float(np.nanmin(dwl_pre.values)):.0f}..{float(np.nanmax(dwl_pre.values)):.0f} ft),
     screen_VI_vs_dWL_decline.png bins long-term VI vs observed decline
Selected map dates: {[str(pd.to_datetime(times[i]) )[:10] for i in sel]}
""")
log("DONE v2")

ds_drought.close(); ds_hsg.close(); ds_k.close(); ds_mask.close()
ds_dwl_pre.close(); ds_dwl_1719.close()
