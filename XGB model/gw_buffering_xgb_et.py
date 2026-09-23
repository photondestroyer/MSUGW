"""
Quantify groundwater buffering of dry-season evaporation (ET) — pilot
implementation of Research Plan Sec. 6 (Phase 0), ET leg only.

Core equation (Plan Sec. 6.1):
    Bi(t) = Fi(t) - F_hat_i(P, SMrz, VPD, Rn, T, PFT, season)
F_hat is an XGBoost regressor trained on merged_datasets; the residual during
detected dry-downs is the candidate buffering signal. Attribution to
groundwater follows only via WTD covariation + well/storage + irrigation
controls (Plan pitfalls table). SIF and VOD legs are NOT implemented here
(no SIF/VOD sensor in merged_datasets) - see assumptions A2, A3.

Datasets used (all read-only, lazy; outputs are NEW files only):
  MODIS_ET_SSEBop_Merged_Ogallala.nc  ET response (dekadal, ~1 km)
  GRIDMET_Merged_Ogallala.nc          P, VPD, Rn(srad), T (daily, 4 km)
  SPL4SMGP_Ogallala_FULL.nc           SMrz root-zone soil moisture (3-hourly)
  DROUGHT_Merged_Ogallala.nc          spi90d/spei90d dry-down detection (4 km)
  MCD12Q1_Merged_Ogallala.nc          PFT = LC_Type1 (annual, 500 m)
  HRES-WTD_2015_Ogallala.nc           static depth-to-water axis (30 m)
  Global_irrigation_Area_Merged...nc  irrigation history (annual, coarse)
  GFSAD1000_V1_2019_Ogallala.nc       cropland mask (2019, 1 km)
  GRACE_Merged_Ogallala.nc            basin-scale storage constraint (monthly)
  merged_ameriflux.nc + AmeriFlux_NEON_sites.xlsx  tower ET validation
  derived_usgs/: aquifer_mask_4km.nc, PAW_rootznaws_4km.nc,
                 dWL_predev_to_2019_ft_4km.nc, dWL_2017_to_2019_ft_4km.nc
                 (PAW = soil covariate; dWL = sustainability screen)

Outputs -> gw_buffering_xgb/ :
  metrics.json, effects.csv, xgb_et_model.json, Bi_drydown_mean_4km.nc,
  *.png (12 figures), run_log.txt

Memory discipline: sources opened lazily (dask); panel built in dekad batches;
full-array .compute() only on <= ~300 MB intermediates.
Run:  python gw_buffering_xgb_et.py
Assumptions: documented inline as A-codes; full register in
TBI_DOCUMENTATION.md Sec. 11.15.
"""

import json
import os
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import xgboost as xgb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------
# 0. Config
# --------------------------------------------------------------------------
BASE = Path(r"G:\MSU_GWB\datasets")
MERGED = BASE / "merged_datasets"
DERIVED = BASE / "derived_usgs"
OUT = BASE / "gw_buffering_xgb"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 42
PIX_STRIDE = 2            # sample every 2nd y/x pixel in aquifer (A8)
BATCH_DEKADS = 32         # dekads per panel-building batch (memory)
FULL_YEARS = list(range(2016, 2022))   # full calendar years in SMAP+SSEBop overlap (A5)
GROW_MONTHS = [5, 6, 7, 8, 9]          # primary dry-season window (A9)
DRYDOWN_SPI = -0.8                     # spi90d threshold for dry-down (A10)
SY = 0.15                              # specific yield, user constant (Sec. 11.3)
XGB_PARAMS = dict(n_estimators=800, max_depth=6, learning_rate=0.05,
                  subsample=0.8, colsample_bytree=0.7, reg_lambda=1.0,
                  min_child_weight=10, tree_method="hist",
                  random_state=SEED, n_jobs=min(8, os.cpu_count() or 4),
                  callbacks=[xgb.callback.EarlyStopping(rounds=50,
                                                        save_best=True)])

LOG_LINES = []


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    LOG_LINES.append(line)


rng = np.random.default_rng(SEED)

# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def grid_extent(da):
    x = np.asarray(da["x"].values, dtype=float)
    y = np.asarray(da["y"].values, dtype=float)
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    if dy < 0:
        return [x[0] - dx / 2, x[-1] + dx / 2, y[-1] + dy / 2, y[0] - dy / 2], "upper"
    return [x[0] - dx / 2, x[-1] + dx / 2, y[0] - dy / 2, y[-1] + dy / 2], "lower"


def decorate(ax, title, bnd, xlabel="Longitude (deg E)", ylabel="Latitude (deg N)"):
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    if bnd is not None:
        bnd.boundary.plot(ax=ax, edgecolor="black", linewidth=0.8)
    ax.grid(True, which="major", color="gray", alpha=0.35,
            linewidth=0.5, linestyle=":")
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=7)


def plot_grid(da, title, cmap="viridis", vmin=None, vmax=None, norm=None,
              cbar_label="", bnd=None, fname=None):
    extent, origin = grid_extent(da)
    fig, ax = plt.subplots(figsize=(7.2, 6.8))
    img = np.ma.masked_invalid(np.asarray(da.values, dtype="float64"))
    im = ax.imshow(img, extent=extent, origin=origin, cmap=cmap,
                   vmin=vmin, vmax=vmax, norm=norm, interpolation="nearest")
    decorate(ax, title, bnd)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label(cbar_label, fontsize=8)
    fig.tight_layout()
    if fname:
        fig.savefig(OUT / fname, dpi=170)
    plt.close(fig)


def rmse(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    return float(np.sqrt(np.mean((a[m] - b[m]) ** 2))) if m.sum() else float("nan")


def mean_ci(a):
    """mean +/- 1.96*SE (normal approx; documented A18)."""
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan"), float("nan"), 0
    se = a.std(ddof=1) / np.sqrt(a.size) if a.size > 1 else 0.0
    return float(a.mean()), float(1.96 * se), int(a.size)


# --------------------------------------------------------------------------
# 1. Statics -> 4 km numpy (aquifer mask, PFT, WTD, irrigation, PAW, dWL)
# --------------------------------------------------------------------------
log("== stage 1: statics ==")
import geopandas as gpd

bnd = gpd.read_file(BASE / "high_plains_quifer" / "hp_bound2010.shp").to_crs("EPSG:4326")

ds_ref = xr.open_dataset(MERGED / "GRIDMET_Merged_Ogallala.nc", chunks={"time": 1})
GY, GX = ds_ref["y"].values, ds_ref["x"].values
H, W = len(GY), len(GX)
log(f"harmonized grid = GRIDMET 4km: {W}x{H} (A4)")

mask = xr.open_dataset(DERIVED / "aquifer_mask_4km.nc", chunks="auto")["aquifer_mask"]
mask_np = (mask.values == 1)
log(f"aquifer pixels: {mask_np.sum()}")

# PFT: LC_Type1 2019 nearest (A6)
mcd = xr.open_dataset(MERGED / "MCD12Q1_Merged_Ogallala.nc", chunks={"time": 1})
pft4 = mcd["LC_Type1"].sel(time="2019-01-01", method="nearest")
pft4 = pft4.interp(y=ds_ref["y"], x=ds_ref["x"], method="nearest",
                   kwargs={"fill_value": np.nan}).compute().values
mcd.close()
log(f"PFT classes present: {np.unique(pft4[np.isfinite(pft4)]).astype(int).tolist()}")

# WTD static: stripe block-average 30 m -> 4 km (exact area mean, bounded RAM).
# Source y may be ascending or descending; handled explicitly below.
wtd_ds = xr.open_dataset(MERGED / "HRES-WTD_2015_Ogallala.nc",
                         chunks={"y": 2048, "x": 2048}, mask_and_scale=False)
wy, wx = wtd_ds["y"].values, wtd_ds["x"].values
Ny, Nx = len(wy), len(wx)
assert wx[0] < wx[-1], "unexpected WTD x order"
y_desc = wy[0] > wy[-1]
wy_asc = wy[::-1] if y_desc else wy
y_edges = np.concatenate(([GY[0] + (GY[0] - GY[1]) / 2],
                          (GY[:-1] + GY[1:]) / 2,
                          [GY[-1] - (GY[-1] - GY[-2]) / 2]))
x_edges = np.concatenate(([GX[0] - (GX[1] - GX[0]) / 2],
                          (GX[:-1] + GX[1:]) / 2,
                          [GX[-1] + (GX[-1] - GX[-2]) / 2]))
row_ranges, col_ranges = [], []
for i in range(H):
    north, south = y_edges[i], y_edges[i + 1]
    a = np.searchsorted(wy_asc, south)
    b = np.searchsorted(wy_asc, north, side="right")
    row_ranges.append((Ny - b, Ny - a) if y_desc else (a, b))
for j in range(W):
    west, east = x_edges[j], x_edges[j + 1]
    col_ranges.append((np.searchsorted(wx, west),
                       np.searchsorted(wx, east, side="right")))
wtd4 = np.full((H, W), np.nan, dtype="float32")
STRIDE = 8  # target-row stripes per read
for s in range(0, H, STRIDE):
    e = min(s + STRIDE, H)
    rs = min(r[0] for r in row_ranges[s:e])
    re = max(r[1] for r in row_ranges[s:e])
    blk = wtd_ds["b1"].isel(time=0, y=slice(rs, re)).compute().values
    for i in range(s, e):
        acc = np.zeros(W, dtype="float64")
        cnt = np.zeros(W, dtype="float64")
        for j in range(W):
            seg = blk[row_ranges[i][0] - rs:row_ranges[i][1] - rs,
                      col_ranges[j][0]:col_ranges[j][1]].ravel()
            seg = seg[np.isfinite(seg)]
            if seg.size:
                acc[j], cnt[j] = seg.sum(), seg.size
        wtd4[i, cnt > 0] = (acc[cnt > 0] / cnt[cnt > 0]).astype("float32")
    log(f"  WTD stripe rows {s}-{e} done")
wtd_ds.close()
vv = wtd4[np.isfinite(wtd4)]
log(f"WTD 4km: valid {vv.size}, median {np.median(vv):.1f} m, p95 {np.percentile(vv,95):.1f} m")
wtd_tert = np.nanpercentile(wtd4[mask_np], [33.33, 66.67])
log(f"WTD tercile cuts (m): {wtd_tert.round(1).tolist()}")

# irrigation: GIR frac(class==1) 2001-2015 + GFSAD cropland (A11)
gir = xr.open_dataset(MERGED / "Global_irrigation_Area_Merged_Ogallala.nc", chunks="auto")
gir4 = gir["classification"].interp(y=ds_ref["y"], x=ds_ref["x"], method="nearest",
                                    kwargs={"fill_value": np.nan}).compute().values
gir.close()
irr_frac = np.nanmean((gir4 == 1).astype(float), axis=0)
# diagnostic: native-resolution GIR shares inside ROI bbox (no interp involved)
vn = gir["classification"].values.reshape(15, -1)
vn = vn[:, np.isfinite(vn).all(axis=0)]
log(f"  GIR native valid cells: {vn.shape[1]}; "
    f"frac(class==1)>=0.5: {np.mean((vn == 1).mean(axis=0) >= 0.5):.2f}")
gfs = xr.open_dataset(MERGED / "GFSAD1000_V1_2019_Ogallala.nc", chunks="auto")
gfs4 = gfs["landcover"].isel(time=0).interp(
    y=ds_ref["y"], x=ds_ref["x"], method="nearest",
    kwargs={"fill_value": np.nan}).compute().values
gfs.close()
crop = np.isin(gfs4, [1, 2])
mgmt = np.full((H, W), 2, dtype=np.int8)          # 0 irrig / 1 rainfed-crop / 2 natural-other
mgmt[irr_frac >= 0.5] = 0
mgmt[(irr_frac < 0.5) & crop] = 1
log(f"mgmt shares in aquifer: irrig {(mgmt[mask_np]==0).mean():.2f}, "
    f"rainfed-crop {(mgmt[mask_np]==1).mean():.2f}, "
    f"natural-other {(mgmt[mask_np]==2).mean():.2f}")

# PAW soil covariate + dWL screens (already 4 km)
paw4 = xr.open_dataset(DERIVED / "PAW_rootznaws_4km.nc", chunks="auto")[
    "root_zone_available_water_storage"].isel(time=0).compute().values
dwl_long = xr.open_dataset(DERIVED / "dWL_predev_to_2019_ft_4km.nc", chunks="auto")[
    "dWL_predev_to_2019"].compute().values
dwl_recent = xr.open_dataset(DERIVED / "dWL_2017_to_2019_ft_4km.nc", chunks="auto")[
    "dWL_2017_to_2019"].compute().values

# sampled pixel index (stride 2, valid statics) (A8)
yy, xx = np.where(mask_np)
keep = (np.isfinite(pft4[yy, xx]) & np.isfinite(wtd4[yy, xx])
        & np.isfinite(paw4[yy, xx]))
yy, xx = yy[keep], xx[keep]
yy_s = yy[::PIX_STRIDE]
xx_s = xx[::PIX_STRIDE]
NS = yy_s.size
log(f"sampled pixels: {NS} (stride {PIX_STRIDE})")

PFT_S = pft4[yy_s, xx_s].astype(int)
WTD_S = wtd4[yy_s, xx_s]
MGMT_S = mgmt[yy_s, xx_s]
PAW_S = paw4[yy_s, xx_s]
wtd_bin = np.digitize(WTD_S, wtd_tert)  # 0 shallow .. 2 deep

# --------------------------------------------------------------------------
# 2. Dekadal panel (batched, lazy sources)
# --------------------------------------------------------------------------
log("== stage 2: dekadal panel ==")
sseb = xr.open_dataset(MERGED / "MODIS_ET_SSEBop_Merged_Ogallala.nc", chunks={"time": 8})
ST = pd.to_datetime(sseb["time"].values)
ok = (ST >= pd.Timestamp(f"{FULL_YEARS[0]}-01-01")) & \
     (ST < pd.Timestamp(f"{FULL_YEARS[-1] + 1}-01-01"))
sseb_idx = np.where(ok)[0]
NT = len(sseb_idx)
log(f"analysis dekads: {NT} ({ST[sseb_idx[0]].date()}..{ST[sseb_idx[-1]].date()}) (A5)")
T0 = ST[sseb_idx].values.astype("datetime64[D]")
T1 = np.append(ST[sseb_idx[1:]].values.astype("datetime64[D]"),
               (ST[sseb_idx[-1]] + pd.Timedelta(days=11)).date()).astype("datetime64[D]")

grid = xr.open_dataset(MERGED / "GRIDMET_Merged_Ogallala.nc",
                       chunks={"time": 500, "y": 287, "x": 181}, mask_and_scale=False)
GT = pd.to_datetime(grid["time"].values).values.astype("datetime64[D]")
smap = xr.open_dataset(MERGED / "SPL4SMGP_Ogallala_FULL.nc", chunks={"time": 500})
SMT = pd.to_datetime(smap["time"].values)
# SMAP lives on a 2-D EASE grid: precompute nearest SMAP cell per sampled
# 4-km pixel once (A13b: nearest-neighbour SMAP->4 km; ~9 km footprint)
from sklearn.neighbors import NearestNeighbors
_slat = smap["lat"].values
_slon = smap["lon"].values
_nn = NearestNeighbors(n_neighbors=1).fit(
    np.stack([_slat.ravel(), _slon.ravel()], axis=1))
_nn_dist, _nn_ind = _nn.kneighbors(np.stack([GY[yy_s], GX[xx_s]], axis=1))
iy_nn, ix_nn = np.unravel_index(_nn_ind[:, 0], _slat.shape)
log(f"  SMAP nearest mapping: max collocation distance {_nn_dist.max():.3f} deg")
del _slat, _slon, _nn, _nn_dist, _nn_ind
dr = xr.open_dataset(MERGED / "DROUGHT_Merged_Ogallala.nc",
                     chunks={"time": 50, "y": 287, "x": 181}, mask_and_scale=False)
DT = pd.to_datetime(dr["time"].values).values.astype("datetime64[D]")

rows = []
for b in range(0, NT, BATCH_DEKADS):
    t0b = time.time()
    bi = sseb_idx[b:b + BATCH_DEKADS]
    nb = len(bi)
    d0, d1 = T0[b], T1[b + nb - 1]
    # --- SSEBop ET batch -> 4 km (A7: ET mm/dekad as stored) ---
    et_b = sseb["et"].isel(time=bi).interp(
        y=ds_ref["y"], x=ds_ref["x"], method="linear",
        kwargs={"fill_value": np.nan}).compute().values.astype("float32")
    # --- GRIDMET daily window (with 95 d antecedent pad) ---
    g0 = np.searchsorted(GT, d0 - np.timedelta64(95, "D"))
    g1 = np.searchsorted(GT, d1, side="right")
    gsel = slice(g0, g1)
    gd = (GT[g0:g1] - np.datetime64("1970-01-01")).astype(int)  # days
    pr = grid["pr"].isel(time=gsel).compute().values.astype("float32")
    tmn = grid["tmmn"].isel(time=gsel).compute().values.astype("float32")
    tmx = grid["tmmx"].isel(time=gsel).compute().values.astype("float32")
    vpd = grid["vpd"].isel(time=gsel).compute().values.astype("float32")
    srad = grid["srad"].isel(time=gsel).compute().values.astype("float32")
    if np.nanmedian(tmx) > 200:  # Kelvin -> degC (A12)
        tmn, tmx = tmn - 273.15, tmx - 273.15
        log("  GRIDMET T detected as Kelvin; converted to degC")
    tmean = (tmn + tmx) / 2
    # dekad bins for days
    edges = np.append(T0[b:b + nb].astype("datetime64[D]").astype(int),
                      T1[b + nb - 1].astype("datetime64[D]").astype(int))
    bin_id = np.searchsorted(edges, gd, side="right") - 1
    inb = (bin_id >= 0) & (bin_id < nb)
    cnt = np.bincount(bin_id[inb], minlength=nb).astype(float)
    # vectorized dekadal sums/means via index tricks
    Pj = pr[:, yy_s, xx_s]                       # (nday, NS)
    Vj = vpd[:, yy_s, xx_s]
    Sj = srad[:, yy_s, xx_s]
    Tj = tmean[:, yy_s, xx_s]
    Pj = np.where(inb[:, None], Pj, np.nan)
    P_dek = np.stack([np.nansum(np.where(bin_id[:, None] == k, Pj, 0), axis=0)
                      for k in range(nb)], axis=0)
    V_dek = np.stack([np.nanmean(np.where(bin_id[:, None] == k, Vj, np.nan), axis=0)
                      for k in range(nb)], axis=0)
    S_dek = np.stack([np.nanmean(np.where(bin_id[:, None] == k, Sj, np.nan), axis=0)
                      for k in range(nb)], axis=0)
    T_dek = np.stack([np.nanmean(np.where(bin_id[:, None] == k, Tj, np.nan), axis=0)
                      for k in range(nb)], axis=0)
    # antecedent P30/P90 ending at dekad start (exclusive) via cumsum
    cs = np.nancumsum(np.where(np.isfinite(Pj), Pj, 0), axis=0)
    di = (T0[b:b + nb].astype("datetime64[D]").astype(int)
          - gd[0])
    P30 = np.stack([cs[di[k]] - (cs[di[k] - 30] if di[k] - 30 >= 0 else 0)
                    for k in range(nb)], axis=0)
    P90 = np.stack([cs[di[k]] - (cs[di[k] - 90] if di[k] - 90 >= 0 else 0)
                    for k in range(nb)], axis=0)
    del pr, tmn, tmx, vpd, srad, tmean, Pj, Vj, Sj, Tj
    # --- SMAP dekadal means + lag1 ---
    s0 = np.searchsorted(SMT.values, np.datetime64(T0[b]))
    s1 = np.searchsorted(SMT.values, np.datetime64(T1[b + nb - 1]), side="right")
    sm = smap["sm_rootzone"].isel(time=slice(s0, s1)).compute().values
    smt = SMT.values[s0:s1]
    sm_b = np.clip(np.searchsorted(
        np.append(T0[b:b + nb].astype("datetime64[D]").astype("datetime64[ns]"),
                  np.datetime64(T1[b + nb - 1])),
        smt, side="right") - 1, 0, nb - 1)
    SM_dek = np.stack([np.nanmean(sm[sm_b == k], axis=0) for k in range(nb)], axis=0)
    SM_dek = SM_dek[:, iy_nn, ix_nn]          # (nb, NS) via precomputed nearest cells
    SM_lag = np.vstack([np.full(NS, np.nan), SM_dek[:-1]])
    del sm
    # --- drought indices: nearest step at/under dekad start (A13) ---
    TT = T0[b:b + nb].astype("datetime64[D]")
    ii = np.array([np.abs(DT - t).argmin() for t in TT])
    spi = dr["spi90d"].isel(time=xr.DataArray(ii, dims="t")).compute().values[:, yy_s, xx_s]
    spei = dr["spei90d"].isel(time=xr.DataArray(ii, dims="t")).compute().values[:, yy_s, xx_s]
    # --- assemble rows ---
    doy = pd.to_datetime(T0[b:b + nb]).dayofyear.values
    for k in range(nb):
        df = pd.DataFrame({
            "dekad": np.full(NS, int(b + k), np.int32),
            "date": np.full(NS, T0[b + k]),
            "year": T0[b + k].astype("datetime64[Y]").astype(int) + 1970,
            "month": pd.to_datetime(T0[b + k]).month,
            "doy_sin": np.sin(2 * np.pi * doy[k] / 365.25),
            "doy_cos": np.cos(2 * np.pi * doy[k] / 365.25),
            "ET": et_b[k][yy_s, xx_s].astype("float32"),
            "P_dekad": P_dek[k].astype("float32"),
            "P30": P30[k].astype("float32"),
            "P90": P90[k].astype("float32"),
            "SMrz": SM_dek[k].astype("float32"),
            "SMrz_lag1": SM_lag[k].astype("float32"),
            "VPD": V_dek[k].astype("float32"),
            "Rn": S_dek[k].astype("float32"),
            "Tmean": T_dek[k].astype("float32"),
            "PFT": PFT_S,
            "PAW": PAW_S,
            "WTD": WTD_S,
            "wtd_bin": wtd_bin,
            "mgmt": MGMT_S,
            "spi90": spi[k].astype("float32"),
            "spei90": spei[k].astype("float32"),
            "pix": np.arange(NS, dtype=np.int32),
        })
        rows.append(df)
    log(f"  batch dekads {b}-{b+nb} done in {time.time()-t0b:.0f}s "
        f"(rows so far {sum(map(len, rows)):,})")

panel = pd.concat(rows, ignore_index=True)
n0 = len(panel)
panel = panel.dropna(subset=["ET", "P_dekad", "P30", "P90", "SMrz", "SMrz_lag1",
                             "VPD", "Rn", "Tmean"])
log(f"panel: {len(panel):,} rows ({100*len(panel)/n0:.1f}% complete)")

grow = panel["month"].isin(GROW_MONTHS).values
log(f"growing-season rows: {grow.sum():,}")

# value sanity (A7)
log(f"ET mm/dekad: med {panel['ET'].median():.1f}, p99 {panel['ET'].quantile(.99):.1f}")
log(f"P_dekad mm: med {panel['P_dekad'].median():.1f}, p99 {panel['P_dekad'].quantile(.99):.1f}")
log(f"Tmean degC: med {panel['Tmean'].median():.1f}; VPD kPa med {panel['VPD'].median():.2f}; "
    f"srad W/m2 med {panel['Rn'].median():.0f}")

# --------------------------------------------------------------------------
# 3. Train XGBoost (hold out 2 driest growing seasons)
# --------------------------------------------------------------------------
log("== stage 3: XGBoost expected-ET model ==")
gs = panel[panel["month"].isin(GROW_MONTHS)]
yr_dry = gs.groupby("year")["spi90"].mean().sort_values()
log("growing-season dryness ranking (mean spi90d):\n" + yr_dry.round(2).to_string())
HOLDOUT = yr_dry.index[:2].tolist()
log(f"HOLDOUT drought years: {HOLDOUT} (transferability test, plan 6.7)")
VAL_YEAR = yr_dry.index[2:3].tolist()
log(f"early-stopping validation year: {VAL_YEAR}")

FEATS = ["P_dekad", "P30", "P90", "SMrz", "SMrz_lag1", "VPD", "Rn", "Tmean",
         "PAW", "doy_sin", "doy_cos"]
pft_dum = pd.get_dummies(panel["PFT"].astype(int), prefix="PFT", dtype="float32")
X_all = pd.concat([panel[FEATS].astype("float32"), pft_dum], axis=1)
FEATURES = X_all.columns.tolist()
log(f"features ({len(FEATURES)}): {FEATURES}")
y_all = panel["ET"].astype("float32").values

is_hold = panel["year"].isin(HOLDOUT).values
is_val = panel["year"].isin(VAL_YEAR).values
is_tr = ~(is_hold | is_val)

model = xgb.XGBRegressor(**XGB_PARAMS)
model.fit(X_all.values[is_tr], y_all[is_tr],
          eval_set=[(X_all.values[is_val], y_all[is_val])],
          verbose=False)
log(f"best iteration: {model.best_iteration}")
model.save_model(str(OUT / "xgb_et_model.json"))

pred = model.predict(X_all.values).astype("float32")
panel["ET_hat"] = pred
panel["Bi"] = panel["ET"] - panel["ET_hat"]   # Plan Eq. 6.1 (ET leg)

mets = {}
for tag, m_ in [("train", is_tr), ("valid", is_val), ("holdout", is_hold)]:
    aa, bb = panel.loc[m_, "ET"].values, panel.loc[m_, "ET_hat"].values
    mets[tag] = {"rmse": rmse(aa, bb), "n": int(m_.sum())}
r2 = lambda a, b: float(1 - np.sum((a - b) ** 2) / np.sum((a - a.mean()) ** 2))
for tag, m_ in [("train", is_tr), ("valid", is_val), ("holdout", is_hold)]:
    a, b = panel.loc[m_, "ET"].values, panel.loc[m_, "ET_hat"].values
    mets[tag].update(mae=float(np.mean(np.abs(a - b))), r2=r2(a, b),
                     n=int(m_.sum()))
log("metrics:\n" + json.dumps(mets, indent=1))

# feature importance (gain)
imp = model.get_booster().get_score(importance_type="gain")
imp = pd.Series({FEATURES[int(k[1:])]: v for k, v in imp.items()}).sort_values()
log("top gain features:\n" + imp.tail(8).round(1).to_string())

# --------------------------------------------------------------------------
# 4. Dry-downs + buffering quantification
# --------------------------------------------------------------------------
log("== stage 4: buffering ==")
# short-baseline dekad-of-year climatology (A14)
panel["doyk"] = pd.to_datetime(panel["date"]).dt.strftime("%m-%d")
clim = panel.groupby("doyk")[["SMrz", "VPD", "P_dekad"]].transform("mean")
panel["SMrz_anom"] = panel["SMrz"] - clim["SMrz"]
panel["VPD_anom"] = panel["VPD"] - clim["VPD"]
dry = ((panel["spi90"] <= DRYDOWN_SPI) & (panel["SMrz_anom"] < 0)
       & (panel["VPD_anom"] > 0)).values
panel["dry"] = dry
log(f"dry-down rows: {dry.sum():,} ({100*dry.mean():.1f}%); "
    f"growing-season dry: {(dry & grow).sum():,}")

dryg = dry & grow
B_all, B_ci, Bn = mean_ci(panel.loc[dryg, "Bi"])
log(f"mean Bi | dry-down, growing season: {B_all:+.2f} +/- {B_ci:.2f} mm/dekad (n={Bn:,})")
m_ = is_hold & dry
aa, bb = panel.loc[m_, "ET"].values, panel.loc[m_, "ET_hat"].values
mets["holdout_drydown"] = {"rmse": rmse(aa, bb), "n": int(m_.sum())}

effects = []


def add_effect(group, label, sub):
    m_, ci_, n_ = mean_ci(panel.loc[sub, "Bi"])
    effects.append(dict(group=group, label=label, mean_Bi_mm_per_dekad=m_,
                        ci95=ci_, n=n_))


for b_, name in [(0, "shallow-T1"), (1, "mid-T2"), (2, "deep-T3")]:
    add_effect("WTD_tercile", name, dryg & (panel["wtd_bin"].values == b_))
for g_, name in [(0, "irrigated"), (1, "rainfed-crop"), (2, "natural-other")]:
    add_effect("mgmt", name, dryg & (panel["mgmt"].values == g_))
for p in sorted(panel["PFT"].unique()):
    add_effect("PFT", f"IGBP_{int(p)}", dryg & (panel["PFT"].values == p))

# GW-access effect: shallow-minus-deep within strata (A16 lite matching)
sh = (panel["wtd_bin"].values == 0) & dryg
dp = (panel["wtd_bin"].values == 2) & dryg
eff_rows = []
for (pft, mg, p3t), idx in panel[dryg].groupby(
        [panel.loc[dryg, "PFT"],
         panel.loc[dryg, "mgmt"],
         pd.qcut(panel.loc[dryg, "P30"], 3, labels=False, duplicates="drop")]).groups.items():
    a = panel.loc[idx, "Bi"].values[panel.loc[idx, "wtd_bin"].values == 0]
    b_ = panel.loc[idx, "Bi"].values[panel.loc[idx, "wtd_bin"].values == 2]
    if len(a) >= 30 and len(b_) >= 30:
        eff_rows.append((float(np.mean(a) - np.mean(b_)), len(a) + len(b_)))
eff_rows = np.array(eff_rows)
tau = float(np.average(eff_rows[:, 0], weights=eff_rows[:, 1])) if len(eff_rows) else float("nan")
log(f"GW-access effect tau (shallow-deep Bi, stratified, {len(eff_rows)} strata): "
    f"{tau:+.2f} mm/dekad")

# wet-year placebo (A17): same contrast in 2 wettest growing seasons
yr_wet = gs.groupby("year")["spi90"].mean().sort_values(ascending=False).index[:2].tolist()
wet = panel["year"].isin(yr_wet).values & grow
pw_rows = []
for (pft, mg, p3t), idx in panel[wet].groupby(
        [panel.loc[wet, "PFT"],
         panel.loc[wet, "mgmt"],
         pd.qcut(panel.loc[wet, "P30"], 3, labels=False, duplicates="drop")]).groups.items():
    a = panel.loc[idx, "Bi"].values[panel.loc[idx, "wtd_bin"].values == 0]
    b_ = panel.loc[idx, "Bi"].values[panel.loc[idx, "wtd_bin"].values == 2]
    if len(a) >= 30 and len(b_) >= 30:
        pw_rows.append((float(np.mean(a) - np.mean(b_)), len(a) + len(b_)))
pw_rows = np.array(pw_rows)
tau_wet = float(np.average(pw_rows[:, 0], weights=pw_rows[:, 1])) if len(pw_rows) else float("nan")
log(f"placebo tau in wet years {yr_wet}: {tau_wet:+.2f} mm/dekad (expect ~0)")

eff_df = pd.DataFrame(effects)
eff_df.to_csv(OUT / "effects.csv", index=False)

# sustainability: pixel-mean dry Bi vs observed decline bins (A19)
pix_mean = panel[dryg].groupby("pix")["Bi"].mean()
sdf = pd.DataFrame({
    "Bi": pix_mean,
    "decl": pd.Series(dwl_long[yy_s, xx_s], index=np.arange(NS))})
sdf["bin"] = pd.cut(sdf["decl"], [-np.inf, -50, -10, 10, np.inf],
                    labels=["decline>50ft", "decline 10-50ft",
                            "stable +-10ft", "rise>10ft"])
sust = sdf.dropna(subset=["Bi", "bin"])
sust_tab = sust.groupby("bin", observed=True)["Bi"].agg(["mean", "count"])

# dry-down mean map on sampled lattice (A8: maps use sampled pixels)
dmap = np.full((H, W), np.nan, dtype="float32")
dmap[yy_s, xx_s] = pix_mean.reindex(pd.Index(np.arange(NS))).values
log("sustainability screen (mean dry Bi by observed decline bin):\n"
    + sust_tab.round(2).to_string())

# volume plausibility (A20): seasonal extra ET vs dWL*Sy.
# Restrict to pixels with >=3 growing-season dry dekads so the seasonal sum
# is meaningful (pixels with 0 dry dekads contribute exactly 0 by construction).
seas = panel[dryg].groupby("pix")["Bi"].agg(["sum", "size"])
seas = seas[seas["size"] >= 3]["sum"]      # mm per growing-season dry dekads
pix_dwl = pd.Series(dwl_long[yy_s, xx_s], index=np.arange(NS))
vv = pd.DataFrame({"extraET_mm": seas, "dWL_ft": pix_dwl}).dropna()
vv["implied_depletion_mm"] = -vv["dWL_ft"] * SY * 304.8
vv = vv[np.isfinite(vv["implied_depletion_mm"])]
r_vol = float(vv[["extraET_mm", "implied_depletion_mm"]].corr().iloc[0, 1]) \
    if len(vv) > 10 else float("nan")
log(f"volume check (n={len(vv):,} pix, >=3 dry dekads): "
    f"mean extraET {vv['extraET_mm'].mean():+.1f} mm/season vs "
    f"mean implied depletion {vv['implied_depletion_mm'].mean():+.1f} mm; "
    f"median {vv['extraET_mm'].median():+.1f} vs "
    f"{vv['implied_depletion_mm'].median():+.1f}; r={r_vol:+.2f}")

with open(OUT / "metrics.json", "w") as f:
    json.dump(dict(metrics=mets, holdout_years=HOLDOUT, valid_year=VAL_YEAR,
                   mean_Bi_dry=dict(mean=B_all, ci95=B_ci, n=Bn),
                   tau_shallow_minus_deep_mm_per_dekad=tau,
                   placebo_wet_years=yr_wet, placebo_tau=tau_wet,
                   volume_check=dict(n_pixels=int(len(vv)),
                                     mean_extraET_mm=float(vv["extraET_mm"].mean()),
                                     mean_implied_depletion_mm=float(
                                         vv["implied_depletion_mm"].mean()),
                                     pearson_r=r_vol),
                   sustainability=sust_tab.reset_index().to_dict("records"),
                   xgb_params={k: str(v) for k, v in XGB_PARAMS.items()},
                   n_rows=len(panel), n_pixels=NS, n_dekads=NT,
                   specific_yield=SY),
              f, indent=1, default=str)

# Bi dry-down mean map (CF-lite, 4 km)
yc = GY.copy()
out_ds = xr.Dataset(
    {"Bi_drydown_mean": (("y", "x"), dmap,
                         {"long_name": "mean XGBoost ET residual during growing-season "
                                       "dry-downs (candidate buffering signal)",
                          "units": "mm/formula dekad", "grid_mapping": "crs"})},
    coords={"y": GY, "x": GX},
    attrs={"Conventions": "CF-1.8",
           "title": "ET buffering signal Bi(t) mean, 2016-2021 growing-season dry-downs",
           "history": "gw_buffering_xgb_et.py; sources read-only"})
out_ds["crs"] = xr.DataArray(np.int32(4326), attrs={
    "grid_mapping_name": "latitude_longitude",
    "longitude_of_prime_meridian": 0.0, "semi_major_axis": 6378137.0,
    "inverse_flattening": 298.257223563, "epsg_code": "EPSG:4326"})
out_ds.to_netcdf(OUT / "Bi_drydown_mean_4km.nc", engine="netcdf4",
                 encoding={"Bi_drydown_mean": {"dtype": "float32", "zlib": True,
                                               "complevel": 3,
                                               "_FillValue": np.float32(9.96921e36)}})
log("saved Bi_drydown_mean_4km.nc")

# --------------------------------------------------------------------------
# 5. Figures
# --------------------------------------------------------------------------
log("== stage 5: figures ==")
a_tr, b_tr = panel.loc[is_tr, "ET"].values, panel.loc[is_tr, "ET_hat"].values
a_ho, b_ho = panel.loc[is_hold, "ET"].values, panel.loc[is_hold, "ET_hat"].values

fig, ax = plt.subplots(figsize=(6, 5.5))
ax.hexbin(a_tr[::7], b_tr[::7], gridsize=60, cmap="Blues", mincnt=1)
ax.plot([0, 70], [0, 70], "r--", lw=1)
ax.set_xlabel("observed SSEBop ET [mm/dekad]")
ax.set_ylabel("XGBoost expected ET [mm/dekad]")
ax.set_title(f"Train fit (R2={mets['train']['r2']:.3f}, RMSE={mets['train']['rmse']:.2f})")
fig.tight_layout()
fig.savefig(OUT / "pred_vs_obs_train.png", dpi=170)
plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 5.5))
ax.hexbin(a_ho, b_ho, gridsize=60, cmap="Greens", mincnt=1)
ax.plot([0, 70], [0, 70], "r--", lw=1)
ax.set_xlabel("observed SSEBop ET [mm/dekad]")
ax.set_ylabel("XGBoost expected ET [mm/dekad]")
ax.set_title(f"Holdout {HOLDOUT} (R2={mets['holdout']['r2']:.3f}, "
             f"RMSE={mets['holdout']['rmse']:.2f})")
fig.tight_layout()
fig.savefig(OUT / "pred_vs_obs_holdout.png", dpi=170)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(panel.loc[is_tr, "Bi"].values[::10], bins=80, alpha=0.6, label="train",
        density=True)
ax.hist(panel.loc[is_hold, "Bi"].values, bins=80, alpha=0.6, label="holdout",
        density=True)
ax.axvline(0, color="k", lw=1)
ax.set_xlabel("Bi = ET_obs - ET_hat [mm/dekad]")
ax.set_title("Residual distribution")
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "residual_hist.png", dpi=170)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.5))
imp.tail(12).plot.barh(ax=ax, color="#2ca25f")
ax.set_xlabel("gain")
ax.set_title("XGBoost feature importance (gain)")
fig.tight_layout()
fig.savefig(OUT / "feat_importance.png", dpi=170)
plt.close(fig)

plot_grid(xr.DataArray(dmap, coords={"y": GY, "x": GX}, dims=("y", "x")),
          "Mean ET residual during growing-season dry-downs\nBi = obs - XGB "
          f"(positive = persisting ET; n={Bn:,} pixel-dekads)",
          cmap="BrBG", vmin=-8, vmax=8, cbar_label="Bi [mm/dekad]",
          bnd=bnd, fname="map_Bi_drydown_mean.png")

w = eff_df[eff_df.group == "WTD_tercile"].set_index("label").loc[
    ["shallow-T1", "mid-T2", "deep-T3"]]
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.bar(w.index, w["mean_Bi_mm_per_dekad"], yerr=w["ci95"], capsize=4,
       color=["#1a9850", "#fee08b", "#d73027"])
ax.axhline(0, color="k", lw=0.8)
ax.set_ylabel("mean Bi | dry-down [mm/dekad]")
ax.set_title("Buffering vs static depth-to-water (HRES-WTD terciles)")
fig.tight_layout()
fig.savefig(OUT / "buffering_by_wtd.png", dpi=170)
plt.close(fig)

mg = eff_df[eff_df.group == "mgmt"].set_index("label").loc[
    ["irrigated", "rainfed-crop", "natural-other"]]
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.bar(mg.index, mg["mean_Bi_mm_per_dekad"], yerr=mg["ci95"], capsize=4,
       color=["#4575b4", "#fee090", "#5aa469"])
ax.axhline(0, color="k", lw=0.8)
ax.set_ylabel("mean Bi | dry-down [mm/dekad]")
ax.set_title("Buffering vs management (GIR+GFSAD stratification)")
fig.tight_layout()
fig.savefig(OUT / "buffering_by_mgmt.png", dpi=170)
plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(["GW-access effect\n(shallow-deep, stratified)",
        f"wet-year placebo\n({', '.join(map(str, yr_wet))})"],
       [tau, tau_wet], color=["#313695", "#a8a8a8"])
ax.axhline(0, color="k", lw=0.8)
ax.set_ylabel("tau [mm/dekad]")
ax.set_title("Causal-lite contrast + required placebo check (Plan 6.2)")
fig.tight_layout()
fig.savefig(OUT / "effect_vs_placebo.png", dpi=170)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.5))
sust.boxplot(column="Bi", by="bin", ax=ax, grid=False)
ax.set_xlabel("observed predev->2019 decline bin")
ax.set_ylabel("pixel-mean dry Bi [mm/dekad]")
ax.set_title("Sustainability screen: persistence vs observed depletion")
fig.suptitle("")
fig.tight_layout()
fig.savefig(OUT / "sustainability_screen.png", dpi=170)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.5))
sub = vv.sample(min(20000, len(vv)), random_state=SEED)
ax.scatter(sub["implied_depletion_mm"], sub["extraET_mm"], s=2, alpha=0.3)
mx = float(max(sub["implied_depletion_mm"].quantile(0.99),
               sub["extraET_mm"].quantile(0.99), 1))
ax.plot([0, mx], [0, mx], "r--", lw=1, label="1:1")
ax.set_xlabel("implied depletion |dWL|*Sy [mm water]")
ax.set_ylabel("seasonal extra ET during dry-downs [mm]")
ax.set_title(f"Volume plausibility (Sy={SY}; medians "
             f"{vv['extraET_mm'].median():+.0f} vs {vv['implied_depletion_mm'].median():+.0f} mm)")
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "volume_plausibility.png", dpi=170)
plt.close(fig)

# GRACE basin constraint (z-scored spatial mean; coarse grid, A21)
gr = xr.open_dataset(MERGED / "GRACE_Merged_Ogallala.nc", chunks="auto")
ser = gr["lwe_thickness"].mean(dim=("y", "x"), skipna=True).compute()
gt = pd.to_datetime(gr["time"].values)
gr.close()
z = (ser.values - np.nanmean(ser.values)) / np.nanstd(ser.values)
fig, ax = plt.subplots(figsize=(9, 3.5))
ax.plot(gt, z, color="#313695", lw=1)
ax.axhline(0, color="k", lw=0.8)
for y0 in HOLDOUT:
    ax.axvspan(pd.Timestamp(f"{y0}-01-01"), pd.Timestamp(f"{y0}-12-31"),
               color="orange", alpha=0.2)
ax.set_ylabel("basin-mean TWS anomaly [z]")
ax.set_title("GRACE basin constraint (holdout drought years shaded)")
fig.tight_layout()
fig.savefig(OUT / "grace_basin.png", dpi=170)
plt.close(fig)

# --------------------------------------------------------------------------
# 6. AmeriFlux point validation (A22)
# --------------------------------------------------------------------------
log("== stage 6: AmeriFlux validation ==")
try:
    sites = pd.read_excel(BASE / "AmeriFlux_NEON_sites.xlsx")
    sites.columns = [c.strip() for c in sites.columns]
    roi_s = sites[(sites["Long"] >= -106.2) & (sites["Long"] <= -96.0)
                  & (sites["Lat"] >= 31.4) & (sites["Lat"] <= 44.0)].copy()
    am = xr.open_dataset(MERGED / "merged_ameriflux.nc", chunks={"record": 1_000_000})
    log(f"ameriflux time attrs: {dict(am['time'].attrs)}; "
        f"sample values: {np.asarray(am['time'].values[:5]).tolist()}")
    lut = {bytes(s).decode().strip(): i
           for i, s in enumerate(am["site_lookup"].values)}
    sseb1 = xr.open_dataset(MERGED / "MODIS_ET_SSEBop_Merged_Ogallala.nc", chunks={"time": 8})
    val_rows = []
    for _, r in roi_s.iterrows():
        sid = str(r["Site_Id"]).strip()
        if sid not in lut:
            continue
        k = lut[sid]
        si = am["site_index"].values
        m_ = si == k
        if m_.sum() < 100:
            continue
        tm = am["time"].values[m_]
        # time units unknown -> infer from magnitude (documented A22):
        # ~1.45e9 with 1800 s steps = Unix seconds, half-hourly records
        tmax = float(np.nanmax(tm[np.isfinite(tm)])) if np.isfinite(tm).any() else 0.0
        unit = "s" if tmax > 1e9 else ("h" if tmax > 1e6 else "D")
        origin = pd.Timestamp("1970-01-01")
        try:
            tdt = pd.to_datetime(tm, unit=unit, origin=origin)
        except Exception:
            continue
        le = am["LE"].values[m_].astype(float)
        ok_ = np.isfinite(le) & (tdt >= pd.Timestamp("2016-01-01")) \
            & (tdt < pd.Timestamp("2022-01-01"))
        if ok_.sum() < 500:
            continue
        df_t = pd.DataFrame({"t": tdt[ok_], "LE": le[ok_]})
        df_t["et_mm"] = df_t["LE"] * 1800 / 2.45e6  # half-hourly W/m2 -> mm;
        # 1800 s * LE[W/m2=J/s/m2] / 2.45e6[J/kg] = kg/m2 = mm water (A22)
        df_t["dek"] = df_t["t"].values.astype("datetime64[D]")
        # map each record to its SSEBop dekad start
        di_ = np.searchsorted(ST.values.astype("datetime64[D]"),
                              df_t["dek"].values, side="right") - 1
        df_t = df_t[(di_ >= sseb_idx[0]) & (di_ < sseb_idx[-1])]
        if len(df_t) < 200:
            continue
        twr = df_t.groupby(
            ST.values.astype("datetime64[D]")[np.clip(
                np.searchsorted(ST.values.astype("datetime64[D]"),
                                df_t["dek"].values, side="right") - 1,
                0, len(ST) - 1)])["et_mm"].sum()
        # nearest SSEBop pixel
        ex = sseb1["x"].values
        ey = sseb1["y"].values
        jx = int(np.abs(ex - r["Long"]).argmin())
        jy = int(np.abs(ey - r["Lat"]).argmin())
        dates = pd.to_datetime(twr.index)
        sse = sseb1["et"].sel(
            x=ex[jx], y=ey[jy], method="nearest").sel(
            time=slice(str(dates.min().date()), str(dates.max().date()))).compute()
        st2 = pd.to_datetime(sse["time"].values)
        common = twr.index.intersection(st2)
        if len(common) < 20:
            continue
        val_rows.append(pd.DataFrame(
            {"site": sid, "date": common,
             "tower_ET": twr.loc[common].values,
             "ssebop_ET": sse.sel(time=common).values}))
    if val_rows:
        V = pd.concat(val_rows, ignore_index=True)
        V.to_csv(OUT / "ameriflux_validation.csv", index=False)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        axes[0].scatter(V["tower_ET"], V["ssebop_ET"], s=4, alpha=0.4)
        axes[0].plot([0, 70], [0, 70], "r--", lw=1)
        axes[0].set_xlabel("tower ET [mm/dekad]")
        axes[0].set_ylabel("SSEBop ET [mm/dekad]")
        axes[0].set_title(f"SSEBop vs towers (RMSE={rmse(V['tower_ET'], V['ssebop_ET']):.2f}, "
                          f"n={len(V):,}, {V['site'].nunique()} sites)")
        axes[1].hist((V["ssebop_ET"] - V["tower_ET"]).values, bins=60, color="#737373")
        axes[1].axvline(0, color="k", lw=1)
        axes[1].set_xlabel("SSEBop - tower [mm/dekad]")
        axes[1].set_title("SSEBop bias distribution")
        fig.tight_layout()
        fig.savefig(OUT / "ameriflux_validation.png", dpi=170)
        plt.close(fig)
        log(f"AmeriFlux validation: {len(V):,} dekads, {V['site'].nunique()} sites, "
            f"SSEBop RMSE={rmse(V['tower_ET'], V['ssebop_ET']):.2f} mm/dekad")
    else:
        log("AmeriFlux validation SKIPPED: no ROI tower with >=20 overlapping dekads")
    am.close()
    sseb1.close()
except Exception as e:
    log(f"AmeriFlux validation SKIPPED ({type(e).__name__}: {str(e)[:200]})")

# --------------------------------------------------------------------------
# wrap-up
# --------------------------------------------------------------------------
sseb.close()
grid.close()
smap.close()
dr.close()
ds_ref.close()
with open(OUT / "run_log.txt", "w") as f:
    f.write("\n".join(LOG_LINES))
log(f"DONE. outputs in {OUT}")
log("files: " + ", ".join(sorted(p.name for p in OUT.iterdir())))
