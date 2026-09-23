"""
Monthly buffering maps + ROI-averaged buffering time series.

Reuses the verified `gw_buffering_xgb_et.py` pipeline (same features, same
sampling, same dry-down definition) but loads the saved XGBoost model instead
of retraining, then:
  1. predicts ET_hat for every sampled pixel-dekad (2016-2021),
  2. Bi = ET_obs - ET_hat  (Plan Sec. 6.1, ET leg),
  3. aggregates Bi to MONTHLY means per pixel -> CF-lite NetCDF (all months),
  4. saves a spatial map PNG for a requested --year/--month (default: 2020-08,
     a holdout drought-year growing-season month),
  5. saves ROI-(aquifer-)averaged monthly Bi time series -> CSV + PNG
     (with monthly dry-down fraction on a twin axis).

Run:  python buffering_monthly_maps_timeseries.py [--year 2020] [--month 8]
Sources read-only; outputs are NEW files in gw_buffering_xgb/.
Assumptions: same A-register as gw_buffering_xgb_et.py (TBI_DOCUMENTATION.md
Sec. 11.15), plus B1-B4 below.
  B1 Monthly Bi = mean of dekadal Bi whose dekad STARTS in that month
     (usually 3 dekads; edge months may have 2-3).
  B2 ROI average = mean over sampled aquifer pixels (stride 2, A8); dry
     fraction = share of sampled pixel-dekads flagged dry-down in that month.
  B3 Requested month must lie in 2016-2021 (SMAP x SSEBop overlap, A5).
  B4 Feature columns are forced to the exact training order (hard-coded from
     the verified run log); any missing PFT dummy is filled with 0.
"""

import argparse
import json
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

warnings.filterwarnings("ignore")

BASE = Path(r"G:\MSU_GWB\datasets")
MERGED = BASE / "merged_datasets"
DERIVED = BASE / "derived_usgs"
OUT = BASE / "gw_buffering_xgb"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 42
PIX_STRIDE = 2
BATCH_DEKADS = 32
FULL_YEARS = list(range(2016, 2022))
GROW_MONTHS = [5, 6, 7, 8, 9]
DRYDOWN_SPI = -0.8

# Exact training feature order (from verified run log, stage 3)
FEATURES = ["P_dekad", "P30", "P90", "SMrz", "SMrz_lag1", "VPD", "Rn",
            "Tmean", "PAW", "doy_sin", "doy_cos",
            "PFT_7", "PFT_10", "PFT_11", "PFT_12", "PFT_13", "PFT_16",
            "PFT_17"]
FEATS = FEATURES[:11]

LOG_LINES = []


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    LOG_LINES.append(line)


def grid_extent(da):
    x = np.asarray(da["x"].values, dtype=float)
    y = np.asarray(da["y"].values, dtype=float)
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    if dy < 0:
        return ([x[0] - dx / 2, x[-1] + dx / 2,
                 y[-1] + dy / 2, y[0] - dy / 2], "upper")
    return ([x[0] - dx / 2, x[-1] + dx / 2,
             y[0] - dy / 2, y[-1] + dy / 2], "lower")


def main(year, month):
    assert year in FULL_YEARS, f"--year must be in {FULL_YEARS}"
    assert 1 <= month <= 12, "--month must be 1..12"
    log(f"requested monthly map: {year}-{month:02d}")

    import geopandas as gpd
    bnd = gpd.read_file(
        BASE / "high_plains_quifer" / "hp_bound2010.shp").to_crs("EPSG:4326")

    # ---------------- 1. statics (subset needed for features) ----------------
    log("== statics ==")
    ds_ref = xr.open_dataset(MERGED / "GRIDMET_Merged_Ogallala.nc",
                             chunks={"time": 1})
    GY, GX = ds_ref["y"].values, ds_ref["x"].values
    H, W = len(GY), len(GX)
    mask = xr.open_dataset(DERIVED / "aquifer_mask_4km.nc",
                           chunks="auto")["aquifer_mask"]
    mask_np = (mask.values == 1)

    mcd = xr.open_dataset(MERGED / "MCD12Q1_Merged_Ogallala.nc", chunks={"time": 1})
    pft4 = mcd["LC_Type1"].sel(time="2019-01-01", method="nearest").interp(
        y=ds_ref["y"], x=ds_ref["x"], method="nearest",
        kwargs={"fill_value": np.nan}).compute().values
    mcd.close()

    paw4 = xr.open_dataset(DERIVED / "PAW_rootznaws_4km.nc", chunks="auto")[
        "root_zone_available_water_storage"].isel(time=0).compute().values

    yy, xx = np.where(mask_np)
    keep = (np.isfinite(pft4[yy, xx]) & np.isfinite(paw4[yy, xx]))
    yy, xx = yy[keep], xx[keep]
    yy_s = yy[::PIX_STRIDE]
    xx_s = xx[::PIX_STRIDE]
    NS = yy_s.size
    log(f"sampled pixels: {NS}")
    PFT_S = pft4[yy_s, xx_s].astype(int)
    PAW_S = paw4[yy_s, xx_s]

    # ---------------- 2. dekadal panel (same build as training) ----------------
    log("== dekadal panel ==")
    sseb = xr.open_dataset(MERGED / "MODIS_ET_SSEBop_Merged_Ogallala.nc",
                           chunks={"time": 8})
    ST = pd.to_datetime(sseb["time"].values)
    ok = (ST >= pd.Timestamp(f"{FULL_YEARS[0]}-01-01")) & \
         (ST < pd.Timestamp(f"{FULL_YEARS[-1] + 1}-01-01"))
    sseb_idx = np.where(ok)[0]
    NT = len(sseb_idx)
    T0 = ST[sseb_idx].values.astype("datetime64[D]")
    T1 = np.append(ST[sseb_idx[1:]].values.astype("datetime64[D]"),
                   (ST[sseb_idx[-1]] + pd.Timedelta(days=11)).date()
                   ).astype("datetime64[D]")

    grid = xr.open_dataset(MERGED / "GRIDMET_Merged_Ogallala.nc",
                           chunks={"time": 500, "y": 287, "x": 181},
                           mask_and_scale=False)
    GT = pd.to_datetime(grid["time"].values).values.astype("datetime64[D]")
    smap = xr.open_dataset(MERGED / "SPL4SMGP_Ogallala_FULL.nc", chunks={"time": 500})
    SMT = pd.to_datetime(smap["time"].values)
    from sklearn.neighbors import NearestNeighbors
    _slat = smap["lat"].values
    _slon = smap["lon"].values
    _nn = NearestNeighbors(n_neighbors=1).fit(
        np.stack([_slat.ravel(), _slon.ravel()], axis=1))
    _, _nn_ind = _nn.kneighbors(np.stack([GY[yy_s], GX[xx_s]], axis=1))
    iy_nn, ix_nn = np.unravel_index(_nn_ind[:, 0], _slat.shape)
    del _slat, _slon, _nn, _nn_ind
    dr = xr.open_dataset(MERGED / "DROUGHT_Merged_Ogallala.nc",
                         chunks={"time": 50, "y": 287, "x": 181},
                         mask_and_scale=False)
    DT = pd.to_datetime(dr["time"].values).values.astype("datetime64[D]")

    rows = []
    for b in range(0, NT, BATCH_DEKADS):
        bi = sseb_idx[b:b + BATCH_DEKADS]
        nb = len(bi)
        d0, d1 = T0[b], T1[b + nb - 1]
        et_b = sseb["et"].isel(time=bi).interp(
            y=ds_ref["y"], x=ds_ref["x"], method="linear",
            kwargs={"fill_value": np.nan}).compute().values.astype("float32")
        g0 = np.searchsorted(GT, d0 - np.timedelta64(95, "D"))
        g1 = np.searchsorted(GT, d1, side="right")
        gd = (GT[g0:g1] - np.datetime64("1970-01-01")).astype(int)
        pr = grid["pr"].isel(time=slice(g0, g1)).compute().values.astype("float32")
        tmn = grid["tmmn"].isel(time=slice(g0, g1)).compute().values.astype("float32")
        tmx = grid["tmmx"].isel(time=slice(g0, g1)).compute().values.astype("float32")
        vpd = grid["vpd"].isel(time=slice(g0, g1)).compute().values.astype("float32")
        srad = grid["srad"].isel(time=slice(g0, g1)).compute().values.astype("float32")
        if np.nanmedian(tmx) > 200:
            tmn, tmx = tmn - 273.15, tmx - 273.15
        tmean = (tmn + tmx) / 2
        edges = np.append(T0[b:b + nb].astype("datetime64[D]").astype(int),
                          T1[b + nb - 1].astype("datetime64[D]").astype(int))
        bin_id = np.searchsorted(edges, gd, side="right") - 1
        inb = (bin_id >= 0) & (bin_id < nb)
        Pj = np.where(inb[:, None], pr[:, yy_s, xx_s], np.nan)
        Vj = vpd[:, yy_s, xx_s]
        Sj = srad[:, yy_s, xx_s]
        Tj = tmean[:, yy_s, xx_s]
        P_dek = np.stack([np.nansum(np.where(bin_id[:, None] == k, Pj, 0), axis=0)
                          for k in range(nb)], axis=0)
        V_dek = np.stack([np.nanmean(np.where(bin_id[:, None] == k, Vj, np.nan), axis=0)
                          for k in range(nb)], axis=0)
        S_dek = np.stack([np.nanmean(np.where(bin_id[:, None] == k, Sj, np.nan), axis=0)
                          for k in range(nb)], axis=0)
        T_dek = np.stack([np.nanmean(np.where(bin_id[:, None] == k, Tj, np.nan), axis=0)
                          for k in range(nb)], axis=0)
        cs = np.nancumsum(np.where(np.isfinite(Pj), Pj, 0), axis=0)
        di = (T0[b:b + nb].astype("datetime64[D]").astype(int) - gd[0])
        P30 = np.stack([cs[di[k]] - (cs[di[k] - 30] if di[k] - 30 >= 0 else 0)
                        for k in range(nb)], axis=0)
        P90 = np.stack([cs[di[k]] - (cs[di[k] - 90] if di[k] - 90 >= 0 else 0)
                        for k in range(nb)], axis=0)
        del pr, tmn, tmx, vpd, srad, tmean, Pj, Vj, Sj, Tj
        s0 = np.searchsorted(SMT.values, np.datetime64(T0[b]))
        s1 = np.searchsorted(SMT.values, np.datetime64(T1[b + nb - 1]), side="right")
        sm = smap["sm_rootzone"].isel(time=slice(s0, s1)).compute().values
        smt = SMT.values[s0:s1]
        sm_b = np.clip(np.searchsorted(
            np.append(T0[b:b + nb].astype("datetime64[D]").astype("datetime64[ns]"),
                      np.datetime64(T1[b + nb - 1])),
            smt, side="right") - 1, 0, nb - 1)
        SM_dek = np.stack([np.nanmean(sm[sm_b == k], axis=0)
                           for k in range(nb)], axis=0)[:, iy_nn, ix_nn]
        SM_lag = np.vstack([np.full(NS, np.nan), SM_dek[:-1]])
        del sm
        TT = T0[b:b + nb].astype("datetime64[D]")
        ii = np.array([np.abs(DT - t).argmin() for t in TT])
        spi = dr["spi90d"].isel(
            time=xr.DataArray(ii, dims="t")).compute().values[:, yy_s, xx_s]
        doy = pd.to_datetime(T0[b:b + nb]).dayofyear.values
        for k in range(nb):
            rows.append(pd.DataFrame({
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
                "spi90": spi[k].astype("float32"),
                "pix": np.arange(NS, dtype=np.int32),
            }))
        log(f"  batch {b}-{b + nb} ({len(rows[-1]) * nb:,} rows this batch)")
    panel = pd.concat(rows, ignore_index=True)
    panel = panel.dropna(subset=["ET", "P_dekad", "P30", "P90", "SMrz",
                                 "SMrz_lag1", "VPD", "Rn", "Tmean"])
    log(f"panel: {len(panel):,} complete rows")

    # ---------------- 3. predict with saved model -> Bi ----------------
    log("== predict Bi with saved model ==")
    model = xgb.XGBRegressor()
    model.load_model(str(OUT / "xgb_et_model.json"))
    pft_dum = pd.get_dummies(panel["PFT"].astype(int), prefix="PFT",
                             dtype="float32")
    X_all = pd.concat([panel[FEATS].astype("float32"), pft_dum],
                      axis=1).reindex(columns=FEATURES, fill_value=0)
    panel["ET_hat"] = model.predict(X_all.values).astype("float32")
    panel["Bi"] = panel["ET"] - panel["ET_hat"]
    log(f"Bi: mean {panel['Bi'].mean():+.3f}, std {panel['Bi'].std():.3f} mm/dekad")

    # dry-down flag (identical definition, A10/A14)
    panel["doyk"] = pd.to_datetime(panel["date"]).dt.strftime("%m-%d")
    clim = panel.groupby("doyk")[["SMrz", "VPD"]].transform("mean")
    panel["dry"] = ((panel["spi90"] <= DRYDOWN_SPI)
                    & ((panel["SMrz"] - clim["SMrz"]) < 0)
                    & ((panel["VPD"] - clim["VPD"]) > 0)).values

    # ---------------- 4. monthly aggregation ----------------
    log("== monthly aggregation ==")
    panel["ym"] = (pd.to_datetime(panel["date"]).dt.strftime("%Y-%m"))
    months = pd.period_range(f"{FULL_YEARS[0]}-01", f"{FULL_YEARS[-1]}-12",
                             freq="M").strftime("%Y-%m")
    NM = len(months)
    mi_of = {m: i for i, m in enumerate(months)}
    Bi_m = np.full((NM, H, W), np.nan, dtype="float32")
    n_m = np.zeros((NM, H, W), dtype="int16")
    g = panel.groupby("ym")
    for m, sub in g:
        if m not in mi_of:
            continue
        i = mi_of[m]
        means = sub.groupby("pix")["Bi"].mean()
        idx = means.index.values
        Bi_m[i, yy_s[idx], xx_s[idx]] = means.values.astype("float32")
        cnts = sub.groupby("pix").size()
        n_m[i, yy_s[cnts.index.values], xx_s[cnts.index.values]] = \
            cnts.values.astype("int16")
    total_months = sum(1 for m in months
                       if np.isfinite(Bi_m[mi_of[m]][mask_np]).any())
    log(f"months with data: {total_months}/{NM}")

    ds_out = xr.Dataset(
        {"Bi_monthly": (("time", "y", "x"), Bi_m,
                        {"long_name": "monthly-mean ET residual Bi = ET_obs - ET_hat "
                                      "(candidate buffering signal)",
                         "units": "mm/dekad", "grid_mapping": "crs"}),
         "n_dekads": (("time", "y", "x"), n_m,
                      {"long_name": "dekads averaged in the month"})},
        coords={"time": pd.to_datetime([m + "-01" for m in months]),
                "y": GY, "x": GX},
        attrs={"Conventions": "CF-1.8",
               "title": "Monthly ET buffering signal, Ogallala ROI, 2016-2021",
               "history": "buffering_monthly_maps_timeseries.py; sources read-only"})
    ds_out["crs"] = xr.DataArray(np.int32(4326), attrs={
        "grid_mapping_name": "latitude_longitude",
        "longitude_of_prime_meridian": 0.0, "semi_major_axis": 6378137.0,
        "inverse_flattening": 298.257223563, "epsg_code": "EPSG:4326"})
    ds_out.to_netcdf(OUT / "Bi_monthly_4km.nc", engine="netcdf4",
                     encoding={"Bi_monthly": {"dtype": "float32", "zlib": True,
                                              "complevel": 3, "_FillValue":
                                              np.float32(9.96921e36)}})
    log("saved Bi_monthly_4km.nc")

    # ---------------- 5. ROI-averaged time series ----------------
    log("== ROI time series ==")
    ts = (panel.groupby("ym")
          .agg(Bi_mean=("Bi", "mean"), Bi_std=("Bi", "std"),
               n=("Bi", "size"),
               dry_frac=("dry", "mean"))
          .reindex(months))
    ts.index.name = "ym"
    ts.to_csv(OUT / "roi_Bi_monthly.csv")
    log(f"saved roi_Bi_monthly.csv ({len(ts)} months)")
    log("ROI Bi monthly: mean "
        f"{ts['Bi_mean'].mean():+.3f}, min {ts['Bi_mean'].min():+.2f} "
        f"({ts['Bi_mean'].idxmin()}), max {ts['Bi_mean'].max():+.2f} "
        f"({ts['Bi_mean'].idxmax()}) mm/dekad")

    t_axis = pd.to_datetime(ts.index)
    fig, ax1 = plt.subplots(figsize=(11, 4.2))
    ax1.plot(t_axis, ts["Bi_mean"].values, color="#2166ac", lw=1.2,
             label="ROI-mean Bi")
    ax1.fill_between(t_axis,
                     (ts["Bi_mean"] - ts["Bi_std"]).values,
                     (ts["Bi_mean"] + ts["Bi_std"]).values,
                     color="#2166ac", alpha=0.15, label="+/-1 SD (pixels)")
    ax1.axhline(0, color="k", lw=0.8)
    ax1.set_ylabel("ROI-mean Bi [mm/dekad]", color="#2166ac")
    ax1.set_title("ROI-averaged ET buffering signal (monthly), Ogallala 2016-2021\n"
                  "Bi = SSEBop ET - XGBoost expected ET; + = persisting ET")
    ax2 = ax1.twinx()
    ax2.bar(t_axis, ts["dry_frac"].values, width=22, color="#d73027",
            alpha=0.35, label="dry-down fraction")
    ax2.set_ylabel("dry-down pixel fraction", color="#d73027")
    ax2.set_ylim(0, 1)
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "roi_Bi_timeseries.png", dpi=170)
    plt.close(fig)
    log("saved roi_Bi_timeseries.png")

    # ---------------- 6. requested monthly map ----------------
    log("== requested map ==")
    key = f"{year}-{month:02d}"
    arr = Bi_m[mi_of[key]]
    v = arr[mask_np & np.isfinite(arr)]
    da = xr.DataArray(arr, coords={"y": GY, "x": GX}, dims=("y", "x"))
    x = GX
    y = GY
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    extent = [x[0] - dx / 2, x[-1] + dx / 2, y[-1] + dy / 2, y[0] - dy / 2]
    fig, ax = plt.subplots(figsize=(7.2, 6.8))
    img = np.ma.masked_invalid(np.asarray(da.values, dtype="float64"))
    vm = float(max(abs(np.nanpercentile(v, 2)), abs(np.nanpercentile(v, 98)), 1))
    im = ax.imshow(img, extent=extent, origin="upper", cmap="BrBG",
                   vmin=-vm, vmax=vm, interpolation="nearest")
    ax.set_title(f"ET buffering signal Bi, {key} (monthly mean)\n"
                 f"positive (green) = ET persisting beyond climate+soil "
                 f"explanation [mm/dekad]", fontsize=9)
    ax.set_xlabel("Longitude (deg E)", fontsize=8)
    ax.set_ylabel("Latitude (deg N)", fontsize=8)
    bnd.boundary.plot(ax=ax, edgecolor="black", linewidth=0.8)
    ax.grid(True, which="major", color="gray", alpha=0.35,
            linewidth=0.5, linestyle=":")
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=7)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label("Bi [mm/dekad]", fontsize=8)
    ndek = int(n_m[mi_of[key]][mask_np].max())
    ax.text(0.01, 0.01, f"valid px: {v.size:,} | dekads in month: {ndek}",
            transform=ax.transAxes, fontsize=7,
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))
    fig.tight_layout()
    fig.savefig(OUT / f"map_Bi_{key}.png", dpi=170)
    plt.close(fig)
    log(f"saved map_Bi_{key}.png (median {np.median(v):+.2f}, "
        f"frac>0 {np.mean(v > 0):.1%})")

    sseb.close()
    grid.close()
    smap.close()
    dr.close()
    ds_ref.close()
    with open(OUT / "run_log_monthly.txt", "w") as f:
        f.write("\n".join(LOG_LINES))
    log("DONE")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, default=2020)
    ap.add_argument("--month", type=int, default=8)
    a = ap.parse_args()
    main(a.year, a.month)
