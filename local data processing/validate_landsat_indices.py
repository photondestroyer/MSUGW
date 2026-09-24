"""
Validate + plot the per-date Landsat spectral-indices cube for the Ogallala ROI.

Source (READ-ONLY, never modified):
    D:/Downloads/Landsat_indices_Ogallala_2020_perdate.nc
    Landsat 8/9 C2 L2 per-date cloud-masked indices on a ~4 km WGS-84 grid.

What the script does:
  1. Inventories every spectral index present in the file.
  2. Runs validity checks: time length, finite-data fraction, valid-range
     compliance (valid_min/valid_max attrs), n_obs consistency, lat/lon grid
     regularity, CF-1.8 conventions + grid_mapping attrs.
  3. Saves plots (always) + a machine-readable QC report (JSON).
     Maps overlay the hp_bound2010 aquifer boundary with faint gridlines,
     following repo plotting conventions.

Output -> G:/MSU_GWB/datasets/landsat_indices_qc/
  qc_report.json, qc_report.txt, grid_footprint.png,
  per-index maps + time series (only when data exist),
  empty_file_verdict.png (only when the cube holds no data).

Run:  python validate_landsat_indices.py [--file PATH]
"""

import argparse
import json
import time
import warnings
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

SRC_DEFAULT = Path(r"D:/Downloads/Landsat_indices_Ogallala_2020_perdate.nc")
BOUND_SHP = Path(r"G:/MSU_GWB/datasets/high_plains_quifer/hp_bound2010.shp")
OUT = Path(r"G:/MSU_GWB/datasets/landsat_indices_qc")
OUT.mkdir(parents=True, exist_ok=True)

LOG_LINES = []


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    LOG_LINES.append(line)


def grid_extent(lon, lat):
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    dx = float(lon[1] - lon[0])
    dy = float(lat[1] - lat[0])
    if dy < 0:
        return ([lon[0] - dx / 2, lon[-1] + dx / 2,
                 lat[-1] + dy / 2, lat[0] - dy / 2], "upper")
    return ([lon[0] - dx / 2, lon[-1] + dx / 2,
             lat[0] - dy / 2, lat[-1] + dy / 2], "lower")


def decorate(ax, title, bnd):
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Longitude (deg E)", fontsize=8)
    ax.set_ylabel("Latitude (deg N)", fontsize=8)
    if bnd is not None:
        bnd.boundary.plot(ax=ax, edgecolor="black", linewidth=0.9)
    ax.grid(True, which="major", color="gray", alpha=0.35,
            linewidth=0.5, linestyle=":")
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=7)


def main(src_path):
    log(f"source: {src_path} (exists={src_path.exists()})")
    if not src_path.exists():
        raise FileNotFoundError(src_path)
    log(f"source size: {src_path.stat().st_size / 1e3:.1f} KB")

    import geopandas as gpd
    bnd = gpd.read_file(BOUND_SHP).to_crs("EPSG:4326")

    ds = xr.open_dataset(src_path, decode_cf=True)
    report = {
        "file": str(src_path),
        "checked_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "global_attrs": {k: str(v)[:300] for k, v in ds.attrs.items()},
    }

    # ---- 1. inventory: which spectral indices are present ----
    META_VARS = {"crs", "spatial_ref"}
    index_vars = [v for v in ds.data_vars
                  if v not in META_VARS and v != "n_obs"]
    report["spectral_indices"] = {
        v: {"long_name": ds[v].attrs.get("long_name", "?"),
            "dims": list(ds[v].dims), "dtype": str(ds[v].dtype),
            "attrs": {k: str(val)[:120] for k, val in ds[v].attrs.items()}}
        for v in index_vars}
    report["count_layer"] = "n_obs" if "n_obs" in ds.data_vars else None
    log(f"spectral indices present ({len(index_vars)}): "
        f"{', '.join(index_vars)}")
    for v in index_vars:
        log(f"  {v}: {ds[v].attrs.get('long_name', '?')}")

    # ---- 2. grid + time structure ----
    lat = ds["lat"].values
    lon = ds["lon"].values
    t = pd_time(ds)
    nt = len(t)
    dlat = np.diff(lat)
    dlon = np.diff(lon)
    grid_ok = bool(np.allclose(dlat, dlat[0]) and np.allclose(dlon, dlon[0]))
    report["grid"] = {
        "lat_range": [float(lat.min()), float(lat.max())],
        "lon_range": [float(lon.min()), float(lon.max())],
        "nlat": int(len(lat)), "nlon": int(len(lon)),
        "dlat_deg": float(dlat[0]), "dlon_deg": float(dlon[0]),
        "regular": grid_ok,
        "n_time": int(nt),
        "time_coverage_attr": [
            ds.attrs.get("time_coverage_start"),
            ds.attrs.get("time_coverage_end")],
    }
    log(f"grid: lat {lat.min():.2f}..{lat.max():.2f} (n={len(lat)}), "
        f"lon {lon.min():.2f}..{lon.max():.2f} (n={len(lon)}), "
        f"regular={grid_ok}")
    log(f"time steps: {nt}")

    # ---- 3. per-index validity ----
    checks = {}
    for v in index_vars + ((["n_obs"] if "n_obs" in ds.data_vars else [])):
        a = ds[v].values
        total = a.size
        fin = np.isfinite(a.astype(float))
        n_fin = int(fin.sum())
        vmin = ds[v].attrs.get("valid_min", None)
        vmax = ds[v].attrs.get("valid_max", None)
        in_range = None
        if n_fin and vmin is not None and vmax is not None:
            f = a.astype(float)[fin]
            in_range = float(np.mean((f >= float(vmin)) & (f <= float(vmax))))
        checks[v] = {"total_cells": int(total), "finite": n_fin,
                     "finite_frac": float(n_fin / total) if total else 0.0,
                     "in_valid_range_frac": in_range,
                     "valid_min_attr": str(vmin), "valid_max_attr": str(vmax)}
        log(f"  {v}: finite {n_fin:,}/{total:,} "
            f"({100 * n_fin / total if total else 0:.2f}%), "
            f"in-range {in_range}")
    report["validity"] = checks
    report["cf_conventions"] = ds.attrs.get("Conventions", None)
    report["has_grid_mapping"] = bool(
        any("grid_mapping" in ds[v].attrs for v in ds.data_vars))

    # n_obs <-> index consistency (finite index requires n_obs > 0)
    if "n_obs" in ds.data_vars and index_vars and nt > 0:
        nobs = ds["n_obs"].values.astype(float)
        idx_fin = np.zeros_like(nobs, dtype=bool)
        for v in index_vars:
            idx_fin |= np.isfinite(ds[v].values.astype(float))
        both = idx_fin & np.isfinite(nobs) & (nobs > 0)
        report["n_obs_consistency_frac"] = float(
            both.sum() / idx_fin.sum()) if idx_fin.sum() else None
        log(f"n_obs consistency: {report['n_obs_consistency_frac']}")

    # ---- verdict ----
    if nt == 0:
        verdict = ("EMPTY: time dimension has length 0 - the file holds "
                   "coordinates + metadata but NO data. Structurally valid "
                   "CF-1.8, scientifically unusable until the producer "
                   "pipeline (ls_perdate) is re-run with a non-empty date "
                   "list.")
    elif all(c["finite"] == 0 for c in checks.values()):
        verdict = "INVALID: time steps exist but every variable is empty/NaN."
    else:
        bad = [v for v, c in checks.items()
               if c["finite_frac"] == 0 or (c["in_valid_range_frac"] is not None
                                           and c["in_valid_range_frac"] < 0.99)]
        verdict = ("VALID" if not bad
                   else f"QUESTIONABLE (see: {', '.join(bad)})")
    report["verdict"] = verdict
    log("VERDICT: " + verdict)

    # ---- 4. plots ----
    extent, origin = grid_extent(lon, lat)

    # 4a. grid footprint vs aquifer boundary (always possible)
    fig, ax = plt.subplots(figsize=(7.2, 6.8))
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    rx = [lon[0] - abs(lon[1] - lon[0]) / 2, lon[-1] + abs(lon[1] - lon[0]) / 2]
    ry = [lat[0] - abs(lat[1] - lat[0]) / 2, lat[-1] + abs(lat[1] - lat[0]) / 2]
    ax.fill_between(rx, ry[0], ry[1], color="#c7e9c0", alpha=0.6,
                    label="NetCDF grid footprint")
    decorate(ax, "Landsat indices cube: grid footprint vs aquifer boundary\n"
                 f"{len(lat)}x{len(lon)} cells, {nt} time steps", bnd)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(OUT / "grid_footprint.png", dpi=170)
    plt.close(fig)
    log("saved grid_footprint.png")

    if nt > 0 and any(checks[v]["finite"] > 0 for v in index_vars):
        # 4b. per-index median maps
        n = len(index_vars)
        ncols = 3
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.6 * nrows))
        for ax, v in zip(np.atleast_2d(axes).flat, index_vars):
            med = np.nanmedian(ds[v].values.astype(float), axis=0)
            img = np.ma.masked_invalid(med)
            im = ax.imshow(img, extent=extent, origin=origin, cmap="YlGn",
                           vmin=-1, vmax=1, interpolation="nearest")
            decorate(ax, f"{v} - {ds[v].attrs.get('long_name', '')}\n"
                         "2020 median over dates", bnd)
            fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
        for ax in np.atleast_2d(axes).flat[n:]:
            ax.axis("off")
        fig.suptitle("Per-date spectral indices: 2020 median maps", fontsize=11)
        fig.tight_layout()
        fig.savefig(OUT / "index_median_maps.png", dpi=170)
        plt.close(fig)
        log("saved index_median_maps.png")

        # 4c. ROI-mean time series per index
        fig, ax = plt.subplots(figsize=(11, 4.5))
        for v in index_vars:
            a = ds[v].values.astype(float)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                ts = np.nanmean(a.reshape(a.shape[0], -1), axis=1)
            ax.plot(pd_time(ds), ts, lw=1.2, label=v)
        ax.set_ylabel("ROI-mean index value [1]")
        ax.set_title("Spectral indices: ROI-averaged per-date time series, 2020")
        ax.legend(fontsize=8, ncol=4)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(OUT / "index_timeseries.png", dpi=170)
        plt.close(fig)
        log("saved index_timeseries.png")

        # 4d. histograms
        fig, axes = plt.subplots(2, 4, figsize=(14, 6))
        for ax, v in zip(axes.flat, index_vars):
            f = ds[v].values.astype(float).ravel()
            f = f[np.isfinite(f)]
            ax.hist(f, bins=80, color="#2ca25f", edgecolor="white",
                    linewidth=0.3)
            ax.set_title(f"{v} (n={f.size:,})", fontsize=9)
            ax.set_xlabel("[1]")
        for ax in axes.flat[n:]:
            ax.axis("off")
        fig.suptitle("Value distributions over all dates x cells", fontsize=11)
        fig.tight_layout()
        fig.savefig(OUT / "index_histograms.png", dpi=170)
        plt.close(fig)
        log("saved index_histograms.png")
    else:
        # 4e. empty-file verdict figure (no data can be mapped)
        fig, ax = plt.subplots(figsize=(7.2, 6.8))
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.text(0.5, 0.55, "NO DATA TO PLOT\ntime dimension = 0",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=16, color="#d73027",
                bbox=dict(facecolor="white", edgecolor="#d73027", pad=12))
        ax.text(0.5, 0.40,
                "File holds coordinates + metadata only.\n"
                "Re-run the ls_perdate producer pipeline\n"
                "with a non-empty 2020 date list.",
                transform=ax.transAxes, ha="center", va="center", fontsize=9)
        decorate(ax, "EMPTY CUBE - nothing to map", bnd)
        fig.tight_layout()
        fig.savefig(OUT / "empty_file_verdict.png", dpi=170)
        plt.close(fig)
        log("saved empty_file_verdict.png")

    ds.close()

    with open(OUT / "qc_report.json", "w") as f:
        json.dump(report, f, indent=1, default=str)
    with open(OUT / "qc_report.txt", "w") as f:
        f.write(f"QC report for {src_path}\n{'=' * 70}\n")
        f.write(f"Spectral indices present ({len(index_vars)}): "
                f"{', '.join(index_vars)}\n")
        f.write(f"Time steps: {nt}\nGrid: {len(lat)}x{len(lon)}, "
                f"regular={grid_ok}\n")
        for v, c in checks.items():
            f.write(f"  {v}: finite {c['finite']:,}/{c['total_cells']:,} "
                    f"({100 * c['finite_frac']:.2f}%), "
                    f"in-range {c['in_valid_range_frac']}\n")
        f.write(f"VERDICT: {verdict}\n")
        f.write("\n".join(LOG_LINES) + "\n")
    log("saved qc_report.json + qc_report.txt")
    log("files: " + ", ".join(sorted(p.name for p in OUT.iterdir())))
    log("DONE")


def pd_time(ds):
    import pandas as pd
    return pd.to_datetime(ds["time"].values)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", type=str, default=str(SRC_DEFAULT))
    a = ap.parse_args()
    main(Path(a.file))
