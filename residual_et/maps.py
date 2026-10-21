"""Maps of the test-phase residual ET over the High Plains aquifer.

    python -m residual_et.maps --profile full --model xgb

Only out-of-fold residuals are drawn: Bi = observed - expected ET anomaly on
dekads the model was tested on, never on dekads it was fitted to. For each
year the map shows the growing-season (May-September) mean of Bi per cell.
One PNG per year plus one panel of all years, in <store>/<model>/maps/.
"""
import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from . import config as C
from . import grid as G
from .store import Store

MODEL_NAME = {"xgb": "XGBoost", "dlstm": "LSTM-HBV"}
MIN_DEKADS = 3
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 0.6})


def yearly_means(bi, tested, starts):
    """{year: (mean Bi per cell, number of test dekads)} for the growing season."""
    year = starts.astype("datetime64[Y]").astype(int) + 1970
    month = (starts.astype("datetime64[M]") - starts.astype("datetime64[Y]").astype("datetime64[M]")).astype(int) + 1
    grow = np.isin(month, C.ANALYSIS_MONTHS) & (tested >= 0)
    out = {}
    for y in np.unique(year[grow]):
        sel = grow & (year == y)
        blk = bi[sel]
        n = np.isfinite(blk).sum(axis=0)
        with np.errstate(invalid="ignore"):
            mean = np.where(n >= MIN_DEKADS, np.nansum(blk, axis=0) / np.maximum(n, 1), np.nan)
        if np.isfinite(mean).any():
            out[int(y)] = (mean, int(sel.sum()))
    return out


def _boundary():
    import geopandas as gpd
    shp = C.ROOT / "HPA_polygon" / "hp_bound2010.shp"
    return gpd.read_file(shp).to_crs("EPSG:4326").dissolve()


def _draw(ax, grid, values, bnd, norm, bbox):
    le, lo = grid.lat_edges, grid.lon_edges
    img = np.ma.masked_invalid(grid.to_map(values))
    im = ax.imshow(img, extent=[lo[0], lo[-1], le[-1], le[0]], origin="upper", cmap="BrBG", norm=norm,
                   interpolation="nearest")
    bnd.boundary.plot(ax=ax, color="0.15", linewidth=0.7)
    ax.set_xlim(bbox[0], bbox[1])
    ax.set_ylim(bbox[2], bbox[3])
    ax.set_aspect(1.0 / np.cos(np.deg2rad(np.mean(bbox[2:]))))
    ax.set_facecolor("white")
    ax.grid(True, color="0.8", linewidth=0.4, linestyle=":")
    ax.set_axisbelow(True)
    ax.tick_params(length=2.5, width=0.5, labelsize=7)
    ax.xaxis.set_major_formatter(lambda v, _: f"{abs(v):g}°W")
    ax.yaxis.set_major_formatter(lambda v, _: f"{v:g}°N")
    return im


def run(profile_name, model):
    profile = C.get_profile(profile_name)
    store, res = Store(profile.store), Store(profile.store / model)
    grid = G.load_grid(profile)
    means = yearly_means(np.asarray(res.arr("Bi")), np.asarray(res.arr("oof_fold")), store.dekad_start)
    if not means:
        raise RuntimeError("no year has three test dekads in the growing season")
    out = res.path / "maps"
    out.mkdir(parents=True, exist_ok=True)
    bnd = _boundary()
    pad = 0.3
    bbox = ((grid.cell_lon.min() - pad, grid.cell_lon.max() + pad, grid.cell_lat.min() - pad,
             grid.cell_lat.max() + pad) if profile.bbox else
            tuple(np.array(bnd.total_bounds)[[0, 2, 1, 3]] + [-pad, pad, -pad, pad]))
    allv = np.concatenate([m[np.isfinite(m)] for m, _ in means.values()])
    lim = float(np.ceil(np.percentile(np.abs(allv), 98)))            # one scale for every year
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
    label = "Residual ET, observed − expected (mm per dekad)"
    note = (f"{MODEL_NAME[model]}; test-phase (out-of-fold) residuals only; May–September mean per "
            f"~4 km cell; cells with fewer than {MIN_DEKADS} test dekads are blank.")
    files = []
    for y, (mean, nd) in means.items():
        fig, ax = plt.subplots(figsize=(5.2, 7.4))
        im = _draw(ax, grid, mean, bnd, norm, bbox)
        ax.set_title(f"Residual ET over the High Plains aquifer, {y}", fontsize=11, fontweight="bold", loc="left")
        v = mean[np.isfinite(mean)]
        ax.text(0.03, 0.03, f"{nd} test dekads\nmedian {np.median(v):+.1f}\n{100 * np.mean(v > 0):.0f}% of cells positive",
                transform=ax.transAxes, fontsize=7.5, va="bottom",
                bbox=dict(facecolor="white", edgecolor="0.6", linewidth=0.5, boxstyle="round,pad=0.3"))
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.06, fraction=0.035, extend="both")
        cb.set_label(label, fontsize=8)
        cb.ax.tick_params(labelsize=7, length=2)
        fig.text(0.5, 0.012, note, ha="center", fontsize=6.3, color="0.3", wrap=True)
        p = out / f"residual_ET_{model}_{y}.png"
        fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        files.append(p)
    n = len(means)
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.9 * ncol, 4.3 * nrow + 0.9), squeeze=False)
    for ax, (y, (mean, nd)) in zip(axes.ravel(), means.items()):
        im = _draw(ax, grid, mean, bnd, norm, bbox)
        ax.set_title(f"{y}  ({nd} dekads)", fontsize=9, fontweight="bold")
        ax.label_outer()
    for ax in axes.ravel()[n:]:
        ax.set_visible(False)
    fig.suptitle(f"Growing-season residual ET over the High Plains aquifer — {MODEL_NAME[model]}",
                 fontsize=12, fontweight="bold", y=0.995)
    cb = fig.colorbar(im, ax=axes, orientation="horizontal", pad=0.05, fraction=0.03, aspect=45, extend="both")
    cb.set_label(label, fontsize=9)
    fig.text(0.5, 0.005, note, ha="center", fontsize=7, color="0.3")
    p = out / f"residual_ET_{model}_all_years.png"
    fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    files.append(p)
    print(f"{model}: years {list(means)}, colour scale ±{lim:.0f} mm/dekad, {len(files)} maps in {out}", flush=True)
    return files


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Maps of test-phase residual ET.")
    ap.add_argument("--profile", default="full", choices=sorted(C.PROFILES))
    ap.add_argument("--model", default="xgb", choices=["xgb", "dlstm"])
    a = ap.parse_args()
    run(a.profile, a.model)
