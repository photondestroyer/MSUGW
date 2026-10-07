"""
Shallow and deep groundwater clusters of the High Plains aquifer for every year from 2000 on.

Same method as plot_gw_shallow_deep_clusters.py (Local Moran statistic, Anselin 1995, Geographical Analysis 27,
93-115), repeated per year on the whole-record water-level file HPA_wells_water_level_timeseries.nc.

Per year
  * measurements of the pre-irrigation season, December (of the year before) to April   [USGS's own convention]
  * records flagged as pumping, recently pumped, dry, obstructed, flowing, destroyed ... are left out, as are
    records that repeat another source and Texas records that are not 'Publishable'
  * one value per well (group): the median of its winter measurements; all wells measured that year are used
  * cell value = median of its wells (at least 2), cells without wells stay empty
  * Local Moran on standardised log10 depth, queen contiguity, conditional permutation, FDR 0.05

Because the set of wells changes from year to year, the tables also report every year's well count and observed
area, and a second depth series restricted to the cells observed in every year.
"""
import re
import warnings
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
import pandas as pd
import shapely
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
from matplotlib.patches import Patch
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = Path(r"G:/MSU_GWB/datasets")
TS = ROOT / "merged_datasets" / "HPA_wells_water_level_timeseries.nc"
GRIDMET = ROOT / "merged_datasets" / "GRIDMET_Ogallala_lab_drive.nc"
BOUND = "zip://G:/MSU_GWB/datasets/high_plains_quifer.zip!high_plains_quifer/hp_bound2010.shp"
OUT = ROOT / "plots" / "gw_shallow_deep_clusters_yearly"
OUT.mkdir(parents=True, exist_ok=True)

DLAT, DLON = 0.25, 0.30
MIN_WELLS = 2
N_PERM = 2999
ALPHA = 0.05
FIRST_YEAR = 2000
MIN_WELLS_YEAR = 1000              # a year with fewer wells than this is reported but not mapped
WINTER = (12, 1, 2, 3, 4)
CLIM_YEARS = (1981, 2010)
FT = 0.3048
BAD = re.compile(r"pump|dry|obstruct|flow|destroy|plug|affect|nearby|foreign|recent|injure|caved|questionable|"
                 r"not.*static|unreliable", re.I)
RNG = np.random.default_rng(20261007)
C = dict(HH="#b2182b", LL="#2166ac", HL="#f4a582", LH="#92c5de", ns="#d9d9d9")


def log(m):
    print(f"[gw-yearly] {m}", flush=True)


# ════════════════════ aquifer and analysis grid (identical to the single-year script) ════════════════════
roi = gpd.read_file(BOUND).to_crs(4326).geometry.union_all()
roi_ea = gpd.GeoSeries([roi], crs=4326).to_crs(5070)
AQ_AREA = float(roi_ea.area.iloc[0]) / 1e6
W, S, E, N = roi.bounds
lon_edges = np.arange(np.floor(W / DLON) * DLON, E + DLON, DLON)
lat_edges = np.arange(np.floor(S / DLAT) * DLAT, N + DLAT, DLAT)
NY, NX = lat_edges.size - 1, lon_edges.size - 1
cells = gpd.GeoDataFrame(dict(iy=np.repeat(np.arange(NY), NX), ix=np.tile(np.arange(NX), NY)),
                         geometry=[shapely.box(lon_edges[j], lat_edges[i], lon_edges[j + 1], lat_edges[i + 1])
                                   for i in range(NY) for j in range(NX)], crs=4326)
cells["aq_km2"] = cells.to_crs(5070).geometry.intersection(roi_ea.iloc[0]).area.values / 1e6
cells = cells[cells.aq_km2 > 0].reset_index(drop=True)
cells["lat"] = lat_edges[cells.iy] + DLAT / 2
cells["lon"] = lon_edges[cells.ix] + DLON / 2
NC = len(cells)
key = {(a, b): k for k, (a, b) in enumerate(zip(cells.iy, cells.ix))}


def to_cell_index(lat, lon):
    iy = np.floor((np.asarray(lat) - lat_edges[0]) / DLAT).astype(int)
    ix = np.floor((np.asarray(lon) - lon_edges[0]) / DLON).astype(int)
    return np.array([key.get((a, b), -1) for a, b in zip(iy, ix)])


# ════════════════════ measurements ════════════════════
with nc4.Dataset(TS) as ds:
    ds.set_auto_maskandscale(False)
    st = pd.DataFrame(dict(lat=ds["lat"][:], lon=ds["lon"][:], group=ds["well_group"][:], inpoly=ds["in_polygon"][:],
                           source=np.array([str(v) for v in ds["source"][:]]),
                           state=np.array([str(v) for v in ds["state"][:]])))
    look = lambda n: np.array([str(v) for v in ds[n + "_text"][:]])[ds[n][:]]
    ob = pd.DataFrame(dict(si=ds["station_index"][:], t=ds["time"][:], d=ds["depth_to_water"][:],
                           year=ds["year"][:].astype(int), dup=ds["duplicate_of_other_source"][:],
                           status=look("status"), src=look("record_source")))
dated = ob.t > -9000
ob["date"] = pd.Timestamp("1900-01-01") + pd.to_timedelta(ob.t.where(dated), unit="D")
ob["month"] = ob.date.dt.month
# water year: December belongs to the following year. Records with a year only (most of Nebraska's) count
# when the source labels them as a spring measurement.
ob["wy"] = np.where(dated, ob.year + (ob.month == 12), ob.year)
ob["winter"] = np.where(dated, ob.month.isin(WINTER), ob.status.str.lower().eq("spring"))
n_all = len(ob)
ob = ob[(ob.wy >= FIRST_YEAR) & (ob.wy <= 2100) & ob.winter]
n_win = len(ob)
bad_status = ob.status.str.contains(BAD) | ((ob.src == "TWDB") & ~ob.status.str.startswith("Publishable"))
keep = (ob.dup == 0) & ~bad_status & (ob.d > -50) & (ob.d < 2000) & (st.inpoly.values[ob.si.values] == 1)
log(f"records {n_all:,}; winter records since {FIRST_YEAR}: {n_win:,}; dropped: repeats of another source "
    f"{int((ob.dup == 1).sum()):,}, flagged status {int(bad_status.sum()):,}, outside polygon "
    f"{int((st.inpoly.values[ob.si.values] == 0).sum()):,}; used {int(keep.sum()):,}")
ob = ob[keep].copy()
ob["group"] = st.group.values[ob.si.values]
ob["dm"] = ob.d * FT
gpos = st.groupby("group")[["lat", "lon"]].first()
wy = ob.groupby(["wy", "group"]).agg(dm=("dm", "median"), src=("src", "first")).reset_index()
wy["cell"] = to_cell_index(gpos.lat.reindex(wy.group).values, gpos.lon.reindex(wy.group).values)
wy = wy[wy.cell >= 0]
wy["state"] = st.groupby("group").state.first().reindex(wy.group).values
YEARS_ALL = sorted(wy.wy.unique())
counts = wy.groupby("wy").size()
ncell = wy.groupby("wy").cell.nunique()
# a year is mapped only if its network is comparable: enough wells, and at least 75 % of the typical number of cells
YEARS = [y for y in YEARS_ALL if counts[y] >= MIN_WELLS_YEAR and ncell[y] >= 0.75 * ncell.median()]
log(f"wells per year: {counts.to_dict()}")
log(f"years analysed: {YEARS[0]} to {YEARS[-1]} ({len(YEARS)}); skipped (network too small that year): "
    f"{[y for y in YEARS_ALL if y not in YEARS]}")
pd.crosstab(wy.wy, wy.state).assign(total=counts).to_csv(OUT / "table_wells_per_year_by_state.csv")
pd.crosstab(wy.wy, wy.src).to_csv(OUT / "table_wells_per_year_by_source.csv")


# ════════════════════ Local Moran for one year ════════════════════
def lisa(depth):
    """depth: cell medians (NaN = no data). Returns cluster labels, global I, its pseudo p, number of tests."""
    obs = np.isfinite(depth)
    idx = np.flatnonzero(obs)
    pos = {(a, b): k for k, (a, b) in enumerate(zip(cells.iy.values[idx], cells.ix.values[idx]))}
    nbrs = [[pos[(a + da, b + db)] for da in (-1, 0, 1) for db in (-1, 0, 1) if (da or db) and (a + da, b + db) in pos]
            for a, b in zip(cells.iy.values[idx], cells.ix.values[idx])]
    has = np.array([len(q) > 0 for q in nbrs])
    x = np.log10(np.clip(depth[idx], 0.3, None))
    z = (x - x[has].mean()) / x[has].std()
    n = len(z)
    lag = np.array([z[q].mean() if q else np.nan for q in nbrs])
    I = float(np.nansum(z[has] * lag[has]) / np.sum(z[has] ** 2))
    p = np.full(n, np.nan)
    for i in np.flatnonzero(has):
        k = len(nbrs[i])
        others = np.delete(z, i)
        pick = np.argpartition(RNG.random((N_PERM, n - 1)), k - 1, axis=1)[:, :k]
        sim = others[pick].mean(axis=1)
        larger = int((sim >= lag[i]).sum())
        p[i] = (min(larger, N_PERM - larger) + 1) / (N_PERM + 1)
    hz = np.flatnonzero(has)
    W_ = np.zeros((hz.size, n))
    for r, i in enumerate(hz):
        W_[r, nbrs[i]] = 1 / len(nbrs[i])
    sims = np.empty(499)
    for r in range(499):
        zp = z.copy()
        zp[hz] = RNG.permutation(z[hz])
        sims[r] = np.sum(zp[hz] * (W_ @ zp)) / np.sum(zp[hz] ** 2)
    pI = (int((sims >= I).sum()) + 1) / 500
    pv = p[has]
    order = np.argsort(pv)
    passed = pv[order] <= ALPHA * np.arange(1, pv.size + 1) / pv.size
    cut = pv[order][np.flatnonzero(passed).max()] if passed.any() else -1.0
    sig = has & (p <= cut)
    quad = np.where(z >= 0, np.where(lag >= 0, "HH", "HL"), np.where(lag >= 0, "LH", "LL"))
    lab = np.full(NC, "no data", dtype=object)
    lab[idx] = np.where(sig, quad, "ns")
    return lab, I, pI, int(has.sum())


# ════════════════════ climate normals per cell (gridMET 1981-2010) ════════════════════
cache = OUT / "_cache_gridmet_normals.npz"
if cache.exists():
    cn = dict(np.load(cache))
else:
    MD = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    with nc4.Dataset(GRIDMET) as ds:
        glat, glon = ds["lat"][:].data, ds["lon"][:].data
        t = pd.to_datetime("1900-01-01") + pd.to_timedelta(ds["time"][:].data, unit="D")
        inside = ds["roi_mask"][:].data == 1
        sel = np.flatnonzero((t.year >= CLIM_YEARS[0]) & (t.year <= CLIM_YEARS[1]))
        mean = {}
        for v in ("pr", "pet", "tmmx", "tmmn", "vpd"):
            tot = np.zeros(inside.shape)
            cnt = np.zeros(inside.shape)
            for y in range(CLIM_YEARS[0], CLIM_YEARS[1] + 1):
                k = sel[t.year.values[sel] == y]
                a = np.ma.filled(ds[v][k[0]:k[-1] + 1].astype("f8"), np.nan)
                tot += np.nansum(a, axis=0)
                cnt += np.isfinite(a).sum(axis=0)
            mean[v] = tot / np.maximum(cnt, 1)
            log(f"gridMET {v} normal done")
    GLON, GLAT = np.meshgrid(glon, glat)
    gc = to_cell_index(GLAT[inside], GLON[inside])
    ok = gc >= 0

    def tc(f):
        s = np.bincount(gc[ok], weights=f[inside][ok], minlength=NC)
        c = np.bincount(gc[ok], minlength=NC)
        return np.where(c > 0, s / np.maximum(c, 1), np.nan)

    cn = dict(P_mm=tc(mean["pr"]) * 365.25, PET_mm=tc(mean["pet"]) * 365.25,
              T_C=tc((mean["tmmx"] + mean["tmmn"]) / 2) - 273.15, VPD_kPa=tc(mean["vpd"]))
    np.savez(cache, **cn)
cn["aridity"] = cn["P_mm"] / cn["PET_mm"]
area = cells.aq_km2.values


def wmean(v, m):
    k = m & np.isfinite(v)
    return float(np.sum(v[k] * area[k]) / np.sum(area[k])) if k.any() else np.nan


# ════════════════════ yearly loop ════════════════════
D = np.full((len(YEARS), NC), np.nan)            # cell median depth (m)
NW = np.zeros((len(YEARS), NC), int)
L = np.empty((len(YEARS), NC), dtype=object)
rows = []
for r, y in enumerate(YEARS):
    sub = wy[wy.wy == y]
    g = sub.groupby("cell").dm
    n = g.size()
    med = g.median()
    NW[r, n.index] = n.values
    good = n.index[n.values >= MIN_WELLS]
    D[r, good] = med[good].values
    lab, I, pI, ntest = lisa(D[r])
    L[r] = lab
    obs = np.isfinite(D[r])
    row = dict(year=y, wells=len(sub), cells_observed=int(obs.sum()), observed_area_km2=area[obs].sum(),
               observed_pct_of_aquifer=100 * area[obs].sum() / AQ_AREA, morans_I=I, morans_p=pI,
               median_depth_wells_m=sub.dm.median(), median_depth_cells_m=float(np.nanmedian(D[r])),
               wells_le_5m_pct=100 * (sub.dm <= 5).mean(), wells_gt_60m_pct=100 * (sub.dm > 60).mean(),
               cells_le_5m_pct_of_observed=100 * area[obs & (D[r] <= 5)].sum() / area[obs].sum(),
               cells_gt_60m_pct_of_observed=100 * area[obs & (D[r] > 60)].sum() / area[obs].sum())
    for c, nm in (("LL", "shallow"), ("HH", "deep"), ("HL", "deep_outlier"), ("LH", "shallow_outlier")):
        m = lab == c
        row[f"{nm}_cells"] = int(m.sum())
        row[f"{nm}_area_km2"] = area[m].sum()
        row[f"{nm}_pct_of_aquifer"] = 100 * area[m].sum() / AQ_AREA
        row[f"{nm}_pct_of_observed"] = 100 * area[m].sum() / area[obs].sum()
        if c in ("LL", "HH"):
            row[f"{nm}_median_depth_m"] = float(np.nanmedian(D[r][m])) if m.any() else np.nan
            row[f"{nm}_lat"] = wmean(cells.lat.values, m)
            row[f"{nm}_lon"] = wmean(cells.lon.values, m)
            for k in ("P_mm", "PET_mm", "aridity", "T_C", "VPD_kPa"):
                row[f"{nm}_{k}"] = wmean(cn[k], m)
    rows.append(row)
    log(f"{y}: wells {len(sub):5d} | observed {row['observed_pct_of_aquifer']:.0f} % | I = {I:.2f} | shallow "
        f"{row['shallow_pct_of_aquifer']:.1f} % | deep {row['deep_pct_of_aquifer']:.1f} % of aquifer")
tab = pd.DataFrame(rows)

# cells observed in every analysed year: a fixed footprint, to separate real change from network change
always = np.isfinite(D).all(axis=0)
tab["fixed_cells_median_depth_m"] = np.nanmedian(D[:, always], axis=1)
tab["fixed_cells_mean_depth_m"] = [wmean(D[r], always) for r in range(len(YEARS))]
for c, nm in (("LL", "shallow"), ("HH", "deep")):
    tab[f"fixed_cells_{nm}_pct"] = [100 * area[always & (L[r] == c)].sum() / area[always].sum() for r in range(len(YEARS))]
tab.to_csv(OUT / "table_yearly_summary.csv", index=False, float_format="%.5g")
log(f"cells observed in every year: {int(always.sum())} = {100 * area[always].sum() / AQ_AREA:.1f} % of the aquifer")

# persistence and trend per cell
yrs = np.array(YEARS, float)
n_obs = np.isfinite(D).sum(axis=0)
n_sh, n_dp = (L == "LL").sum(axis=0), (L == "HH").sum(axis=0)
slope = np.full(NC, np.nan)
ptrend = np.full(NC, np.nan)
for j in np.flatnonzero(n_obs >= max(10, int(0.6 * len(YEARS)))):
    k = np.isfinite(D[:, j])
    slope[j] = stats.theilslopes(D[k, j], yrs[k])[0]
    ptrend[j] = stats.kendalltau(yrs[k], D[k, j]).pvalue
cellt = cells.drop(columns="geometry").assign(years_observed=n_obs, years_shallow=n_sh, years_deep=n_dp,
                                               pct_years_shallow=100 * n_sh / np.maximum(n_obs, 1),
                                               pct_years_deep=100 * n_dp / np.maximum(n_obs, 1),
                                               depth_trend_m_per_yr=slope, trend_p=ptrend,
                                               observed_every_year=always.astype(int))
for r, y in enumerate(YEARS):
    cellt[f"depth_m_{y}"] = D[r]
    cellt[f"cluster_{y}"] = L[r]
cellt.to_csv(OUT / "table_cells_by_year.csv", index=False, float_format="%.5g")

# ════════════════════ figures ════════════════════
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "savefig.dpi": 200, "figure.dpi": 100})
ASPECT = 1 / np.cos(np.radians(37.7))
outline = gpd.GeoSeries([roi], crs=4326)
clipped = cells.copy()
clipped["geometry"] = cells.geometry.intersection(roi)


def base(ax, title=None, labels=True):
    outline.boundary.plot(ax=ax, color="k", linewidth=0.5)
    ax.set_xlim(W - 0.3, E + 0.3)
    ax.set_ylim(S - 0.3, N + 0.3)
    ax.set_aspect(ASPECT)
    if labels:
        ax.set_xlabel("Longitude (deg E)")
        ax.set_ylabel("Latitude (deg N)")
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    if title:
        ax.set_title(title, loc="left")


# ---- A: cluster map per year ----
ncol = 7
nrow = int(np.ceil(len(YEARS) / ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(2.35 * ncol, 3.9 * nrow), constrained_layout=True)
for ax in axs.flat:
    ax.axis("off")
for r, y in enumerate(YEARS):
    ax = axs.flat[r]
    ax.axis("on")
    for c in ("ns", "LL", "HH", "LH", "HL"):
        m = L[r] == c
        if m.any():
            clipped[m].plot(ax=ax, color=C[c], edgecolor="none")
    base(ax, labels=False)
    ax.set_title(f"{y}  n={tab.wells[r]:,}\nshallow {tab.shallow_pct_of_aquifer[r]:.0f}%  deep "
                 f"{tab.deep_pct_of_aquifer[r]:.0f}%", fontsize=8, loc="left")
fig.legend(handles=[Patch(facecolor=C["LL"], label="Shallow cluster"), Patch(facecolor=C["HH"], label="Deep cluster"),
                    Patch(facecolor=C["HL"], label="Deep outlier"), Patch(facecolor=C["LH"], label="Shallow outlier"),
                    Patch(facecolor=C["ns"], label="Not significant"),
                    Patch(facecolor="white", edgecolor="0.6", label="No wells")],
           loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.5, -0.03))
fig.suptitle("Shallow and deep groundwater clusters by year, High Plains aquifer (winter depth to water; Local Moran, "
             "FDR 0.05; percentages of aquifer area)", fontsize=11)
fig.savefig(OUT / "figA_cluster_maps_by_year.png", bbox_inches="tight")
plt.close(fig)

# ---- B: yearly metrics ----
fig, axs = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
ax = axs[0, 0]
ax.bar(tab.year, tab.wells, color="0.6")
ax.set_ylabel("Wells used")
ax2 = ax.twinx()
ax2.plot(tab.year, tab.observed_pct_of_aquifer, "k-o", ms=3)
ax2.set_ylabel("Observed area (% of aquifer)")
ax2.set_ylim(0, 100)
ax.set_title("(a) Wells and coverage", loc="left")
ax = axs[0, 1]
ax.plot(tab.year, tab.shallow_pct_of_aquifer, "-o", color=C["LL"], ms=3, label="Shallow, % of aquifer")
ax.plot(tab.year, tab.deep_pct_of_aquifer, "-o", color=C["HH"], ms=3, label="Deep, % of aquifer")
ax.plot(tab.year, tab.shallow_pct_of_observed, "--", color=C["LL"], label="Shallow, % of observed area")
ax.plot(tab.year, tab.deep_pct_of_observed, "--", color=C["HH"], label="Deep, % of observed area")
ax.set_ylabel("Cluster area (%)")
ax.legend(fontsize=7)
ax.set_title("(b) Cluster area", loc="left")
ax = axs[0, 2]
ax.plot(tab.year, tab.morans_I, "k-o", ms=3)
ax.set_ylabel("Global Moran's I of log depth")
ax.set_title("(c) Strength of spatial clustering", loc="left")
ax = axs[1, 0]
ax.plot(tab.year, tab.shallow_median_depth_m, "-o", color=C["LL"], ms=3, label="Shallow cluster")
ax.set_ylabel("Median depth, shallow cluster (m)", color=C["LL"])
ax2 = ax.twinx()
ax2.plot(tab.year, tab.deep_median_depth_m, "-o", color=C["HH"], ms=3)
ax2.set_ylabel("Median depth, deep cluster (m)", color=C["HH"])
ax.set_title("(d) Depth to water inside the clusters", loc="left")
ax = axs[1, 1]
ax.plot(tab.year, tab.median_depth_cells_m, "-o", color="0.5", ms=3, label="All observed cells (network changes)")
ax.plot(tab.year, tab.fixed_cells_median_depth_m, "k-o", ms=3,
        label=f"Cells observed every year ({int(always.sum())}), median")
ax.plot(tab.year, tab.fixed_cells_mean_depth_m, "k--", label="Cells observed every year, area mean")
ax.set_ylabel("Depth to water (m)")
ax.invert_yaxis()
ax.legend(fontsize=7)
ax.set_title("(e) Aquifer-wide depth to water", loc="left")
ax = axs[1, 2]
ax.plot(tab.year, tab.cells_le_5m_pct_of_observed, "-o", color=C["LL"], ms=3, label="0-5 m (Kollet & Maxwell range)")
ax.plot(tab.year, tab.cells_gt_60m_pct_of_observed, "-o", color=C["HH"], ms=3, label="> 60 m")
ax.set_ylabel("Share of observed area (%)")
ax.legend(fontsize=7)
ax.set_title("(f) Fixed depth classes", loc="left")
for ax in axs.flat:
    ax.grid(color="0.9", linewidth=0.4)
    ax.set_xlabel("Year (winter measurement)")
fig.suptitle("Yearly metrics of shallow and deep groundwater, High Plains aquifer", fontsize=12)
fig.savefig(OUT / "figB_yearly_metrics.png")
plt.close(fig)

# ---- C: persistence ----
fig, axs = plt.subplots(1, 3, figsize=(15, 7.8), constrained_layout=True)
seen = n_obs > 0
for ax, v, cmap, ttl in ((axs[0], cellt.pct_years_shallow, "Blues", "(a) Share of observed years in the SHALLOW cluster (%)"),
                         (axs[1], cellt.pct_years_deep, "Reds", "(b) Share of observed years in the DEEP cluster (%)")):
    clipped[seen].assign(v=v[seen].values).plot(ax=ax, column="v", cmap=cmap, vmin=0, vmax=100, legend=True,
                                                edgecolor="white", linewidth=0.15,
                                                legend_kwds=dict(orientation="horizontal", shrink=0.8, pad=0.06))
    clipped[~seen].plot(ax=ax, color="none", edgecolor="0.6", linewidth=0.2, hatch="///")
    base(ax, ttl)
clipped[seen].assign(v=n_obs[seen]).plot(ax=axs[2], column="v", cmap="viridis", vmin=0, vmax=len(YEARS), legend=True,
                                         edgecolor="white", linewidth=0.15,
                                         legend_kwds=dict(orientation="horizontal", shrink=0.8, pad=0.06))
clipped[~seen].plot(ax=axs[2], color="none", edgecolor="0.6", linewidth=0.2, hatch="///")
base(axs[2], f"(c) Years with data (of {len(YEARS)})")
core_sh = 100 * area[(cellt.pct_years_shallow >= 80) & (n_obs >= 10)].sum() / AQ_AREA
core_dp = 100 * area[(cellt.pct_years_deep >= 80) & (n_obs >= 10)].sum() / AQ_AREA
fig.suptitle(f"Persistence of the clusters, {YEARS[0]}-{YEARS[-1]}: cells in the shallow cluster in at least 80 % of "
             f"their years cover {core_sh:.1f} % of the aquifer, in the deep cluster {core_dp:.1f} %", fontsize=11)
fig.savefig(OUT / "figC_cluster_persistence.png")
plt.close(fig)

# ---- D: trend ----
fig, axs = plt.subplots(1, 2, figsize=(12.5, 7.8), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.2, 1]))
has_t = np.isfinite(slope)
lim = float(np.nanpercentile(np.abs(slope), 97))
clipped[has_t].assign(v=slope[has_t]).plot(ax=axs[0], column="v", cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim),
                                           legend=True, edgecolor="white", linewidth=0.15,
                                           legend_kwds=dict(orientation="horizontal", shrink=0.8, pad=0.06,
                                                            label="Trend of depth to water (m per year); red = deepening"))
sigc = has_t & (ptrend < 0.05)
clipped[sigc].centroid.plot(ax=axs[0], color="k", markersize=2)
clipped[~has_t].plot(ax=axs[0], color="none", edgecolor="0.6", linewidth=0.2, hatch="///")
base(axs[0], f"(a) Theil-Sen trend {YEARS[0]}-{YEARS[-1]} (dots: Mann-Kendall p < 0.05)")
ever_sh, ever_dp = (cellt.pct_years_shallow >= 50).values & has_t, (cellt.pct_years_deep >= 50).values & has_t
other = has_t & ~ever_sh & ~ever_dp
axs[1].hist([slope[ever_sh], slope[other], slope[ever_dp]], bins=np.linspace(-lim, lim, 31), stacked=True,
            color=[C["LL"], C["ns"], C["HH"]], label=[f"Mostly shallow cluster (median {np.median(slope[ever_sh]):+.3f})",
                                                      f"Other (median {np.median(slope[other]):+.3f})",
                                                      f"Mostly deep cluster (median {np.median(slope[ever_dp]):+.3f})"])
axs[1].axvline(0, color="k", linewidth=0.6)
axs[1].set_xlabel("Trend of depth to water (m per year)")
axs[1].set_ylabel("Cells")
axs[1].legend(fontsize=8)
axs[1].set_title("(b) Trends by cluster membership", loc="left")
fig.suptitle(f"Change of depth to water per cell: {int(has_t.sum())} cells with enough years; deepening "
             f"{100 * area[has_t & (slope > 0)].sum() / area[has_t].sum():.0f} % of their area, significantly "
             f"{100 * area[sigc & (slope > 0)].sum() / area[has_t].sum():.0f} %", fontsize=11)
fig.savefig(OUT / "figD_depth_trend.png")
plt.close(fig)

# ---- E: climate envelope of each year's clusters ----
fig, axs = plt.subplots(1, 4, figsize=(16, 4.2), constrained_layout=True)
for ax, (k, lab) in zip(axs, (("P_mm", "Precipitation normal (mm/yr)"), ("aridity", "Aridity index P / reference ET"),
                              ("T_C", "Mean temperature normal (deg C)"), ("VPD_kPa", "VPD normal (kPa)"))):
    ax.plot(tab.year, tab[f"shallow_{k}"], "-o", color=C["LL"], ms=3, label="Shallow cluster")
    ax.plot(tab.year, tab[f"deep_{k}"], "-o", color=C["HH"], ms=3, label="Deep cluster")
    ax.axhline(wmean(cn[k], np.ones(NC, bool)), color="k", linestyle="--", linewidth=0.8, label="Whole aquifer")
    ax.set_ylabel(lab)
    ax.set_xlabel("Year")
    ax.grid(color="0.9", linewidth=0.4)
axs[0].legend(fontsize=8)
fig.suptitle(f"Climate normals ({CLIM_YEARS[0]}-{CLIM_YEARS[1]}, gridMET) of the area each year's clusters occupy",
             fontsize=11)
fig.savefig(OUT / "figE_climate_of_clusters_by_year.png")
plt.close(fig)

# ════════════════════ console ════════════════════
pd.set_option("display.width", 260)
print(tab[["year", "wells", "observed_pct_of_aquifer", "morans_I", "shallow_pct_of_aquifer", "deep_pct_of_aquifer",
           "shallow_pct_of_observed", "deep_pct_of_observed", "shallow_median_depth_m", "deep_median_depth_m",
           "fixed_cells_median_depth_m", "fixed_cells_mean_depth_m", "cells_le_5m_pct_of_observed",
           "cells_gt_60m_pct_of_observed"]].round(2).to_string(index=False))
print(tab[["year", "shallow_P_mm", "deep_P_mm", "shallow_aridity", "deep_aridity", "shallow_T_C", "deep_T_C",
           "shallow_lat", "deep_lat"]].round(2).to_string(index=False))
a = tab.fixed_cells_mean_depth_m.values
log(f"fixed cells area-mean depth: {a[0]:.2f} m in {YEARS[0]} -> {a[-1]:.2f} m in {YEARS[-1]}; Theil-Sen "
    f"{stats.theilslopes(a, yrs)[0]:+.3f} m/yr (Mann-Kendall p = {stats.kendalltau(yrs, a).pvalue:.2g})")
log(f"core areas (>= 80 % of years): shallow {core_sh:.1f} %, deep {core_dp:.1f} % of the aquifer")
log(f"written to {OUT}")
