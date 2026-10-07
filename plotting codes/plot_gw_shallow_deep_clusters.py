"""
Shallow and deep groundwater areas of the High Plains (Ogallala) aquifer as spatial clusters, with their climate.

Method
------
Spatial clusters are found with the Local Moran statistic of
    Anselin, L., 1995, Local Indicators of Spatial Association - LISA. Geographical Analysis 27, 93-115,
    doi:10.1111/j.1538-4632.1995.tb00338.x
applied to depth to water (2019, USGS wells), exactly as defined there:
    z_i  = standardised log10(median depth to water of cell i)
    I_i  = z_i * sum_j w_ij z_j           (w = row-standardised queen contiguity between observed cells)
    significance by conditional permutation (9,999), one-sided pseudo p, false-discovery-rate control at 0.05
    High-High  = deep cell among deep cells      -> DEEP groundwater cluster
    Low-Low    = shallow cell among shallow cells -> SHALLOW groundwater cluster
    High-Low / Low-High = spatial outliers, everything else not significant
A second, fixed classification uses the 1-5 m critical water-table depth of
    Kollet, S.J., and Maxwell, R.M., 2008, Water Resources Research 44, W02402, doi:10.1029/2007WR006004
(the range in which the land-surface energy budget is most sensitive to groundwater; Little Washita, Oklahoma).

Data
----
Groundwater (folder derived_usgs): GW_wells_HPA_CF.nc (depth to water, 2019 campaign), dWL_*_4km.nc,
K_hydraulic_conductivity_mday_4km.nc, PAW_rootznaws_4km.nc, GDE_30arcsec_Ogallala_CF.nc.
Climate: GRIDMET_Ogallala_lab_drive.nc (gridMET from the Climatology Lab), 1981-2010 normals.
Known defects of the derived_usgs files are handled here, not hidden: K = 0 cells and K cells outside the
aquifer mask are dropped; PAW is treated as mm (its 'cm' label is wrong).

Nothing is interpolated: a cell without wells stays empty.
"""
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
import xarray as xr
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm
from matplotlib.patches import Patch
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = Path(r"G:/MSU_GWB/datasets")
USGS = ROOT / "merged_datasets" / "USGS data" / "derived_usgs"
GRIDMET = ROOT / "merged_datasets" / "GRIDMET_Ogallala_lab_drive.nc"
BOUND = "zip://G:/MSU_GWB/datasets/high_plains_quifer.zip!high_plains_quifer/hp_bound2010.shp"
OUT = ROOT / "plots" / "gw_shallow_deep_clusters"
OUT.mkdir(parents=True, exist_ok=True)

DLAT, DLON = 0.25, 0.30            # analysis cell, about 28 km x 26 km at 38 N
MIN_WELLS = 2                      # wells needed for a cell value
N_PERM = 9999
ALPHA = 0.05
CLIM_YEARS = (1981, 2010)
FT = 0.3048
DEPTH_BINS = [0, 5, 15, 30, 60, 1e9]
DEPTH_LABELS = ["0-5 m (land-surface coupled)", "5-15 m", "15-30 m", "30-60 m", "> 60 m"]
RNG = np.random.default_rng(20261007)
C = dict(HH="#b2182b", LL="#2166ac", HL="#f4a582", LH="#92c5de", ns="#d9d9d9")
NAME = dict(HH="Deep cluster (High-High)", LL="Shallow cluster (Low-Low)", HL="Deep outlier (High-Low)",
            LH="Shallow outlier (Low-High)", ns="Not significant")


def log(m):
    print(f"[gw-clusters] {m}")


# ════════════════════════════════════════════════════════════════════════════
# 1. Aquifer, analysis grid, wells
# ════════════════════════════════════════════════════════════════════════════
roi = gpd.read_file(BOUND).to_crs(4326).geometry.union_all()
roi_ea = gpd.GeoSeries([roi], crs=4326).to_crs(5070)
AQ_AREA = float(roi_ea.area.iloc[0]) / 1e6
W, S, E, N = roi.bounds
lon_edges = np.arange(np.floor(W / DLON) * DLON, E + DLON, DLON)
lat_edges = np.arange(np.floor(S / DLAT) * DLAT, N + DLAT, DLAT)
NY, NX = lat_edges.size - 1, lon_edges.size - 1
boxes = [shapely.box(lon_edges[j], lat_edges[i], lon_edges[j + 1], lat_edges[i + 1])
         for i in range(NY) for j in range(NX)]
cells = gpd.GeoDataFrame(dict(iy=np.repeat(np.arange(NY), NX), ix=np.tile(np.arange(NX), NY)), geometry=boxes, crs=4326)
inter = cells.to_crs(5070).geometry.intersection(roi_ea.iloc[0])
cells["aq_km2"] = inter.area.values / 1e6
cells["cell_km2"] = cells.to_crs(5070).area.values / 1e6
cells = cells[cells.aq_km2 > 0].reset_index(drop=True)
cells["lat"] = lat_edges[cells.iy] + DLAT / 2
cells["lon"] = lon_edges[cells.ix] + DLON / 2
log(f"aquifer {AQ_AREA:,.0f} km2; analysis grid {DLAT} x {DLON} deg; {len(cells)} cells touch the aquifer")


def cell_of(lat, lon):
    iy = np.floor((np.asarray(lat) - lat_edges[0]) / DLAT).astype(int)
    ix = np.floor((np.asarray(lon) - lon_edges[0]) / DLON).astype(int)
    return iy, ix


key = {(a, b): k for k, (a, b) in enumerate(zip(cells.iy, cells.ix))}


def to_cell_index(lat, lon):
    iy, ix = cell_of(lat, lon)
    return np.array([key.get((a, b), -1) for a, b in zip(iy, ix)])


with nc4.Dataset(USGS / "GW_wells_HPA_CF.nc") as w:
    w.set_auto_maskandscale(False)
    txt = lambda n: np.array([b"".join(r).decode().strip() for r in w[n][:]])
    st = pd.DataFrame(dict(lat=w["lat"][:], lon=w["lon"][:], well_depth_ft=w["well_depth_ft"][:],
                           state=txt("state")))
    ob = pd.DataFrame(dict(si=w["station_index"][:], camp=txt("campaign"), wl_ft=w["water_level_ft"][:]))
w19 = ob[ob.camp == "2019"].groupby("si").wl_ft.first()
wells = st.loc[w19.index].copy()
wells["dtw_m"] = w19.values * FT
wells["well_depth_m"] = wells.well_depth_ft * FT
wells["cell"] = to_cell_index(wells.lat, wells.lon)
wells = wells[wells.cell >= 0]
log(f"wells with a 2019 depth to water: {len(wells)}; median {wells.dtw_m.median():.1f} m, "
    f"range {wells.dtw_m.min():.1f} to {wells.dtw_m.max():.1f} m")

g = wells.groupby("cell")
cells["n_wells"] = g.size().reindex(cells.index).fillna(0).astype(int).values
cells["dtw_m"] = g.dtw_m.median().reindex(cells.index).values
cells["dtw_iqr_m"] = (g.dtw_m.quantile(0.75) - g.dtw_m.quantile(0.25)).reindex(cells.index).values
cells["well_depth_m"] = g.well_depth_m.median().reindex(cells.index).values
cells.loc[cells.n_wells < MIN_WELLS, ["dtw_m", "dtw_iqr_m", "well_depth_m"]] = np.nan
obs = cells.dtw_m.notna().values
OBS_AREA = cells.aq_km2[obs].sum()
log(f"cells with at least {MIN_WELLS} wells: {obs.sum()} of {len(cells)} = {100 * OBS_AREA / AQ_AREA:.1f} % of the "
    f"aquifer area; wells per observed cell: median {int(cells.n_wells[obs].median())}")

# ════════════════════════════════════════════════════════════════════════════
# 2. Local Moran clusters (Anselin 1995)
# ════════════════════════════════════════════════════════════════════════════
idx = np.flatnonzero(obs)
pos = {(a, b): k for k, (a, b) in enumerate(zip(cells.iy[idx], cells.ix[idx]))}
nbrs = []
for a, b in zip(cells.iy[idx], cells.ix[idx]):
    nbrs.append([pos[(a + da, b + db)] for da in (-1, 0, 1) for db in (-1, 0, 1)
                 if (da or db) and (a + da, b + db) in pos])
has = np.array([len(n) > 0 for n in nbrs])
x = np.log10(np.clip(cells.dtw_m.values[idx], 0.3, None))
z = (x - x[has].mean()) / x[has].std()
n = len(z)
lag = np.array([z[nb].mean() if nb else np.nan for nb in nbrs])
I_loc = z * lag
I_glob = float(np.nansum(z[has] * lag[has]) / np.sum(z[has] ** 2))

p_loc = np.full(n, np.nan)
for i in np.flatnonzero(has):
    k = len(nbrs[i])
    others = np.delete(z, i)
    pick = np.argpartition(RNG.random((N_PERM, n - 1)), k - 1, axis=1)[:, :k]      # k distinct cells, no replacement
    sim = others[pick].mean(axis=1)
    larger = int((sim >= lag[i]).sum())
    p_loc[i] = (min(larger, N_PERM - larger) + 1) / (N_PERM + 1)
# global Moran's I: permute the whole map
sim_g = np.empty(N_PERM)
hz = np.flatnonzero(has)
for r in range(N_PERM):
    zp = z.copy()
    zp[hz] = RNG.permutation(z[hz])
    sim_g[r] = np.sum(zp[hz] * np.array([zp[nbrs[i]].mean() for i in hz])) / np.sum(zp[hz] ** 2)
p_glob = (int((sim_g >= I_glob).sum()) + 1) / (N_PERM + 1)
# Benjamini-Hochberg false discovery rate over the local tests
pv = p_loc[has]
order = np.argsort(pv)
thresh = ALPHA * np.arange(1, pv.size + 1) / pv.size
passed = pv[order] <= thresh
cut = pv[order][np.flatnonzero(passed).max()] if passed.any() else -1.0
sig = has & (p_loc <= cut)
quad = np.where(z >= 0, np.where(lag >= 0, "HH", "HL"), np.where(lag >= 0, "LH", "LL"))
cells["cluster"] = "no data"
cells.loc[idx, "cluster"] = np.where(sig, quad, "ns")
cells.loc[idx, "z"] = z
cells.loc[idx, "lag"] = lag
cells.loc[idx, "local_I"] = I_loc
cells.loc[idx, "p_local"] = p_loc
log(f"global Moran's I of log depth = {I_glob:.3f} (pseudo p = {p_glob:.4f}, {N_PERM} permutations); "
    f"FDR p cut-off {cut:.4f}; significant cells {int(sig.sum())} of {int(has.sum())}")

cells["depth_class"] = pd.cut(cells.dtw_m, DEPTH_BINS, labels=DEPTH_LABELS, include_lowest=True)
wells["depth_class"] = pd.cut(wells.dtw_m.clip(lower=0), DEPTH_BINS, labels=DEPTH_LABELS, include_lowest=True)
wells["cluster"] = cells.cluster.values[wells.cell.values]

# ════════════════════════════════════════════════════════════════════════════
# 3. Climate normals from gridMET (1981-2010) and aquifer properties per cell
# ════════════════════════════════════════════════════════════════════════════
MONTH_DAYS = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
with nc4.Dataset(GRIDMET) as ds:
    glat, glon = ds["lat"][:].data, ds["lon"][:].data
    t = pd.to_datetime("1900-01-01") + pd.to_timedelta(ds["time"][:].data, unit="D")
    inside = ds["roi_mask"][:].data == 1
    sel = np.flatnonzero((t.year >= CLIM_YEARS[0]) & (t.year <= CLIM_YEARS[1]))
    months, years = t.month.values, t.year.values
    clim = {}
    for v in ("pr", "pet", "tmmx", "tmmn", "vpd"):
        tot = np.zeros((12,) + inside.shape)
        cnt = np.zeros((12,) + inside.shape)
        for y in range(CLIM_YEARS[0], CLIM_YEARS[1] + 1):
            k = sel[years[sel] == y]
            a = np.ma.filled(ds[v][k[0]:k[-1] + 1].astype("f8"), np.nan)
            mm = months[k]
            for m in range(12):
                part = a[mm == m + 1]
                tot[m] += np.nansum(part, axis=0)
                cnt[m] += np.isfinite(part).sum(axis=0)
        clim[v] = tot / np.maximum(cnt, 1)                    # mean daily value per calendar month
        log(f"gridMET {v}: {CLIM_YEARS[0]}-{CLIM_YEARS[1]} monthly normals done")
GLON, GLAT = np.meshgrid(glon, glat)
gcell = to_cell_index(GLAT[inside], GLON[inside])
ok = gcell >= 0
monthly = dict(P=clim["pr"] * MONTH_DAYS[:, None, None], PET=clim["pet"] * MONTH_DAYS[:, None, None],
               T=(clim["tmmx"] + clim["tmmn"]) / 2 - 273.15, VPD=clim["vpd"])


def to_cells(field):
    """Mean of the gridMET cells (inside the aquifer) that fall in each analysis cell."""
    s = np.bincount(gcell[ok], weights=field[inside][ok], minlength=len(cells))
    c = np.bincount(gcell[ok], minlength=len(cells))
    return np.where(c > 0, s / np.maximum(c, 1), np.nan)


cm = {k: np.stack([to_cells(v[m]) for m in range(12)]) for k, v in monthly.items()}      # (12, n cells)
cells["P_mm"] = cm["P"].sum(0)
cells["PET_mm"] = cm["PET"].sum(0)
cells["T_C"] = (cm["T"] * MONTH_DAYS[:, None]).sum(0) / MONTH_DAYS.sum()
cells["VPD_kPa"] = (cm["VPD"] * MONTH_DAYS[:, None]).sum(0) / MONTH_DAYS.sum()
cells["aridity_P_PET"] = cells.P_mm / cells.PET_mm
cells["deficit_mm"] = cells.P_mm - cells.PET_mm
cells["P_warm_frac"] = cm["P"][3:9].sum(0) / cm["P"].sum(0)         # April-September share of precipitation
cells["T_Jan_C"], cells["T_Jul_C"] = cm["T"][0], cm["T"][6]


def grid_to_cells(path, var, clean=None):
    d = xr.open_dataset(path)
    a = d[var].values.squeeze().astype("f8")
    if clean is not None:
        a = clean(a)
    yy, xx = np.meshgrid(d.y.values, d.x.values, indexing="ij")
    good = np.isfinite(a)
    ci = to_cell_index(yy[good], xx[good])
    k = ci >= 0
    s = np.bincount(ci[k], weights=a[good][k], minlength=len(cells))
    c = np.bincount(ci[k], minlength=len(cells))
    return np.where(c > 0, s / np.maximum(c, 1), np.nan)


mask4 = xr.open_dataset(USGS / "aquifer_mask_4km.nc").aquifer_mask.values == 1
cells["dWL_predev_2019_m"] = grid_to_cells(USGS / "dWL_predev_to_2019_ft_4km.nc", "dWL_predev_to_2019") * FT
cells["dWL_2017_2019_m"] = grid_to_cells(USGS / "dWL_2017_to_2019_ft_4km.nc", "dWL_2017_to_2019") * FT
cells["K_m_day"] = grid_to_cells(USGS / "K_hydraulic_conductivity_mday_4km.nc", "K_mday",
                                 lambda a: np.where(mask4 & (a > 0), a, np.nan))
cells["PAW_mm"] = grid_to_cells(USGS / "PAW_rootznaws_4km.nc", "root_zone_available_water_storage")
gde = xr.open_dataset(USGS / "GDE_30arcsec_Ogallala_CF.nc", mask_and_scale=True)
cells["GDE_frac_of_analysed"] = grid_to_cells(USGS / "GDE_30arcsec_Ogallala_CF.nc", "GDE_frac_AA")
gde.close()

# ════════════════════════════════════════════════════════════════════════════
# 4. Tables
# ════════════════════════════════════════════════════════════════════════════
ORDER = ["LL", "HH", "LH", "HL", "ns"]
METRICS = ["dtw_m", "well_depth_m", "P_mm", "PET_mm", "aridity_P_PET", "deficit_mm", "T_C", "T_Jan_C", "T_Jul_C",
           "VPD_kPa", "P_warm_frac", "dWL_predev_2019_m", "dWL_2017_2019_m", "K_m_day", "PAW_mm",
           "GDE_frac_of_analysed"]


def wmean(v, w):
    k = np.isfinite(v)
    return float(np.sum(v[k] * w[k]) / np.sum(w[k])) if k.any() else np.nan


rows = []
for c in ORDER + ["observed", "no data", "aquifer"]:
    m = (cells.cluster == c) if c in ORDER + ["no data"] else (obs if c == "observed" else np.ones(len(cells), bool))
    m = np.asarray(m)
    sub = cells[m]
    r = dict(group=NAME.get(c, c), cells=int(m.sum()), area_km2=sub.aq_km2.sum(),
             pct_of_aquifer=100 * sub.aq_km2.sum() / AQ_AREA,
             pct_of_observed_area=100 * sub.aq_km2.sum() / OBS_AREA if c in ORDER + ["observed"] else np.nan,
             wells=int(sub.n_wells.sum()), pct_of_wells=100 * sub.n_wells.sum() / cells.n_wells.sum(),
             lat_centroid=wmean(sub.lat.values, sub.aq_km2.values), lon_centroid=wmean(sub.lon.values, sub.aq_km2.values))
    for k in METRICS:
        r[k + "_mean"] = wmean(sub[k].values, sub.aq_km2.values)
        r[k + "_median"] = float(np.nanmedian(sub[k])) if sub[k].notna().any() else np.nan
        r[k + "_p25"] = float(np.nanpercentile(sub[k], 25)) if sub[k].notna().any() else np.nan
        r[k + "_p75"] = float(np.nanpercentile(sub[k], 75)) if sub[k].notna().any() else np.nan
    rows.append(r)
summary = pd.DataFrame(rows)
summary.to_csv(OUT / "table1_cluster_summary.csv", index=False, float_format="%.4g")

tests = []
LL, HH = cells[cells.cluster == "LL"], cells[cells.cluster == "HH"]
for k in METRICS:
    a, b = LL[k].dropna(), HH[k].dropna()
    if len(a) > 2 and len(b) > 2:
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rho, prho = stats.spearmanr(cells.dtw_m[obs], cells[k][obs], nan_policy="omit")
        tests.append(dict(variable=k, shallow_median=a.median(), deep_median=b.median(),
                          difference_deep_minus_shallow=b.median() - a.median(), mann_whitney_p=p,
                          rank_biserial=1 - 2 * u / (len(a) * len(b)),      # + : deep cluster has larger values
                          spearman_rho_with_depth_all_cells=rho, spearman_p=prho, n_shallow=len(a), n_deep=len(b)))
tests = pd.DataFrame(tests)
tests.to_csv(OUT / "table2_shallow_vs_deep_tests.csv", index=False, float_format="%.4g")

dc = cells[obs].groupby("depth_class", observed=False).agg(cells=("aq_km2", "size"), area_km2=("aq_km2", "sum"))
dc["pct_of_observed_area"] = 100 * dc.area_km2 / OBS_AREA
dc["pct_of_aquifer"] = 100 * dc.area_km2 / AQ_AREA
wc = wells.groupby("depth_class", observed=False).size()
dc["wells"] = wc
dc["pct_of_wells"] = 100 * wc / len(wells)
dc.to_csv(OUT / "table3_fixed_depth_classes.csv", float_format="%.4g")

mon = []
for c in ("LL", "HH", "observed"):
    m = np.asarray(cells.cluster == c) if c != "observed" else obs
    wt = cells.aq_km2.values[m]
    for k in cm:
        vals = [wmean(cm[k][mth][m], wt) for mth in range(12)]
        mon.append(dict(group=NAME.get(c, "All observed cells"), variable=k, **{f"m{q + 1:02d}": v for q, v in enumerate(vals)}))
mon = pd.DataFrame(mon)
mon.to_csv(OUT / "table4_monthly_climatology.csv", index=False, float_format="%.4g")
cells.drop(columns="geometry").to_csv(OUT / "table5_cells.csv", index=False, float_format="%.5g")

# ════════════════════════════════════════════════════════════════════════════
# 5. Figures
# ════════════════════════════════════════════════════════════════════════════
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "savefig.dpi": 220, "figure.dpi": 110})
ASPECT = 1 / np.cos(np.radians(37.7))
outline = gpd.GeoSeries([roi], crs=4326)


def base(ax, title):
    outline.boundary.plot(ax=ax, color="k", linewidth=0.7)
    ax.set_xlim(W - 0.3, E + 0.3)
    ax.set_ylim(S - 0.3, N + 0.3)
    ax.set_aspect(ASPECT)
    ax.set_xlabel("Longitude (deg E)")
    ax.set_ylabel("Latitude (deg N)")
    ax.grid(color="0.85", linewidth=0.4)
    ax.set_title(title, loc="left")


clipped = cells.copy()
clipped["geometry"] = cells.geometry.intersection(roi)

# ---- Figure 1: depth to water ----
fig, axs = plt.subplots(1, 2, figsize=(11.5, 8.2), constrained_layout=True)
norm = LogNorm(1, 150)
sc = axs[0].scatter(wells.lon, wells.lat, c=wells.dtw_m.clip(lower=1), s=3, cmap="viridis_r", norm=norm, linewidths=0)
base(axs[0], f"(a) Depth to water at {len(wells):,} wells, 2019")
clipped[obs].plot(ax=axs[1], column="dtw_m", cmap="viridis_r", norm=norm, edgecolor="white", linewidth=0.2)
clipped[~obs].plot(ax=axs[1], color="none", edgecolor="0.6", linewidth=0.2, hatch="///")
base(axs[1], f"(b) Cell median ({DLAT} x {DLON} deg, at least {MIN_WELLS} wells)")
cb = fig.colorbar(sc, ax=axs, orientation="horizontal", shrink=0.6, pad=0.02, extend="both")
cb.set_label("Depth to water below land surface (m)")
axs[1].legend(handles=[Patch(facecolor="none", edgecolor="0.6", hatch="///", label="no wells (left empty)")],
              loc="lower left", frameon=True)
fig.suptitle("High Plains aquifer: depth to water, 2019 (USGS wells)", fontsize=12)
fig.savefig(OUT / "fig1_depth_to_water_2019.png")
plt.close(fig)

# ---- Figure 2: LISA clusters ----
fig, axs = plt.subplots(1, 2, figsize=(12.5, 8.2), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.15, 1]))
handles = []
for c in ORDER:
    m = cells.cluster == c
    if m.any():
        clipped[m].plot(ax=axs[0], color=C[c], edgecolor="white", linewidth=0.2)
    handles.append(Patch(facecolor=C[c], label=f"{NAME[c]}: {100 * cells.aq_km2[m].sum() / OBS_AREA:.1f} % of observed "
                                                f"area, {100 * cells.aq_km2[m].sum() / AQ_AREA:.1f} % of aquifer"))
clipped[~obs].plot(ax=axs[0], color="none", edgecolor="0.6", linewidth=0.2, hatch="///")
handles.append(Patch(facecolor="none", edgecolor="0.6", hatch="///",
                     label=f"No wells: {100 * (1 - OBS_AREA / AQ_AREA):.1f} % of aquifer"))
base(axs[0], "(a) Shallow and deep groundwater clusters (Local Moran, FDR 0.05)")
axs[0].legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.07), frameon=False, fontsize=8)
cl = cells.cluster.values[idx]
for c in ORDER:
    m = cl == c
    axs[1].scatter(z[m], lag[m], s=14, color=C[c], edgecolor="k" if c != "ns" else "none", linewidths=0.3, label=NAME[c])
xx = np.linspace(np.nanmin(z), np.nanmax(z), 10)
axs[1].plot(xx, I_glob * xx, "k--", linewidth=1, label=f"slope = global Moran's I = {I_glob:.2f}")
axs[1].axhline(0, color="0.5", linewidth=0.5)
axs[1].axvline(0, color="0.5", linewidth=0.5)
axs[1].set_xlabel("Standardised log10 depth to water of the cell")
axs[1].set_ylabel("Mean of its neighbours (spatial lag)")
axs[1].set_title(f"(b) Moran scatterplot (pseudo p = {p_glob:.4f}, {N_PERM:,} permutations)", loc="left")
axs[1].legend(fontsize=8, loc="upper left")
fig.suptitle("Spatial clusters of depth to water, High Plains aquifer (method: Anselin 1995)", fontsize=12)
fig.savefig(OUT / "fig2_shallow_deep_clusters.png")
plt.close(fig)

# ---- Figure 3: fixed depth classes ----
cols = ["#2166ac", "#67a9cf", "#fddbc7", "#ef8a62", "#b2182b"]
fig, axs = plt.subplots(1, 2, figsize=(12, 8.2), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.3, 1]))
for lab, col in zip(DEPTH_LABELS, cols):
    m = (cells.depth_class == lab).values
    if m.any():
        clipped[m].plot(ax=axs[0], color=col, edgecolor="white", linewidth=0.2)
clipped[~obs].plot(ax=axs[0], color="none", edgecolor="0.6", linewidth=0.2, hatch="///")
base(axs[0], "(a) Depth-to-water class of the cell median")
axs[0].legend(handles=[Patch(facecolor=col, label=f"{lab}: {dc.pct_of_observed_area[lab]:.1f} % of observed area")
                       for lab, col in zip(DEPTH_LABELS, cols)], loc="upper center", bbox_to_anchor=(0.5, -0.07),
              frameon=False, fontsize=8)
yy = np.arange(len(DEPTH_LABELS))
axs[1].barh(yy - 0.2, dc.pct_of_observed_area.values, 0.4, color=cols, edgecolor="k", linewidth=0.4, label="share of observed area")
axs[1].barh(yy + 0.2, dc.pct_of_wells.values, 0.4, color=cols, edgecolor="k", linewidth=0.4, hatch="..", alpha=0.6,
            label="share of wells")
for q in range(len(yy)):
    axs[1].text(dc.pct_of_observed_area.values[q] + 0.5, q - 0.2, f"{dc.pct_of_observed_area.values[q]:.1f} %", va="center", fontsize=8)
    axs[1].text(dc.pct_of_wells.values[q] + 0.5, q + 0.2, f"{dc.pct_of_wells.values[q]:.1f} %", va="center", fontsize=8)
axs[1].set_yticks(yy, DEPTH_LABELS)
axs[1].invert_yaxis()
axs[1].set_xlabel("Per cent")
axs[1].set_title("(b) Shares by class", loc="left")
axs[1].legend(fontsize=8, loc="lower right")
fig.suptitle("Fixed depth classes (0-5 m: critical range of Kollet & Maxwell 2008)", fontsize=12)
fig.savefig(OUT / "fig3_fixed_depth_classes.png")
plt.close(fig)

# ---- Figure 4: monthly climatology of the clusters ----
fig, axs = plt.subplots(2, 3, figsize=(13, 7.2), constrained_layout=True)
mx = np.arange(1, 13)
panels = [("P", "Precipitation (mm per month)"), ("PET", "Reference ET, grass (mm per month)"),
          ("T", "Mean air temperature (deg C)"), ("VPD", "Vapour pressure deficit (kPa)")]
grp = [("LL", C["LL"], "Shallow cluster"), ("HH", C["HH"], "Deep cluster"), ("observed", "k", "All observed cells")]
series = {}
for c, col, lab in grp:
    m = np.asarray(cells.cluster == c) if c != "observed" else obs
    wt = cells.aq_km2.values[m]
    series[c] = {k: np.array([wmean(cm[k][q][m], wt) for q in range(12)]) for k in cm}
    q25 = {k: np.nanpercentile(cm[k][:, m], 25, axis=1) for k in cm}
    q75 = {k: np.nanpercentile(cm[k][:, m], 75, axis=1) for k in cm}
    for ax, (k, yl) in zip(axs.flat, panels):
        ax.plot(mx, series[c][k], color=col, marker="o", ms=3, linewidth=1.5 if c != "observed" else 1,
                linestyle="-" if c != "observed" else "--", label=lab)
        if c != "observed":
            ax.fill_between(mx, q25[k], q75[k], color=col, alpha=0.15, linewidth=0)
        ax.set_ylabel(yl)
    axs[1, 1].plot(mx, series[c]["P"] - series[c]["PET"], color=col, marker="o", ms=3,
                   linestyle="-" if c != "observed" else "--", label=lab)
axs[1, 1].axhline(0, color="0.5", linewidth=0.5)
axs[1, 1].set_ylabel("P minus reference ET (mm per month)")
for ax, ttl in zip(axs.flat[:5], "abcde"):
    ax.set_xticks(mx, list("JFMAMJJASOND"))
    ax.set_title(f"({ttl})", loc="left")
    ax.grid(color="0.9", linewidth=0.4)
axs[0, 0].legend(fontsize=8)
ax = axs[1, 2]
ax.axis("off")
lines = [f"{'':24s}{'Shallow':>9s}{'Deep':>9s}{'All':>9s}"]
for k, lab, f in (("P_mm", "P (mm/yr)", "{:.0f}"), ("PET_mm", "Ref. ET (mm/yr)", "{:.0f}"),
                  ("aridity_P_PET", "P / ref. ET", "{:.2f}"), ("deficit_mm", "P - ref. ET (mm/yr)", "{:.0f}"),
                  ("T_C", "T mean (C)", "{:.1f}"), ("T_Jan_C", "T January (C)", "{:.1f}"),
                  ("T_Jul_C", "T July (C)", "{:.1f}"), ("VPD_kPa", "VPD (kPa)", "{:.2f}"),
                  ("P_warm_frac", "Apr-Sep share of P", "{:.2f}"), ("dtw_m", "Depth to water (m)", "{:.1f}")):
    v = [summary.loc[summary.group == NAME[c] if c in NAME else summary.group == c, k + "_mean"].iloc[0]
         for c in ("LL", "HH", "observed")]
    lines.append(f"{lab:24s}" + "".join(f"{f.format(q):>9s}" for q in v))
ax.text(0, 1, f"(f) Area-weighted means, {CLIM_YEARS[0]}-{CLIM_YEARS[1]}\n\n" + "\n".join(lines), family="monospace",
        fontsize=8.5, va="top")
fig.suptitle(f"Climate of the shallow and deep groundwater clusters (gridMET normals {CLIM_YEARS[0]}-{CLIM_YEARS[1]}; "
             f"bands = interquartile range of cells)", fontsize=11)
fig.savefig(OUT / "fig4_climatology_by_cluster.png")
plt.close(fig)

# ---- Figure 5: climate maps with the cluster outlines ----
fig, axs = plt.subplots(1, 3, figsize=(15, 7.8), constrained_layout=True)
for ax, (k, lab, cmap) in zip(axs, (("P_mm", "Precipitation (mm per year)", "YlGnBu"),
                                    ("aridity_P_PET", "Aridity index P / reference ET", "BrBG"),
                                    ("T_C", "Mean air temperature (deg C)", "RdYlBu_r"))):
    clipped.plot(ax=ax, column=k, cmap=cmap, legend=True, edgecolor="none",
                 legend_kwds=dict(orientation="horizontal", shrink=0.8, pad=0.06, label=lab))
    for c, ls in (("LL", C["LL"]), ("HH", C["HH"])):
        m = (cells.cluster == c).values
        if m.any():
            gpd.GeoSeries([clipped[m].geometry.union_all()], crs=4326).boundary.plot(ax=ax, color=ls, linewidth=1.6)
    base(ax, lab)
axs[0].legend(handles=[Patch(facecolor="none", edgecolor=C["LL"], linewidth=1.6, label="Shallow cluster"),
                       Patch(facecolor="none", edgecolor=C["HH"], linewidth=1.6, label="Deep cluster")], loc="lower left")
fig.suptitle(f"Climate normals {CLIM_YEARS[0]}-{CLIM_YEARS[1]} (gridMET) with the groundwater clusters outlined", fontsize=12)
fig.savefig(OUT / "fig5_climate_maps_with_clusters.png")
plt.close(fig)

# ---- Figure 6: aquifer and climate properties by cluster ----
props = [("dtw_m", "Depth to water (m)"), ("well_depth_m", "Well depth (m)"),
         ("dWL_predev_2019_m", "Level change, predev. to 2019 (m)"), ("dWL_2017_2019_m", "Level change, 2017-19 (m)"),
         ("K_m_day", "Hydraulic conductivity (m/day)"), ("PAW_mm", "Root-zone available water (mm)"),
         ("P_mm", "Precipitation (mm/yr)"), ("aridity_P_PET", "P / reference ET"),
         ("GDE_frac_of_analysed", "GDE share of analysed area")]
fig, axs = plt.subplots(3, 3, figsize=(12.5, 10), constrained_layout=True)
for ax, (k, lab) in zip(axs.flat, props):
    data = [cells.loc[cells.cluster == c, k].dropna().values for c in ("LL", "ns", "HH")]
    bp = ax.boxplot(data, tick_labels=["Shallow", "Not sig.", "Deep"], patch_artist=True, showfliers=False, widths=0.6)
    for patch, c in zip(bp["boxes"], ("LL", "ns", "HH")):
        patch.set_facecolor(C[c])
    for med in bp["medians"]:
        med.set_color("k")
    row = tests[tests.variable == k]
    extra = f"  (shallow vs deep p = {row.mann_whitney_p.iloc[0]:.1e})" if len(row) else ""
    ax.set_title(lab + extra, loc="left", fontsize=8.5)
    ax.grid(axis="y", color="0.9", linewidth=0.4)
    ax.text(0.99, 0.02, "n = " + " / ".join(str(len(d)) for d in data), transform=ax.transAxes, ha="right", fontsize=7)
fig.suptitle("Aquifer and climate properties of the clusters (boxes: quartiles of cells; whiskers 1.5 IQR; "
             "Mann-Whitney test)", fontsize=11)
fig.savefig(OUT / "fig6_properties_by_cluster.png")
plt.close(fig)

# ════════════════════════════════════════════════════════════════════════════
# 6. Console summary
# ════════════════════════════════════════════════════════════════════════════
pd.set_option("display.width", 250)
print(summary[["group", "cells", "area_km2", "pct_of_aquifer", "pct_of_observed_area", "wells", "pct_of_wells",
               "lat_centroid", "lon_centroid", "dtw_m_median", "dtw_m_p25", "dtw_m_p75", "P_mm_mean", "PET_mm_mean",
               "aridity_P_PET_mean", "T_C_mean", "VPD_kPa_mean"]].round(2).to_string(index=False))
print(tests.round(4).to_string(index=False))
print(dc.round(2).to_string())
print("wells by cluster:", wells.cluster.value_counts().to_dict())
print("well-level: share <= 5 m %.1f %%, <= 10 m %.1f %%, > 30 m %.1f %%, > 60 m %.1f %%" % (
    100 * (wells.dtw_m <= 5).mean(), 100 * (wells.dtw_m <= 10).mean(), 100 * (wells.dtw_m > 30).mean(),
    100 * (wells.dtw_m > 60).mean()))
log(f"figures and tables written to {OUT}")
