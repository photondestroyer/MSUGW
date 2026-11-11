"""Planning checks: (a) flux towers inside the aquifer polygon, (b) where water levels rose, (c) wells usable for a two-drought panel."""
import json, sys
import numpy as np, netCDF4 as nc
import geopandas as gpd
from shapely.geometry import Point
from pyproj import Transformer
from scipy import ndimage as ndi

ROOT = "G:/MSU_GWB/datasets"
out = {}

# (a) towers
poly = gpd.read_file(f"{ROOT}/HPA_polygon/hp_bound2010.shp").to_crs(4326)
geom = poly.union_all() if hasattr(poly, "union_all") else poly.unary_union
towers = {"US-Ne1 (Mead irrigated maize)": (41.1651, -96.4766), "US-Ne2 (Mead irrigated rotation)": (41.1649, -96.4701),
          "US-Ne3 (Mead rainfed rotation)": (41.1797, -96.4397)}
proj = poly.to_crs(5070); g5 = proj.union_all() if hasattr(proj, "union_all") else proj.unary_union
tr = Transformer.from_crs(4326, 5070, always_xy=True)
out["towers"] = {}
for k, (la, lo) in towers.items():
    p = Point(*tr.transform(lo, la))
    out["towers"][k] = dict(inside=bool(g5.contains(p)), km_to_boundary=round(g5.boundary.distance(p) / 1000, 1))
print("TOWERS", out["towers"])
print("polygon bounds", [round(v, 3) for v in geom.bounds], "area km2", round(g5.area / 1e6))

# (b) water-level change map
d = nc.Dataset(f"{ROOT}/merged_datasets/USGS_HPA_native/USGS_HPA_water_level_change_500m.nc")
x, y = d["x"][:], d["y"][:]
m = d["aquifer_mask"][:] == 1
tb = Transformer.from_crs(5070, 4326, always_xy=True)
for name in ("dWL_predev_to_2019", "dWL_2017_to_2019"):
    a = np.ma.filled(d[name][:].astype("float64"), np.nan)
    a[(a < -1e30) | ~m] = np.nan
    ft = a; n = np.isfinite(ft).sum()
    rep = dict(cells=int(n), area_km2=round(n * 0.25), zero_share=float((ft == 0).sum() / n),
               rise_ge10ft_km2=round(float((ft >= 10).sum() * 0.25)), rise_ge25ft_km2=round(float((ft >= 25).sum() * 0.25)),
               decline_ge25ft_km2=round(float((ft <= -25).sum() * 0.25)), decline_ge100ft_km2=round(float((ft <= -100).sum() * 0.25)),
               values=np.unique(ft[np.isfinite(ft)])[:40].round(1).tolist())
    if name == "dWL_predev_to_2019":
        lab, k = ndi.label(np.nan_to_num(ft, nan=0) >= 10)
        sizes = ndi.sum(np.ones_like(lab), lab, range(1, k + 1)) * 0.25
        top = np.argsort(sizes)[::-1][:8]
        patches = []
        for t in top:
            yy, xx = np.nonzero(lab == t + 1)
            lo, la = tb.transform(float(x[xx].mean()), float(y[yy].mean()))
            patches.append(dict(km2=round(float(sizes[t])), lat=round(la, 2), lon=round(lo, 2), max_ft=float(np.nanmax(ft[lab == t + 1]))))
        rep["largest_rise_patches"] = patches
    out[name] = rep
    print(name, rep)
d.close()

# (c) wells
w = nc.Dataset(f"{ROOT}/merged_datasets/HPA_wells_water_level_timeseries.nc")
st = w["station_index"][:]; yr = w["year"][:]; dtw = np.ma.filled(w["depth_to_water"][:], np.nan) * 0.3048
dup = w["duplicate_of_other_source"][:]
tm = np.ma.filled(w["time"][:], np.nan)
state = np.array(w["state"][:], dtype=object); inpoly = w["in_polygon"][:] == 1
lat = w["lat"][:]; lon = w["lon"][:]
ok = np.isfinite(dtw) & (dtw >= 0) & (dtw < 300) & (dup == 0) & inpoly[st] & (yr >= 2000) & (yr <= 2026)
print("usable obs 2000-2026:", int(ok.sum()), "of", ok.size)
states = sorted(set(state[inpoly]))
print("states:", states)
tab = {}
for s in states:
    ms = ok & (state[st] == s)
    tab[s] = {int(y): int(np.unique(st[ms & (yr == y)]).size) for y in (2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024, 2025)}
    print(" wells measured per year", s, tab[s])
out["wells_per_year"] = tab

def period_median(lo, hi):
    mm = ok & (yr >= lo) & (yr <= hi)
    order = np.argsort(st[mm], kind="stable")
    s_, v_ = st[mm][order], dtw[mm][order]
    u, first = np.unique(s_, return_index=True)
    med = np.array([np.median(c) for c in np.split(v_, first[1:])])
    return dict(zip(u.tolist(), med.tolist()))
A, B = period_median(2010, 2013), period_median(2020, 2023)
both = sorted(set(A) & set(B))
a = np.array([A[k] for k in both]); b = np.array([B[k] for k in both]); ids = np.array(both)
print("wells measured in both 2010-13 and 2020-23:", len(both))
panel = {}
for lab_, lo_, hi_ in (("<=3 m", 0, 3), ("3-6 m", 3, 6), ("6-10 m", 6, 10), ("10-30 m", 10, 30), (">30 m", 30, 1e9)):
    mm = (a > lo_) & (a <= hi_) if lo_ > 0 else (a <= hi_)
    ch = b[mm] - a[mm]
    panel[lab_] = dict(n=int(mm.sum()), median_change_m=float(np.median(ch)) if mm.any() else None,
                       deeper_gt1m=int((ch > 1).sum()), shallower_gt1m=int((ch < -1).sum()), abs_gt2m=int((np.abs(ch) > 2).sum()),
                       by_state={s: int((state[ids[mm]] == s).sum()) for s in states})
    print(" 2010-13 depth", lab_, panel[lab_])
out["two_drought_panel"] = panel
sh = a <= 10
out["shallow_panel_n"] = int(sh.sum())
# shallow wells: how many are within the aquifer by 1-degree latitude band
band = np.floor(lat[ids[sh]]).astype(int)
out["shallow_panel_by_lat"] = {int(k): int(v) for k, v in zip(*np.unique(band, return_counts=True))}
print(" shallow (<=10 m) panel wells by latitude band:", out["shallow_panel_by_lat"])
# month of measurement, 2008+
mo = ((tm[ok & (yr >= 2008)] % 365.25) / 30.44).astype(int) + 1
ok2 = ok & (yr >= 2008) & np.isfinite(tm)
dates = nc.num2date(tm[ok2][:200000], "days since 1900-01-01", only_use_cftime_datetimes=False)
mon = np.bincount([dd.month for dd in dates], minlength=13)[1:]
out["month_share_first200k"] = (mon / mon.sum()).round(3).tolist()
print(" month share of measurements (sample of 200k, 2008+):", out["month_share_first200k"])
w.close()
json.dump(out, open(sys.argv[1], "w"), indent=1, default=float)
