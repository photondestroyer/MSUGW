"""Planning check: how much of the aquifer can serve as rain-only reference, by climate decile.
Read-only, bounded memory. Sources: lab gridMET (aridity), old HRES-WTD file (caveat: partly coarse), LANID."""
import json, sys, time
import numpy as np, netCDF4 as nc
from pyproj import Transformer

ROOT = "G:/MSU_GWB/datasets/merged_datasets"
OUT = sys.argv[1]
t0 = time.time()

# ---- 1. aridity index P / PET, 1985-2014, per lab cell
g = nc.Dataset(f"{ROOT}/GRIDMET_Ogallala_lab_latest.nc")
g.set_auto_maskandscale(False)
tv = g["time"]
dates = nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)
years = np.array([d.year for d in dates])
lat, lon = g["lat"][:].astype("float64"), g["lon"][:].astype("float64")
frac = g["roi_fraction"][:].astype("float64")
tot = {}
for name in ("pr", "pet"):
    v = g[name]
    sf, fill = float(v.scale_factor), int(v._FillValue)
    s = np.zeros(v.shape[1:], "float64"); n = np.zeros(v.shape[1:], "int64")
    for y in range(1985, 2015):
        i = np.flatnonzero(years == y)
        a = v[i[0]:i[-1] + 1]
        ok = a != fill
        s += np.where(ok, a, 0).sum(0) * sf; n += ok.sum(0)
    tot[name] = s / np.maximum(n, 1) * 365.25
    print(name, "dims", v.dimensions, "done", round(time.time() - t0), "s", flush=True)
g.close()
ai = tot["pr"] / np.maximum(tot["pet"], 1e-6)
aq = frac >= 0.5
edges = np.quantile(ai[aq], np.linspace(0, 1, 11))
dec = np.clip(np.searchsorted(edges, ai, side="right"), 1, 10)
print("aquifer cells", int(aq.sum()), "AI decile edges", np.round(edges, 3).tolist())
print("P by decile (mm/yr)", [round(float(tot['pr'][aq & (dec == d)].mean())) for d in range(1, 11)])

# ---- 2. sample the water-table file on a ~1 km lattice, chunk by chunk
h = nc.Dataset(f"{ROOT}/HRES-WTD_2015_Ogallala.nc")
b = h["b1"]; hx, hy = h["x"][:], h["y"][:]
print("HRES x", hx[0], hx[-1], "step", np.median(np.diff(hx)), "| y", hy[0], hy[-1], "step", np.median(np.diff(hy)), flush=True)
K = 33
P_lon, P_lat, P_w = [], [], []
ny, nx = b.shape[1:]
for r0 in range(0, ny, 2048):
    rows = np.arange(r0, min(r0 + 2048, ny)); rows = rows[rows % K == K // 2]
    if rows.size == 0: continue
    for c0 in range(0, nx, 2048):
        cols = np.arange(c0, min(c0 + 2048, nx)); cols = cols[cols % K == K // 2]
        if cols.size == 0: continue
        blk = b[0, r0:min(r0 + 2048, ny), c0:min(c0 + 2048, nx)]
        blk = np.ma.filled(blk, np.nan) if np.ma.isMaskedArray(blk) else blk
        sub = blk[np.ix_(rows - r0, cols - c0)]
        m = np.isfinite(sub)
        if m.any():
            yy, xx = np.meshgrid(hy[rows], hx[cols], indexing="ij")
            P_lon.append(xx[m]); P_lat.append(yy[m]); P_w.append(sub[m])
    print("  rows", r0, round(time.time() - t0), "s", flush=True)
h.close()
P_lon, P_lat, P_w = map(np.concatenate, (P_lon, P_lat, P_w))
print("sample points with a water-table value:", P_w.size)

# ---- 3. lab cell of each point
iy = np.abs(lat[:, None] - P_lat[None, :1]).argmin(0) if False else np.rint((P_lat - lat[0]) / (lat[1] - lat[0])).astype(int)
ix = np.rint((P_lon - lon[0]) / (lon[1] - lon[0])).astype(int)
inside = (iy >= 0) & (iy < lat.size) & (ix >= 0) & (ix < lon.size)
P_lon, P_lat, P_w, iy, ix = P_lon[inside], P_lat[inside], P_w[inside], iy[inside], ix[inside]
keep = aq[iy, ix]
P_lon, P_lat, P_w, iy, ix = P_lon[keep], P_lat[keep], P_w[keep], iy[keep], ix[keep]
P_dec = dec[iy, ix]
print("points in aquifer cells (roi_fraction >= 0.5):", P_w.size)

# ---- 4. LANID ever-irrigated and frequency at each point (row bands)
L = nc.Dataset(f"{ROOT}/LANID_v3_30m_1997-2017_Ogallala.nc")
lx, ly = L["x"][:], L["y"][:]
tr = Transformer.from_crs(4326, 5070, always_xy=True)
X, Y = tr.transform(P_lon, P_lat)
cx = np.rint((X - lx[0]) / (lx[1] - lx[0])).astype(int)
cy = np.rint((Y - ly[0]) / (ly[1] - ly[0])).astype(int)
okL = (cx >= 0) & (cx < lx.size) & (cy >= 0) & (cy < ly.size)
freq = np.full(P_w.size, 255, "uint8")
vF = L["derived"]["irrFreq"]; vF.set_auto_maskandscale(False)
order = np.argsort(cy)
for r0 in range(0, ly.size, 512):
    lo, hi = np.searchsorted(cy[order], [r0, r0 + 512])
    idx = order[lo:hi]; idx = idx[okL[idx]]
    if idx.size == 0: continue
    band = vF[r0:min(r0 + 512, ly.size), :]
    freq[idx] = band[cy[idx] - r0, cx[idx]]
L.close()
print("LANID irrFreq values at points:", dict(zip(*[a.tolist() for a in np.unique(freq, return_counts=True)])))
have = freq != 255
never = have & (freq == 0)

# ---- 5. tables
rows = []
for d in range(1, 11):
    m = (P_dec == d) & have
    n = int(m.sum())
    w = P_w[m]; nv = never[m]
    rows.append(dict(decile=d, n=n, ai_lo=float(edges[d - 1]), ai_hi=float(edges[d]),
                     wtd_med=float(np.median(w)),
                     gt15=float((w > 15).mean()), gt30=float((w > 30).mean()), gt60=float((w > 60).mean()),
                     never_irr=float(nv.mean()),
                     ref15=float(((w > 15) & nv).mean()), ref30=float(((w > 30) & nv).mean()), ref60=float(((w > 60) & nv).mean()),
                     nat_le2=float(((w <= 2) & nv).mean()), nat_2_5=float(((w > 2) & (w <= 5) & nv).mean()),
                     nat_5_10=float(((w > 5) & (w <= 10) & nv).mean()), nat_10_30=float(((w > 10) & (w <= 30) & nv).mean())))
    print(rows[-1], flush=True)
ever = have & (freq > 0)
fr = freq[ever]
summ = dict(n_points=int(have.sum()), ever_irrigated_share=float(ever.sum() / have.sum()),
            ever_irr_le2yrs=float((fr <= 2).mean()), ever_irr_le5yrs=float((fr <= 5).mean()), ever_irr_ge15yrs=float((fr >= 15).mean()),
            wtd_quantiles={q: float(np.quantile(P_w[have], q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95)},
            ai_edges=[float(e) for e in edges])
print(summ)
json.dump(dict(rows=rows, summary=summ), open(OUT, "w"), indent=1)
np.savez_compressed(OUT.replace(".json", ".npz"), lon=P_lon.astype("float32"), lat=P_lat.astype("float32"), wtd=P_w.astype("float32"), dec=P_dec.astype("int8"), freq=freq)
print("done", round(time.time() - t0), "s")
