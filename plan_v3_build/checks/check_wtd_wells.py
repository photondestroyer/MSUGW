"""Planning check: modelled water-table depth (old HRES-WTD file) against measured wells inside the aquifer."""
import json, sys
import numpy as np, netCDF4 as nc
ROOT = "G:/MSU_GWB/datasets/merged_datasets"
w = nc.Dataset(f"{ROOT}/HPA_wells_water_level_timeseries.nc")
st = w["station_index"][:]; yr = w["year"][:]; dtw = np.ma.filled(w["depth_to_water"][:], np.nan) * 0.3048
dup = w["duplicate_of_other_source"][:]; inpoly = w["in_polygon"][:] == 1
lat = np.ma.filled(w["lat"][:], np.nan); lon = np.ma.filled(w["lon"][:], np.nan)
state = np.array(w["state"][:], dtype=object); w.close()
ok = np.isfinite(dtw) & (dtw >= 0) & (dtw < 300) & (dup == 0) & inpoly[st] & (yr >= 2000) & (yr <= 2023)
order = np.argsort(st[ok], kind="stable"); s_, v_ = st[ok][order], dtw[ok][order]
u, first = np.unique(s_, return_index=True)
obs = np.array([np.median(c) for c in np.split(v_, first[1:])])
h = nc.Dataset(f"{ROOT}/HRES-WTD_2015_Ogallala.nc"); b = h["b1"]; hx, hy = h["x"][:], h["y"][:]
cx = np.rint((lon[u] - hx[0]) / (hx[1] - hx[0])).astype(int); cy = np.rint((lat[u] - hy[0]) / (hy[1] - hy[0])).astype(int)
good = (cx >= 0) & (cx < hx.size) & (cy >= 0) & (cy < hy.size)
mod = np.full(u.size, np.nan)
key = (cy // 2048) * 1000 + (cx // 2048)
for k in np.unique(key[good]):
    m = good & (key == k); r0, c0 = (k // 1000) * 2048, (k % 1000) * 2048
    blk = b[0, r0:r0 + 2048, c0:c0 + 2048]; blk = np.ma.filled(blk, np.nan) if np.ma.isMaskedArray(blk) else blk
    mod[m] = blk[cy[m] - r0, cx[m] - c0]
h.close()
f = np.isfinite(mod)
o, m_, stt = obs[f], mod[f], state[u][f]
print("wells with a 2000-2023 measurement and a map value:", o.size)
out = dict(n=int(o.size), pearson=float(np.corrcoef(o, m_)[0, 1]), spearman=float(np.corrcoef(np.argsort(np.argsort(o)), np.argsort(np.argsort(m_)))[0, 1]),
           bias=float(np.mean(m_ - o)), rmse=float(np.sqrt(np.mean((m_ - o) ** 2))), mae=float(np.mean(np.abs(m_ - o))))
print("all:", out)
cls = [("<=5", 0, 5), ("5-10", 5, 10), ("10-30", 10, 30), (">30", 30, 1e9)]
c_o = np.digitize(o, [5, 10, 30]); c_m = np.digitize(m_, [5, 10, 30])
conf = np.zeros((4, 4), int)
for a, bb in zip(c_o, c_m): conf[a, bb] += 1
print("confusion rows = measured class, cols = map class (<=5, 5-10, 10-30, >30):"); print(conf)
out["confusion"] = conf.tolist()
out["by_measured_class"] = {}
for i, (lab, lo, hi) in enumerate(cls):
    mm = c_o == i
    out["by_measured_class"][lab] = dict(n=int(mm.sum()), bias=float(np.mean(m_[mm] - o[mm])), rmse=float(np.sqrt(np.mean((m_[mm] - o[mm]) ** 2))), mae=float(np.mean(np.abs(m_[mm] - o[mm]))))
    print(" measured", lab, out["by_measured_class"][lab])
# decision-relevant rates
out["map_gt30_but_measured_le10"] = float(((c_m == 3) & (c_o <= 1)).sum() / max((c_m == 3).sum(), 1))
out["map_gt30_but_measured_le15"] = float(((m_ > 30) & (o <= 15)).sum() / max((m_ > 30).sum(), 1))
out["map_gt15_but_measured_le5"] = float(((m_ > 15) & (o <= 5)).sum() / max((m_ > 15).sum(), 1))
out["map_gt15_but_measured_le10"] = float(((m_ > 15) & (o <= 10)).sum() / max((m_ > 15).sum(), 1))
out["map_gt10_but_measured_le5"] = float(((m_ > 10) & (o <= 5)).sum() / max((m_ > 10).sum(), 1))
out["map_le5_and_measured_le5"] = float(((m_ <= 5) & (o <= 5)).sum() / max((m_ <= 5).sum(), 1))
out["map_le5_but_measured_gt10"] = float(((m_ <= 5) & (o > 10)).sum() / max((m_ <= 5).sum(), 1))
out["measured_le5_found_by_map_le5"] = float(((m_ <= 5) & (o <= 5)).sum() / max((o <= 5).sum(), 1))
out["measured_le5_found_by_map_le10"] = float(((m_ <= 10) & (o <= 5)).sum() / max((o <= 5).sum(), 1))
for k in list(out)[-9:]: print(" ", k, round(out[k], 4))
out["by_state"] = {}
for s in sorted(set(stt)):
    mm = stt == s
    if mm.sum() >= 50:
        out["by_state"][s] = dict(n=int(mm.sum()), bias=float(np.mean(m_[mm] - o[mm])), rmse=float(np.sqrt(np.mean((m_[mm] - o[mm]) ** 2))))
        print(" ", s, out["by_state"][s])
json.dump(out, open(sys.argv[1], "w"), indent=1)
