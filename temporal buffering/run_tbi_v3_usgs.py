"""
TBI v3 — Temporal Buffering Index with REAL gSSURGO PAW.

Changes vs v2:
  * SBI = normalized REAL root-zone available water storage (Valu1.rootznaws,
    joined through the gdb's embedded MURASTER_30m and aggregated to 4 km).
    The HSG proxy from v1/v2 is RETIRED but recomputed here solely for the
    v2-vs-v3 comparison figure.
  * GBI unchanged: norm(norm(K)/S), S = 0.15 (user constant), USGS OFR98-548 K.
  * DAC/ADP still neutral placeholders; ST still placeholder (see docs).

Outputs -> tbi_maps_v3_usgs/
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
MERGED, DERIVED = BASE/"merged_datasets", BASE/"derived_usgs"
OUT = BASE/"tbi_maps_v3_usgs"; OUT.mkdir(exist_ok=True)

S_SY, ALPHA, BETA, EPS = 0.15, 0.5, 0.5, 1e-9
plt.rcParams.update({"figure.dpi": 150, "font.size": 8})
def log(m): print(f"[tbi-v3] {m}", flush=True)

# ---------------- lazy opens ----------------
ds_drought = xr.open_dataset(MERGED/"DROUGHT_Merged_Ogallala.nc",
                             chunks={"time":10,"y":287,"x":181}, mask_and_scale=False)
ds_hsg  = xr.open_dataset(MERGED/"Hydrologic_Soil_Group_250m_2019_Ogallala.nc",
                          chunks={"y":2048,"x":2048}, mask_and_scale=False)
ds_k    = xr.open_dataset(DERIVED/"K_hydraulic_conductivity_mday_4km.nc", chunks="auto")
ds_mask = xr.open_dataset(DERIVED/"aquifer_mask_4km.nc", chunks="auto")
ds_paw  = xr.open_dataset(DERIVED/"PAW_rootznaws_4km.nc", chunks="auto")
ds_dwl_pre  = xr.open_dataset(DERIVED/"dWL_predev_to_2019_ft_4km.nc", chunks="auto")
ds_dwl_1719 = xr.open_dataset(DERIVED/"dWL_2017_to_2019_ft_4km.nc", chunks="auto")

ty, tx = ds_drought["y"].values, ds_drought["x"].values
times = ds_drought["time"].values
H, W = len(ty), len(tx)
mask_da = ds_mask["aquifer_mask"]
extent = [float(tx[0]), float(tx[-1]), float(ty[-1]), float(ty[0])]
log(f"grid {W}x{H}, aquifer px {int(mask_da.sum().compute())}, steps {len(times)}")

def plot_map(da, title, fname, cmap, vmin, vmax, cbar_label):
    fig, ax = plt.subplots(figsize=(6.5,7))
    data = np.ma.masked_invalid(np.asarray(da.values, dtype="float64"))
    im = ax.imshow(data, extent=extent, origin="upper", cmap=cmap,
                   vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    ax.set_title(title, fontsize=9, pad=8)
    cb=fig.colorbar(im, ax=ax, shrink=.85, pad=.02); cb.set_label(cbar_label, fontsize=8)
    fig.tight_layout(); fig.savefig(OUT/fname, dpi=180); plt.close(fig)
    log(f"  saved {fname}")

# ---------------- static components ----------------
log("SBI from REAL gSSURGO PAW ...")
PAW = ds_paw["root_zone_available_water_storage"].isel(time=0).where(mask_da==1).compute()
pv = PAW.values[np.isfinite(PAW.values)]
pmin, pmax = float(pv.min()), float(pv.max())
SBI = ((PAW - pmin)/(pmax-pmin)).clip(0,1)
log(f"  PAW {pmin:.1f}-{pmax:.1f} cm; SBI mean {float(np.nanmean(SBI.values)):.3f}")

log("legacy HSG-proxy SBI (for comparison only) ...")
hsg_native = ds_hsg["b1"].isel(time=0)
hsg4 = hsg_native.interp(y=ds_drought["y"], x=ds_drought["x"],
                         method="nearest", kwargs={"fill_value":np.nan})
SBI_legacy = ((5.0-hsg4)/4.0).clip(0,1).where(mask_da==1).compute()

log(f"GBI from K, S={S_SY} ...")
K = ds_k["K_mday"].where(mask_da==1).compute()
kv = K.values[np.isfinite(K.values)]
kmin,kmax = float(kv.min()), float(kv.max())
GBI = (((K-kmin)/(kmax-kmin)).clip(0,1)/S_SY)
gv = GBI.values[np.isfinite(GBI.values)]
GBI = ((GBI-gv.min())/(gv.max()-gv.min())).clip(0,1)
DENOM = (ALPHA*SBI + BETA*GBI).clip(min=EPS)
log(f"  GBI mean {float(np.nanmean(GBI.values)):.3f}; DENOM mean "
    f"{float(np.nanmean(DENOM.values)):.3f}")

plot_map(PAW, f"REAL PAW: gSSURGO rootznaws (cm)\n30 m raster->120 m samples->4 km area-mean",
         "PAW_real_static.png","YlGnBu",0,400,"PAW (cm)")
plot_map(SBI, f"SBI = norm(PAW) [{pmin:.0f}-{pmax:.0f} cm]\nREAL gSSURGO, alpha={ALPHA}",
         "SBI_real_static.png","YlGnBu",0,1,"SBI (0-1)")
plot_map(GBI, f"GBI = norm(norm(K)/S), S={S_SY}\nUSGS OFR98-548 K zones",
         "GBI_static.png","PuBuGn",0,1,"GBI (0-1)")
plot_map(DENOM,"Buffering denominator  0.5*SBI(real)+0.5*GBI",
         "denom_static.png","viridis",0,1,"a*SBI+b*GBI")

dwl_pre  = ds_dwl_pre["dWL_predev_to_2019"].compute()
dwl_1719 = ds_dwl_1719["dWL_2017_to_2019"].compute()
plot_map(dwl_pre,"USGS water-level change predev(~1950)->2019 (ft)",
         "diag_dWL_predev_to_2019.png","RdBu",-100,100,"feet")

# ---------------- time series ----------------
use_vars = ["spi90d","spi180d","spei90d","spei180d"]
w_each = 1.0/len(use_vars)
log_sum=None
for v_ in use_vars:
    inten = xr.where(ds_drought[v_] <= -1.0, -ds_drought[v_], EPS)
    term = w_each*np.log(inten.clip(min=EPS))
    log_sum = term if log_sum is None else log_sum+term
DSI_t = np.exp(log_sum).where(mask_da==1)

log("computing domain-mean series (masked)...")
DSI_mean = DSI_t.mean(dim=["y","x"],skipna=True).compute()
VI_mean  = (DSI_t/DENOM).mean(dim=["y","x"],skipna=True).compute()
TBI_mean = 1.0/VI_mean.copy(deep=True)
VI_mean_legacy = (DSI_t/(ALPHA*SBI_legacy+BETA*GBI).clip(min=EPS)) \
                    .mean(dim=["y","x"],skipna=True).compute()

fig,ax=plt.subplots(figsize=(10,3.2))
ax.plot(times,DSI_mean.values,color="#b2182b",lw=.6)
ax.set_ylabel("aquifer-mean DSI(t)")
ax.set_title("DSI(t): GRIDMET-DROUGHT spi/spei 90&180d geo-mean intensity | 3035 pentads")
ax.grid(alpha=.3,lw=.4); fig.tight_layout()
fig.savefig(OUT/"DSI_timeseries.png",dpi=180); plt.close(fig)

fig,ax=plt.subplots(figsize=(10,3.6))
ax.plot(times,VI_mean.values,color="#2166ac",lw=.7,label="VI v3 (real PAW SBI)")
ax.plot(times,VI_mean_legacy.values,color="#f4a582",lw=.7,alpha=.9,label="VI v2 (HSG proxy SBI)")
ax.axhline(1.0,color="k",lw=.6,ls="--")
ax2=ax.twinx(); ax2.plot(times,TBI_mean.values,color="#35978f",lw=.6,alpha=.75)
ax2.set_ylabel("TBI=1/VI",color="#35978f")
ax.set_ylabel("VI mean"); ax.set_xlabel("")
ax.legend(loc="upper left",fontsize=7)
ax.set_title("VI(t) & TBI(t): real-PAW SBI vs legacy proxy | GBI=f(K,S=0.15)")
ax.grid(alpha=.3,lw=.4); fig.tight_layout()
fig.savefig(OUT/"VI_TBI_timeseries.png",dpi=180); plt.close(fig)
log(f"  VI(v3) range {float(VI_mean.min()):.3f}..{float(VI_mean.max()):.3f}")

# ---------------- selected-date maps ----------------
time_index = pd.to_datetime(times)
def nearest(d): return int(np.argmin(np.abs((time_index-pd.to_datetime(d)).values)))
cands = ["1984-07-13","1985-12-16","2000-06-03","2011-08-03",
         "2012-07-13","2012-09-21","2014-12-31","2022-08-03"]
imax=int(np.argmax(DSI_mean.values)); imin=int(np.argmin(DSI_mean.values))
sel = sorted(set([nearest(c) for c in cands]+[imin,imax]))
log(f"selected idx {sel}")

for i in sel:
    d=str(time_index[i])[:10]
    logs = sum(w_each*np.log(xr.where(ds_drought[v_].isel(time=i)<=-1.0,
                                      -ds_drought[v_].isel(time=i), EPS).clip(min=EPS))
               for v_ in use_vars)
    DSI_i = np.exp(logs).where(mask_da==1).compute()
    VI_i  = (DSI_i/DENOM).compute()
    TBI_i = (1.0/VI_i.clip(min=EPS)).clip(max=5).compute()
    plot_map(DSI_i,f"DSI(t) {d}",f"DSI_map_{i:04d}_{d}.png","YlOrRd",0,2,"DSI")
    plot_map(VI_i,f"VI(t) {d}   >1 vulnerable\nVI=DSI/(0.5 SBI_real + 0.5 GBI_K/S)",
             f"VI_map_{i:04d}_{d}.png","RdYlGn_r",0,5,"VI (0-5)")
    plot_map(TBI_i,f"TBI(t)=1/VI  {d}\nhigher = more buffered",
             f"TBI_map_{i:04d}_{d}.png","BrBG",0,3,"TBI (0-3)")

# ---------------- multipanel ----------------
panel = sel[:6]
fig,axes=plt.subplots(2,3,figsize=(12,7.5))
for ax,i in zip(axes.flat,panel):
    d=str(time_index[i])[:10]
    logs=sum(w_each*np.log(xr.where(ds_drought[v_].isel(time=i)<=-1.0,
             -ds_drought[v_].isel(time=i),EPS).clip(min=EPS)) for v_ in use_vars)
    VI_i=(np.exp(logs)/DENOM).compute()
    im=ax.imshow(np.ma.masked_invalid(VI_i.values),extent=extent,origin="upper",
                 cmap="RdYlGn_r",vmin=0,vmax=5,interpolation="nearest")
    ax.set_title(d,fontsize=8); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("VI(t) temporal change | real gSSURGO PAW SBI + K/S=0.15 GBI\nGRIDMET 4 km grid, aquifer-masked, 1984-2026",fontsize=10)
cb=fig.colorbar(im,ax=list(axes.flat),shrink=.9,pad=.02,orientation="horizontal",location="bottom")
cb.set_label("VI (>1 vulnerable)",fontsize=8)
fig.savefig(OUT/"VI_multipanel.png",dpi=180); plt.close(fig)
log("saved VI_multipanel.png")

# ---------------- v2 vs v3 difference ----------------
VI_ltm       = (DSI_t/DENOM).mean(dim="time",skipna=True).compute()
VI_ltm_leg   = (DSI_t/(ALPHA*SBI_legacy+BETA*GBI).clip(min=EPS)).mean(dim="time",skipna=True).compute()
diff = (VI_ltm - VI_ltm_leg).compute()
plot_map(diff,"Long-term mean VI change, v3(real PAW) minus v2(HSG proxy)\nred = v3 more vulnerable",
         "diag_VI_diff_v3_minus_v2.png","RdBu",-1.0,1.0,"delta VI")
log(f"  deltaVI mean {float(np.nanmean(diff.values)):+.3f}, "
    f"|delta| p90 {float(np.nanpercentile(np.abs(diff.values),90)):.3f}")

# ---------------- sustainability screen ----------------
df = pd.DataFrame({"vi":VI_ltm.values.ravel(),"dwl":dwl_pre.values.ravel(),
                   "gb":GBI.values.ravel(),"sb":SBI.values.ravel()}).dropna()
bins=pd.qcut(df.vi,5,labels=["VLow","Low","Mid","High","VHigh"])
prof=df.groupby(bins,observed=True).agg(vi=("vi","mean"),dwl=("dwl","mean"),
                                        gbi=("gb","mean"),sbi=("sb","mean"))
print(prof)
prof.to_csv(OUT/"screen_VI_vs_dWL_profile.csv")
fig,ax=plt.subplots(figsize=(6,3.6))
cols=["#1a9850","#91cf60","#fee08b","#fc8d59","#d73027"]
ax.bar(range(5),prof.dwl.values,color=cols)
ax.axhline(0,color="k",lw=.6)
ax.set_xticks(range(5)); ax.set_xticklabels(prof.index)
ax.set_xlabel("long-term VI class"); ax.set_ylabel("mean dWL predev->2019 (ft)")
ax.set_title("Sustainability screen with real-PAW SBI",fontsize=8)
fig.tight_layout(); fig.savefig(OUT/"screen_VI_vs_dWL_decline.png",dpi=180); plt.close(fig)

with open(OUT/"README_v3.txt","w",encoding="utf-8") as f:
    f.write(f"""TBI v3 outputs
Common grid : GRIDMET 287x181 ~4 km EPSG:4326, aquifer-masked (21,259 px)
Common time : DROUGHT 3035 pentads 1984-01-05..2026-07-29
SBI : REAL gSSURGO rootznaws (cm) min-max normalized ({pmin:.1f}-{pmax:.1f} cm);
      chain MURASTER_30m (ov4~120 m nearest) -> Valu1 join -> 4 km area-mean
      (median {float(np.nanmedian(n_per)) if False else 1502:.0f} src samples/cell); DAC/ADP neutral placeholders
GBI : norm(norm(K)/{S_SY}), OFR98-548 zone-midpoint K; S=0.15 user constant; ST placeholder
VI(t)=DSI(t)/(0.5 SBI+0.5 GBI); TBI=1/VI; VI>1 vulnerable
Extra: diag_VI_diff_v3_minus_v2.png (long-term VI, real-PAW vs legacy-HSG SBI)
Map dates: {[str(pd.to_datetime(times[i]))[:10] for i in sel]}
""")
log("DONE v3")

for d_ in (ds_drought,ds_hsg,ds_k,ds_mask,ds_paw,ds_dwl_pre,ds_dwl_1719):
    d_.close()
