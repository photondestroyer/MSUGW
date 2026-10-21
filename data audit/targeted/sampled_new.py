import xarray as xr, numpy as np, dask
dask.config.set(scheduler="synchronous")
R="G:/MSU_GWB/datasets/merged_datasets/"
def summ(name, a, fill=None):
    a=np.asarray(a); m=np.isfinite(a) if a.dtype.kind=="f" else np.ones(a.shape,bool)
    if fill is not None: m&=(a!=fill)
    v=a[m]
    if not v.size: print(f"   {name:18s} no valid values"); return
    q=np.percentile(v,[0.1,1,50,99,99.9])
    print(f"   {name:18s} valid={m.mean():.4f} min={v.min():.4g} p0.1={q[0]:.4g} p1={q[1]:.4g} med={q[2]:.4g} p99={q[3]:.4g} p99.9={q[4]:.4g} max={v.max():.4g} n<0={int((v<0).sum())}")
e=xr.open_dataset(R+"ECO_L3T_JET_Ogallala_ET (1).nc", chunks={}, mask_and_scale=False)
idx=list(range(0,e.sizes["time"],10)); print("## ECO_L3T_JET (1): sampled", len(idx), "of", e.sizes["time"], "records")
roi=e["roi_mask"].values==1
for k in ("ETdaily","PTJPLSMinst","STICinst","MOD16inst","BESSinst","ETinstUncertainty"):
    a=e[k].isel(time=idx).values; summ(k,a)
    if k=="ETdaily":
        nv=np.isfinite(a).reshape(len(idx),-1).sum(1); print("   records with no valid ETdaily:", int((nv==0).sum()), "of", len(idx), "| valid outside roi_mask:", int((np.isfinite(a)&~roi).sum()))
L=xr.open_dataset(R+"LANDSAT_C2L2_indices_500m_p031r031_Ogallala.nc", chunks={}, mask_and_scale=False)
idx=list(range(0,L.sizes["time"],10)); print("## LANDSAT p031r031: sampled", len(idx), "of", L.sizes["time"])
for k in ("NDVI","EVI","EVI2","SAVI","kNDVI","NDMI","NMDI","MSI","NDWI","MNDWI","AWEInsh","AWEIsh","WI2015"):
    summ(k, L[k].isel(time=idx).values)
nc=L["n_clear"].isel(time=idx).values; nv=L["n_valid"].isel(time=idx).values; m=(nc!=65535)&(nv!=65535)
print("   n_clear > n_valid cells:", int((nc[m]>nv[m]).sum()), "| max n_valid", int(nv[m].max()) if m.any() else None, "(a 500 m cell holds about 278 30 m pixels)")
nd=L["NDVI"].isel(time=idx).values; print("   NDVI finite where n_clear==0 or fill:", int((np.isfinite(nd)&~((nc!=65535)&(nc>0))).sum()))
W=xr.open_dataset(R+"USGS_HPA_native/USGS_HPA_water_level_change_500m.nc", mask_and_scale=False)
am=W["aquifer_mask"].values==1
for k in ("dWL_predev_to_2019","dWL_2017_to_2019"):
    a=W[k].values; ok=a>-3e38; summ(k,a,fill=a.min() if (~ok).any() else None)
    print(f"      exactly 0 inside aquifer: {(a[am&ok]==0).mean():.3f} | valid outside aquifer: {int((ok&~am).sum())}")
for k in ("K_lower","K_upper"):
    a=W[k].values; summ(k,a,fill=-1)
print("   K_limits_conflict=1 cells:", int((W["K_limits_conflict"].values==1).sum()), "| K_lower>K_upper:", int(((W["K_lower"].values>W["K_upper"].values)&(W["K_lower"].values!=-1)).sum()))
