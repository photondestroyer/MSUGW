import netCDF4, numpy as np
R="G:/MSU_GWB/datasets/merged_datasets/"
def raw(ds,name):
    v=ds[name]; v.set_auto_maskandscale(False); return v
def compare(small, big, label):
    a=netCDF4.Dataset(R+small); b=netCDF4.Dataset(R+big)
    ta=raw(a,"time")[:]; tb=raw(b,"time")[:]
    ua, ub = a["time"].units, b["time"].units
    da=netCDF4.num2date(ta,ua,a["time"].calendar,only_use_cftime_datetimes=False,only_use_python_datetimes=True)
    db=netCDF4.num2date(tb,ub,b["time"].calendar,only_use_cftime_datetimes=False,only_use_python_datetimes=True)
    pos={d:i for i,d in enumerate(db)}
    idx=[pos.get(d,-1) for d in da]
    print(f"== {label}: {len(da)} steps in small file, matched in big: {sum(i>=0 for i in idx)}")
    print("   x equal:", np.array_equal(raw(a,'x')[:],raw(b,'x')[:]), " y equal:", np.array_equal(raw(a,'y')[:],raw(b,'y')[:]))
    for name,v in a.variables.items():
        if v.ndim!=3 or name not in b.variables: continue
        va=raw(a,name); vb=raw(b,name)
        ndiff=0; nsteps_diff=0; maxabs=0.0; first=None
        for k,i in enumerate(idx):
            if i<0: continue
            x=va[k]; y=vb[i]
            neq = ~((x==y)|(np.isnan(x)&np.isnan(y))) if x.dtype.kind=='f' else (x!=y)
            n=int(neq.sum())
            if n:
                ndiff+=n; nsteps_diff+=1; first=first or (k,str(da[k]))
                if x.dtype.kind=='f': maxabs=max(maxabs,float(np.nanmax(np.abs(x[neq]-y[neq]))))
        print(f"   {name:14s} differing cells: {ndiff:>10,}  steps with any diff: {nsteps_diff:>4}  max|diff|={maxabs:.4g} first={first}")
compare("GRIDMET_2016_Ogallala.nc","GRIDMET_Merged_Ogallala.nc","gridMET 2016 vs merged")
compare("modis_et_v5_dekadal_2020_Ogallala.nc","MODIS_ET_SSEBop_Merged_Ogallala.nc","SSEBop 2020 vs merged")
