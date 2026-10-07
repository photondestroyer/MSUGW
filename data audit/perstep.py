import glob, os, sys, numpy as np, netCDF4
R="G:/MSU_GWB/datasets/merged_datasets/"
cache={}
def dates(fn):
    if fn in cache: return cache[fn]
    ds=netCDF4.Dataset(R+fn)
    tv=ds["time"]; t=tv[:]
    d=netCDF4.num2date(t, tv.units, getattr(tv,"calendar","standard"), only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    cache[fn]=d; return d
only=sys.argv[1:]
for p in sorted(glob.glob("results/perstep/*.npz")):
    base=os.path.basename(p)[:-4]; fn,var=base.split("__")
    if only and not any(o in fn for o in only): continue
    z=np.load(p); n=z["n"]; mean=z["mean"]; crc=z["crc"]
    if len(n)<3: continue
    try: d=dates(fn)
    except Exception: d=None
    npx=int(np.nanmax(n)) 
    lab=lambda i: str(d[i])[:10] if d is not None and len(d)==len(n) else str(i)
    out=[]
    empty=np.where(n==0)[0]
    if len(empty): out.append(f"allEmptySteps={len(empty)} first={lab(empty[0])} last={lab(empty[-1])}")
    med=np.median(n[n>0]) if (n>0).any() else 0
    part=np.where((n>0)&(n<0.5*med))[0]
    if len(part): out.append(f"partial(<50%median)={len(part)} e.g.{[lab(i) for i in part[:4]]}")
    # exact duplicate slabs (non-empty)
    nz=n>0
    _,first_idx,counts=np.unique(crc[nz],return_index=True,return_counts=True)
    dup_groups=int((counts>1).sum())
    if dup_groups:
        idx=np.where(nz)[0]
        c2=crc[nz]; vals,cnt=np.unique(c2,return_counts=True); dv=vals[cnt>1][:4]
        ex=[[lab(i) for i in idx[c2==v]][:4] for v in dv]
        out.append(f"duplicateSlabGroups={dup_groups} e.g.{ex}")
    consec=int(((crc[1:]==crc[:-1])&nz[1:]).sum())
    if consec: out.append(f"consecutiveIdentical={consec}")
    # level shift: z-score of step-to-step mean change
    m=mean[np.isfinite(mean)]
    if len(m)>10:
        dm=np.diff(mean); s=np.nanstd(dm)
        big=np.where(np.abs(dm)>8*s)[0] if s>0 else []
        if len(big) and len(big)<40: out.append(f"meanJumps>8sd={len(big)} e.g.{[lab(i+1) for i in big[:4]]}")
    vm=z["vmap"]
    if vm.ndim==2 and vm.size>1:
        full=(vm==vm.max()).mean(); zero=(vm==0).mean()
        out.append(f"pixelValidMap: allSteps={full:.3f} never={zero:.3f}")
    print(f"{fn[:24]:24s} {var[:26]:26s} steps={len(n)} px/step={npx} "+(" | ".join(out) if out else "ok"))
