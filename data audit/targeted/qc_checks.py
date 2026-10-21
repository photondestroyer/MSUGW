import netCDF4, numpy as np
from collections import Counter
R="G:/MSU_GWB/datasets/merged_datasets/"
for fn in ("GRIDMET_Ogallala_lab_latest.nc","GRIDMET_DROUGHT_Ogallala_lab_colab.nc","GRIDMET_DROUGHT_Ogallala_lab.nc","GRIDMET_Ogallala_lab_drive.nc"):
    ds=netCDF4.Dataset(R+fn); ds.set_auto_maskandscale(False)
    t=netCDF4.num2date(ds["time"][:], ds["time"].units, ds["time"].calendar)
    st=ds["status"][:]; names=ds["status"].flag_meanings.split()
    print("##",fn, "| status:", {names[k]:int(v) for k,v in Counter(st.tolist()).items()})
    for k in sorted(set(st.tolist())):
        idx=np.where(st==k)[0]; print("    ",names[k], str(t[idx[0]])[:10], "->", str(t[idx[-1]])[:10])
    qcs=[k for k in ds.variables if k.endswith("_qc")]
    if qcs:
        v=ds[qcs[0]]; print("   qc attrs:", {a:str(v.getncattr(a))[:250] for a in v.ncattrs()})
    flagged={}
    for q in qcs:
        a=ds[q][:]; nz=np.nonzero(a)[0]
        if nz.size: flagged[q]=(int(nz.size), dict(Counter(a[nz].tolist())), [str(t[i])[:10] for i in nz[:6]])
    for q,v in flagged.items(): print("   ",q, v)
    extra=[a for a in ds.ncattrs() if a in ("qc_summary","known_source_issues","flag_notes","variables_missing","history")]
    for a in extra: print("   :",a,"=",str(ds.getncattr(a))[:400])
    print("   n data vars:", len([k for k,v in ds.variables.items() if v.ndim==3]))
    ds.close()
