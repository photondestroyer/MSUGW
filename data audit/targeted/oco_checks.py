import netCDF4, numpy as np, datetime as dt
from collections import Counter
R="G:/MSU_GWB/datasets/merged_datasets/"
for fn in ("OCO2_L2_Lite_SIF_Ogallala.nc","OCO3_L2_Lite_SIF_Ogallala.nc"):
    ds=netCDF4.Dataset(R+fn); ds.set_auto_maskandscale(False)
    print("##",fn, "status:", ds.getncattr("processing_status") if "processing_status" in ds.ncattrs() else None)
    for a in ("time_coverage_start","time_coverage_end","date_modified","committed_days","source_missing_values"):
        if a in ds.ncattrs(): print("  ",a,"=",str(ds.getncattr(a))[:200])
    d=ds["day_time"][:].astype(int); days=np.array([dt.date(1970,1,1)+dt.timedelta(int(x)) for x in d])
    gaps=[(days[i],days[i+1],(days[i+1]-days[i]).days) for i in range(len(days)-1) if (days[i+1]-days[i]).days>1]
    print("  days",len(days),days[0],days[-1],"missing-day runs:",len(gaps),"missing days total:",sum(g[2]-1 for g in gaps))
    for g in sorted(gaps,key=lambda g:-g[2])[:10]: print("     gap",g[0],"->",g[1],g[2]-1,"days missing")
    ver=[str(x) for x in ds["day_collection_version"][:]]; print("  versions per day:",Counter(ver))
    nf=ds["day_n_soundings"][:]; first=ds["day_first_sounding"][:]
    ns=len(ds.dimensions["sounding"])
    print("  sum day_n_soundings =",int(nf.sum()),"sounding dim =",ns, "contiguous:", bool(np.all(first[1:]==first[:-1]+nf[:-1])))
    print("  days with 0 soundings kept:",int((nf==0).sum()), " soundings/day median", float(np.median(nf)))
    di=ds["day_index"][:]
    print("  day_index monotone:", bool(np.all(np.diff(di)>=0)))
    inp=ds["in_polygon"][:]; print("  in_polygon:",Counter(inp.tolist()))
    for q in ("Quality_Flag","SimplyGoodOrBadQualityFlag","Science_sounding_qual_flag","Cloud_cloud_flag_abp","Metadata_MeasurementMode","collection_version"):
        if q in ds.variables: print("  ",q,dict(sorted(Counter(ds[q][:].tolist()).items())))
    # per-year sounding counts
    yrs=Counter([days[i].year for i in di[::1]])
    print("  soundings per year:",dict(sorted(yrs.items())))
    dtm=ds["Delta_Time"][:]; t93=ds["Geolocation_time_tai93"][:]
    print("  Delta_Time sorted:", bool(np.all(np.diff(dtm)>=0)))
    sid=ds["Metadata_SoundingId"][:]; print("  SoundingId unique:", len(np.unique(sid))==len(sid), " sorted:", bool(np.all(np.diff(sid)>0)))
    ds.close()
