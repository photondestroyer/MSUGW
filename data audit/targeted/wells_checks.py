import netCDF4, numpy as np, datetime as dt
from collections import Counter
R="G:/MSU_GWB/datasets/merged_datasets/"
ds=netCDF4.Dataset(R+"HPA_wells_water_level_timeseries.nc"); ds.set_auto_maskandscale(False)
for a in ds.ncattrs():
    if a not in ("title","summary"): print(" :",a,"=",str(ds.getncattr(a))[:300])
src=np.array(ds["record_source"][:]); rst=[str(x) for x in ds["record_source_text"][:]]
t=ds["time"][:]; yr=ds["year"][:]; dtw=ds["depth_to_water"][:]; si=ds["station_index"][:]
ss=ds["station_source"][:]; ssn=[str(x) for x in ds["source_name"][:]]
print("obs per record_source:", {rst[k]:int(v) for k,v in Counter(src.tolist()).items()})
print("station_index sorted:", bool(np.all(np.diff(si)>=0)), " stations with no obs:", len(ss)-len(np.unique(si)))
nodate=t==-9999
print("obs without date:", int(nodate.sum()), "by source", {rst[k]:int(v) for k,v in Counter(src[nodate].tolist()).items()})
# year from time vs year var
dates=netCDF4.num2date(t[~nodate],"days since 1900-01-01","standard")
ty=np.array([d.year for d in dates]); bad=ty!=yr[~nodate]
print("time-year != year var:", int(bad.sum()))
for name,m in (("year>2026",yr>2026),("year<1900",yr<1900)):
    idx=np.where(m)[0]; print(name, len(idx), [(int(yr[i]), rst[src[i]], float(dtw[i])) for i in idx[:8]])
print("time decoded extremes:", [str(x) for x in netCDF4.num2date([t[~nodate].min(), t[~nodate].max()],"days since 1900-01-01")])
early=np.where((~nodate)&(t<0))[0]; print("dates before 1900-01-01:", len(early), [(str(netCDF4.num2date(t[i],"days since 1900-01-01")), rst[src[i]]) for i in early[:5]])
print("depth fill:", int((dtw==-9999).sum()), " negative depths (above land surface):", int(((dtw<0)&(dtw!=-9999)).sum()),
      " ==-99:", int((dtw==-99).sum()), " >1000 ft:", int((dtw>1000).sum()), " >2000 ft:", int((dtw>2000).sum()))
for name,m in (("-99",dtw==-99),(">2000",dtw>2000)):
    idx=np.where(m)[0]; print("  ",name,[(rst[src[i]], float(dtw[i]), int(yr[i])) for i in idx[:6]])
wd=ds["well_depth"][:]; print("well_depth 99999:", int((wd==99999).sum()), " >5000:", int(((wd>5000)&(wd!=99999)).sum()), "max<99999:", float(wd[(wd<99999)].max()))
dup=ds["duplicate_of_other_source"][:]; print("duplicate_of_other_source=1:", int(dup.sum()))
inp=ds["in_polygon"][:]; print("stations in polygon:", Counter(inp.tolist()), " by source:", {ssn[k]:int(((ss==k)&(inp==1)).sum()) for k in range(3)})
lat=ds["lat"][:]; lon=ds["lon"][:]; print("lat range", lat.min(), lat.max(), " lon range", lon.min(), lon.max(), " lat>44:", int((lat>44).sum()))
# exact duplicate records within a source
key=np.stack([si.astype(np.int64), np.round(np.where(nodate,-1,t)).astype(np.int64), np.round(dtw*100).astype(np.int64), yr.astype(np.int64)],1)
u,c=np.unique(key,axis=0,return_counts=True); print("exact duplicate records (same well, date, depth, year):", int((c-1)[c>1].sum()))
# yearly coverage per source (dated + year-only)
for k in range(3):
    m=src==k; yc=Counter(yr[m].tolist())
    print(" ",rst[k], "years", min(yc), "-", max(k2 for k2 in yc if k2<2100), "| obs 2000-2015 per yr ~", int(np.mean([yc.get(y,0) for y in range(2000,2016)])), "| 2018-2025:", [yc.get(y,0) for y in range(2018,2026)])
