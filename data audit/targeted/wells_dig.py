import netCDF4, numpy as np
from collections import Counter
ds=netCDF4.Dataset("G:/MSU_GWB/datasets/merged_datasets/HPA_wells_water_level_timeseries.nc"); ds.set_auto_maskandscale(False)
si=ds["station_index"][:]; dtw=ds["depth_to_water"][:]; t=ds["time"][:]; yr=ds["year"][:]
wid=ds["well_id"]; st=ds["state"]; wd=ds["well_depth"][:]
stt=[str(x) for x in ds["status_text"][:]]; status=ds["status"][:]
big=np.where(dtw>1000)[0]
print("depth>1000 ft wells:", Counter([wid[int(si[i])] for i in big]))
for w in sorted(set(int(si[i]) for i in big)):
    m=si==w; print("  ",wid[w], st[w], "well_depth", wd[w], "n obs", int(m.sum()), "depth range", float(dtw[m].min()), float(dtw[m].max()))
odd=np.where((t!=-9999)&(t<0))[0]
print("pre-1900 dates:", [(wid[int(si[i])], float(t[i]), int(yr[i]), float(dtw[i])) for i in odd])
y=np.where(yr>2026)[0]; print("year 3785:", [(wid[int(si[i])], int(yr[i]), float(t[i]), float(dtw[i])) for i in y])
neg=(dtw<0); print("negative depth status:", Counter([stt[k] for k in status[neg]]).most_common(6))
print("negative depth by source:", Counter([wid[int(s)].split(':')[0] for s in si[neg]]))
print("outside polygon stations by state:", Counter([st[i] for i in np.where(ds['in_polygon'][:]==0)[0]]).most_common(8))
