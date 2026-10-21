import xarray as xr, numpy as np
ds=xr.open_dataset("G:/MSU_GWB/datasets/merged_datasets/GRIDMET_Ogallala_lab_latest.nc", chunks={}, mask_and_scale=False)
def day(v,d): return ds[v].sel(time=d).values
def same(v,d1,d2): a,b=day(v,d1),day(v,d2); return bool(np.array_equal(a,b))
print("old §1.3 duplicated days, now:")
print("  tmmn 2017-12-31 == 2018-12-31 == 2019-12-31:", same("tmmn","2017-12-31","2018-12-31"), same("tmmn","2018-12-31","2019-12-31"))
print("  tmmx 2017-12-31 == 2018-12-31:", same("tmmx","2017-12-31","2018-12-31"))
print("  tmmn 2010-12-31 == 2011-01-01:", same("tmmn","2010-12-31","2011-01-01"))
print("  srad 2013-12-30 == 2013-12-31:", same("srad","2013-12-30","2013-12-31"), " pr:", same("pr","2013-12-30","2013-12-31"))
print("  th 2011-03-01 == 2012-02-29:", same("th","2011-03-01","2012-02-29"), " vs:", same("vs","2011-03-01","2012-02-29"))
fill=32767
for v,d in (("srad","2014-12-22"),("th","2025-05-10"),("vs","2025-05-10"),("th","1981-12-31")):
    a=day(v,d); a=a[a!=fill]; print(f"  {v} {d}: unique raw values {np.unique(a)[:5]} (n={a.size})")
th=ds["th"].sel(time=slice("2011-01-01","2014-12-31")).values
print("  th raw min/max 2011-2014:", int(th[th!=fill].min()), int(th[th!=fill].max()))
pr=ds["pr"].sel(time="2015-10-23").values; m=pr!=fill; print("  pr max 2015-10-23 mm:", float(pr[m].max())*0.1)
