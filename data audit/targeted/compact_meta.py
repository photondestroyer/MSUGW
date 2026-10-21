import sys, os, netCDF4
sys.path.insert(0, "G:/MSU_GWB/datasets/data audit")
from audit_paths import nc_files, rel, walk
KEEP = ("units","scale_factor","add_offset","_FillValue","missing_value","valid_range","valid_min","valid_max","flag_values","flag_meanings","calendar","long_name")
for f in nc_files(sys.argv[1:]):
    ds = netCDF4.Dataset(f)
    print(f"\n######## {rel(f)}  ({os.path.getsize(f)/1e6:,.1f} MB)")
    for gp, g in walk(ds):
        ga = {k: str(g.getncattr(k))[:160] for k in g.ncattrs()}
        if gp == "/":
            for k in ("Conventions","title","history","processing_status","time_coverage_start","time_coverage_end","committed_records","source","comment"):
                if k in ga: print(f"  :{k} = {ga[k]}")
            print("  dims:", {k: len(d) for k, d in g.dimensions.items()})
        else:
            print(f"  [group {gp}] dims={ {k: len(d) for k, d in g.dimensions.items()} } nvars={len(g.variables)}")
            if len(g.variables) > 12:
                print("     vars:", ", ".join(list(g.variables)[:40]))
                continue
        for k, v in g.variables.items():
            a = {kk: str(v.getncattr(kk))[:70] for kk in KEEP if kk in v.ncattrs()}
            print(f"   {str(v.dtype):8s} {k}{v.dimensions} {a}")
    ds.close()
