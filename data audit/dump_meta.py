"""Full metadata dump (no data reads) of every .nc under merged_datasets -> one text file each."""
import glob
import os
import sys
import netCDF4

ROOT = r"G:/MSU_GWB/datasets/merged_datasets"
OUT = os.path.dirname(os.path.abspath(__file__)) + "/meta"
os.makedirs(OUT, exist_ok=True)

files = sorted(glob.glob(ROOT + "/*.nc") + glob.glob(ROOT + "/USGS data/derived_usgs/**/*.nc", recursive=True))
for f in files:
    name = os.path.basename(f)
    lines = [f"FILE {f}  size={os.path.getsize(f):,}"]
    try:
        ds = netCDF4.Dataset(f)
    except Exception as e:  # noqa: BLE001
        lines.append(f"OPEN FAILED: {type(e).__name__}: {e}")
        open(f"{OUT}/{name}.txt", "w", encoding="utf-8").write("\n".join(lines))
        print(name, "OPEN FAILED", e)
        continue
    lines.append(f"format={ds.data_model} disk_format={ds.disk_format}")
    lines.append("== global attrs")
    for k in ds.ncattrs():
        v = str(ds.getncattr(k))
        lines.append(f"  :{k} = {v[:300]}")
    lines.append("== dims")
    for k, d in ds.dimensions.items():
        lines.append(f"  {k} = {len(d)}{' (unlimited)' if d.isunlimited() else ''}")
    lines.append("== variables")
    for k, v in ds.variables.items():
        try:
            ch = v.chunking()
            flt = v.filters()
        except Exception:  # noqa: BLE001
            ch, flt = "?", "?"
        lines.append(f"  {v.dtype} {k}{v.dimensions} shape={v.shape} chunks={ch} filters={flt}")
        for a in v.ncattrs():
            lines.append(f"      {k}:{a} = {str(v.getncattr(a))[:200]}")
    for g in ds.groups:
        lines.append(f"GROUP {g}")
    ds.close()
    open(f"{OUT}/{name}.txt", "w", encoding="utf-8").write("\n".join(lines))
    print(name, "ok", len(lines), "lines")
