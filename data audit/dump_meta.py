"""Full metadata dump (no data reads) of every .nc under merged_datasets -> one text file each."""
import os
import sys

import netCDF4

from audit_paths import OUT, nc_files, rel, walk

META = OUT + "/meta"
os.makedirs(META, exist_ok=True)

for f in nc_files(sys.argv[1:]):
    name = rel(f).replace("/", "__")
    lines = [f"FILE {f}  size={os.path.getsize(f):,}"]
    try:
        ds = netCDF4.Dataset(f)
    except Exception as e:  # noqa: BLE001
        lines.append(f"OPEN FAILED: {type(e).__name__}: {e}")
        open(f"{META}/{name}.txt", "w", encoding="utf-8").write("\n".join(lines))
        print(name, "OPEN FAILED", e)
        continue
    lines.append(f"format={ds.data_model} disk_format={ds.disk_format}")
    for gpath, g in walk(ds):
        lines.append(f"######## GROUP {gpath}")
        lines.append("== attrs")
        for k in g.ncattrs():
            v = str(g.getncattr(k))
            lines.append(f"  :{k} = {v[:400]}")
        lines.append("== dims")
        for k, d in g.dimensions.items():
            lines.append(f"  {k} = {len(d)}{' (unlimited)' if d.isunlimited() else ''}")
        lines.append("== variables")
        for k, v in g.variables.items():
            try:
                ch = v.chunking()
                flt = v.filters()
            except Exception:  # noqa: BLE001
                ch, flt = "?", "?"
            lines.append(f"  {v.dtype} {k}{v.dimensions} shape={v.shape} chunks={ch} filters={flt}")
            for a in v.ncattrs():
                lines.append(f"      {k}:{a} = {str(v.getncattr(a))[:300]}")
    ds.close()
    open(f"{META}/{name}.txt", "w", encoding="utf-8").write("\n".join(lines))
    print(name, "ok", len(lines), "lines")
