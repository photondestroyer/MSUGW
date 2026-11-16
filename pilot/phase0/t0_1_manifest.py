"""Task 0.1  Manifest.

"List every existing data file with its units, period and product version. Archive the outputs of the earlier
residual-ET index."
Physical idea: none of the later checks mean anything unless every input is pinned to a version; Table 1 of
Paper 1 ("Data sources, versions and periods") is built from this manifest.
Decisions: none.   Deliverables: manifest.csv (one row per input), manifest.json (with variables and units),
archive_inventory.csv (earlier outputs, inventoried in place, not moved).
"""
import json
import os
import subprocess

TEXT_ATTRS = ("title", "summary", "source", "version", "product_version", "references", "time_coverage_start",
              "time_coverage_end", "processing_status", "institution")


def _nc_info(p):
    import netCDF4 as nc
    out = {}
    with nc.Dataset(p) as ds:
        out["attrs"] = {a: str(ds.getncattr(a))[:300] for a in ds.ncattrs() if a in TEXT_ATTRS}
        out["dims"] = {k: len(v) for k, v in ds.dimensions.items()}
        out["variables"] = {n: dict(units=str(getattr(v, "units", "")), long_name=str(getattr(v, "long_name", ""))[:120],
                                    dtype=str(v.dtype), dims=list(v.dimensions))
                            for n, v in list(ds.variables.items())[:60]}
    return out


def _tif_info(p):
    import rasterio
    with rasterio.open(p) as r:
        return dict(crs=str(r.crs), res=list(r.res), shape=[r.height, r.width], nodata=r.nodata,
                    dtype=r.dtypes[0], bounds=list(r.bounds))


def _folder_info(p):
    """Tiles / year folders: files per year (top-level digit folders) and total size."""
    years, total, n = {}, 0, 0
    for dp, _, files in os.walk(p):
        for f in files:
            if f.endswith((".nc", ".tif", ".png", ".txt", ".csv")):
                s = os.path.getsize(os.path.join(dp, f))
                total += s
                n += 1
                rel = os.path.relpath(dp, p).split(os.sep)[0]
                if rel.isdigit():
                    years[rel] = years.get(rel, 0) + 1
    return dict(files=n, bytes=total, files_per_year=dict(sorted(years.items())))


def _git_dates(path):
    try:
        out = subprocess.run(["git", "-C", os.path.dirname(path) or ".", "log", "--format=%ad", "--date=short", "--",
                              path], capture_output=True, text=True, timeout=60).stdout.split()
        return (out[-1], out[0]) if out else ("", "")
    except (OSError, subprocess.SubprocessError):
        return ("", "")


def run(ctx):
    import pandas as pd
    from .. import config as C
    from ..core import provenance as P
    from ..common import report

    # ══════ STEP 0.1-1  every source dataset: existence, identity, version, period, variables and units ══════
    rows, full = [], {}
    for key in C.DATASETS:
        p = C.src(key)
        row = dict(key=key, path=p, exists=os.path.exists(p))
        info = {}
        if row["exists"]:
            if os.path.isdir(p):
                info = _folder_info(p)
                row.update(kind="folder", bytes=info["bytes"], files=info["files"])
            else:
                fp = P.fingerprint(p)
                row.update(kind="file", bytes=fp["size"], mtime=fp["mtime"], sampled_md5=fp["sampled_md5"],
                           header_md5=fp["header"])
                try:
                    if p.endswith(".nc"):
                        info = _nc_info(p)
                        a = info["attrs"]
                        row.update(title=a.get("title", ""), period=f"{a.get('time_coverage_start', '')} to "
                                   f"{a.get('time_coverage_end', '')}".strip(" to"), status=a.get("processing_status", ""),
                                   units="; ".join(f"{k}: {v['units']}" for k, v in info["variables"].items()
                                                   if v["units"])[:400])
                    elif p.endswith((".tif", ".txt")) and ("bot" in p or p.endswith(".tif")):
                        info = _tif_info(p)
                        row.update(units=f"crs {info['crs']}; res {info['res']}; nodata {info['nodata']}")
                except Exception as e:
                    row["note"] = f"header not read: {type(e).__name__}: {e}"
        rows.append(row)
        full[key] = dict(row=row, info=info)
        ctx.log.info("STEP 0.1-1  %-20s %s", key, "ok" if row["exists"] else "MISSING")
    man = pd.DataFrame(rows)

    # ══════ STEP 0.1-2  archive inventory of earlier outputs (kept in place; listed with sizes and checksums) ══════
    arch = []
    for rel in C.EARLIER_OUTPUTS:
        d = os.path.join(C.REPO_ROOT, rel)
        if not os.path.isdir(d):
            arch.append(dict(folder=rel, file="", exists=False))
            continue
        for dp, _, files in os.walk(d):
            if "__pycache__" in dp:
                continue
            for f in files:
                p = os.path.join(dp, f)
                s = os.path.getsize(p)
                first, last = _git_dates(p)
                arch.append(dict(folder=rel, file=os.path.relpath(p, C.REPO_ROOT).replace("\\", "/"), exists=True,
                                 bytes=s, sampled_md5=P.sampled_md5(p), git_first=first, git_last=last))
    arch = pd.DataFrame(arch)
    ctx.log.info("STEP 0.1-2  %d files of earlier outputs inventoried", int(arch.get("exists", pd.Series()).sum()))

    # ══════ STEP 0.1-3  write ══════
    report.table(ctx, man, "manifest", units={"bytes": "bytes"})
    report.text(ctx, "manifest.json", json.dumps(full, indent=1, default=str))
    report.table(ctx, arch, "archive_inventory", units={"bytes": "bytes"},
                 notes="earlier residual-ET index outputs, inventoried in place (not moved or changed)")
    ctx.complete(dict(inputs=int(len(man)), missing=man.loc[~man.exists, "key"].tolist(),
                      archived_files=int(arch.get("exists", pd.Series(dtype=bool)).sum())))
