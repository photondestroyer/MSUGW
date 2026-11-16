"""python -m pilot sync-inputs [--dry-run]: put what Colab needs on Google Drive (via the desktop mount G:).

Copies (local -> G:/My Drive/MSUGWB/...):
  * every input whose Colab path is under PILOT_inputs/ (KGS tables, PLSS, GMD, AIM-HPA 2011-2017, Haacker grids,
    the USGS surfaces, WWC5, the Mieno repository, the local SMAP fallback, the LANID 1997-2020 file and the
    LANID 2018-2020 CONUS maps); shapefiles with all sidecar files
  * the pilot code -> MSUGWB/code/pilot
Each file is copied to <name>.part, compared by size and sampled MD5, then renamed. A destination that already
matches is skipped. Sources are only read. --dry-run lists what would be copied and the total size.
"""
import glob
import os
import shutil

from . import config as C
from .core import provenance as P

AIM_YEARS = range(2011, 2018)
HAACKER = ["README.txt", "AqBase_m.tif", "AqBase_m.tfw", "AqBase_m.tif.aux.xml"] + \
          [f"WTelev_Clip_{y}.{e}" for y in range(2011, 2017) for e in ("tif", "tfw", "tif.aux.xml")]


def _pairs():
    out = []
    root = f"{C.LOCAL_DRIVE}/PILOT_inputs"
    for key, (local, colab, _) in C.DATASETS.items():
        if not colab.startswith("{inputs}") or not os.path.exists(local):
            continue
        dst = root + colab[len("{inputs}"):]
        if os.path.isdir(local):
            if key == "aim_dir":
                files = [f for f in os.listdir(local) if f[:4].isdigit() and int(f[:4]) in AIM_YEARS]
            elif key == "haacker_dir":
                files = [f for f in HAACKER if os.path.exists(os.path.join(local, f))]
            elif key == "lanid_1820_dir":                      # the yearly maps only (not maxIrrExt)
                files = [f for f in os.listdir(local) if f.startswith("lanid") and f.endswith(".tif")]
            else:
                files = [os.path.relpath(os.path.join(dp, f), local) for dp, ds, fs in os.walk(local)
                         if ".git" not in dp.split(os.sep) for f in fs]
            out += [(os.path.join(local, f), os.path.join(dst, f)) for f in files]
        elif local.endswith(".shp"):
            stem = local[:-4]
            out += [(p, os.path.join(os.path.dirname(dst), os.path.basename(p))) for p in glob.glob(stem + ".*")]
        else:
            out.append((local, dst))
    code_dst = f"{C.LOCAL_DRIVE}/code/pilot"
    for dp, ds, fs in os.walk(C.PKG_DIR):
        if "__pycache__" in dp:
            continue
        for f in fs:
            if f.endswith((".py", ".yaml", ".R", ".md", ".ipynb")):
                p = os.path.join(dp, f)
                out.append((p, os.path.join(code_dst, os.path.relpath(p, C.PKG_DIR))))
    return [(a.replace("\\", "/"), b.replace("\\", "/")) for a, b in out]


def _same(a, b):
    return os.path.exists(b) and os.path.getsize(a) == os.path.getsize(b) and P.sampled_md5(a) == P.sampled_md5(b)


def run(dry_run=True):
    if C.RUN_ENV != "local" or not os.path.isdir(C.LOCAL_DRIVE):
        raise SystemExit(f"sync-inputs runs locally with Google Drive mounted at {C.LOCAL_DRIVE}")
    pairs = _pairs()
    todo = [(a, b) for a, b in pairs if not _same(a, b)]
    gb = sum(os.path.getsize(a) for a, _ in todo) / 1e9
    print(f"{len(pairs)} files, {len(todo)} to copy ({gb:.2f} GB)")
    if dry_run:
        for a, b in todo[:40]:
            print(f"  {a}\n    -> {b}")
        if len(todo) > 40:
            print(f"  ... and {len(todo) - 40} more")
        return
    for i, (a, b) in enumerate(todo, 1):
        os.makedirs(os.path.dirname(b), exist_ok=True)
        part = b + ".part"
        shutil.copyfile(a, part)
        if os.path.getsize(part) != os.path.getsize(a) or P.sampled_md5(part) != P.sampled_md5(a):
            os.remove(part)
            raise SystemExit(f"copy did not verify: {a}")
        os.replace(part, b)
        if i % 20 == 0 or i == len(todo):
            print(f"  {i}/{len(todo)} copied")
