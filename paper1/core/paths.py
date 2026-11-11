"""Where every input and output lives, for the current RUN_ENV / RUN_MODE.

Nothing here reads data. Inputs come from the registry in config.DATASETS; outputs go to
<out_root>/<stage>/<signature>/. A registry path that does not exist is searched for by file name
under the data root (LOCATE_BY_NAME) so a file uploaded to a different Drive folder is still found -
the place it was found is logged and written into provenance.
"""
import os
import shutil
import zipfile

from .. import config as C
from .errors import SourceError


def _join(root, rel):
    return rel if os.path.isabs(rel) or rel.startswith("/") else f"{root}/{rel}"


def _locate_by_name(name, root, depth=3):
    """Breadth-first search for a file name under root (bounded depth; Drive FUSE listings are slow)."""
    level = [root]
    for _ in range(depth + 1):
        nxt = []
        for d in level:
            try:
                entries = list(os.scandir(d))
            except OSError:
                continue
            for e in entries:
                if e.is_file() and e.name == name:
                    return e.path.replace("\\", "/")
                if e.is_dir() and not e.name.startswith((".", "_test", "_focus")) and e.name != C.OUTPUT_FOLDER:
                    nxt.append(e.path)
        level = nxt
    return None


def dataset_path(key, must_exist=True):
    """Absolute path of a single-file dataset of the registry."""
    spec = C.DATASETS[key]
    if spec["kind"] != "file":
        raise ValueError(f"{key} is a tiled dataset; use tile_path()")
    if C.RUN_ENV == "local" and spec.get("local"):
        p = spec["local"]
    else:
        p = _join(C.data_root(), spec["drive"])
    p = p.replace("\\", "/")
    for bad in spec.get("refuse", []):
        if os.path.basename(p) == bad:
            raise SourceError(f"{key}: {bad} is refused (see config.DATASETS['{key}']['refuse'])")
    if os.path.exists(p) or not must_exist:
        return p
    if C.LOCATE_BY_NAME:
        found = _locate_by_name(os.path.basename(p), C.data_root())
        if found and os.path.basename(found) not in spec.get("refuse", []):
            return found
    raise SourceError(f"{key}: not found at {p}" + (" (nor by name under the data root)" if C.LOCATE_BY_NAME else ""))


def tile_dir(key):
    spec = C.DATASETS[key]
    return f"{C.data_root()}/{spec['tag']}{C.SOURCE_SUBFOLDER}"


def tile_path(key, tile, year=None):
    """<TAG>/<year>/<tile>/<TAG>_<tile>_<year>_Ogallala.nc (NWI: <TAG>/<tile>/<TAG>_<tile>_Ogallala.nc)."""
    spec = C.DATASETS[key]
    base = tile_dir(key)
    if spec["per_year"]:
        if year is None:
            raise ValueError(f"{key} needs a year")
        return f"{base}/{year}/{tile}/{spec['tag']}_{tile}_{year}_{C.ROI_LABEL}.nc"
    return f"{base}/{tile}/{spec['tag']}_{tile}_{C.ROI_LABEL}.nc"


def frames_path(key="landsat"):
    return f"{tile_dir(key)}/frames.json"


def aquifer_shapefile():
    """The hp_bound2010 polygon; on Colab the zip on Drive is unpacked to scratch once."""
    p = dataset_path("aquifer")
    if p.endswith(".shp"):
        return p
    dest = f"{C.scratch_root()}/aquifer"
    shp = f"{dest}/{C.DATASETS['aquifer']['member']}"
    if not os.path.exists(shp):
        os.makedirs(dest, exist_ok=True)
        with zipfile.ZipFile(p) as z:
            z.extractall(dest)
    if not os.path.exists(shp):
        raise SourceError(f"aquifer: {C.DATASETS['aquifer']['member']} not inside {p}")
    return shp


def stage_dir(stage, sig):
    return f"{C.out_root()}/{stage}/{sig}"


def shared_dir(*parts):
    d = "/".join([C.shared_root(), *parts])
    os.makedirs(d, exist_ok=True)
    return d


def scratch_dir(*parts):
    d = "/".join([C.scratch_root(), *parts])
    os.makedirs(d, exist_ok=True)
    return d


def stage_to_scratch(path, log=None):
    """Copy a big single file to the VM's local disk once (Drive FUSE random reads are slow).
    Locally the file is used where it is."""
    if C.RUN_ENV != "colab":
        return path
    dst = f"{scratch_dir('staged')}/{os.path.basename(path)}"
    st = os.stat(path)
    if os.path.exists(dst) and os.path.getsize(dst) == st.st_size:
        return dst
    if log:
        log.info("staging %s (%.2f GB) to scratch", path, st.st_size / 1e9)
    tmp = dst + ".part"
    shutil.copyfile(path, tmp)
    if os.path.getsize(tmp) != st.st_size:
        raise SourceError(f"staging {path}: size mismatch")
    os.replace(tmp, dst)
    return dst
