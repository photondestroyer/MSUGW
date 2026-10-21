"""Shared file discovery for the audit scripts: every NetCDF under merged_datasets, groups included."""
import glob
import os

ROOT = "G:/MSU_GWB/datasets/merged_datasets"
HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.environ.get("AUDIT_RUN", "2026-10-20")      # results of each audit go to their own folder
OUT = f"{HERE}/results/{RUN}"
os.makedirs(OUT, exist_ok=True)


def nc_files(only=()):
    files = sorted(glob.glob(ROOT + "/**/*.nc", recursive=True))
    return [f for f in files if not only or any(o in f for o in only)]


def raster_files(only=()):
    files = sorted(glob.glob(ROOT + "/**/*.tif", recursive=True) + glob.glob(ROOT + "/**/*.gpkg", recursive=True))
    return [f for f in files if not only or any(o in f for o in only)]


def walk(group, prefix=""):
    """Yield (path, group) for the root group and every subgroup, depth first."""
    yield prefix or "/", group
    for name, g in group.groups.items():
        yield from walk(g, f"{prefix}/{name}")


def get_var(ds, path):
    """Variable by '/group/sub/name' path (or plain name for the root group)."""
    parts = [p for p in path.split("/") if p]
    g = ds
    for p in parts[:-1]:
        g = g.groups[p]
    return g.variables[parts[-1]]


def rel(f):
    return os.path.relpath(f, ROOT).replace("\\", "/")
