"""CDL class-group pixel counts per block and year (for P1.1 class eligibility and dryland shares).

Per block-year, for every group of cdl_groups.GROUPS (aquifer pixels only):
    cdl_n          pixels of the group
    cdl_n_irr      ... of which LANID ever-irrigated (derived/maxIrrExt == 1)
    cdl_n_assessed ... of which inside LANID's coverage (ever-irrigated known)
The LANID masks come from the `static` stage (same 30 m tile lattice, exact centre lookup).
"""
import os

import numpy as np

from .. import config as C
from .. import cdl_groups as CG
from ..core import paths, ncio
from . import tiles as TL
from . import static as ST


def counts(cdl, roi, ever, assessed, nby, nbx):
    grp = CG.group_of(cdl)
    nb = nby * nbx
    pix_block = np.arange(nb, dtype="int64").reshape(nby, nbx).repeat(C.BLOCK_PX, 0).repeat(C.BLOCK_PX, 1)
    ok = roi & (grp != 255)
    key = pix_block[ok] * CG.N + grp[ok]
    out = {}
    for name, m in (("cdl_n", None), ("cdl_n_irr", ever[ok]), ("cdl_n_assessed", assessed[ok])):
        k = key if m is None else key[m]
        c = np.bincount(k, minlength=nb * CG.N).reshape(nby, nbx, CG.N)
        out[name] = np.moveaxis(c, -1, 0).astype("uint16")
    return out


def run(ctx):
    """Stage `classes`."""
    from ..core import provenance as P
    T, G = TL.load_table(ctx.dep("tiles"))
    sdir = ctx.dep("static")
    missing = []
    for t in sorted(G):
        g = G[t]
        ever, assessed = ST.open_lanid(sdir, t)
        roi, _ = ST.roi_window(g)
        for year in C.years():
            ctx.check_budget()
            p = f"{ctx.out}/{year}/{t}/P1_classes_{t}_{year}.nc"
            if ncio.is_complete(p, ctx.sig):
                continue
            src = paths.tile_path("cdl", t, year)
            if not os.path.exists(src):
                missing.append(f"{t}/{year}")
                ctx.unit(f"{t}/{year}", "missing", note=src)
                continue
            with ncio.open_ro(src) as ds:
                cdl = ncio.read_window(ds["cropland"], g.frame_rows(), g.frame_cols(), 255)
            V = counts(cdl, roi, ever, assessed, g.nby, g.nbx)
            variables = {k: dict(dims=("group", "by", "bx"), data=v, dtype="u2") for k, v in V.items()}
            variables["group_name"] = dict(dims=("group",), data=np.array(CG.NAMES, dtype=object), dtype="str")
            prov = ctx.record({"cdl": P.fingerprint(src)}, dict(tile=t, year=year))
            ncio.write_nc(p, dict(group=CG.N, by=g.nby, bx=g.nbx), variables,
                          dict(tile=t, year=year, title="CDL class-group pixel counts per block"), prov, ctx.scratch)
            ctx.unit(f"{t}/{year}", "done")
    ctx.complete(dict(missing=missing))


def open_tile_year(stage_dir, tile, year):
    p = f"{stage_dir}/{year}/{tile}/P1_classes_{tile}_{year}.nc"
    if not os.path.exists(p):
        return None
    with ncio.open_ro(p) as ds:
        return {k: np.asarray(ds[k][...]).astype("int64") for k in ("cdl_n", "cdl_n_irr", "cdl_n_assessed")}
