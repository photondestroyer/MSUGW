"""P1.1 dryland rule: which crops are mostly grown dryland in each county and year.

A crop group is dryland-eligible in (county, year) when, among its CDL pixels inside LANID's coverage,
the share OUTSIDE LANID's ever-irrigated area (1997-2017, derived/maxIrrExt) is >= DRYLAND_SHARE_MIN,
and there are at least DRYLAND_MIN_PIXELS such pixels. Counts come from primary, valid blocks only, so
the zone overlap is not counted twice.
Deviation from the plan (recorded in docs): LANID stands in for well and point-of-diversion registries;
LANID is derived from Landsat greenness and ends in 2017.
In test mode only the processed window enters the counts, so shares describe the window, not the county.
"""
import numpy as np

from .. import config as C
from .. import cdl_groups as CG
from ..step0 import tiles as TL, static as ST, classes as CL


def table(T, G, static_dir, classes_dir, years):
    tiles = list(T["tiles"])
    acc = {}
    for ti, t in enumerate(tiles):
        sl = TL.tile_slice(T, ti)
        g = G[t]
        use = (T["valid"][sl] & T["primary"][sl]).reshape(g.nby, g.nbx)
        county = ST.open_tile(static_dir, t)["county"]
        for y in years:
            Cn = CL.open_tile_year(classes_dir, t, y)
            if Cn is None:
                continue
            for name in ("cdl_n_irr", "cdl_n_assessed"):
                a = Cn[name][:, use]                 # (group, n_blocks)
                cty = county[use]
                for gi in CG.CROP_GROUPS:
                    s = np.bincount(np.maximum(cty, 0), weights=a[gi] * (cty >= 0),
                                    minlength=max(int(cty.max()) + 1, 1) if cty.size else 1)
                    for c in np.nonzero(s)[0]:
                        key = (int(c), int(y), int(gi))
                        acc.setdefault(key, [0, 0])
                        acc[key][0 if name == "cdl_n_irr" else 1] += int(s[c])
    keys = sorted(acc)
    out = dict(county=np.array([k[0] for k in keys], "int32"), year=np.array([k[1] for k in keys], "int16"),
               group=np.array([k[2] for k in keys], "int16"),
               n_irr=np.array([acc[k][0] for k in keys], "int64"), n_assessed=np.array([acc[k][1] for k in keys], "int64"))
    with np.errstate(invalid="ignore", divide="ignore"):
        out["share_outside"] = np.where(out["n_assessed"] > 0, 1 - out["n_irr"] / out["n_assessed"], np.nan)
    out["eligible"] = (out["share_outside"] >= C.DRYLAND_SHARE_MIN) & (out["n_assessed"] >= C.DRYLAND_MIN_PIXELS)
    return out


def lookup(D):
    """dict (county, year, group) -> eligible."""
    return {(int(c), int(y), int(g)): bool(e) for c, y, g, e in zip(D["county"], D["year"], D["group"], D["eligible"])}


def run(ctx):
    """Stage `dryland`."""
    T, G = TL.load_table(ctx.dep("tiles"))
    D = table(T, G, ctx.dep("static"), ctx.dep("classes"), C.years())
    np.savez_compressed(f"{ctx.out}/dryland.npz", **D)
    summary = dict(rows=int(D["county"].size), eligible=int(D["eligible"].sum()),
                   by_group={CG.NAMES[g]: int(D["eligible"][D["group"] == g].sum()) for g in CG.CROP_GROUPS},
                   note="test mode: window only" if C.RUN_MODE == "test" else "")
    ctx.complete(summary)


def load(stage_dir):
    with np.load(f"{stage_dir}/dryland.npz") as z:
        return {k: z[k] for k in z.files}
