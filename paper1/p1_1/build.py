"""Stage `reference` (P1.1): yearly rain-only reference blocks with one flag bit per rule.

Output per tile: reference/<tile>/P1_reference_<tile>.nc, dims (year, by, bx):
    flags (uint32, rules.BITS), reference, reference_gde, reference_irr3, reference_ring (rules.SETS),
    group (dominant CDL group, 255 none), group_prev (the same of the previous calendar year: a fallow year
    stores water the next crop uses), functional (cdl_groups.FUNC_NAMES), purity, valid_months, woody_px,
    depth_used_m, depth_class_m (rules.depth_rule); aridity_decile (by, bx) of the block's lab cell; county.
"""
import numpy as np

from .. import config as C
from ..core import ncio
from ..step0 import tiles as TL, static as ST, blocks as BK, classes as CL, drivers as DR
from . import rules as R, dryland as DY

PER_YEAR = R.SETS + ["flags", "group", "group_prev", "functional", "purity", "valid_months", "woody_px",
                     "depth_used_m", "depth_class_m"]


def block_aridity(T, sl, g, D):
    lut = DR.cell_lookup(D)
    iy, ix = T["lab_iy"][sl], T["lab_ix"][sl]
    cell = np.where(iy >= 0, lut[np.maximum(iy, 0), np.maximum(ix, 0)], -1)
    dec = np.where(cell >= 0, D["aridity_decile"][np.maximum(cell, 0)], 0)
    return dec.reshape(g.nby, g.nbx).astype("int8"), cell.reshape(g.nby, g.nbx)


def run(ctx):
    T, G = TL.load_table(ctx.dep("tiles"))
    sdir, bdir, cdir = ctx.dep("static"), ctx.dep("blocks"), ctx.dep("classes")
    dry = DY.lookup(DY.load(ctx.dep("dryland")))
    D = DR.load(ctx.dep("drivers"))
    years = C.years()
    totals = dict(block_years=0, **{k: 0 for k in R.SETS})
    for ti, t in enumerate(T["tiles"]):
        ctx.check_budget()
        p = f"{ctx.out}/{t}/P1_reference_{t}.nc"
        g = G[t]
        sl = TL.tile_slice(T, ti)
        if ncio.is_complete(p, ctx.sig):
            continue
        S = ST.open_tile(sdir, t)
        valid = T["valid"][sl].reshape(g.nby, g.nbx)
        primary = T["primary"][sl].reshape(g.nby, g.nbx)
        dec, _ = block_aridity(T, sl, g, D)
        out = {k: [] for k in PER_YEAR}
        prev = CL.open_tile_year(cdir, t, years[0] - 1)
        prev_group = R.class_of(prev["cdl_n"])[0] if prev is not None else np.full((g.nby, g.nbx), 255, "uint8")
        for y in years:
            B = BK.open_tile_year(bdir, t, y)
            if B is None:
                B = dict(water_px=np.zeros((12, g.nby, g.nbx), "uint16"),
                         **{f"{n}_n": np.zeros((12, g.nby, g.nbx), "uint16") for n in ("ET", "NIRv", "NDMI")})
            Cn = CL.open_tile_year(cdir, t, y)
            E = R.evaluate(S, B, Cn, valid, primary, y, dry, decile=dec)
            E["group_prev"] = prev_group
            prev_group = E["group"]
            for k in out:
                out[k].append(E[k])
            totals["block_years"] += int((valid & primary).sum())
            for k in R.SETS:
                totals[k] += int(E[k].sum())
        dims3 = ("year", "by", "bx")
        variables = {k: dict(dims=dims3, data=np.stack(v), dtype=str(np.stack(v).dtype)) for k, v in out.items()}
        variables["year"] = dict(dims=("year",), data=np.array(years, "int16"), dtype="i2")
        variables["aridity_decile"] = dict(dims=("by", "bx"), data=dec, dtype="i1")
        variables["county"] = dict(dims=("by", "bx"), data=S["county"].astype("int32"), dtype="i4")
        variables["flags"]["attrs"] = dict(flag_masks=np.array([1 << b for b in R.BITS.values()], "u4"),
                                           flag_meanings=" ".join(R.BITS), main_rules=" ".join(R.MAIN))
        variables["functional"]["attrs"] = dict(flag_values=np.arange(len(R.CG.FUNC_NAMES), dtype="u1"),
                                                flag_meanings=" ".join(R.CG.FUNC_NAMES))
        prov = ctx.record({}, dict(tile=t, depth_accepted=dict(C.DEPTH_ACCEPTED)))
        ncio.write_nc(p, dict(year=len(years), by=g.nby, bx=g.nbx), variables,
                      dict(tile=t, title="P1.1 rain-only reference blocks"), prov, ctx.scratch)
        ctx.unit(t, "done", reference=int(np.stack(out["reference"]).sum()),
                 ring=int(np.stack(out["reference_ring"]).sum()))
    ctx.complete(totals)


def open_tile(stage_dir, tile):
    with ncio.open_ro(f"{stage_dir}/{tile}/P1_reference_{tile}.nc") as ds:
        return {k: np.asarray(ds[k][...]) for k in ds.variables}
