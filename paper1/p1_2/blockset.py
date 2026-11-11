"""Block sets: the fold-independent arrays of a set of blocks that P1.2 trains on or predicts for.

A block set holds, per block: identity and position, the driver cell, bucket bin, county, aridity decile,
static inputs (STATIC_INPUTS), and per calendar year the monthly targets Y (ET, NIRv, NDMI; NaN when the
block-month fails the TARGET_MIN_PIXEL_FRACTION coverage), the reference flag, the rule flags and the
dominant CDL group. Nothing in it depends on a fold; windows, calendars and scalers are built later.
Plan v3 adds: MEAS (block, year, month, MEAS_INPUTS) observation support from the block files; FUNC (functional
group), GRP_PREV (CDL group of the previous year), RING (reference_ring) per year; DIAG (block, DIAG_NAMES)
terrain / canal / irrigation distances and `state` (FIPS) - DIAG and state serve the verification checks only and
are refused as model inputs (config.EXPOSURE_VARS).

Stage `samples` stores a stratified sample of reference blocks (all their years):
    strata = aridity decile x the block's most frequent dominant group in its reference years
    each block is kept with probability min(1, SAMPLE_BLOCKS_PER_STRATUM / stratum size), decided by a
    hash of its block id (so the sample does not depend on tile order or reruns).
"""
import numpy as np

from .. import config as C
from .. import cdl_groups as CG
from ..core import ncio
from ..step0 import tiles as TL, static as ST, blocks as BK, drivers as DR
from ..p1_1 import build as RB

DIAG_NAMES = ["hand_mean", "hand_frac", "facc_max_km2", "canal_dist_m", "lanid_dist_m"]


def _hash_u(block_id):
    return ((np.asarray(block_id, "uint64") * np.uint64(2654435761)) % np.uint64(2 ** 32)).astype("float64") / 2 ** 32


def tile_blockset(T, G, ti, dirs, D, lut, select=None):
    t = T["tiles"][ti]
    g = G[t]
    sl = TL.tile_slice(T, ti)
    S = ST.open_tile(dirs["static"], t)
    R = RB.open_tile(dirs["reference"], t)
    n = g.nby * g.nbx
    sel = np.arange(n) if select is None else np.asarray(select)
    ya = C.agg_years()
    Y = np.full((sel.size, len(ya), 12, 3), np.nan, "float32")
    MEAS = np.full((sel.size, len(ya), 12, len(C.MEAS_INPUTS)), np.nan, "float32")
    for k, y in enumerate(ya):
        B = BK.open_tile_year(dirs["blocks"], t, y)
        if B is None:
            continue
        need = np.maximum(B["roi_px"].astype("float64") * C.TARGET_MIN_PIXEL_FRACTION, 1)
        for j, nm in enumerate(C.TARGETS):
            v = B[f"{nm}_mean"].astype("float32")
            v = np.where(B[f"{nm}_n"] >= need, v, np.nan)
            Y[:, k, :, j] = v.reshape(12, n)[:, sel].T
        for j, nm in enumerate(C.MEAS_INPUTS):
            if nm in B:
                MEAS[:, k, :, j] = B[nm].astype("float32").reshape(12, n)[:, sel].T
    iy, ix = T["lab_iy"][sl][sel], T["lab_ix"][sl][sel]
    cell = np.where(iy >= 0, lut[np.maximum(iy, 0), np.maximum(ix, 0)], -1)
    stat = np.stack([S[v].ravel()[sel].astype("float32") for v in C.STATIC_INPUTS], -1)
    diag = np.stack([S[v].ravel()[sel].astype("float32") if v in S else np.full(sel.size, np.nan, "float32")
                     for v in DIAG_NAMES], -1)
    geoid = ST.open_county_geoid(dirs["static"], t)
    cty = S["county"].ravel()[sel].astype("int64")
    state = np.array([int(geoid[c][:2]) if 0 <= c < len(geoid) else 0 for c in cty], "int16")
    ny = len(R["year"])
    yr2 = lambda name, fill: (R[name].reshape(ny, n)[:, sel].T.copy() if name in R
                              else np.full((sel.size, ny), fill, "uint8"))
    return dict(block_id=T["block_id"][sl][sel], tile_idx=np.full(sel.size, ti, "int16"), index_in_tile=sel.astype("int32"),
                lon=T["lon"][sl][sel], lat=T["lat"][sl][sel], valid=T["valid"][sl][sel], primary=T["primary"][sl][sel],
                cell=cell.astype("int32"), bucket_bin=S["bucket_bin"].ravel()[sel].astype("int8"),
                county=S["county"].ravel()[sel].astype("int32"),
                aridity=R["aridity_decile"].ravel()[sel].astype("int8"), static=stat, Y=Y,
                REF=R["reference"].reshape(len(R["year"]), n)[:, sel].T.copy(),
                REF_GDE=R["reference_gde"].reshape(len(R["year"]), n)[:, sel].T.copy(),
                FLAGS=R["flags"].reshape(len(R["year"]), n)[:, sel].T.copy(),
                GRP=R["group"].reshape(len(R["year"]), n)[:, sel].T.copy(),
                MEAS=MEAS, DIAG=diag, state=state, RING=yr2("reference_ring", 0), GRP_PREV=yr2("group_prev", 255),
                FUNC=yr2("functional", CG.FUNC_NONE))


def concat(parts):
    if not parts:
        return None
    return {k: np.concatenate([p[k] for p in parts], 0) for k in parts[0]}


def dominant_group(GRP, REF):
    out = np.full(GRP.shape[0], 255, "uint8")
    for i in range(GRP.shape[0]):
        g = GRP[i][REF[i] == 1]
        g = g[g != 255]
        if g.size:
            out[i] = np.bincount(g, minlength=CG.N).argmax()
    return out


def run(ctx):
    """Stage `samples`."""
    T, G = TL.load_table(ctx.dep("tiles"))
    dirs = dict(static=ctx.dep("static"), reference=ctx.dep("reference"), blocks=ctx.dep("blocks"))
    D = DR.load(ctx.dep("drivers"))
    lut = DR.cell_lookup(D)
    # pass 1: stratum sizes from the reference files
    sizes = {}
    cand = {}
    for ti, t in enumerate(T["tiles"]):
        R = RB.open_tile(dirs["reference"], t)
        sl = TL.tile_slice(T, ti)
        n = R["reference"].shape[1] * R["reference"].shape[2]
        ref = R["reference"].reshape(R["reference"].shape[0], n).T
        grp = R["group"].reshape(R["group"].shape[0], n).T
        ok = (ref.sum(1) > 0) & T["valid"][sl] & T["primary"][sl]
        idx = np.nonzero(ok)[0]
        dom = dominant_group(grp[idx], ref[idx])
        dec = R["aridity_decile"].ravel()[idx]
        cand[ti] = (idx, dec, dom)
        for d, gr in zip(dec, dom):
            sizes[(int(d), int(gr))] = sizes.get((int(d), int(gr)), 0) + 1
    cap = C.TEST_SAMPLE_BLOCKS_PER_STRATUM if C.RUN_MODE == "test" else C.SAMPLE_BLOCKS_PER_STRATUM
    prob = {k: min(1.0, cap / v) for k, v in sizes.items()}
    parts = []
    for ti, (idx, dec, dom) in cand.items():
        if idx.size == 0:
            continue
        u = _hash_u(T["block_id"][TL.tile_slice(T, ti)][idx])
        p = np.array([prob[(int(d), int(g))] for d, g in zip(dec, dom)])
        keep = idx[u < p]
        if keep.size:
            parts.append(tile_blockset(T, G, ti, dirs, D, lut, keep))
        ctx.unit(f"{T['tiles'][ti]}", "done", candidates=int(idx.size), kept=int(keep.size))
    BS = concat(parts)
    if BS is None:
        raise RuntimeError("no reference blocks to sample: check the P1.1 verification report")
    np.savez_compressed(f"{ctx.out}/samples.npz", years_agg=np.array(C.agg_years()), years=np.array(C.years()), **BS)
    ctx.complete(dict(blocks=int(BS["block_id"].size), strata={f"{d}|{CG.NAMES[g] if g < CG.N else 'none'}": v
                                                               for (d, g), v in sorted(sizes.items())},
                      reference_block_years=int(BS["REF"].sum())))


def load(stage_dir):
    with np.load(f"{stage_dir}/samples.npz") as z:
        return {k: z[k] for k in z.files}
