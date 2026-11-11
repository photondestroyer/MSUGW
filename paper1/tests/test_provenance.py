"""Provenance tests: do the outputs trace back to the source files in the configured data root?

Run against whatever RUN_ENV / RUN_MODE config.py (or P1_RUN_ENV / P1_RUN_MODE) selects - locally that is
LOCAL_DATA_ROOT with the tile-years you downloaded; in Colab, Drive. Every test SKIPs (saying what is
missing) when the files or stage outputs it needs are not there.
    (a) source identity       inputs are the files the pipeline was written for (headers, types, scales, grids)
    (b) output -> input chain outputs cite inputs whose fingerprints match the files here, and were made by
                              the current code and settings (staleness)
    (c) value tracing         random blocks / cells recomputed by tests/oracle.py (no paper1 imports)
    (d) flags and samples     reference bits recompute from block statistics; training blocks are reference
    (e) determinism           rebuilding a tile-year gives identical arrays
    mutation                  a wrong scale factor in the oracle is caught
"""
import glob
import json
import os

import numpy as np

from .. import config as C
from ..core import ncio, provenance as P
from . import oracle as O
from .util import stage_dir

N_BLOCKS = 12
RNG = np.random.default_rng(11)


def _need(*stages):
    dirs = {s: stage_dir(s) for s in stages}
    missing = [s for s, d in dirs.items() if d is None]
    return dirs, (f"SKIP stage(s) {missing} not complete for RUN_ENV={C.RUN_ENV} RUN_MODE={C.RUN_MODE} "
                  f"under {C.out_root()}" if missing else None)


def _table(tiles_dir):
    from ..step0 import tiles as TL
    return TL.load_table(tiles_dir)


def _sample_blocks(T, n=N_BLOCKS):
    idx = np.nonzero(T["valid"])[0]
    return RNG.choice(idx, min(n, idx.size), replace=False) if idx.size else idx


def _pos(T, G, k):
    t = str(T["tiles"][T["tile_idx"][k]])
    g = G[t]
    return t, g, int(T["br"][k]), int(T["bc"][k]), int(T["br"][k]) - g.br0, int(T["bc"][k]) - g.bc0


# ------------------------------------------------------------------ (a)
def test_sources_are_the_expected_files():
    from .. import checks
    from ..core import paths
    from ..core.errors import SourceError
    problems, found = [], 0
    for key, spec in C.DATASETS.items():
        if spec["kind"] != "file" or key not in checks.FILE_CHECKS:
            continue
        try:
            p = paths.dataset_path(key)
        except SourceError:
            continue
        found += 1
        problems += checks.FILE_CHECKS[key](p)
    try:
        from ..step0 import tiles as TL
        frames = TL.load_frames()
        for r in checks.inventory(TL.selected_tiles(frames), C.agg_years()):
            if r["exists"]:
                found += 1
                problems += checks.check_tile_file(r["dataset"], r["path"], frames.get(r["tile"]))
    except SourceError:
        pass
    if not found:
        return f"SKIP no source file found under {C.data_root()}"
    assert not problems, problems


def test_refused_files_are_refused():
    from ..core import paths
    from ..core.errors import SourceError
    import copy
    old = copy.deepcopy(C.DATASETS)
    try:
        C.DATASETS["hres_wtd"]["drive"] = "x/HRES-WTD_2015_Ogallala.nc"
        C.DATASETS["hres_wtd"]["local"] = "x/HRES-WTD_2015_Ogallala.nc"
        try:
            paths.dataset_path("hres_wtd", must_exist=False)
            raise AssertionError("the stretched HRES-WTD export was accepted")
        except SourceError:
            pass
    finally:
        C.DATASETS = old


# ------------------------------------------------------------------ (b)
def test_outputs_cite_inputs_that_match_the_files_here():
    from ..pipeline import STAGES
    names = ["blocks", "drivers", "static", "classes", "reference", "hand"]
    dirs, skip = _need(*names)
    if skip:
        return skip
    stage = {s["name"]: s for s in STAGES}
    checked = 0
    for n, d in dirs.items():
        files = glob.glob(f"{d}/**/*.nc", recursive=True)
        for f in files[:20]:
            prov = ncio.json_attr(f, "provenance")
            assert prov["signature"] == os.path.basename(d), f"{f}: made under another signature"
            assert prov["code_hash"] == P.code_hash(stage[n]["modules"]), f"{f}: made by different code (stale)"
            for key, fp in prov["inputs"].items():
                if isinstance(fp, dict) and "sampled_md5" in fp and os.path.exists(fp["path"]):
                    assert P.same_file(fp, P.fingerprint(fp["path"])), f"{f}: input {key} changed since ({fp['path']})"
                    checked += 1
    assert checked > 0


# ------------------------------------------------------------------ (c)
def test_block_statistics_trace_to_the_tile_files():
    dirs, skip = _need("tiles", "blocks")
    if skip:
        return skip
    from ..core import paths
    from ..step0 import blocks as BK
    T, G = _table(dirs["tiles"])
    n = 0
    for k in _sample_blocks(T):
        t, g, br, bc, i, j = _pos(T, G, k)
        y = int(RNG.choice(C.agg_years()))
        m = int(RNG.integers(1, 13))
        B = BK.open_tile_year(dirs["blocks"], t, y)
        if B is None:
            continue
        lp, op = paths.tile_path("landsat", t, y), paths.tile_path("openet", t, y)
        # the oracle finds the block from the block centre's lon/lat, not from the table's indices
        assert O.block_of_lonlat(t, float(T["lon"][k]), float(T["lat"][k])) == (br, bc)
        for var, path in (("NIRv", lp), ("NDMI", lp), ("et_ensemble_mad", op)):
            if not os.path.exists(path):
                continue
            v = O.tile_block_values(path, t, br, bc, y, m, var)
            if v is None:
                continue
            v = v[np.isfinite(v)]
            name = "ET" if var.startswith("et_") else var
            assert int(B[f"{name}_n"][m - 1, i, j]) == v.size, (t, y, m, name, B[f"{name}_n"][m - 1, i, j], v.size)
            if v.size:
                assert abs(float(B[f"{name}_mean"][m - 1, i, j]) - v.mean()) <= 1e-4 * max(1, abs(v.mean())), (t, y, m, name)
                n += 1
    if n == 0:
        return "SKIP no block-month with data among the sampled blocks"


def test_monthly_drivers_trace_to_daily_gridmet():
    dirs, skip = _need("drivers")
    if skip:
        return skip
    from ..core import paths
    from ..step0 import drivers as DR
    D = DR.load(dirs["drivers"])
    src = paths.dataset_path("gridmet")
    for _ in range(6):
        c = int(RNG.integers(0, D["cell_iy"].size))
        k = int(RNG.integers(12, D["months"].size))
        mo = D["months"][k]
        y, m = int(str(mo)[:4]), int(str(mo)[5:7])
        iy, ix = int(D["cell_iy"][c]), int(D["cell_ix"][c])
        for out, var in (("P", "pr"), ("PET", "pet"), ("TMIN", "tmmn"), ("VPD", "vpd")):
            want = O.lab_month(src, var, iy, ix, y, m)
            got = float(D[out][k, c])
            assert (np.isnan(want) and np.isnan(got)) or abs(got - want) <= 1e-3 * max(1, abs(want)), (out, y, m, got, want)
        want = O.lab_window_sum(src, iy, ix, f"{y:04d}-{m:02d}-01", 90)
        got = float(D["P90d"][k, c])
        assert (np.isnan(want) and np.isnan(got)) or abs(got - want) <= 1e-3 * max(1, want), ("P90d", got, want)


def test_static_layers_trace_to_their_sources():
    dirs, skip = _need("tiles", "static")
    if skip:
        return skip
    from ..core import paths
    from ..step0 import static as ST, tiles as TL
    T, G = _table(dirs["tiles"])
    hres, dem, lanid = paths.dataset_path("hres_wtd"), paths.dataset_path("dem"), paths.dataset_path("lanid")
    n = 0
    for k in _sample_blocks(T, 8):
        t, g, br, bc, i, j = _pos(T, G, k)
        S = ST.open_tile(dirs["static"], t)
        lp = TL.first_source_file("landsat", t, C.agg_years())
        roi = O.roi_of_block(lp, t, br, bc)
        if not roi.any():
            continue
        wmin, wmean, cov = O.wtd_block(hres, t, br, bc, roi)
        if np.isfinite(wmin):
            assert abs(S["wtd_min"][i, j] - wmin) < 1e-4 and abs(S["wtd_mean"][i, j] - wmean) < 1e-3, (t, br, bc)
        assert abs(S["wtd_coverage"][i, j] - cov) < 1e-6
        z = O.dem_block_mean(dem, t, br, bc, roi)
        assert abs(S["elevation"][i, j] - z) < 1e-2, (S["elevation"][i, j], z)
        d = O.lanid_block_distance(lanid, t, br, bc, roi)
        got = float(S["lanid_dist_m"][i, j])
        if got <= C.IRR_BUFFER_M:
            assert abs(got - d) < 1e-3, (got, d)
        else:
            assert d >= got - 1e-3 or not np.isfinite(d)
        n += 1
    if n == 0:
        return "SKIP no sampled block inside the aquifer"


def test_class_counts_trace_to_cdl():
    dirs, skip = _need("tiles", "classes")
    if skip:
        return skip
    from ..core import paths
    from ..step0 import classes as CL, tiles as TL
    from .. import cdl_groups as CG
    T, G = _table(dirs["tiles"])
    for k in _sample_blocks(T, 6):
        t, g, br, bc, i, j = _pos(T, G, k)
        y = int(RNG.choice(C.years()))
        Cn = CL.open_tile_year(dirs["classes"], t, y)
        if Cn is None:
            continue
        (X, Y), _ = O.block_pixel_centres(t, br, bc)
        roi = O.roi_of_block(TL.first_source_file("landsat", t, C.agg_years()), t, br, bc)
        with O.open_raw(paths.tile_path("cdl", t, y)) as ds:
            r, c = O.file_rows_cols(ds, X[0], Y[:, 0])
            codes = np.full((16, 16), 255)
            for a in range(16):
                for b in range(16):
                    if r[a] >= 0 and c[b] >= 0:
                        codes[a, b] = int(ds["cropland"][r[a], c[b]])
        valid = roi & (codes != 255)
        assert Cn["cdl_n"][:, i, j].sum() == valid.sum()
        assert Cn["cdl_n"][CG.NAMES.index("grass_pasture"), i, j] >= ((codes == 176) & roi).sum()
        assert Cn["cdl_n"][CG.NAMES.index("winter_wheat"), i, j] >= ((codes == 24) & roi).sum()


def test_hand_follows_the_stored_d8_directions():
    dirs, skip = _need("hand")
    if skip:
        return skip
    files = sorted(glob.glob(f"{dirs['hand']}/hu8/HAND_*.nc"))
    if not files:
        return "SKIP no HAND files"
    traced = 0
    for f in files[:3]:
        with O.open_raw(f) as ds:
            fd = np.asarray(ds["fdir"][:])
            h = O.decoded(ds["hand"], np.asarray(ds["hand"][:]))
            ch = np.asarray(ds["channel"][:]).astype(bool) & np.asarray(ds["inside"][:]).astype(bool)
        if not (fd > 0).any():
            return "SKIP HAND files were made without D8 directions (pysheds stand-in)"
        hc = h[ch & np.isfinite(h)]
        assert hc.size == 0 or np.mean(np.abs(hc) < 1.0) >= 0.95, "S0.6 verify (1): channel cells must have HAND < 1 m"
        cand = np.argwhere(np.isfinite(h) & (fd > 0))
        for r, c in cand[RNG.choice(len(cand), min(20, len(cand)), replace=False)]:
            got, want = O.hand_trace(f, int(r), int(c))
            if np.isfinite(want):
                assert abs(got - want) <= 0.051 + 1e-3, (f, r, c, got, want)      # stored at 0.1 m
                traced += 1
    assert traced > 0


# ------------------------------------------------------------------ (d)
def test_reference_flags_recompute_from_block_statistics():
    dirs, skip = _need("tiles", "static", "reference")
    if skip:
        return skip
    from ..step0 import static as ST
    from ..p1_1 import build as RB, rules as RU
    T, G = _table(dirs["tiles"])
    for t in T["tiles"]:
        S = ST.open_tile(dirs["static"], str(t))
        R = RB.open_tile(dirs["reference"], str(t))
        f = R["flags"].astype("int64")
        with np.errstate(invalid="ignore"):
            cov = S["wtd_coverage"] >= C.WTD_MIN_COVERAGE
            wtd = (S["wtd_min"][None] >= R["depth_used_m"]) & cov[None]           # the depth rule used per block-year
            cls = (S["wtd_min"][None] >= R["depth_class_m"]) & cov[None]
            wtd30 = (S["wtd_min"] >= C.WTD_MIN_M) & cov
            hand = (S["hand_frac"] >= C.HAND_PASS_FRACTION) & (S["hand_coverage"] >= C.HAND_PASS_FRACTION)
            assessed = S["lanid_unassessed_dist_m"] > C.IRR_BUFFER_M
            irr = (S["lanid_dist_m"] >= C.IRR_BUFFER_M) & assessed
            irr3 = (S["lanid3_dist_m"] >= C.IRR_BUFFER_M) & assessed
            canal = S["canal_dist_m"] >= C.CANAL_BUFFER_M
        for name, want in (("WTD", wtd), ("CLASSDEPTH", cls), ("WTD30", wtd30[None]), ("HAND", hand[None]),
                           ("IRR", irr[None]), ("IRR3", irr3[None]), ("CANAL", canal[None])):
            got = ((f >> RU.BITS[name]) & 1).astype(bool)
            assert (got == np.broadcast_to(want, got.shape)).all(), name
        assert (R["depth_used_m"] >= R["depth_class_m"]).all() and (R["depth_used_m"] >= min(C.DEPTH_RUNGS_M)).all()
        ref = (f & RU.MAIN_MASK) == RU.MAIN_MASK
        assert (ref == R["reference"].astype(bool)).all()
        ring_mask = RU.mask_of([b for b in RU.MAIN if b != "WTD"] + ["CLASSDEPTH"])
        ring = ((f & ring_mask) == ring_mask) & ~((f >> RU.BITS["WTD"]) & 1).astype(bool)
        assert (ring == R["reference_ring"].astype(bool)).all() and not (ring & ref).any()


def test_training_samples_are_reference_blocks():
    dirs, skip = _need("tiles", "reference", "samples")
    if skip:
        return skip
    from ..p1_2 import blockset as BKS
    from ..p1_1 import build as RB
    T, G = _table(dirs["tiles"])
    BS = BKS.load(dirs["samples"])
    assert (BS["REF"].sum(1) > 0).all() and BS["primary"].all() and BS["valid"].all()
    for ti in np.unique(BS["tile_idx"]):
        t = str(T["tiles"][ti])
        R = RB.open_tile(dirs["reference"], t)
        sel = BS["tile_idx"] == ti
        n = R["reference"].shape[1] * R["reference"].shape[2]
        ref = R["reference"].reshape(R["reference"].shape[0], n).T[BS["index_in_tile"][sel]]
        assert (ref == BS["REF"][sel]).all()


# ------------------------------------------------------------------ (e) and mutation
def test_rebuilding_a_tile_year_gives_identical_arrays():
    dirs, skip = _need("tiles", "blocks")
    if skip:
        return skip
    from ..step0 import blocks as BK
    T, G = _table(dirs["tiles"])
    t = str(T["tiles"][0])
    y = C.agg_years()[-1]
    stored = BK.open_tile_year(dirs["blocks"], t, y)
    if stored is None:
        return f"SKIP no block file for {t} {y}"
    V, _, _ = BK.build_tile_year(G[t], y)
    for k, v in V.items():
        assert np.array_equal(np.asarray(stored[k]), v, equal_nan=v.dtype.kind == "f"), k


def test_a_wrong_scale_factor_is_caught():
    """The oracle comparison must fail when the decoding is wrong (guards against vacuous tests)."""
    dirs, skip = _need("tiles", "blocks")
    if skip:
        return skip
    from ..core import paths
    from ..step0 import blocks as BK
    T, G = _table(dirs["tiles"])
    for k in _sample_blocks(T, 30):
        t, g, br, bc, i, j = _pos(T, G, k)
        y = C.agg_years()[-1]
        B = BK.open_tile_year(dirs["blocks"], t, y)
        lp = paths.tile_path("landsat", t, y)
        if B is None or not os.path.exists(lp):
            continue
        v = O.tile_block_values(lp, t, br, bc, y, 7, "NIRv")
        if v is None or not np.isfinite(v).any():
            continue
        wrong = np.nanmean(v) * 2.0                     # as if scale_factor were 2e-4 instead of 1e-4
        assert abs(float(B["NIRv_mean"][6, i, j]) - wrong) > 1e-4 * max(1, abs(wrong))
        return None
    return "SKIP no block with July data"
