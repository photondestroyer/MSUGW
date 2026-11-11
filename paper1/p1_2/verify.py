"""Stage `verify_p12`: the P1.2 verification tests of the plan (v3) and gate G1b.

From the out-of-fold results of `train` (test items of every fold):
(1) skill: RMSE of M1's median against M0 (the block's own climatology) per aridity decile and target - M1 must
    beat M0 in every decile with >= 100 test block-months, and beat M2 overall.
(2) calibration overall, and the mass bound: share of held-out reference block-water-years whose predicted median
    ET summed over the water year exceeds (P - storage change) x (1 + V2_TOLERANCE) - at most
    PASS['mass_violation_max'].
(3) planted subsidy (per fold, two variants) and the class-flip test (per fold).
Variants on fold (0, 0): without pre-season inputs, without measurement inputs, with the ET / PET ratio target.

From the deployed predictor (`predict`: every block-month predicted by the model of its space fold and year
group, so out of fold in space and in time - the 2012 drought is predicted by models trained on 2015-2025):
residual r = observed - calibrated median; the "null SD" of a target is the SD of r over reference block-months
of the same functional group and aridity tercile. Sums of r, r^2 and of exceedances are streamed per stratum
(see Acc), so the stage never holds all block-months.
(4) coverage of the 90 % interval per stratum (target x functional x tercile x drought class x regime) with
    >= stratum_min_n() reference block-months: within PASS['coverage'].
(5) G1b drought bias: mean r of reference blocks in months with SPI-90 <= DROUGHT_DRY_SPI, per target x functional
    x tercile, separately for the years of each DROUGHT_TEST_YEARS entry: within +-PASS['drought_bias_sd'] null
    SD, and for ET within +-PASS['drought_bias_et_mm'] mm; the same for the wettest months (SPI >= wet).
(6) steps at the observation-regime boundaries: mean r over the 12 months after minus the 12 months before each
    boundary of OBS_REGIMES, per target: below PASS['regime_step_sd'] null SD. Residual against each measurement
    input (binned): range of bin means reported against the same limit.
(7) known-zero sets: share of block-months above the calibrated 95 % quantile among held-out reference blocks
    and the accepted ring - at most PASS['fpr_max']; the Kansas set (no water right, no LANID label) is not
    assessed (WIMAS not on disk). The run-on set (dryland crops failing only HAND) is diagnostic: its excess by
    drought class is reported (the coupling slope of P1.4 decides).
(8) ring test, per functional group x aridity decile: mean r of ring blocks minus that of reference blocks, as a
    90 % interval (two one-sided tests) inside +-PASS['ring_mean_sd'] null SD for every target, and ring
    exceedance <= PASS['fpr_max']. The standard error counts block-years (block-months / 12), not block-months,
    because the months of a block-year share their errors. The coupling-slope criterion waits for P1.4. Strata that
    pass are written to reports/depth_accepted_proposal.json - for review, then copy into config.DEPTH_ACCEPTED.
(9) terrain and canals: binned mean r of reference blocks against HAND, flow accumulation and canal distance:
    lowest-to-highest bin difference below PASS['terrain_sd'] null SD.
(10) 2008-2009 crop layers: residual SD of reference blocks in 2008-2009 over 2010-2012, below
    PASS['crop_layer_sd_ratio']. Winter cereal after fallow and after winter cereal: mean r within
    +-PASS['rotation_sd'] null SD. LANID drift: mean r of reference blocks 2018-2025 minus 2008-2017 per state,
    from yearly means - flagged when the 95 % interval lies above zero.
Gate G1b (decided here; recorded for the steps after P1.2) fails on (4), (5), (6) regime steps or (7).
"""
import json
import os
import warnings

import numpy as np

from .. import config as C
from .. import cdl_groups as CG
from ..core.drive import write_text
from ..step0 import tiles as TL, blocks as BK, static as ST, drivers as DR
from ..p1_1 import build as RB, rules as RU
from . import folds as FO, train as TR, predict as PR, features as FE, baselines as BL

SET_NAMES = ["reference", "ring", "run_on"]
BIN_EDGES = dict(hand_mean=[10, 15, 20, 30, 50, np.inf], facc_max_km2=[0, 0.01, 0.1, 1, 10, np.inf],
                 canal_dist_m=[800, 1600, 3200, 6400, np.inf], oli_share=[0, 0.25, 0.5, 0.75, 1.0001],
                 nobs_mean=[0, 1, 2, 3, 5, np.inf], obs_day_mean=[1, 8, 15, 22, 32],
                 ET_scenes_mean=[0, 1, 2, 3, 5, np.inf])


# ------------------------------------------------------------------ streamed sufficient statistics
class Acc:
    """Per key (tuple of ints over `fields`): n, sum r, sum r^2, count below q05, count above q95."""

    def __init__(self, fields):
        self.fields = list(fields)
        self.d = {}
        self._kv = None

    def _arrays(self):
        if self._kv is None:
            self._kv = (np.array(list(self.d.keys()), "int64"), np.array(list(self.d.values()), "float64"))
        return self._kv

    def add(self, cols, r, lo, hi):
        if r.size == 0:
            return
        self._kv = None
        K = np.stack([np.asarray(c, "int64").ravel() for c in cols], 1)
        u, inv = np.unique(K, axis=0, return_inverse=True)
        inv = np.asarray(inv).ravel()
        m = len(u)
        st = np.stack([np.bincount(inv, minlength=m), np.bincount(inv, r, m), np.bincount(inv, r * r, m),
                       np.bincount(inv, lo.astype(float), m), np.bincount(inv, hi.astype(float), m)], 1)
        for key, row in zip(map(tuple, u.tolist()), st):
            a = self.d.get(key)
            if a is None:
                self.d[key] = row.astype("float64")
            else:
                a += row

    def group(self, by, **where):
        """{tuple of `by` values: stats (5,)} summed over the keys whose fields equal / are in `where`."""
        if not self.d:
            return {}
        K, V = self._arrays()
        keep = np.ones(len(K), bool)
        for f, v in where.items():
            col = K[:, self.fields.index(f)]
            keep &= np.isin(col, np.atleast_1d(v))
        K, V = K[keep], V[keep]
        if not len(K):
            return {}
        idx = [self.fields.index(f) for f in by]
        u, inv = np.unique(K[:, idx], axis=0, return_inverse=True)
        out = np.zeros((len(u), 5))
        np.add.at(out, np.asarray(inv).ravel(), V)
        return {tuple(k): v for k, v in zip(u.tolist(), out)}


def mean_sd(st):
    n = st[0]
    if n <= 0:
        return float("nan"), float("nan"), 0
    mu = st[1] / n
    var = max(st[2] / n - mu * mu, 0.0)
    return float(mu), float(np.sqrt(var)), int(n)


# ------------------------------------------------------------------ out-of-fold results (train stage)
def skill(oofs):
    rep = {}
    for t, name in enumerate(C.TARGETS):
        y = np.concatenate([o["y"][..., t][o["ymask"][..., t]] for o in oofs])
        m1 = np.concatenate([o["m1"][..., t, 1][o["ymask"][..., t]] for o in oofs])
        m2 = np.concatenate([o["m2"][..., t, 1][o["ymask"][..., t]] for o in oofs])
        m0 = np.concatenate([o["m0"][..., t, 1][o["ymask"][..., t]] for o in oofs])
        lo = np.concatenate([o["m1"][..., t, 0][o["ymask"][..., t]] for o in oofs])
        hi = np.concatenate([o["m1"][..., t, 2][o["ymask"][..., t]] for o in oofs])
        dec = np.concatenate([np.broadcast_to(o["aridity"][:, None], o["ymask"].shape[:2])[o["ymask"][..., t]]
                              for o in oofs])
        rmse = lambda a, b, m: float(np.sqrt(np.mean((a[m] - b[m]) ** 2))) if m.sum() else float("nan")
        both = np.isfinite(m0) & np.isfinite(m1)
        by = {}
        for d in range(1, 11):
            m = both & (dec == d)
            if m.sum() >= 100:
                by[d] = dict(n=int(m.sum()), rmse_m1=rmse(y, m1, m), rmse_m0=rmse(y, m0, m),
                             skill=1 - rmse(y, m1, m) / rmse(y, m0, m))
        m2ok = np.isfinite(m2) & np.isfinite(m1)
        cov = float(np.mean((y >= lo) & (y <= hi))) if y.size else float("nan")
        r1, r2 = rmse(y, m1, m2ok), rmse(y, m2, m2ok)
        rep[name] = dict(n=int(y.size), by_decile=by,
                         # None = not testable here (no M0 / M2 value for the held-out block-months)
                         beats_m0_everywhere=(all(v["skill"] > 0 for v in by.values()) if by else None),
                         rmse_m1=r1, rmse_m2=r2, beats_m2=(r1 < r2 if np.isfinite(r2) else None),
                         coverage=cov, coverage_passes=C.PASS["coverage"][0] <= cov <= C.PASS["coverage"][1])
    return rep


def mass_violation(oofs):
    n = bad = 0
    for o in oofs:
        if "bound" not in o:
            continue
        ref = o["ymask"][..., 0].any(1)
        et = o["m1"][..., 0, 1].sum(1)
        ok = ref & np.isfinite(et) & np.isfinite(o["bound"])
        n += int(ok.sum())
        bad += int((et[ok] > o["bound"][ok] * (1 + C.V2_TOLERANCE)).sum())
    share = bad / n if n else float("nan")
    return dict(block_water_years=n, violating=bad, share=share,
                passes=(share <= C.PASS["mass_violation_max"]) if n else None)


def plants(metas):
    out = {}
    for variant in ("season_only", "with_preseason"):
        out[variant] = {}
        for name in C.TARGETS:
            vals = [m["plant"][variant][name] for m in metas if m.get("plant") and m["plant"][variant].get(name)]
            if not vals:
                continue
            w = np.array([v["n"] for v in vals], "float64")
            est = float(np.average([v["estimated"] for v in vals], weights=w))
            out[variant][name] = dict(planted=vals[0]["planted"], estimated=est, bias=est / vals[0]["planted"] - 1,
                                      passes=abs(est / vals[0]["planted"] - 1) < C.PASS["plant_bias"], folds=len(vals))
    return out


def class_flips(metas):
    out = {}
    for name in C.TARGETS:
        vals = [m["class_flip"][name] for m in metas if m.get("class_flip") and m["class_flip"].get(name)]
        if not vals:
            continue
        n = sum(v["n"] for v in vals)
        rf = sum(v["rate_flipped"] * v["n"] for v in vals) / n
        ro = sum(v["rate_original"] * v["n"] for v in vals) / n
        dq = [v.get("mean_abs_change_q50") for v in vals if v.get("mean_abs_change_q50") is not None]
        out[name] = dict(n=n, rate_flipped=rf, rate_original=ro, passes=bool(rf <= C.PASS["fpr_max"]),
                         mean_abs_change_q50=float(np.mean(dq)) if dq else None)
    return out


def variants(tdir, main00):
    """Fold (0, 0) with and without pre-season inputs / measurement inputs, and the ET / PET ratio model."""
    out = {}
    for v in ("nopre", "nomeas", "ratio"):
        p = f"{tdir}/k0_g0_{v}/oof.npz"
        if not os.path.exists(p):
            continue
        with np.load(p) as z:
            o = {k: z[k] for k in z.files}
        row = {}
        for t, name in enumerate(C.TARGETS):
            m = main00["ymask"][..., t] & o["ymask"][..., t]
            if v == "ratio" and t == 0:
                y_v = o["y"][..., 0] * o["pet"]
                q_v = o["m1"][..., 0, 1] * o["pet"]
            else:
                y_v, q_v = o["y"][..., t], o["m1"][..., t, 1]
            y0, q0 = main00["y"][..., t], main00["m1"][..., t, 1]
            ok = m & np.isfinite(q_v) & np.isfinite(q0) & np.isfinite(y_v)
            if not ok.any():
                continue
            row[name] = dict(n=int(ok.sum()), rmse_main=float(np.sqrt(np.mean((y0[ok] - q0[ok]) ** 2))),
                             rmse_variant=float(np.sqrt(np.mean((y_v[ok] - q_v[ok]) ** 2))),
                             coverage_main=BL.coverage(y0[ok], main00["m1"][..., t, 0][ok], main00["m1"][..., t, 2][ok]))
        out[v] = row
    return out


# ------------------------------------------------------------------ deployed predictor: streamed residuals
def residual_pass(T, G, dirs, D):
    """Stream every tile-year: residuals of the reference, ring and run-on sets into accumulators."""
    lut = DR.cell_lookup(D)
    months = D["months"]
    base = Acc(["set", "target", "func", "decile", "drought", "regime", "year", "month"])
    by_state = Acc(["target", "state", "year"])
    by_bin = Acc(["target", "var", "bin"])
    rotation = Acc(["target", "func", "func_prev"])
    bin_vars = [v for v in BIN_EDGES if v in ("hand_mean", "facc_max_km2", "canal_dist_m") or v in C.MEAS_INPUTS]
    starts = [np.datetime64(s, "M") for _, s in C.OBS_REGIMES]
    for ti, t in enumerate(T["tiles"]):
        g = G[t]
        sl = TL.tile_slice(T, ti)
        R = RB.open_tile(dirs["reference"], t)
        S = ST.open_tile(dirs["static"], t)
        geoid = ST.open_county_geoid(dirs["static"], t)
        state = np.array([[int(geoid[c][:2]) if 0 <= c < len(geoid) else 0 for c in row] for row in S["county"]])
        iy, ix = T["lab_iy"][sl], T["lab_ix"][sl]
        cell = np.where(iy >= 0, lut[np.maximum(iy, 0), np.maximum(ix, 0)], -1).reshape(g.nby, g.nbx)
        use = (T["valid"][sl] & T["primary"][sl]).reshape(g.nby, g.nbx)
        dec = R["aridity_decile"].astype("int64")
        fn_all = R["functional"] if "functional" in R else CG.functional_of(R["group"])
        for yi, y in enumerate(R["year"]):
            y = int(y)
            P = PR.open_tile_year(dirs["predict"], t, y)
            B = BK.open_tile_year(dirs["blocks"], t, y)
            if P is None or B is None:
                continue
            sets = np.stack([R["reference"][yi].astype(bool), R["reference_ring"][yi].astype(bool)
                             if "reference_ring" in R else np.zeros(use.shape, bool),
                             RU.run_on_set(R["flags"][yi], R["group"][yi])]) & use[None]
            if not sets.any():
                continue
            fn = fn_all[yi].astype("int64")
            fprev = CG.functional_of(R["group_prev"][yi]) if "group_prev" in R else np.full(use.shape, CG.FUNC_NONE)
            mi = np.searchsorted(months, np.datetime64(f"{y}-01", "M") + np.arange(12))
            mi_ok = (mi < months.size)
            spi = np.full((12,) + use.shape, np.nan, "float32")
            if "SPI90" in D:
                for m in range(12):
                    if mi_ok[m] and months[mi[m]] == np.datetime64(f"{y}-{m + 1:02d}", "M"):
                        spi[m] = np.where(cell >= 0, D["SPI90"][mi[m]][np.maximum(cell, 0)], np.nan)
            drought = DR.drought_class(spi)
            regime = np.array([int(np.searchsorted(starts, np.datetime64(f"{y}-{m + 1:02d}", "M"), side="right") - 1)
                               for m in range(12)])
            need = np.maximum(B["roi_px"].astype("float64") * C.TARGET_MIN_PIXEL_FRACTION, 1)
            month = np.broadcast_to(np.arange(12)[:, None, None], (12,) + use.shape)
            for ti_, name in enumerate(C.TARGETS):
                obs = np.where(B[f"{name}_n"] >= need, B[f"{name}_mean"], np.nan)
                q05, q50, q95 = P[f"{name}_q05"], P[f"{name}_q50"], P[f"{name}_q95"]
                ok = np.isfinite(obs) & np.isfinite(q50)
                r = obs - q50
                lo, hi = obs < q05, obs > q95
                for si in range(3):
                    m = ok & sets[si][None]
                    if not m.any():
                        continue
                    b3 = lambda a: np.broadcast_to(a, (12,) + use.shape)[m]
                    base.add([np.full(int(m.sum()), si), np.full(int(m.sum()), ti_), b3(fn), b3(dec), drought[m],
                              regime[month[m]], np.full(int(m.sum()), y), month[m]], r[m], lo[m], hi[m])
                    if si != 0:
                        continue
                    by_state.add([np.full(int(m.sum()), ti_), b3(state), np.full(int(m.sum()), y)], r[m], lo[m], hi[m])
                    rotation.add([np.full(int(m.sum()), ti_), b3(fn), b3(fprev)], r[m], lo[m], hi[m])
                    for vi, v in enumerate(bin_vars):
                        src = B[v] if v in B else (np.broadcast_to(S[v], (12,) + use.shape) if v in S else None)
                        if src is None:
                            continue
                        x = np.broadcast_to(src, (12,) + use.shape)[m]
                        e = np.asarray(BIN_EDGES[v], "float64")
                        with np.errstate(invalid="ignore"):
                            bi = np.searchsorted(e, x, side="right") - 1
                        okb = np.isfinite(x) & (bi >= 0) & (bi < len(e) - 1)
                        by_bin.add([np.full(int(okb.sum()), ti_), np.full(int(okb.sum()), vi), bi[okb]],
                                   r[m][okb], lo[m][okb], hi[m][okb])
    return dict(base=base, by_state=by_state, by_bin=by_bin, rotation=rotation, bin_vars=bin_vars)


def null_sd(base):
    """{(target, func, tercile): SD of r over reference block-months}, and {target: pooled SD}."""
    out = {}
    for (t, f, d), st in base.group(["target", "func", "decile"], set=0).items():
        key = (t, f, int(BL.tercile(d)))
        out[key] = out.get(key, np.zeros(5)) + st
    sd = {k: mean_sd(v)[1] for k, v in out.items()}
    pooled = {t: mean_sd(st)[1] for (t,), st in base.group(["target"], set=0).items()}
    return sd, pooled


def _sd_of(sd, pooled, t, f, terc):
    v = sd.get((t, f, terc))
    return v if v and np.isfinite(v) and v > 0 else pooled.get(t, float("nan"))


def coverage_strata(base):
    rows, failing = {}, []
    nmin = C.stratum_min_n()
    acc = {}
    for (t, f, d, dr, rg), st in base.group(["target", "func", "decile", "drought", "regime"], set=0).items():
        key = (t, f, int(BL.tercile(d)), dr, rg)
        acc[key] = acc.get(key, np.zeros(5)) + st
    for (t, f, terc, dr, rg), st in sorted(acc.items()):
        if st[0] < nmin:
            continue
        cov = 1 - (st[3] + st[4]) / st[0]
        name = f"{C.TARGETS[t]}|{CG.FUNC_NAMES[f]}|T{terc + 1}|{['dry', 'mid', 'wet'][dr]}|{C.OBS_REGIMES[rg][0]}"
        ok = C.PASS["coverage"][0] <= cov <= C.PASS["coverage"][1]
        rows[name] = dict(n=int(st[0]), coverage=float(cov), passes=bool(ok))
        if not ok:
            failing.append(name)
    return dict(strata=rows, failing=failing, assessed=len(rows))


def drought_bias(base, sd, pooled):
    out, failing = {}, []
    nmin = C.stratum_min_n()
    for label, years in C.DROUGHT_TEST_YEARS.items():
        for cls, cname in ((0, "dry"), (2, "wet")):
            acc = {}
            for (t, f, d), st in base.group(["target", "func", "decile"], set=0, drought=cls, year=years).items():
                key = (t, f, int(BL.tercile(d)))
                acc[key] = acc.get(key, np.zeros(5)) + st
            for (t, f, terc), st in sorted(acc.items()):
                mu, _, n = mean_sd(st)
                if n < nmin:
                    continue
                s = _sd_of(sd, pooled, t, f, terc)
                z = mu / s if s > 0 else float("nan")
                ok = abs(z) <= C.PASS["drought_bias_sd"] and (t != 0 or abs(mu) <= C.PASS["drought_bias_et_mm"])
                name = f"{label}|{cname}|{C.TARGETS[t]}|{CG.FUNC_NAMES[f]}|T{terc + 1}"
                out[name] = dict(n=n, mean=mu, mean_null_sd=z, passes=bool(ok))
                if not ok:
                    failing.append(name)
    return dict(strata=out, failing=failing, assessed=len(out))


def regime_steps(base, pooled):
    out, failing = {}, []
    for k in range(1, len(C.OBS_REGIMES)):
        b = np.datetime64(C.OBS_REGIMES[k][1], "M")
        before = [b - 12 + i for i in range(12)]
        after = [b + i for i in range(12)]
        for t in range(len(C.TARGETS)):
            sums = {}
            for side, ms in (("before", before), ("after", after)):
                st = np.zeros(5)
                for mth in ms:
                    yy, mm = int(str(mth)[:4]), int(str(mth)[5:7]) - 1
                    for _, v in base.group(["target"], set=0, target=t, year=yy, month=mm).items():
                        st += v
                sums[side] = st
            if sums["before"][0] < C.stratum_min_n() or sums["after"][0] < C.stratum_min_n():
                continue
            step = mean_sd(sums["after"])[0] - mean_sd(sums["before"])[0]
            z = step / pooled[t] if pooled.get(t) else float("nan")
            ok = abs(z) < C.PASS["regime_step_sd"]
            name = f"{C.OBS_REGIMES[k][1]}|{C.TARGETS[t]}"
            out[name] = dict(step=step, step_null_sd=z, passes=bool(ok))
            if not ok:
                failing.append(name)
    return dict(boundaries=out, failing=failing, assessed=len(out))


def binned(by_bin, bin_vars, pooled):
    out, flagged = {}, []
    for vi, v in enumerate(bin_vars):
        for t in range(len(C.TARGETS)):
            g = by_bin.group(["bin"], target=t, var=vi)
            rows = {int(b): mean_sd(st) for (b,), st in sorted(g.items()) if st[0] >= max(10, C.stratum_min_n() // 10)}
            if len(rows) < 2:
                continue
            means = [m for m, _, _ in rows.values()]
            diff = (means[-1] - means[0]) / pooled[t] if pooled.get(t) else float("nan")
            rng = (max(means) - min(means)) / pooled[t] if pooled.get(t) else float("nan")
            e = BIN_EDGES[v]
            name = f"{v}|{C.TARGETS[t]}"
            out[name] = dict(bins={f"{e[b]}-{e[b + 1]}": dict(mean=m, n=n) for b, (m, _, n) in rows.items()},
                             first_to_last_null_sd=diff, range_null_sd=rng)
            lim = C.PASS["terrain_sd"] if v in ("hand_mean", "facc_max_km2", "canal_dist_m") else C.PASS["regime_step_sd"]
            if abs(diff) >= lim:
                flagged.append(name)
    return dict(variables=out, flagged=flagged)


def known_zero(base):
    out = {}
    for si, sname in enumerate(SET_NAMES):
        for (t,), st in base.group(["target"], set=si).items():
            rate = st[4] / st[0] if st[0] else float("nan")
            out[f"{sname}|{C.TARGETS[t]}"] = dict(n=int(st[0]), above_q95=float(rate),
                                                    passes=bool(rate <= C.PASS["fpr_max"]) if sname != "run_on" else None)
    run_on = {}
    for (t, dr), st in base.group(["target", "drought"], set=2).items():
        mu, _, n = mean_sd(st)
        run_on[f"{C.TARGETS[t]}|{['dry', 'mid', 'wet'][dr]}"] = dict(n=n, mean_residual=mu,
                                                                    above_q95=float(st[4] / st[0]) if st[0] else None)
    failing = [k for k, v in out.items() if v["passes"] is False]
    return dict(sets=out, run_on_by_drought=run_on, failing=failing,
                not_assessed=["Kansas fields with no water right and no LANID label (WIMAS not on disk)"])


def ring_test(base, sd, pooled):
    """Per functional group x decile: TOST of mean(ring) - mean(reference) within +-ring_mean_sd null SD."""
    z90 = 1.645
    rows, proposal = {}, {}
    ring = base.group(["func", "decile", "target"], set=1)
    ref = base.group(["func", "decile", "target"], set=0)
    strata = sorted({(f, d) for f, d, _ in ring})
    for f, d in strata:
        per_t, all_ok = {}, True
        for t in range(len(C.TARGETS)):
            a, b = ring.get((f, d, t)), ref.get((f, d, t))
            if a is None or b is None or a[0] < C.stratum_min_n() or b[0] < C.stratum_min_n():
                per_t[C.TARGETS[t]] = None
                all_ok = False
                continue
            ma, sa, na = mean_sd(a)
            mb, sb, nb = mean_sd(b)
            s = _sd_of(sd, pooled, t, f, int(BL.tercile(d)))
            se = np.sqrt(sa ** 2 / max(na / 12, 1) + sb ** 2 / max(nb / 12, 1))
            diff = (ma - mb) / s
            lo, hi = (ma - mb - z90 * se) / s, (ma - mb + z90 * se) / s
            exc = a[4] / a[0]
            ok = (lo > -C.PASS["ring_mean_sd"]) and (hi < C.PASS["ring_mean_sd"]) and exc <= C.PASS["fpr_max"]
            per_t[C.TARGETS[t]] = dict(n_ring=na, n_reference=nb, diff_null_sd=diff, ci90=[lo, hi], ring_above_q95=exc,
                                       passes=bool(ok))
            all_ok &= bool(ok)
        name = f"{CG.FUNC_NAMES[f]}|{d}"
        rows[name] = dict(targets=per_t, passes=all_ok, coupling_slope="pending (P1.4)")
        rung = C.DEPTH_BY_FUNCTIONAL.get(CG.FUNC_NAMES[f])
        if all_ok and rung is not None and rung < C.WTD_MIN_M:
            proposal[name] = rung
    return dict(strata=rows, proposal=proposal)


def crop_layer_test(base):
    out = {}
    for t in range(len(C.TARGETS)):
        a = sum(base.group(["target"], set=0, target=t, year=[2008, 2009]).values(), np.zeros(5))
        b = sum(base.group(["target"], set=0, target=t, year=[2010, 2011, 2012]).values(), np.zeros(5))
        if a[0] < C.stratum_min_n() or b[0] < C.stratum_min_n():
            continue
        ratio = mean_sd(a)[1] / mean_sd(b)[1]
        out[C.TARGETS[t]] = dict(sd_ratio=ratio, passes=bool(ratio < C.PASS["crop_layer_sd_ratio"]))
    return out


def rotation_test(rotation, pooled):
    wc, fa = CG.FUNC_NAMES.index("winter_cereal"), CG.FUNC_NAMES.index("fallow")
    out = {}
    for t in range(len(C.TARGETS)):
        for prev, label in ((fa, "after_fallow"), (wc, "after_winter_cereal")):
            st = rotation.d.get((t, wc, prev))
            if st is None or st[0] < C.stratum_min_n():
                continue
            mu, _, n = mean_sd(st)
            z = mu / pooled[t] if pooled.get(t) else float("nan")
            out[f"{C.TARGETS[t]}|{label}"] = dict(n=n, mean_null_sd=z, passes=bool(abs(z) <= C.PASS["rotation_sd"]))
    return out


def drift_test(by_state, pooled):
    (b0, b1), (a0, a1) = C.DRIFT_YEARS
    out, flagged = {}, []
    for t in range(len(C.TARGETS)):
        states = sorted({key[1] for key in by_state.d if key[0] == t})
        for s in states:
            yearly = {y: mean_sd(st)[0] for (y,), st in by_state.group(["year"], target=t, state=s).items()
                      if st[0] >= max(10, C.stratum_min_n() // 10)}
            before = [v for y, v in yearly.items() if b0 <= y <= b1]
            after = [v for y, v in yearly.items() if a0 <= y <= a1]
            if len(before) < 2 or len(after) < 2:
                continue
            d = np.mean(after) - np.mean(before)
            se = np.sqrt(np.var(after, ddof=1) / len(after) + np.var(before, ddof=1) / len(before))
            lo = d - 1.96 * se
            name = f"{s:02d}|{C.TARGETS[t]}"
            out[name] = dict(diff=float(d), diff_null_sd=float(d / pooled[t]) if pooled.get(t) else None,
                             ci95_low=float(lo), years_before=len(before), years_after=len(after))
            if lo > 0:
                flagged.append(name)
    return dict(states=out, upward_drift=flagged)


def gate_g1b(rep):
    reasons = []
    for key, label in (("coverage_strata", "coverage"), ("drought_bias", "drought bias"),
                       ("regime_steps", "regime step"), ("known_zero", "known-zero exceedance")):
        f = rep[key]["failing"]
        if f:
            reasons.append(f"{label}: {', '.join(f[:6])}" + (f" (+{len(f) - 6} more)" if len(f) > 6 else ""))
    not_assessed = [label for key, label in (("coverage_strata", "coverage per stratum"),
                                             ("drought_bias", "drought bias on unseen droughts"),
                                             ("regime_steps", "regime steps")) if not rep[key]["assessed"]]
    return dict(name="G1b", passes=not reasons, reasons=reasons,
                not_assessed=not_assessed + rep["known_zero"]["not_assessed"])


def run(ctx):
    T, G = TL.load_table(ctx.dep("tiles"))
    F = FO.load(ctx.dep("folds"))
    D = DR.load(ctx.dep("drivers"))
    tdir = ctx.dep("train")
    metas, oofs = [], []
    for k in range(F["folds"]):
        for g in range(len(F["year_groups"])):
            m, o = TR.load_fold(tdir, k, g)
            metas.append(m)
            oofs.append(o)
    rep = dict(skill=skill(oofs), mass_bound=mass_violation(oofs), plant=plants(metas), class_flip=class_flips(metas),
               variants=variants(tdir, oofs[0]))
    dirs = dict(reference=ctx.dep("reference"), blocks=ctx.dep("blocks"), predict=ctx.dep("predict"),
                static=ctx.dep("static"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        A = residual_pass(T, G, dirs, D)
    sd, pooled = null_sd(A["base"])
    rep["null_sd_pooled"] = {C.TARGETS[t]: v for t, v in pooled.items()}
    rep["coverage_strata"] = coverage_strata(A["base"])
    rep["drought_bias"] = drought_bias(A["base"], sd, pooled)
    rep["regime_steps"] = regime_steps(A["base"], pooled)
    rep["binned"] = binned(A["by_bin"], A["bin_vars"], pooled)
    rep["known_zero"] = known_zero(A["base"])
    rep["ring_test"] = ring_test(A["base"], sd, pooled)
    rep["crop_layer_2008_2009"] = crop_layer_test(A["base"])
    rep["rotation"] = rotation_test(A["rotation"], pooled)
    rep["lanid_drift"] = drift_test(A["by_state"], pooled)
    rep["gate"] = gate_g1b(rep)
    outdir = f"{ctx.out}/reports"
    os.makedirs(outdir, exist_ok=True)
    write_text(f"{outdir}/verify_p12.json", json.dumps(rep, indent=1, default=float))
    write_text(f"{outdir}/verify_p12.md", markdown(rep))
    write_text(f"{outdir}/depth_accepted_proposal.json", json.dumps(dict(
        note="Strata whose ring passed the ring test. Review verify_p12.md (ring test) and copy the entries you "
             "accept into config.DEPTH_ACCEPTED, then rerun `reference` and `train`. The coupling-slope criterion "
             "(P1.4) is not yet part of this proposal.",
        DEPTH_ACCEPTED=rep["ring_test"]["proposal"]), indent=1))
    summ = {n: dict(beats_m0=rep["skill"][n]["beats_m0_everywhere"], beats_m2=rep["skill"][n]["beats_m2"],
                    coverage=rep["skill"][n]["coverage"]) for n in C.TARGETS}
    ctx.complete(dict(gate=rep["gate"], targets=summ, ring_proposal=rep["ring_test"]["proposal"]))
    ctx.log.info("P1.2 verification written to %s; gate G1b %s", outdir, "passes" if rep["gate"]["passes"] else
                 "FAILS: " + "; ".join(rep["gate"]["reasons"]))


def _yn(v):
    return "n/a" if v is None else ("PASS" if v else "FAIL") if isinstance(v, bool) else v


def markdown(rep):
    gt = rep["gate"]
    L = ["# P1.2 verification", "",
         f"**Gate G1b: {'PASS' if gt['passes'] else 'FAIL'}**" + ("" if gt["passes"] else " - " + "; ".join(gt["reasons"])),
         "", "Not assessed: " + ("; ".join(gt["not_assessed"]) or "none"), ""]
    L += ["## (1) skill and calibration (out-of-fold test items)", "",
          "| target | test block-months | M1 beats M0 in every decile | M1 RMSE | M2 RMSE | M1 beats M2 | 90 % coverage |",
          "|---|---|---|---|---|---|---|"]
    for n, s in rep["skill"].items():
        L.append(f"| {n} | {s['n']:,} | {_yn(s['beats_m0_everywhere'])} | {s['rmse_m1']:.4g} | {s['rmse_m2']:.4g} | "
                 f"{_yn(s['beats_m2'])} | {s['coverage']:.3f} |")
    mb = rep["mass_bound"]
    L += ["", f"Mass bound: {mb['violating']:,} of {mb['block_water_years']:,} held-out reference block-water-years "
              f"exceed it ({100 * mb['share']:.2f} %; at most {100 * C.PASS['mass_violation_max']:.0f} %): {_yn(mb['passes'])}", ""]
    L += ["## (3) planted subsidy and class flip", "", "| variant | target | planted | estimated | bias | passes |",
          "|---|---|---|---|---|---|"]
    for var, d in rep["plant"].items():
        for n, v in d.items():
            L.append(f"| {var} | {n} | {v['planted']:.4g} | {v['estimated']:.4g} | {100 * v['bias']:+.1f} % | {_yn(v['passes'])} |")
    L += ["", "| class flip | block-months | above q95 flipped | above q95 original | mean change of median | passes |",
          "|---|---|---|---|---|---|"]
    L += [f"| {n} | {v['n']:,} | {v['rate_flipped']:.3f} | {v['rate_original']:.3f} | "
          f"{'n/a' if v.get('mean_abs_change_q50') is None else format(v['mean_abs_change_q50'], '.4g')} | {_yn(v['passes'])} |"
          for n, v in rep["class_flip"].items()]
    L += ["", "## Variants on fold (0, 0)", "", "| variant | target | RMSE main | RMSE variant |", "|---|---|---|---|"]
    for v, d in rep["variants"].items():
        L += [f"| {v} | {n} | {r['rmse_main']:.4g} | {r['rmse_variant']:.4g} |" for n, r in d.items()]
    for key, title in (("coverage_strata", "(4) coverage per stratum"), ("drought_bias", "(5) drought bias (G1b)")):
        r = rep[key]
        L += ["", f"## {title}", "", f"{r['assessed']} strata assessed; failing: {', '.join(r['failing']) or 'none'}"]
    rs = rep["regime_steps"]
    L += ["", "## (6) regime steps", "", "| boundary / target | step (null SD) | passes |", "|---|---|---|"]
    L += [f"| {k} | {v['step_null_sd']:+.3f} | {_yn(v['passes'])} |" for k, v in rs["boundaries"].items()]
    L += ["", "Binned residuals beyond their limit: " + (", ".join(rep["binned"]["flagged"]) or "none")]
    kz = rep["known_zero"]
    L += ["", "## (7) known-zero sets", "", "| set / target | block-months | above q95 | passes |", "|---|---|---|---|"]
    L += [f"| {k} | {v['n']:,} | {v['above_q95']:.3f} | {_yn(v['passes'])} |" for k, v in kz["sets"].items()]
    L += ["", "Run-on set by drought class (diagnostic): " +
          ", ".join(f"{k}: mean {v['mean_residual']:+.3g} (n={v['n']:,})" for k, v in kz["run_on_by_drought"].items())]
    rt = rep["ring_test"]
    L += ["", "## (8) ring test", "", "| stratum | passes | ET / NIRv / NDMI difference (null SD, 90 % interval) |", "|---|---|---|"]
    for k, v in rt["strata"].items():
        cells = []
        for n in C.TARGETS:
            x = v["targets"].get(n)
            cells.append("n/a" if x is None else f"{x['diff_null_sd']:+.2f} [{x['ci90'][0]:+.2f}, {x['ci90'][1]:+.2f}]")
        L.append(f"| {k} | {_yn(v['passes'])} | {' / '.join(cells)} |")
    L += ["", f"Proposed DEPTH_ACCEPTED (review before copying into config): {json.dumps(rt['proposal'])}", ""]
    L += ["## (10) crop layers, rotation, LANID drift", "",
          "2008-2009 SD ratio: " + (", ".join(f"{k}: {v['sd_ratio']:.2f} ({_yn(v['passes'])})"
                                              for k, v in rep["crop_layer_2008_2009"].items()) or "not assessed"),
          "", "Winter cereal: " + (", ".join(f"{k}: {v['mean_null_sd']:+.2f} ({_yn(v['passes'])})"
                                             for k, v in rep["rotation"].items()) or "not assessed"),
          "", "Upward drift after 2017 (state|target): " + (", ".join(rep["lanid_drift"]["upward_drift"]) or "none"), ""]
    return "\n".join(L)
