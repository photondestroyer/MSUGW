"""Stage `verify_p11`: the P1.1 verification tests of the plan (v3) and gate G1a.

V1  reference blocks per aridity decile x FUNCTIONAL group (distinct blocks and block-years); strata with fewer
    than MIN_REF_BLOCKS_PER_STRATUM distinct blocks are listed as out of domain (they go to P1.3).
V2  plausibility (diagnostic, not a selection rule): water-year ET <= (P - storage change) * (1 + V2_TOLERANCE)
    in >= V2_PASS_SHARE of reference block-years with 12 valid ET months. Storage change = bucket storage at
    the end of September minus a year earlier (bucket bin of the block). Failures are diagnosed: a failure
    rate that rises near irrigation means contamination; a rate that is uniform within a class points to the
    ET or rain product (plan v3).
V3  yearly maps of the reference set (PNG) - the reference NetCDF files hold the same data.
V4  representativeness: standardised mean difference of V4_VARS between reference blocks and all non-irrigated
    blocks (valid, counted copy, no LANID pixel within IRR_BUFFER_M), per aridity decile; |SMD| <= V4_SMD_MAX.
    Clay and slope are named by the plan but are not among the static layers; they are listed as not available.
P1.0 valid share: share of valid block-months (ET, NIRv and NDMI) in April-September 2012 per functional group x
    decile, against 2014-2016; strata below P10_VALID_SHARE_MIN are listed (2012 results there are marked).
Ring (information): blocks that only the looser class rung admits, per functional group x decile.
Gate G1a (halts `samples` in focus/full mode unless "G1a" is in GATE_OVERRIDE) fails when
    the main set fails V2, or no stratum is in domain, or DEPTH_ACCEPTED is malformed.
    Out-of-domain strata do not fail the gate: their fallback (P1.3) is part of the design.
    Not assessed (data not on disk): the Kansas registry comparison (WIMAS) and the ring test, which needs
    the P1.2 models and is part of gate G1b's report.
Every number is reported for the main set and the sensitivity sets.
"""
import json
import os

import numpy as np

from .. import config as C
from .. import cdl_groups as CG
from ..core.drive import write_text
from ..step0 import tiles as TL, static as ST, blocks as BK, drivers as DR
from . import build as RB, rules as RU

KEYS = ["reference", "reference_gde", "reference_irr3"]


def _func(R):
    return R["functional"] if "functional" in R else CG.functional_of(R["group"])


def v1_counts(T, G, rdir, key):
    blocks, byears = {}, {}
    for t in T["tiles"]:
        R = RB.open_tile(rdir, t)
        ref = R[key].astype(bool)
        dec = np.broadcast_to(R["aridity_decile"], ref.shape)
        fn = _func(R)
        for d in range(1, 11):
            for fi, name in enumerate(CG.FUNC_NAMES):
                m = ref & (dec == d) & (fn == fi)
                if m.any():
                    k = f"{d}|{name}"
                    byears[k] = byears.get(k, 0) + int(m.sum())
                    blocks[k] = blocks.get(k, 0) + int(m.any(0).sum())
    out_of_domain = sorted(k for k, v in blocks.items() if v < C.MIN_REF_BLOCKS_PER_STRATUM)
    return dict(distinct_blocks=blocks, block_years=byears, out_of_domain=out_of_domain,
                strata_in_domain=len(blocks) - len(out_of_domain), total_distinct=sum(blocks.values()))


def v2_water_balance(T, G, rdir, bdir, sdir, D, key):
    lut = DR.cell_lookup(D)
    months = D["months"]
    caps = np.asarray(C.BUCKET_CAPACITY_BINS_MM, "float64")
    n_ok = n_all = 0
    ratios, fails, dists, funcs = [], [], [], []
    for ti, t in enumerate(T["tiles"]):
        g = G[t]
        sl = TL.tile_slice(T, ti)
        R = RB.open_tile(rdir, t)
        S = ST.open_tile(sdir, t)
        fn = _func(R)
        iy, ix = T["lab_iy"][sl], T["lab_ix"][sl]
        cell = np.where(iy >= 0, lut[np.maximum(iy, 0), np.maximum(ix, 0)], -1).reshape(g.nby, g.nbx)
        for k, y in enumerate(R["year"]):
            prev, cur = BK.open_tile_year(bdir, t, int(y) - 1), BK.open_tile_year(bdir, t, int(y))
            if prev is None or cur is None:
                continue
            et = np.concatenate([prev["ET_mean"][9:12], cur["ET_mean"][0:9]], 0)
            ref = R[key][k].astype(bool) & np.isfinite(et).all(0) & (cell >= 0) & (S["bucket_bin"] >= 0)
            if not ref.any():
                continue
            m0 = np.datetime64(f"{int(y) - 1}-10", "M")
            mi = np.searchsorted(months, m0 + np.arange(12))
            if mi[-1] >= months.size or (months[mi] != m0 + np.arange(12)).any():
                continue
            c = cell[ref]
            P = D["P"][mi][:, c].sum(0)
            b = S["bucket_bin"][ref].astype(int)
            s_end = D["BUCKET"][mi[-1], b, c]
            prev_sep = np.searchsorted(months, m0 - 1)
            if months[prev_sep] != m0 - 1:
                continue
            s_start = D["BUCKET"][prev_sep, b, c]
            dS = (s_end - s_start) * caps[b]
            bound = P - dS
            e = et[:, ref].sum(0)
            ok = np.isfinite(bound) & np.isfinite(e)
            passed = e[ok] <= bound[ok] * (1 + C.V2_TOLERANCE)
            n_all += int(ok.sum())
            n_ok += int(passed.sum())
            with np.errstate(invalid="ignore", divide="ignore"):
                ratios.append((e[ok] / bound[ok]).astype("float32"))
            fails.append(~passed)
            dists.append(S["lanid_dist_m"][ref][ok])
            funcs.append(fn[k][ref][ok])
    share = n_ok / n_all if n_all else float("nan")
    r = np.concatenate(ratios) if ratios else np.array([], "float32")
    q = np.nanquantile(r, [0.05, 0.5, 0.95]).tolist() if r.size else []
    return dict(block_years=n_all, passing=n_ok, share=share, passes=bool(n_all and share >= C.V2_PASS_SHARE),
                et_over_available_q05_q50_q95=q, diagnosis=v2_diagnosis(fails, dists, funcs))


def v2_diagnosis(fails, dists, funcs):
    """Failure rate by distance to the nearest LANID pixel and by functional group."""
    if not fails:
        return {}
    f, d, fn = np.concatenate(fails), np.concatenate(dists), np.concatenate(funcs)
    edges = [C.IRR_BUFFER_M, 2 * C.IRR_BUFFER_M, 4 * C.IRR_BUFFER_M, np.inf]
    by_dist, lo = {}, 0.0
    for hi in edges:
        m = (d >= lo) & (d < hi)
        if m.any():
            by_dist[f"{lo:.0f}-{hi:.0f} m"] = dict(n=int(m.sum()), failure_rate=float(f[m].mean()))
        lo = hi
    by_func = {CG.FUNC_NAMES[i]: dict(n=int((fn == i).sum()), failure_rate=float(f[fn == i].mean()))
               for i in np.unique(fn)}
    rates = [v["failure_rate"] for v in by_dist.values()]
    near_contamination = len(rates) >= 2 and rates[0] > 2 * max(rates[-1], 0.01)
    return dict(by_distance_to_lanid=by_dist, by_functional=by_func,
                reading=("failures concentrate near irrigation: contamination" if near_contamination else
                         "no concentration near irrigation: if V2 fails, suspect the ET or rain product "
                         "(compare with rain-fed flux towers)"))


def block_climate(D, cells):
    """Mean annual P and PET (mm) over the run's years for driver cells (-1 -> NaN)."""
    ys = C.years()
    m = (D["months"] >= np.datetime64(f"{ys[0]}-01")) & (D["months"] <= np.datetime64(f"{ys[-1]}-12"))
    out = {}
    for k in ("P", "PET"):
        with np.errstate(invalid="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                a = np.nanmean(D[k][m], 0) * 12
        out[f"{k}_annual"] = np.where(cells >= 0, a[np.maximum(cells, 0)], np.nan)
    return out


def v4_representativeness(T, G, rdir, sdir, D, key):
    lut = DR.cell_lookup(D)
    vals = {v: dict(ref=[], tgt=[]) for v in C.V4_VARS}
    decs = dict(ref=[], tgt=[])
    for ti, t in enumerate(T["tiles"]):
        g = G[t]
        sl = TL.tile_slice(T, ti)
        R = RB.open_tile(rdir, t)
        S = ST.open_tile(sdir, t)
        iy, ix = T["lab_iy"][sl], T["lab_ix"][sl]
        cell = np.where(iy >= 0, lut[np.maximum(iy, 0), np.maximum(ix, 0)], -1).reshape(g.nby, g.nbx)
        clim = block_climate(D, cell)
        use = (T["valid"][sl] & T["primary"][sl]).reshape(g.nby, g.nbx)
        ref = R[key].astype(bool).any(0) & use
        with np.errstate(invalid="ignore"):
            tgt = use & (S["lanid_dist_m"] >= C.IRR_BUFFER_M)
        for name, m in (("ref", ref), ("tgt", tgt)):
            decs[name].append(R["aridity_decile"][m])
            for v in C.V4_VARS:
                src = clim[v] if v in clim else S.get(v)
                vals[v][name].append(np.asarray(src, "float64")[m] if src is not None else np.full(int(m.sum()), np.nan))
    dr, dt = np.concatenate(decs["ref"]), np.concatenate(decs["tgt"])
    table, flagged = {}, []
    for d in range(1, 11):
        row = {}
        for v in C.V4_VARS:
            a, b = np.concatenate(vals[v]["ref"])[dr == d], np.concatenate(vals[v]["tgt"])[dt == d]
            a, b = a[np.isfinite(a)], b[np.isfinite(b)]
            if a.size < 2 or b.size < 2:
                continue
            sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
            smd = float((a.mean() - b.mean()) / sd) if sd > 0 else 0.0
            row[v] = dict(smd=smd, n_ref=int(a.size), n_target=int(b.size))
            if abs(smd) > C.V4_SMD_MAX:
                flagged.append(f"{d}|{v}")
        if row:
            table[str(d)] = row
    return dict(by_decile=table, beyond_limit=flagged, not_available=["clay", "slope"],
                passes=not flagged)


def p10_valid_share(T, G, rdir, bdir):
    """Share of valid block-months, April-September, 2012 against 2014-2016, per functional group x decile."""
    acc = {}
    for ti, t in enumerate(T["tiles"]):
        g = G[t]
        sl = TL.tile_slice(T, ti)
        R = RB.open_tile(rdir, t)
        use = (T["valid"][sl] & T["primary"][sl]).reshape(g.nby, g.nbx)
        fn = _func(R)
        for k, y in enumerate(R["year"]):
            period = "2012" if int(y) == 2012 else ("2014-2016" if 2014 <= int(y) <= 2016 else None)
            if period is None:
                continue
            B = BK.open_tile_year(bdir, t, int(y))
            if B is None:
                continue
            _, ok = RU.valid_months(B, B["roi_px"].astype("float64"))
            ok = ok[3:9]                                                  # April - September
            for d in range(1, 11):
                for fi in np.unique(fn[k][use]):
                    m = use & (R["aridity_decile"] == d) & (fn[k] == fi)
                    if not m.any():
                        continue
                    key = f"{d}|{CG.FUNC_NAMES[fi]}"
                    a = acc.setdefault(key, {"2012": [0, 0], "2014-2016": [0, 0]})
                    a[period][0] += int(ok[:, m].sum())
                    a[period][1] += int(m.sum()) * 6
    out, low = {}, []
    for key, a in sorted(acc.items()):
        s12 = a["2012"][0] / a["2012"][1] if a["2012"][1] else float("nan")
        s14 = a["2014-2016"][0] / a["2014-2016"][1] if a["2014-2016"][1] else float("nan")
        out[key] = dict(share_2012=s12, share_2014_2016=s14, block_months_2012=a["2012"][1])
        if np.isfinite(s12) and s12 < C.P10_VALID_SHARE_MIN:
            low.append(key)
    return dict(by_stratum=out, below_minimum_2012=low, assessed=bool(any(a["2012"][1] for a in acc.values())))


def ring_counts(T, G, rdir):
    out = {}
    for t in T["tiles"]:
        R = RB.open_tile(rdir, t)
        ring = R["reference_ring"].astype(bool)
        dec = np.broadcast_to(R["aridity_decile"], ring.shape)
        fn = _func(R)
        for d in range(1, 11):
            for fi, name in enumerate(CG.FUNC_NAMES):
                m = ring & (dec == d) & (fn == fi)
                if m.any():
                    out[f"{d}|{name}"] = out.get(f"{d}|{name}", 0) + int(m.any(0).sum())
    return out


def v3_maps(T, G, rdir, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    paths_ = []
    per_year = {}
    for ti, t in enumerate(T["tiles"]):
        sl = TL.tile_slice(T, ti)
        R = RB.open_tile(rdir, t)
        use = (T["valid"][sl] & T["primary"][sl])
        for k, y in enumerate(R["year"]):
            ref = R["reference"][k].ravel().astype(bool)
            ring = R["reference_ring"][k].ravel().astype(bool)
            per_year.setdefault(int(y), []).append((T["lon"][sl][use], T["lat"][sl][use], ref[use], ring[use]))
    for y, parts in sorted(per_year.items()):
        lo = np.concatenate([p[0] for p in parts])
        la = np.concatenate([p[1] for p in parts])
        rf = np.concatenate([p[2] for p in parts])
        rg = np.concatenate([p[3] for p in parts])
        other = ~rf & ~rg
        fig, ax = plt.subplots(figsize=(6, 7), dpi=130)
        ax.scatter(lo[other], la[other], s=0.4, c="#d9d8d3", linewidths=0, label="other aquifer blocks")
        ax.scatter(lo[rg], la[rg], s=0.4, c="#e3a33b", linewidths=0, label=f"ring, class rung only ({int(rg.sum()):,})")
        ax.scatter(lo[rf], la[rf], s=0.4, c="#2a78d6", linewidths=0, label=f"reference ({int(rf.sum()):,})")
        ax.set_title(f"P1.1 rain-only reference blocks, {y}", fontsize=10)
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
        ax.set_aspect(1 / np.cos(np.radians(np.nanmean(la))) if la.size else 1)
        ax.legend(loc="lower left", fontsize=8, markerscale=12, frameon=False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        p = f"{outdir}/reference_{y}.png"
        fig.savefig(p, bbox_inches="tight")
        plt.close(fig)
        paths_.append(os.path.basename(p))
    return paths_


def gate_g1a(rep, nwi_missing=()):
    reasons = []
    main = rep["reference"]
    if not main["V2"]["passes"]:
        reasons.append(f"V2: {100 * main['V2']['share']:.1f} % of reference block-years within the water bound "
                       f"(needed {100 * C.V2_PASS_SHARE:.0f} %); see the V2 diagnosis")
    if main["V1"]["strata_in_domain"] == 0:
        reasons.append(f"V1: no functional group x decile stratum has {C.MIN_REF_BLOCKS_PER_STRATUM} reference blocks")
    reasons += RU.check_depth_accepted()
    not_assessed = ["Kansas registry comparison (WIMAS not on disk)",
                    "ring test (needs the P1.2 models: reported by verify_p12)"]
    if nwi_missing:
        how = "passed every block" if C.nwi_missing_policy() == "pass" else "failed every block"
        not_assessed.append(f"NWI rule: no NWI tile for {', '.join(nwi_missing)} - the rule {how} there")
    return dict(name="G1a", passes=not reasons, reasons=reasons, out_of_domain=main["V1"]["out_of_domain"],
                not_assessed=not_assessed)


def run(ctx):
    T, G = TL.load_table(ctx.dep("tiles"))
    rdir, bdir, sdir = ctx.dep("reference"), ctx.dep("blocks"), ctx.dep("static")
    D = DR.load(ctx.dep("drivers"))
    rep = {}
    for key in KEYS:
        rep[key] = dict(V1=v1_counts(T, G, rdir, key), V2=v2_water_balance(T, G, rdir, bdir, sdir, D, key),
                        V4=v4_representativeness(T, G, rdir, sdir, D, key))
    rep["ring_distinct_blocks"] = ring_counts(T, G, rdir)
    rep["P10_valid_share"] = p10_valid_share(T, G, rdir, bdir)
    outdir = f"{ctx.out}/reports"
    os.makedirs(outdir, exist_ok=True)
    rep["V3_maps"] = v3_maps(T, G, rdir, outdir)
    nwi_missing = [str(t) for t in T["tiles"] if int(ST.open_tile(sdir, t)["nwi_available"].max()) == 0]
    rep["gate"] = gate_g1a(rep, nwi_missing)
    write_text(f"{outdir}/verify_p11.json", json.dumps(rep, indent=1, default=float))
    write_text(f"{outdir}/verify_p11.md", markdown(rep))
    ctx.complete(dict(gate=rep["gate"], V1_strata_in_domain=rep["reference"]["V1"]["strata_in_domain"],
                      V1_out_of_domain=len(rep["reference"]["V1"]["out_of_domain"]),
                      V2_share=rep["reference"]["V2"]["share"], V2_passes=rep["reference"]["V2"]["passes"],
                      V4_passes=rep["reference"]["V4"]["passes"]))
    ctx.log.info("P1.1 verification written to %s; gate G1a %s", outdir, "passes" if rep["gate"]["passes"] else
                 "FAILS: " + "; ".join(rep["gate"]["reasons"]))


def markdown(rep):
    L = ["# P1.1 verification", ""]
    gt = rep["gate"]
    L += [f"**Gate G1a: {'PASS' if gt['passes'] else 'FAIL'}**" + ("" if gt["passes"] else " - " + "; ".join(gt["reasons"])),
          "", "Not assessed: " + "; ".join(gt["not_assessed"]), ""]
    for key in KEYS:
        v1, v2, v4 = rep[key]["V1"], rep[key]["V2"], rep[key]["V4"]
        L += [f"## {key}", "",
              f"V1: {v1['total_distinct']:,} distinct reference blocks; {v1['strata_in_domain']} strata (aridity decile x "
              f"functional group) with >= {C.MIN_REF_BLOCKS_PER_STRATUM} blocks; {len(v1['out_of_domain'])} out of domain.", "",
              "| decile | functional group | distinct blocks | block-years |", "|---|---|---|---|"]
        for k in sorted(v1["distinct_blocks"], key=lambda s: (int(s.split('|')[0]), s)):
            d, c = k.split("|")
            L.append(f"| {d} | {c} | {v1['distinct_blocks'][k]:,} | {v1['block_years'][k]:,} |")
        q = ", ".join(f"{x:.2f}" for x in v2["et_over_available_q05_q50_q95"])
        L += ["", f"V2: {v2['passing']:,} of {v2['block_years']:,} reference block-years "
                  f"({100 * v2['share']:.1f} %) have water-year ET <= (P - dS) x {1 + C.V2_TOLERANCE:.2f}; "
                  f"needed {100 * C.V2_PASS_SHARE:.0f} %: **{'PASS' if v2['passes'] else 'FAIL'}**. "
                  f"ET / (P - dS) quantiles 5/50/95 %: {q}", ""]
        dg = v2.get("diagnosis") or {}
        if dg:
            L += [f"V2 diagnosis: {dg['reading']}", "", "| distance to LANID | block-years | failure rate |", "|---|---|---|"]
            L += [f"| {k} | {v['n']:,} | {100 * v['failure_rate']:.1f} % |" for k, v in dg["by_distance_to_lanid"].items()]
            L += [""]
        L += [f"V4 representativeness (|SMD| <= {C.V4_SMD_MAX}): **{'PASS' if v4['passes'] else 'FAIL'}**"
              + (f"; beyond the limit: {', '.join(v4['beyond_limit'])}" if v4["beyond_limit"] else "")
              + f". Not available: {', '.join(v4['not_available'])}.", ""]
    rc = rep["ring_distinct_blocks"]
    L += ["## Ring (blocks only the class rung admits)", "",
          ", ".join(f"{k}: {v:,}" for k, v in sorted(rc.items())) or "none", ""]
    pv = rep["P10_valid_share"]
    L += ["## P1.0 valid share, April-September", ""]
    if pv["assessed"]:
        L += [f"Strata below {100 * C.P10_VALID_SHARE_MIN:.0f} % in 2012: {', '.join(pv['below_minimum_2012']) or 'none'}", "",
              "| stratum | 2012 | 2014-2016 |", "|---|---|---|"]
        L += [f"| {k} | {v['share_2012']:.2f} | {v['share_2014_2016']:.2f} |" for k, v in pv["by_stratum"].items()]
    else:
        L += ["2012 is not among the processed years: not assessed."]
    L += ["", "V3 maps: " + ", ".join(rep["V3_maps"]), ""]
    return "\n".join(L)
