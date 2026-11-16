"""Task 1.3  SMAP dates and the soil-moisture state S.

"For each reference footprint and month, apply the quality flags and list the valid retrieval dates. Flag
supported months. Then compute county monthly means from supported footprints only. ... Do not move the peak
month to avoid a gap, and do not fill a gap with Level 4."
Physical idea: S is the drought as the land around the cohort feels it, without irrigation. It only means that
when the satellite actually saw the footprint often enough in the month - hence the supported-month rule.
Decisions: D11 (footprint value), D15 (weights), D16 (overpass), D17 (quality), D18 (climatology), D19 (valid
county month).   Deliverables: footprint_month.csv, county_month_S.csv, smap_gaps.csv, support_by_month.png.
"""
import numpy as np


def _valid(sm, flag, d17):
    ok = np.isfinite(sm) & (sm >= 0.02) & (sm <= 0.5)
    if d17.get("retrieval_qual_flag_bit0", 0) == 0:
        f = np.ma.filled(flag, 65534).astype("int64")
        ok &= (f != 65534) & ((f & 1) == 0)
    return ok


def run(ctx):
    import pandas as pd
    from .. import config as C
    from .. import decisions as D
    from .. import definitions as DEF
    from ..common import io, report, grids

    d16, d17, d18 = D.value("D16"), D.value("D17"), D.value("D18")
    w = report.read(ctx.dep("t1_2"), "county_reference_weights")
    fps = sorted(w.fp.unique())

    # ══════ STEP 1.3-1  the 16 nested 9 km cells of every reference footprint, located in the file ══════
    p, smap_label = io.smap_source()
    ds = io.smap_open(p)
    r9, c9 = grids.check_ease_axes(ds["x"][:], ds["y"][:])
    rowpos = {r: i for i, r in enumerate(r9)}
    colpos = {c: j for j, c in enumerate(c9)}
    nest = {}
    for f in fps:
        r36, c36 = int(f[1:f.index("C")]), int(f[f.index("C") + 1:])
        cells = [(rowpos[r], colpos[c]) for r in range(4 * r36, 4 * r36 + 4) for c in range(4 * c36, 4 * c36 + 4)
                 if r in rowpos and c in colpos]
        nest[f] = (np.array([a for a, _ in cells], int), np.array([b for _, b in cells], int))
    dates = io.smap_dates(ds)
    keep = np.where((dates.year >= d18["years"][0]) & (dates.year <= d18["years"][1])
                    & dates.month.isin(C.SEASON_MONTHS))[0]
    ctx.log.info("STEP 1.3-1  %d reference footprints; file %s (%s); %d season days in %d-%d", len(fps), p,
                 smap_label, len(keep), *d18["years"])

    # ══════ STEP 1.3-2  daily footprint values: overpass (D16), quality (D17), >= 8 of 16 cells (D11) ══════
    passes = {"am": ["am"], "pm": ["pm"], "both": ["am", "pm"]}[d16["overpass"]]
    vals = {f: np.full(len(keep), np.nan) for f in fps}
    for i0 in range(0, len(keep), 256):
        idx = keep[i0:i0 + 256]
        per_pass = []
        for ps in passes:
            sm = np.ma.filled(ds[f"soil_moisture_{ps}"][idx].astype("float64"), np.nan)
            fl = ds[f"retrieval_qual_flag_{ps}"][idx]
            sm[~_valid(sm, fl, d17)] = np.nan
            per_pass.append(sm)
        for f, (rr, cc) in nest.items():
            if len(rr) == 0:
                continue
            fvals = [DEF.footprint_value(a[:, rr, cc]) for a in per_pass]
            v = fvals[0] if len(fvals) == 1 else np.nanmean(np.vstack(fvals), axis=0)  # a date counts once
            vals[f][i0:i0 + len(idx)] = v
    ds.close()
    day = dates[keep]

    # ══════ STEP 1.3-3  supported months (8 dates, one in each third) and monthly footprint means ══════
    rows = []
    for f in fps:
        v = vals[f]
        ok = np.isfinite(v)
        for (y, m), k in pd.Series(np.arange(len(day))).groupby([day.year, day.month]):
            k = k.to_numpy()
            good = k[ok[k]]
            sup, n, thirds = DEF.is_supported_month(day[good].day.to_numpy(), grids.days_in_month(y, m))
            rows.append(dict(fp=f, year=int(y), month=int(m), n_dates=n, thirds_ok=thirds, supported=sup,
                             sm_mean=float(np.mean(v[good])) if len(good) else np.nan))
    fm = pd.DataFrame(rows)

    # ══════ STEP 1.3-4  climatology by footprint and calendar month (D18) and standardised anomalies ══════
    sup = fm[fm.supported]
    clim = sup.groupby(["fp", "month"]).sm_mean.agg(["mean", "std", "count"]).reset_index()
    clim.loc[clim["count"] < d18["min_years"], ["mean", "std"]] = np.nan
    fm = fm.merge(clim.rename(columns={"mean": "clim_mean", "std": "clim_sd", "count": "clim_years"}),
                  on=["fp", "month"], how="left")
    fm["anomaly"] = np.where(fm.supported, DEF.anomaly(fm.sm_mean, fm.clim_mean, fm.clim_sd), np.nan)
    ctx.log.info("STEP 1.3-4  %.1f%% of footprint-months supported", 100 * fm.supported.mean())

    # ══════ STEP 1.3-5  county S: equal weights over supported reference footprints, >= D19 of them ══════
    out = []
    j = w.merge(fm[fm.year.isin(C.PANEL_YEARS)], on="fp", how="left")
    for (c, y, m), g in j.groupby(["county", "year", "month"]):
        S, n = DEF.soil_moisture_state(g.anomaly.to_numpy(), g.supported.fillna(False).to_numpy())
        out.append(dict(county=c, year=int(y), month=int(m), S=S, n_footprints=n,
                        n_reference=int(len(g)), S_valid=bool(np.isfinite(S))))
    cm = pd.DataFrame(out)
    allc = report.read(ctx.dep("counties"), "counties")
    allc = allc[allc.study.astype(bool)].county.astype(str)
    grid = pd.MultiIndex.from_product([sorted(allc), C.PANEL_YEARS, C.SEASON_MONTHS],
                                      names=["county", "year", "month"]).to_frame(index=False)
    cm = grid.merge(cm, on=["county", "year", "month"], how="left")
    cm["S_valid"] = cm["S_valid"].astype("boolean").fillna(False).astype(bool)
    cm["n_footprints"] = cm["n_footprints"].fillna(0).astype(int)
    cm["S_pct"] = cm["S"].rank(pct=True) * 100

    # ══════ STEP 1.3-6  gaps: season days of the panel with no valid retrieval on any reference footprint ══════
    all_days = pd.date_range(f"{C.PANEL_YEARS[0]}-05-01", f"{C.PANEL_YEARS[-1]}-09-30", freq="D")
    all_days = all_days[all_days.month.isin(C.SEASON_MONTHS)]
    seen = pd.DatetimeIndex(day[np.any(np.vstack([np.isfinite(vals[f]) for f in fps]), axis=0)]) if fps else day[:0]
    missing = all_days.difference(seen)
    gaps = []
    if len(missing):
        run_start = prev = missing[0]
        for t in list(missing[1:]) + [None]:
            if t is None or (t - prev).days > 1:
                gaps.append(dict(start=run_start.date(), end=prev.date(), days=(prev - run_start).days + 1))
                if t is not None:
                    run_start = t
            if t is not None:
                prev = t
    gaps = pd.DataFrame(gaps)
    long_gaps = gaps[gaps.days >= 7] if len(gaps) else gaps
    file_end = str(dates.max().date())
    ctx.log.info("STEP 1.3-6  %d gaps of 7+ days; file ends %s", len(long_gaps), file_end)

    # ══════ STEP 1.3-7  write ══════
    report.table(ctx, fm, "footprint_month", parquet=True, units=dict(sm_mean="m3/m3", anomaly="z"),
                 notes=f"SMAP {smap_label}; overpass {d16['overpass']}")
    report.table(ctx, cm, "county_month_S", parquet=True, units=dict(S="z (equal-weight mean of footprint anomalies)"),
                 notes=f"SMAP {smap_label}")
    report.table(ctx, gaps, "smap_gaps", units={"days": "days"},
                 notes=f"expected from the pilot text: {C.SMAP_GAPS_EXPECTED}; file ends {file_end}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    s = fm[fm.year.isin(C.PANEL_YEARS)].groupby(["year", "month"]).supported.mean().unstack()
    fig, ax = plt.subplots(figsize=(6, 3.5))
    im = ax.imshow(s.to_numpy(), cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_yticks(range(len(s.index)), s.index)
    ax.set_xticks(range(len(s.columns)), s.columns)
    ax.set_xlabel("month")
    fig.colorbar(im, ax=ax, label="share of reference footprints with a supported month")
    report.figure(ctx, fig, "support_by_month", units={"share": "fraction"}, mask="reference footprints",
                  caption=f"SMAP {smap_label}")
    ctx.complete(dict(smap_source=smap_label, file_ends=file_end, footprints=len(fps),
                      supported_share=float(fm.supported.mean()), county_months_valid=float(cm.S_valid.mean()),
                      gaps_7d=long_gaps.astype(str).to_dict("records") if len(long_gaps) else []))
