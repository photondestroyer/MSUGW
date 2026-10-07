"""SSEBop dekadal ET against flux towers, screened and gap-aware.

    python -m residual_et.towers

Gives the observation-error term of the uncertainty budget. The old comparison
(finding U3) summed whatever half-hours existed, so a dekad with gaps looked
like a dekad with low ET, and it kept LE values up to 14,980 W m-2 and leaked
-9999 sentinels. Here a dekad counts only with enough valid records, and its
ET is the MEAN flux over those records times the length of the dekad.

Independent of the feature store: towers are compared with the mean of the
3 x 3 SSEBop pixels around them.
"""
import json

import numpy as np

from . import config as C
from . import climate as K
from . import sources as S

LE_RANGE = (-100.0, 800.0)       # W m-2; outside this a half-hourly latent heat flux is not physical
LATENT_HEAT = 2.45e6             # J kg-1
MIN_COVERAGE = 0.8               # share of the dekad that must hold valid records


def screen_le(le):
    """NaN for fill, sentinels and non-physical values."""
    le = np.asarray(le, "float64").copy()
    le[~np.isfinite(le) | (le <= -9000) | (le < LE_RANGE[0]) | (le > LE_RANGE[1])] = np.nan
    return le


def dekad_et(seconds, le, starts):
    """Tower ET (mm per dekad) and coverage for each dekad in `starts`.
    seconds: record times as Unix seconds; le: screened latent heat flux (W m-2)."""
    t = np.asarray(seconds, "float64")
    step = float(np.median(np.diff(t))) if t.size > 1 else 1800.0
    s0 = (np.asarray(starts, "datetime64[D]") - np.datetime64("1970-01-01")).astype("float64") * 86400.0
    s1 = (K.dekad_end(starts) - np.datetime64("1970-01-01")).astype("float64") * 86400.0
    et = np.full(len(starts), np.nan)
    cov = np.zeros(len(starts))
    i0, i1 = np.searchsorted(t, s0), np.searchsorted(t, s1)
    for k in range(len(starts)):
        v = le[i0[k]:i1[k]]
        ok = np.isfinite(v)
        cov[k] = ok.sum() * step / (s1[k] - s0[k])
        if cov[k] >= MIN_COVERAGE:
            et[k] = v[ok].mean() * (s1[k] - s0[k]) / LATENT_HEAT
    return et, cov


def site_blocks(ds, chunk=4_000_000):
    """First and one-past-last record of each site (records of a site are contiguous)."""
    n = ds.dimensions["record"].size
    blocks, prev, start = {}, None, 0
    for a in range(0, n, chunk):
        si = np.asarray(ds["site_index"][a:a + chunk])
        cut = np.where(np.diff(si) != 0)[0] + 1
        if prev is not None and si[0] != prev:
            cut = np.concatenate([[0], cut])
        for c in cut:
            sid = int(si[c - 1]) if c > 0 else int(prev)
            if sid in blocks:
                raise S.SourceError(f"ameriflux: records of site index {sid} are not contiguous")
            blocks[sid] = (start, a + int(c))
            start = a + int(c)
        prev = si[-1]
    blocks[int(prev)] = (start, n)
    return blocks


def _scores(tower, sat):
    ok = np.isfinite(tower) & np.isfinite(sat)
    if ok.sum() < 10:
        return dict(n=int(ok.sum()))
    d = sat[ok] - tower[ok]
    return dict(n=int(ok.sum()), bias_ssebop_minus_tower=float(d.mean()), rmse=float(np.sqrt(np.mean(d ** 2))),
                unbiased_rmse=float(d.std()), correlation=float(np.corrcoef(sat[ok], tower[ok])[0, 1]))


def run():
    import pandas as pd
    sites = pd.read_excel(C.SRC["sites"])
    sites.columns = [c.strip() for c in sites.columns]
    r = S.Raster(C.SRC["ssebop"], "et")
    lat0, lat1, lon0, lon1 = r.lat.min(), r.lat.max(), r.lon.min(), r.lon.max()
    sites = sites[(sites["Lat"] > lat0) & (sites["Lat"] < lat1) & (sites["Long"] > lon0) & (sites["Long"] < lon1)]
    am = S.open_ro(C.SRC["ameriflux"])
    try:
        names = [bytes(x).split(b"\x00")[0].decode().strip() for x in np.asarray(am["site_lookup"][:])]
        blocks = site_blocks(am)
        use = []
        for _, row in sites.iterrows():
            sid = str(row["Site_Id"]).strip()
            if sid in names and names.index(sid) in blocks:
                iy, ix = int(np.abs(r.lat - row["Lat"]).argmin()), int(np.abs(r.lon - row["Long"]).argmin())
                use.append((sid, names.index(sid), iy, ix))
        sat = np.full((r.days.size, len(use)), np.nan)
        for a in range(0, r.days.size, 36):                           # one pass over SSEBop
            blk = r.read(a, min(a + 36, r.days.size), slice(None), slice(None))
            for j, (_, _, iy, ix) in enumerate(use):
                with np.errstate(invalid="ignore"):
                    sat[a:a + blk.shape[0], j] = np.nanmean(
                        blk[:, max(iy - 1, 0):iy + 2, max(ix - 1, 0):ix + 2].reshape(blk.shape[0], -1), axis=1)
        tower = np.full_like(sat, np.nan)
        info = []
        for j, (sid, k, _, _) in enumerate(use):
            b0, b1 = blocks[k]
            t = np.asarray(am["time"][b0:b1], "float64")
            le_raw = np.asarray(am["LE"][b0:b1], "float64")
            le = screen_le(le_raw)
            tower[:, j], cov = dekad_et(t, le, r.days)
            info.append(dict(site=sid, records=int(b1 - b0), screened_out=int(np.isfinite(le_raw).sum() - np.isfinite(le).sum()),
                             dekads_used=int(np.isfinite(tower[:, j]).sum()),
                             dekads_rejected_for_gaps=int(((cov > 0) & (cov < MIN_COVERAGE)).sum())))
    finally:
        am.close()
        r.close()
    month = (r.days.astype("datetime64[M]") - r.days.astype("datetime64[Y]").astype("datetime64[M]")).astype(int) + 1
    grow = np.isin(month, C.ANALYSIS_MONTHS)
    out = dict(le_range_w_m2=list(LE_RANGE), min_coverage=MIN_COVERAGE, ssebop_window="3 x 3 pixels",
               n_sites_in_box=len(use), n_sites_with_dekads=int(sum(i["dekads_used"] > 0 for i in info)),
               all_dekads=_scores(tower, sat), growing_season=_scores(tower[grow], sat[grow]),
               sites=info)
    d = C.STORE / "towers"
    d.mkdir(parents=True, exist_ok=True)
    (d / "summary.json").write_text(json.dumps(out, indent=1))
    np.save(d / "tower_et.npy", tower)
    np.save(d / "ssebop_et.npy", sat)
    print(json.dumps({k: v for k, v in out.items() if k != "sites"}, indent=1), flush=True)
    return out


if __name__ == "__main__":
    run()
