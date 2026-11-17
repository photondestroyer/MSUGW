"""The pilot's definitions (GW_ET_High_Plains_pilot.docx, Section 4 and Section 8) as code.

One small, pure function per term. Every task module calls these and never re-implements a rule;
tests/test_definitions.py checks that, and `python -m pilot docs` writes docs/pilot/DEFINITIONS.md from TERMS
below (the "written label rules" deliverable of Task 0.2).

Conventions: arrays are numpy, NaN = unknown. A function takes the decision values it depends on as an
argument `d` (a dict from decisions.yaml) defaulting to the current log, so tests can pass their own.
"""
import numpy as np

from . import config as C
from . import decisions as D


def _d(did, d):
    return D.value(did) if d is None else d


# ===================================================================================== Cell (D03)
def cell_regular(area_ha, d=None):
    """Cell: "A quarter-section of about 65 ha." A quarter whose area lies outside D03's range is kept in the
    cell table but flagged irregular (short sections, correction lines) and excluded from cohort and rainfed."""
    d = _d("D03", d)
    a = np.asarray(area_ha, float)
    return (a >= d["min_ha"]) & (a <= d["max_ha"])


def build_cells(sections, d=None):
    """Sections (GeoDataFrame, any projected CRS) -> cells: each section split into 4 quarters along its own
    axes (common/cells.split_quarters). Adds cell_id = <section id><NE|NW|SW|SE>, area_ha, regular."""
    from .common.cells import split_quarters
    cells = split_quarters(sections)
    cells["area_ha"] = cells.geometry.area / 1e4
    cells["regular"] = cell_regular(cells["area_ha"].to_numpy(), d)
    return cells


# ===================================================================================== Cohort (D04, D05, D06, D09)
def is_irrigated_year(share, d=None):
    """Irrigated cell-year: irrigated share of the cell's mapped pixels >= D04 min_share. NaN -> False."""
    s = np.asarray(share, float)
    return np.nan_to_num(s, nan=-1.0) >= _d("D04", d)["min_share"]


def is_corn_year(share, d=None):
    """Corn cell-year: share of the cell's pixels in CDL corn >= D05 min_share. NaN -> False."""
    s = np.asarray(share, float)
    return np.nan_to_num(s, nan=-1.0) >= _d("D05", d)["min_share"]


def map_share(n_irrigated, n_labelled, n_pixels, d=None):
    """Irrigated share of a cell in one map and year: irrigated pixels / labelled pixels. Unknown (NaN) when
    the map labels fewer than D06 min_map_coverage of the cell's pixels, so a cell cut by the edge of a map
    (LANID holds labels only inside the aquifer polygon) is never labelled from a sliver.
    Returns (share, coverage)."""
    hit, lab, px = (np.asarray(x, float) for x in (n_irrigated, n_labelled, n_pixels))
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = np.where(px > 0, lab / px, 0.0)
        share = np.where(lab > 0, hit / lab, np.nan)
    return np.where(cov >= _d("D06", d)["min_map_coverage"], share, np.nan), cov


def map_roles(d=None):
    """D06: which irrigation map is primary and which is the check, as the column prefixes of the cell table
    ('lanid' or 'aim'). Returns (primary, check)."""
    d = _d("D06", d)
    key = {"LANID": "lanid", "AIM-HPA": "aim"}
    if d["primary"] not in key or d["check"] not in key or d["primary"] == d["check"]:
        raise ValueError(f"D06 primary/check must be LANID and AIM-HPA (one each), not {d['primary']!r}/{d['check']!r}")
    return key[d["primary"]], key[d["check"]]


def maps_agree(primary_share, check_share, d04=None, d06=None):
    """D06: per year, do the primary and the check map give the same irrigated / not-irrigated label (D04 on
    each)? Returns (n_years_agreeing, agree_ok). A year with either map unknown counts as not agreeing."""
    a, l = np.asarray(primary_share, float), np.asarray(check_share, float)
    known = np.isfinite(a) & np.isfinite(l)
    same = known & (is_irrigated_year(a, d04) == is_irrigated_year(l, d04))
    n = same.sum(axis=-1)
    return n, n >= _d("D06", d06)["min_agree_years"]


def served_by_surface_water(nearest_sw_km, d=None):
    """"Not served by surface water" (D09): a WIMAS surface-water point of diversion with irrigation use within
    D09 radius_km of the cell centre. nearest_sw_km: distance to the nearest such point (inf when none)."""
    return np.asarray(nearest_sw_km, float) <= _d("D09", d)["radius_km"]


COHORT_REASONS = ["irregular cell", "no irrigation map", "surface water", "irrigated < 4 of 2011-2015",
                  "corn < 2 of 2013-2015", "irrigation maps disagree"]


def is_cohort(primary_share, check_share, corn_share, surface_water, regular, d04=None, d05=None, d06=None):
    """Cohort: "Cells irrigated in at least four of the five years 2011-2015, with corn in at least two of
    2013-2015, and not served by surface water. The cohort is fixed: a cell stays in it if it later stops
    irrigating or changes crop."
    primary_share, check_share: (n, 5) yearly irrigated shares 2011-2015 of the primary map (the one that
    decides "irrigated", D06: LANID) and of the check map (AIM-HPA; used only for the agreement test);
    corn_share: (n, 3) for 2013-2015; surface_water, regular: (n,) bool.
    Returns (cohort bool (n,), reason str array; '' for cohort cells). The first failing reason is given, in
    the order of COHORT_REASONS. Later years are never read, which is what keeps the cohort fixed."""
    prim = np.asarray(primary_share, float)
    irr_n = is_irrigated_year(prim, d04).sum(axis=1)
    corn_n = is_corn_year(corn_share, d05).sum(axis=1)
    _, agree = maps_agree(prim, check_share, d04, d06)
    tests = [~np.asarray(regular, bool), ~np.isfinite(prim).all(axis=1), np.asarray(surface_water, bool),
             irr_n < C.COHORT_MIN_IRR_YEARS, corn_n < C.COHORT_MIN_CORN_YEARS, ~agree]
    reason = np.full(prim.shape[0], "", dtype=object)
    for why, fails in zip(COHORT_REASONS, tests):
        reason[(reason == "") & fails] = why
    return reason == "", reason


# ===================================================================================== Rainfed cells (D07, D08)
RAINFED_REASONS = ["irregular cell", "no irrigation map", "irrigated in 2011-2015", "not cropland"]


def is_rainfed(primary_share, check_share, cropland_share, regular, d=None):
    """Rainfed cells: "Cropland cells in the same county with no irrigation in 2011-2015" (matching is
    match_rainfed). No irrigation = at most D07 max_irr_share in EVERY year in BOTH maps (both must be known);
    cropland = at least D07 min_cropland_share in every year. Returns (rainfed bool, reason)."""
    d = _d("D07", d)
    aim, lan, crop = (np.asarray(x, float) for x in (primary_share, check_share, cropland_share))
    known = np.isfinite(aim).all(axis=1) & np.isfinite(lan).all(axis=1)
    dry = (np.nan_to_num(aim, nan=1) <= d["max_irr_share"]).all(axis=1) & \
          (np.nan_to_num(lan, nan=1) <= d["max_irr_share"]).all(axis=1)
    cropland = (np.nan_to_num(crop, nan=0) >= d["min_cropland_share"]).all(axis=1)
    tests = [~np.asarray(regular, bool), ~known, ~dry, ~cropland]
    reason = np.full(aim.shape[0], "", dtype=object)
    for why, fails in zip(RAINFED_REASONS, tests):
        reason[(reason == "") & fails] = why
    return reason == "", reason


def irrigated_later(later_shares, d=None):
    """"Check whether any were irrigated later": a cell with more than D07 max_irr_share irrigated in any later
    year of either map. later_shares: (n, years, maps); a map without a year (AIM-HPA ends in 2017) is NaN
    there and does not count. The flag is reported; it never moves a cell between groups (the groups are fixed
    from 2011-2015)."""
    return (np.nan_to_num(np.asarray(later_shares, float), nan=0.0) > _d("D07", d)["max_irr_share"]).any(axis=(1, 2))


def dominant_class(class_shares, names, d=None):
    """Per cell-year: the crop class covering >= D08 dominant_min_share of the cell, else ''.
    class_shares: (n, years, k) shares of the k classes in `names`."""
    d = _d("D08", d)
    s = np.nan_to_num(np.asarray(class_shares, float), nan=0.0)
    k = s.argmax(axis=-1)
    top = np.take_along_axis(s, k[..., None], axis=-1)[..., 0]
    out = np.asarray(names, dtype=object)[k]
    out[top < d["dominant_min_share"]] = ""
    return out


def crop_class(dominant_by_year, d=None):
    """Crop history of 2013-2015: the class dominant in at least D08 dominant_min_years of the years, else
    'rotation' (or 'none' when no year has a dominant crop)."""
    d = _d("D08", d)
    dom = np.asarray(dominant_by_year, dtype=object)
    out = np.full(dom.shape[0], "rotation", dtype=object)
    for i, row in enumerate(dom):
        vals = [v for v in row if v]
        if not vals:
            out[i] = "none"
            continue
        u, c = np.unique(np.asarray(vals, dtype=str), return_counts=True)
        if c.max() >= d["dominant_min_years"]:
            out[i] = u[c.argmax()]
    return out


def match_rainfed(county, crop, rzws_mm, cohort, rainfed, d=None):
    """Rainfed matching: a rainfed cell is matched when a cohort cell of the SAME county and crop class has
    root-zone storage within D08 rzws_tol_mm. Returns (rainfed_matched bool, cohort_n_matches int) per cell."""
    d = _d("D08", d)
    import pandas as pd
    df = pd.DataFrame(dict(county=np.asarray(county), crop=np.asarray(crop), rz=np.asarray(rzws_mm, float),
                           co=np.asarray(cohort, bool), rf=np.asarray(rainfed, bool)))
    matched = np.zeros(len(df), bool)
    n_match = np.zeros(len(df), int)
    tol = d["rzws_tol_mm"]
    for _, g in df[df.co | df.rf].groupby(["county", "crop"]):
        co = g[g.co & np.isfinite(g.rz)]
        rf = g[g.rf & np.isfinite(g.rz)]
        if co.empty or rf.empty:
            continue
        a = np.sort(co.rz.to_numpy())
        lo = np.searchsorted(a, rf.rz.to_numpy() - tol, "left")
        hi = np.searchsorted(a, rf.rz.to_numpy() + tol, "right")
        matched[rf.index.to_numpy()] = hi > lo
        b = np.sort(rf.rz.to_numpy())
        lo = np.searchsorted(b, co.rz.to_numpy() - tol, "left")
        hi = np.searchsorted(b, co.rz.to_numpy() + tol, "right")
        n_match[co.index.to_numpy()] = hi - lo
    return matched, n_match


# ===================================================================================== Deficit, peak month (D01)
def deficit(ref_et_mm, precip_mm):
    """Deficit: "County reference ET minus precipitation from gridMET" (D01 picks grass ETo), in mm."""
    return np.asarray(ref_et_mm, float) - np.asarray(precip_mm, float)


def peak_month(deficit_jun, deficit_jul, deficit_aug):
    """Peak month: "The month among June, July and August with the largest deficit in that county and year.
    It is chosen from weather only." Ties go to the earlier month. NaN in any month -> 0 (undefined)."""
    m = np.stack([deficit_jun, deficit_jul, deficit_aug], axis=-1).astype(float)
    out = np.array(C.PEAK_CANDIDATES)[np.argmax(np.nan_to_num(m, nan=-np.inf), axis=-1)]
    return np.where(np.isfinite(m).all(axis=-1), out, 0)


# ===================================================================================== Dry conditions, episodes (D20, D10)
def season_dry(season_deficit, d=None):
    """Dry season: "seasons above the 70th percentile of seasonal deficit", the percentile taken over all
    Kansas county-years passed in. Returns (dry bool, threshold)."""
    x = np.asarray(season_deficit, float)
    thr = float(np.nanpercentile(x, _d("D20", d)["season_pct"]))
    return np.nan_to_num(x, nan=-np.inf) > thr, thr


def month_dry(S, d=None):
    """Dry month: "months below the 30th percentile of the soil-moisture state S ... Treat a month with missing
    S as dry." Percentile over all Kansas county-months passed in. Returns (dry bool, threshold)."""
    d = _d("D20", d)
    s = np.asarray(S, float)
    thr = float(np.nanpercentile(s, d["s_pct"])) if np.isfinite(s).any() else np.nan
    dry = s < thr
    if d["missing_s_is_dry"]:
        dry = dry | ~np.isfinite(s)
    return dry, thr


def dry_episodes(year, county, dry, neighbours, d=None):
    """Dry episode: "A drought that covers neighbouring counties in the same period counts as one episode."
    Dry county-years of one year that are linked through shared borders form one episode (D10).
    neighbours: {county: set of bordering counties}. Returns (episode id per row, -1 where not dry;
    multi-year id that also merges consecutive dry years sharing a county, reported alongside)."""
    year, county, dry = np.asarray(year), np.asarray(county), np.asarray(dry, bool)
    n = len(year)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        a, b = find(i), find(j)
        if a != b:
            parent[max(a, b)] = min(a, b)

    idx = {(y, c): i for i, (y, c) in enumerate(zip(year, county)) if dry[i]}
    for (y, c), i in idx.items():
        for nb in neighbours.get(c, ()):
            j = idx.get((y, nb))
            if j is not None:
                union(i, j)
    roots = np.array([find(i) for i in range(n)])
    ep = np.full(n, -1)
    lab = {r: k for k, r in enumerate(sorted({roots[i] for i in idx.values()}))}
    for i in idx.values():
        ep[i] = lab[roots[i]]
    for (y, c), i in idx.items():                      # multi-year: same county dry in consecutive years
        j = idx.get((y + 1, c))
        if j is not None:
            union(i, j)
    roots = np.array([find(i) for i in range(n)])
    ep_multi = np.full(n, -1)
    lab = {r: k for k, r in enumerate(sorted({roots[i] for i in idx.values()}))}
    for i in idx.values():
        ep_multi[i] = lab[roots[i]]
    return ep, ep_multi


# ===================================================================================== Aquifer thickness A (D22-D26)
def thickness(wte_m, wte_se_m, base_m, base_se_m):
    """Aquifer thickness, A: "Pre-season (January-April) water-table elevation, from recovered (static) water
    levels, minus the aquifer base elevation, with its uncertainty." The two errors are independent."""
    a = np.asarray(wte_m, float) - np.asarray(base_m, float)
    se = np.sqrt(np.asarray(wte_se_m, float) ** 2 + np.asarray(base_se_m, float) ** 2)
    return a, se


def thickness_thirds(cohort_mean_A, d=None):
    """"Thirds of A (thin, middle, thick) are cut once over the whole cohort and used for display" (D26).
    Returns the two cut points."""
    q = _d("D26", d)["quantiles"]
    return [float(x) for x in np.nanquantile(np.asarray(cohort_mean_A, float), q)]


def assign_third(A, cuts):
    """thin (A <= first cut), middle, thick (A > second cut); '' where A is unknown. Cuts from thickness_thirds."""
    a = np.asarray(A, float)
    out = np.where(a <= cuts[0], "thin", np.where(a <= cuts[1], "middle", "thick")).astype(object)
    out[~np.isfinite(a)] = ""
    return out


def well_support(well_km, d=None):
    """A cell-year has a measured well close enough to know its water table: nearest well of that pre-season
    within D25 max_well_km (the pilot's 8 km support distance)."""
    return np.asarray(well_km, float) <= _d("D25", d)["max_well_km"]


def thickness_poor(A, A_se, well_km, third, cuts, d=None):
    """Poor thickness estimate (D25), any of:
    absolute test        SE of A above D25 max_se_m;
    third-confidence     a thin cell-year whose interval A + third_z*SE reaches the thick cut, or a thick one whose
                         A - third_z*SE reaches the thin cut (it could belong to the opposite extreme third);
                         middle cells are not tested (the models use A as continuous);
    support              nearest measured well beyond max_well_km; A not positive or unknown.
    third: the cell's third (assign_third), cuts: the two cut points (thickness_thirds)."""
    d = _d("D25", d)
    a, se, km = (np.asarray(x, float) for x in (A, A_se, well_km))
    th = np.asarray(third, dtype=object)
    z = d["third_z"]
    with np.errstate(invalid="ignore"):
        conflict = ((th == "thin") & (a + z * se > cuts[1])) | ((th == "thick") & (a - z * se < cuts[0]))
        return ~np.isfinite(a) | (a <= 0) | (se > d["max_se_m"]) | conflict | ~well_support(km, d)


# ===================================================================================== Footprints, S (D11-D19)
def footprint_value(values, d=None):
    """A footprint's soil moisture on a date: mean of its valid nested 9 km retrievals, NaN unless at least
    D11 min_nested_9km of the 16 are valid. values: (..., 16)."""
    v = np.asarray(values, float)
    n = np.isfinite(v).sum(axis=-1)
    s = np.nansum(v, axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    return np.where(n >= _d("D11", d)["min_nested_9km"], m, np.nan)


def footprint_observed(n_nested_in_file, d=None):
    """The SMAP file holds enough of a footprint to ever give it a value: at least D11 min_nested_9km of its 16
    nested 9 km cells lie on the file's grid with data."""
    return np.asarray(n_nested_in_file) >= _d("D11", d)["min_nested_9km"]


def irrigation_free(irr_max, label_cov, d=None):
    """The irrigation part of the reference-footprint rule (D12): largest yearly irrigated share under
    max_irr_share, with the label map covering at least min_label_coverage of the footprint."""
    d = _d("D12", d)
    return (np.asarray(irr_max, float) < d["max_irr_share"]) & (np.asarray(label_cov, float) >= d["min_label_coverage"])


def month_third(day, days_in_month):
    """1 for days 1-10, 2 for 11-20, 3 for 21-end of month."""
    day = np.asarray(day)
    return np.where(day <= 10, 1, np.where(day <= 20, 2, 3))


def is_supported_month(valid_days, days_in_month):
    """Supported month: "For one footprint: valid retrievals on at least eight different dates, with at least one
    in each third of the month. Count a date once when it has both a morning and an evening retrieval. Do not
    combine footprints to reach eight dates." valid_days: day-of-month numbers of ONE footprint (duplicates
    allowed). Returns (supported, n_dates, thirds_ok)."""
    days = np.unique(np.asarray(valid_days, int))
    n = len(days)
    thirds_ok = set(month_third(days, days_in_month).tolist()) >= {1, 2, 3} if n else False
    return (n >= 8) and thirds_ok, n, thirds_ok


def is_reference_footprint(irr_max, label_cov, dist_km, rain_ratio, rzws_ratio, rzws_known, radius_km,
                           d12=None, d14=None):
    """Reference footprint: "A SMAP footprint of about 36 km with less than 5% irrigated area, within 100 km of
    the county centroid, with similar rainfall and soils." For one county; arrays over footprints.
    irr_max: largest yearly irrigated share 2011-2015 (D12); label_cov: share of the footprint the label map
    covers; rain_ratio / rzws_ratio: footprint / county normal; rzws_known: soils cover enough of it.
    Returns (bool, reason)."""
    d12, d14 = _d("D12", d12), _d("D14", d14)
    irr, cov, km, rr, sr = (np.asarray(x, float) for x in (irr_max, label_cov, dist_km, rain_ratio, rzws_ratio))
    known = np.asarray(rzws_known, bool)
    tests = [("beyond radius", ~(km <= radius_km)),
             ("irrigation map coverage", ~(cov >= d12["min_label_coverage"])),
             ("irrigated >= 5%", ~irrigation_free(irr, np.ones_like(irr), d12)),
             ("rainfall not similar", ~(np.abs(rr - 1) <= d14["rain_tol_frac"])),
             ("soils not similar", known & ~(np.abs(sr - 1) <= d14["rzws_tol_frac"]))]
    reason = np.full(irr.shape, "", dtype=object)
    for why, fails in tests:
        reason[(reason == "") & fails] = why
    return reason == "", reason


def footprint_radius(n_ref_at_radius, d=None):
    """D13: the radius used for every county. n_ref_at_radius: {radius_km: array of reference-footprint counts
    per county}. 100 km unless fewer than 80% of counties reach 3, then 150 km (the one allowed widening)."""
    d = _d("D13", d)
    share = float(np.mean(np.asarray(n_ref_at_radius[d["radius_km"]]) >= d["min_footprints"]))
    return (d["radius_km"] if share >= d["min_share_counties"] else d["fallback_km"]), share


def anomaly(monthly_mean, clim_mean, clim_sd):
    """Standardised by calendar month: (value - mean) / SD of that footprint and calendar month."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return (np.asarray(monthly_mean, float) - clim_mean) / np.where(np.asarray(clim_sd) > 0, clim_sd, np.nan)


def soil_moisture_state(anomalies, supported, d15=None, d19=None):
    """Soil-moisture state, S: "SMAP Level-3 soil moisture (SPL3SMP_E) averaged over reference footprints, as an
    anomaly standardised by calendar month." For one county-month: equal weights (D15) over its SUPPORTED
    reference footprints; NaN unless at least D19 min_supported_footprints. Returns (S, n_footprints)."""
    _d("D15", d15)
    a = np.asarray(anomalies, float)
    ok = np.asarray(supported, bool) & np.isfinite(a)
    n = int(ok.sum())
    if n < _d("D19", d19)["min_supported_footprints"]:
        return np.nan, n
    return float(a[ok].mean()), n


# ===================================================================================== ET support (D21, D29)
def et_model_ok(valid_share, d=None):
    """A model counts for a cell-month when at least D21 min_valid_share of the cell's pixels hold a value with
    a positive observation count."""
    return np.nan_to_num(np.asarray(valid_share, float), nan=0.0) >= _d("D21", d)["min_valid_share"]


def et_supported(n_models, d=None):
    """Supported ET cell-month: at least D21 min_models models count (pilot target: four)."""
    return np.asarray(n_models) >= _d("D21", d)["min_models"]


# ===================================================================================== Sample masks (Task 1.7)
MASKS = ["mask_sm", "mask_peak", "mask_common", "mask_timing"]


def sample_masks(role, et_ok, s_valid, a_ok, is_peak):
    """The four sample masks of Task 1.7, per cell-month row (May-September):
      mask_sm      soil-moisture contrast: peak month, ET supported, county S valid, thickness usable (cohort)
      mask_peak    peak-demand contrast:   peak month, ET supported, thickness usable (cohort)
      mask_common  in both
      mask_timing  monthly timing model:   every season month, ET supported, thickness usable (cohort)
    Rainfed rows follow the same rules without the thickness condition (they give the baseline)."""
    role = np.asarray(role, dtype=object)
    et_ok, s_valid, a_ok, is_peak = (np.asarray(x, bool) for x in (et_ok, s_valid, a_ok, is_peak))
    cohort = role == "cohort"
    rainfed = role == "rainfed"
    a_ok = np.where(cohort, a_ok, True)
    base = (cohort | rainfed) & et_ok & a_ok
    sm = base & is_peak & s_valid
    peak = base & is_peak
    return dict(mask_sm=sm, mask_peak=peak, mask_common=sm & peak, mask_timing=base)


# ===================================================================================== Contrasts and models (Section 8)
MODELS = dict(
    C1=dict(adds=["A:S_peak"], y="peak-month ET (mm)", fe=["cell", "county_year"]),
    P1=dict(adds=["A:D_peak"], y="peak-month ET (mm)", fe=["cell", "county_year"]),
    J1=dict(adds=["A:S_peak", "A:D_peak"], y="peak-month ET (mm)", fe=["cell", "county_year"]),
    T1=dict(adds=["A:X_month", "A:is_peak", "A:X_month:is_peak"], y="monthly ET May-Sep (mm)",
            fe=["cell_month", "county_month_year"]),
)
COMMON_TERMS = ["A", "A:D_season", "rzws:D_season", "crop:D_season"]   # "Terms in every model"


def contrast_sm(beta_AxS, A_thin, A_thick, S_p10, S_p50):
    """Soil-moisture contrast: thin-minus-thick peak-month ET at dry S (10th percentile) minus the same at
    median S. In a model with b * A * S that is b * (A_thin - A_thick) * (S_p10 - S_p50)."""
    return beta_AxS * (A_thin - A_thick) * (S_p10 - S_p50)


def contrast_peak(beta_AxD, A_thin, A_thick, D_p90, D_p50):
    """Peak-demand contrast: thin-minus-thick peak-month ET at a high peak-month deficit (90th percentile)
    minus the same at the median: b * (A_thin - A_thick) * (D_p90 - D_p50)."""
    return beta_AxD * (A_thin - A_thick) * (D_p90 - D_p50)


# ===================================================================================== registry (docs + tests)
TERMS = [
    ("Cell", "A quarter-section of about 65 ha.", ["build_cells", "cell_regular"], ["D03"], "0.2a"),
    ("Cohort", "Cells irrigated in at least four of the five years 2011-2015, with corn in at least two of "
     "2013-2015, and not served by surface water. The cohort is fixed.",
     ["map_share", "map_roles", "is_irrigated_year", "is_corn_year", "maps_agree", "served_by_surface_water",
      "is_cohort"],
     ["D04", "D05", "D06", "D09"], "0.2c"),
    ("Rainfed cells", "Cropland cells in the same county with no irrigation in 2011-2015, matched to the cohort "
     "on crop history and root-zone water storage.",
     ["is_rainfed", "irrigated_later", "dominant_class", "crop_class", "match_rainfed"],
     ["D07", "D08"], "0.2c"),
    ("Deficit", "County reference ET minus precipitation from gridMET, monthly and May-September sums.",
     ["deficit"], ["D01"], "1.1"),
    ("Peak month", "The month among June, July and August with the largest deficit in that county and year.",
     ["peak_month"], [], "1.1"),
    ("Dry episode", "A drought that covers neighbouring counties in the same period counts as one episode.",
     ["dry_episodes"], ["D10"], "1.1"),
    ("Dry conditions", "Seasons above the 70th percentile of seasonal deficit and months below the 30th "
     "percentile of S; a month with missing S is dry.", ["season_dry", "month_dry"], ["D20"], "dry"),
    ("Aquifer thickness, A", "Pre-season water-table elevation from recovered levels minus the aquifer base "
     "elevation, with its uncertainty. Thirds cut once over the cohort, for display.",
     ["thickness", "thickness_thirds", "assign_third", "well_support", "thickness_poor"],
     ["D22", "D23", "D24", "D25", "D26"], "2.4"),
    ("Soil-moisture state, S", "SMAP SPL3SMP_E averaged over reference footprints, as an anomaly standardised by "
     "calendar month.", ["footprint_value", "anomaly", "soil_moisture_state"], ["D11", "D15", "D16", "D17", "D18",
                                                                                 "D19"], "1.3"),
    ("Reference footprint", "A SMAP footprint of about 36 km with less than 5% irrigated area, within 100 km of "
     "the county centroid, with similar rainfall and soils.",
     ["is_reference_footprint", "irrigation_free", "footprint_observed", "footprint_radius"],
     ["D11", "D12", "D13", "D14"], "1.2"),
    ("Supported month", "Valid retrievals on at least eight different dates, at least one in each third of the "
     "month; footprints are not combined.", ["is_supported_month", "month_third"], ["D16", "D17"], "1.3"),
    ("Supported ET cell-month", "At least four models cover the cell.", ["et_model_ok", "et_supported"],
     ["D21", "D29"], "1.4"),
    ("Sample masks", "Soil-moisture contrast, peak-demand contrast, common sample, monthly timing.",
     ["sample_masks"], [], "1.7"),
    ("Contrasts and models", "Thin-minus-thick ET at dry S (or high peak deficit) minus the same at the median; "
     "models C1, P1, J1, T1.", ["contrast_sm", "contrast_peak"], ["D28"], "2.6, 2.9"),
]
