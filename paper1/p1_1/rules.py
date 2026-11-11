"""P1.1 rules, one bit each, evaluated per block-year from exposure / administrative data only.

No rule looks at the vegetation signal it later serves as a reference for, except two data checks:
open water (MNDWI and AWEIsh, water indices) and DATA (enough valid months of the three targets).

bit  name        passes when
 0   ROI         >= BLOCK_MIN_ROI_FRACTION of the block is inside the aquifer
 1   WTD         HRES-WTD min over the block >= the depth rule USED for the block (below) and valid coverage
                 >= WTD_MIN_COVERAGE
 2   WTD15       same with WTD_SENSITIVITY_M[0] (ladder rung, stored)
 3   WTD60       same with WTD_SENSITIVITY_M[1] (ladder rung, stored)
 4   HAND        share of pixels with HAND > HAND_MIN_M >= HAND_PASS_FRACTION (and coverage >= that share)
 5   NWI         wetland + riparian pixels <= WETLAND_NOISE_PX (no NWI tile -> NWI_MISSING_POLICY)
 6   WATER       no month with >= WATER_NOISE_PX open-water pixels and <= WATER_MAX_MONTHS months with any
 7   IRR         no LANID ever-irrigated pixel within IRR_BUFFER_M, and LANID covers that whole distance
 8   REGISTRY    no registered irrigation well / point of diversion within IRR_BUFFER_M (none configured -> pass)
 9   CLASS       dominant CDL group has purity >= CLASS_PURITY_MIN and is natural, or a crop that is
                 dryland-eligible in the block's county and year
10   DATA        >= MIN_VALID_MONTHS months with ET, NIRv and NDMI each covering >= TARGET_MIN_PIXEL_FRACTION
11   GDE         max GDE fraction <= GDE_SENSITIVITY_MAX, or not assessed (sensitivity only; user decision)
12   PRIMARY     the counted copy of a block in the zone overlap
plan v3:
13   WTD10       depth >= DEPTH_RUNGS_M[-1] (ladder rung, stored)
14   WTD30       depth >= WTD_MIN_M (the strict rung, stored)
15   CANAL       no NHD canal / ditch within CANAL_BUFFER_M
16   FIRE        no fire perimeter in the year or the year before (slot: passes while FIRE_PERIMETERS is None)
17   IRR3        as IRR, but only pixels irrigated in >= LANID_FREQ_MIN_YEARS years count (sensitivity)
18   CLASSDEPTH  depth >= the CLASS rung of the block (DEPTH_BY_FUNCTIONAL; WOODY_DEPTH_M with woody cover)

Depth rule (user, 2026-11-11: class-specific ladder; ring test reviewed before a rung is used):
    class rung = DEPTH_BY_FUNCTIONAL[functional group], raised to WOODY_DEPTH_M when the block holds
                 >= WOODY_NOISE_PX shrub / forest pixels; WTD_MIN_M for a block without class
    used rung  = max(class rung, DEPTH_ACCEPTED.get("<functional>|<decile>", WTD_MIN_M))
    So the reference uses the strict 30 m everywhere until a stratum is accepted in config, and an accepted
    rung can never be looser than the class rung.
reference       = all of MAIN
reference_gde   = reference and GDE                          (sensitivity, GDE diagnostic)
reference_irr3  = MAIN with IRR replaced by IRR3             (sensitivity, LANID labels of 1-2 years ignored)
reference_ring  = MAIN without WTD, plus CLASSDEPTH, and WTD fails: the blocks only the looser class rung
                  admits - predicted out of fold and compared with the reference in the ring test (P1.2 verify)
run_on_set      = dryland crop blocks failing only HAND (the run-on diagnostic set; was the known-zero set)
"""
import numpy as np

from .. import config as C
from .. import cdl_groups as CG

BITS = dict(ROI=0, WTD=1, WTD15=2, WTD60=3, HAND=4, NWI=5, WATER=6, IRR=7, REGISTRY=8, CLASS=9, DATA=10,
            GDE=11, PRIMARY=12, WTD10=13, WTD30=14, CANAL=15, FIRE=16, IRR3=17, CLASSDEPTH=18)
MAIN = ["ROI", "WTD", "HAND", "NWI", "WATER", "IRR", "REGISTRY", "CANAL", "FIRE", "CLASS", "DATA", "PRIMARY"]
MAIN_MASK = sum(1 << BITS[b] for b in MAIN)
SETS = ["reference", "reference_gde", "reference_irr3", "reference_ring"]


def bit(name):
    return np.uint32(1 << BITS[name])


def mask_of(names):
    return sum(1 << BITS[b] for b in names)


def wtd_pass(S, depth):
    with np.errstate(invalid="ignore"):
        return (S["wtd_min"] >= depth) & (S["wtd_coverage"] >= C.WTD_MIN_COVERAGE)


def class_of(cdl_n):
    """Dominant group and purity from (group, by, bx) counts; 255 / NaN where no CDL pixel."""
    tot = cdl_n.sum(0)
    dom = cdl_n.argmax(0).astype("uint8")
    with np.errstate(invalid="ignore", divide="ignore"):
        purity = np.where(tot > 0, cdl_n.max(0) / tot, np.nan).astype("float32")
    return np.where(tot > 0, dom, 255).astype("uint8"), purity


def class_eligible(dom, county, year, dry):
    elig = np.zeros(dom.shape, bool)
    for gi in np.unique(dom[dom != 255]):
        sel = dom == gi
        kind = CG.ELIGIBILITY[CG.NAMES[gi]]
        if kind == "natural":
            elig[sel] = True
        elif kind == "crop":
            cs = county[sel]
            elig[sel] = [dry.get((int(c), int(year), int(gi)), False) for c in cs]
    return elig


def valid_months(B, roi_px):
    need = np.maximum(roi_px * C.TARGET_MIN_PIXEL_FRACTION, 1)
    ok = (B["ET_n"] >= need) & (B["NIRv_n"] >= need) & (B["NDMI_n"] >= need)
    return ok.sum(0).astype("uint8"), ok


def woody_pixels(cdl_n):
    return cdl_n[CG.WOODY_GROUPS].sum(0)


def depth_rule(func, woody_px, decile):
    """(class rung, used rung) in m per block - see the module docstring."""
    by_f = np.array([C.DEPTH_BY_FUNCTIONAL.get(n, C.WTD_MIN_M) for n in CG.FUNC_NAMES], "float64")
    cls = by_f[np.asarray(func, "int64")]
    cls = np.where(np.asarray(woody_px) >= C.WOODY_NOISE_PX, np.maximum(cls, C.WOODY_DEPTH_M), cls)
    acc = np.full(cls.shape, C.WTD_MIN_M, "float64")
    for key, rung in C.DEPTH_ACCEPTED.items():
        name, dec = key.split("|")
        acc[(np.asarray(func) == CG.FUNC_NAMES.index(name)) & (np.asarray(decile) == int(dec))] = float(rung)
    return cls, np.maximum(cls, acc)


def check_depth_accepted():
    """Problems with config.DEPTH_ACCEPTED (unknown group, decile out of 1-10, rung not on the ladder)."""
    probs = []
    for key, rung in C.DEPTH_ACCEPTED.items():
        try:
            name, dec = key.split("|")
            ok = name in CG.FUNC_NAMES and name != "none" and 1 <= int(dec) <= 10
        except ValueError:
            ok = False
        if not ok:
            probs.append(f"DEPTH_ACCEPTED key {key!r} is not '<functional group>|<decile 1-10>'")
        if float(rung) not in [float(r) for r in C.DEPTH_RUNGS_M]:
            probs.append(f"DEPTH_ACCEPTED[{key!r}] = {rung} is not a rung of DEPTH_RUNGS_M {C.DEPTH_RUNGS_M}")
    return probs


def evaluate(S, B, Cn, valid, primary, year, dry, decile=None):
    """S static (by,bx) arrays; B block monthly stats; Cn class counts; valid/primary (by,bx) bool;
    decile (by,bx) aridity decile of the block (0 unknown) for DEPTH_ACCEPTED."""
    probs = check_depth_accepted()
    if probs:
        raise ValueError("; ".join(probs))
    if C.IRRIGATION_REGISTRIES:
        raise NotImplementedError("irrigation registries are configured but no loader exists yet "
                                  "(see docs/PIPELINE.md, 'registry slot')")
    if C.FIRE_PERIMETERS:
        raise NotImplementedError("FIRE_PERIMETERS is set but no MTBS loader exists yet (plan v3 slot)")
    decile = np.zeros(valid.shape, "int64") if decile is None else np.asarray(decile)
    cdl_n = Cn["cdl_n"] if Cn is not None else np.zeros((CG.N,) + valid.shape, "int64")
    if Cn is not None:
        dom, purity = class_of(cdl_n)
    else:
        dom, purity = np.full(valid.shape, 255, "uint8"), np.full(valid.shape, np.nan, "float32")
    func = CG.functional_of(dom)
    woody = woody_pixels(cdl_n)
    cls_rung, used_rung = depth_rule(func, woody, decile)

    f = np.zeros(valid.shape, "uint32")

    def put(name, ok):
        nonlocal f
        f |= np.where(ok, bit(name), 0).astype("uint32")
    with np.errstate(invalid="ignore"):
        put("ROI", valid)
        put("WTD", wtd_pass(S, used_rung))
        put("WTD15", wtd_pass(S, C.WTD_SENSITIVITY_M[0]))
        put("WTD60", wtd_pass(S, C.WTD_SENSITIVITY_M[1]))
        put("WTD10", wtd_pass(S, C.DEPTH_RUNGS_M[-1]))
        put("WTD30", wtd_pass(S, C.WTD_MIN_M))
        put("CLASSDEPTH", wtd_pass(S, cls_rung))
        put("HAND", (S["hand_frac"] >= C.HAND_PASS_FRACTION) & (S["hand_coverage"] >= C.HAND_PASS_FRACTION))
        wet = S["nwi_wetland_px"].astype(int) + S["nwi_riparian_px"].astype(int)
        put("NWI", np.where(S["nwi_available"] == 1, wet <= C.WETLAND_NOISE_PX, C.nwi_missing_policy() == "pass"))
        w = B["water_px"].astype(int)
        put("WATER", (w.max(0) < C.WATER_NOISE_PX) & ((w > 0).sum(0) <= C.WATER_MAX_MONTHS))
        assessed = S["lanid_unassessed_dist_m"] > C.IRR_BUFFER_M
        put("IRR", (S["lanid_dist_m"] >= C.IRR_BUFFER_M) & assessed)
        put("IRR3", (S["lanid3_dist_m"] >= C.IRR_BUFFER_M) & assessed)
        put("CANAL", S["canal_dist_m"] >= C.CANAL_BUFFER_M)
        put("CLASS", class_eligible(dom, S["county"], year, dry) & (purity >= C.CLASS_PURITY_MIN))
        gde_ok = (S["gde_max"] <= C.GDE_SENSITIVITY_MAX) | ~(S["gde_assessed"] > 0)
    put("REGISTRY", np.ones(valid.shape, bool))
    put("FIRE", np.ones(valid.shape, bool))
    nvm, _ = valid_months(B, S["roi_px"].astype("float64"))
    put("DATA", nvm >= C.MIN_VALID_MONTHS)
    put("GDE", gde_ok)
    put("PRIMARY", primary & valid)
    ref = (f & MAIN_MASK) == MAIN_MASK
    has = lambda name: (f & bit(name)) > 0
    irr3_mask = mask_of([b for b in MAIN if b != "IRR"] + ["IRR3"])
    ring_mask = mask_of([b for b in MAIN if b != "WTD"] + ["CLASSDEPTH"])
    sets = dict(reference=ref, reference_gde=ref & has("GDE"), reference_irr3=(f & irr3_mask) == irr3_mask,
                reference_ring=((f & ring_mask) == ring_mask) & ~has("WTD"))
    out = {k: v.astype("uint8") for k, v in sets.items()}
    out.update(flags=f, group=dom, functional=func.astype("uint8"), purity=purity, valid_months=nvm,
               woody_px=np.minimum(woody, 65535).astype("uint16"), depth_used_m=used_rung.astype("float32"),
               depth_class_m=cls_rung.astype("float32"))
    return out


def run_on_set(flags, group):
    """Dryland CROP blocks that fail only the HAND rule: the run-on diagnostic set of plan v3 (low ground
    receives run-on, so any excess there must scale with rain). Was the known-zero set of revision 2."""
    need = MAIN_MASK & ~int(bit("HAND"))
    only_hand = ((flags & need) == need) & ((flags & bit("HAND")) == 0)
    crop = np.isin(group, CG.CROP_GROUPS)
    return only_hand & crop


def explain(flags):
    """Names of the main rules a flags value fails (for reports)."""
    return [b for b in MAIN if not (int(flags) >> BITS[b]) & 1]
