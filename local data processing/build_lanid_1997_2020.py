"""LANID annual irrigation 1997-2020 over the Ogallala aquifer: the v3 file plus the 2018-2020 maps, as a NEW file.

    python "local data processing/build_lanid_1997_2020.py" --check     checks only (writes a JSON report)
    python "local data processing/build_lanid_1997_2020.py"             checks, then builds and verifies the file

Inputs (only read, never changed):
  merged_datasets/LANID_v3_30m_1997-2017_Ogallala.nc     LANID-US v3.0 (Zenodo 5548555), cut to the aquifer
  pilot test/LANID extended/lanid2018.tif, lanid2019.tif, lanid2020.tif     annual maps 2018-2020 (CONUS)
  pilot test/LANID extended/maxIrrExt_0020.tif           maximum irrigated extent 2000-2020 (CONUS)
Output (new):  merged_datasets/LANID_30m_1997-2020_Ogallala.nc
  irrigated(time=1997..2020, y, x) on the v3 file's grid, window and aquifer mask; 1997-2017 copied value for value,
  2018-2020 added; group `derived` = the v3 derived layers unchanged + maxIrrExt_0020; group `derived_irrAreaGain`
  copied unchanged.

Grid. The v3 file and maxIrrExt_0020 share one 30 m lattice (pixel edges at 15 m mod 30): maxIrrExt_0020 is read by
a whole-pixel shift. The 2018-2020 maps lie on a lattice offset by about -4.99 m (x) and -6.06 m (y), so no
whole-pixel shift aligns them: each v3 cell takes the value of the 2018-2020 pixel that contains its CENTRE
(nearest neighbour; the two pixels overlap by about 66%). This is the only resampling in the file, it is written
into the file's attributes, and step C-3 checks that this mapping (no extra shift) aligns best.
Values. v3: 1 irrigated, 0 not irrigated, 255 outside the aquifer. The 2018-2020 maps hold 0/1 with 0 also tagged
as nodata; inside the aquifer (all land) 0 is read as "not irrigated", as in v3.
"""
import datetime as _dt
import hashlib
import json
import math
import os
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))).replace("\\", "/")
SRC_NC = f"{REPO}/merged_datasets/LANID_v3_30m_1997-2017_Ogallala.nc"
EXT = f"{REPO}/pilot test/LANID extended"
NEW_YEARS = [2018, 2019, 2020]
OUT = f"{REPO}/merged_datasets/LANID_30m_1997-2020_Ogallala.nc"
SCRATCH = os.environ.get("LANID_SCRATCH", "H:/MSU_GWB/pilot_scratch/lanid_build")
STRIP = 512                      # rows per strip (= the chunk height of the file)
CHECK_EVERY = 8                  # checks read one strip in every CHECK_EVERY (12.5 % of the rows)
FILL = np.uint8(255)


def log(*a):
    print(_dt.datetime.now().strftime("%H:%M:%S"), *a, flush=True)


def md5_full(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(16 << 20), b""):
            h.update(b)
    return h.hexdigest()


def nc_grid(ds):
    gt = [float(v) for v in str(ds["crs"].GeoTransform).split()]
    return gt[0], gt[3], gt[1], gt[5]                       # x0, y0 (upper-left corner), dx, dy


class NewMaps:
    """The 2018-2020 maps (and maxIrrExt_0020) read onto rows r0..r1 of the v3 window."""

    def __init__(self, x0, y0, w):
        import rasterio
        self.src = {y: rasterio.open(f"{EXT}/lanid{y}.tif") for y in NEW_YEARS}
        self.mx = rasterio.open(f"{EXT}/maxIrrExt_0020.tif")
        self.x0, self.y0, self.w = x0, y0, w
        t = self.src[NEW_YEARS[0]].transform
        for r in list(self.src.values()):
            if r.transform != t or r.crs != self.src[NEW_YEARS[0]].crs:
                raise SystemExit("the 2018-2020 maps are not on one grid")
        self.t = t
        xc = x0 + (np.arange(w) + 0.5) * 30.0                # v3 cell centres
        self.cols = np.floor((xc - t.c) / t.a).astype(np.int64)
        mt = self.mx.transform
        dc, dr = (x0 - mt.c) / mt.a, (y0 - mt.f) / mt.e
        if abs(dc - round(dc)) > 1e-6 or abs(dr - round(dr)) > 1e-6:
            raise SystemExit("maxIrrExt_0020 is not on the v3 lattice")
        self.mcol0, self.mrow0 = int(round(dc)), int(round(dr))

    def rows_of(self, r0, r1, shift_r=0):
        yc = self.y0 - (np.arange(r0, r1) + 0.5) * 30.0
        return np.floor((yc - self.t.f) / self.t.e).astype(np.int64) + shift_r

    def read(self, year, r0, r1, shift_r=0, shift_c=0):
        """Values of the map pixel containing each v3 cell centre (shifted by whole pixels for the shift test)."""
        from rasterio.windows import Window
        rows = self.rows_of(r0, r1, shift_r)
        cols = self.cols + shift_c
        ra, rb, ca, cb = rows.min(), rows.max() + 1, cols.min(), cols.max() + 1
        a = self.src[year].read(1, window=Window(ca, ra, cb - ca, rb - ra), boundless=True, fill_value=255)
        return a[np.ix_(rows - ra, cols - ca)]

    def read_max(self, r0, r1):
        from rasterio.windows import Window
        return self.mx.read(1, window=Window(self.mcol0, self.mrow0 + r0, self.w, r1 - r0), boundless=True,
                            fill_value=255)


def jaccard(a, b):
    u = np.count_nonzero(a | b)
    return np.count_nonzero(a & b) / u if u else float("nan")


# ═══════════════════════════════════════════════════════════════════ STEP C  checks
def check(ds, nm):
    """C-1 irrigated share by year; C-2 agreement between consecutive years (v3-v3 baseline vs the 2017/2018 seam
    vs new-new); C-3 shift test of the nearest-centre mapping; C-4 containment in maxIrrExt_0020."""
    H, W = ds["roi_mask"].shape
    years = [int(v) for v in ds["time"][:]]
    acc = {k: 0 for k in ("n",)}
    irr = {y: 0 for y in years + NEW_YEARS}
    inter, union = {}, {}
    pairs = [(2015, 2016), (2016, 2017), (2017, 2018), (2018, 2019), (2019, 2020)]
    shifts = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1)]
    sh_i, sh_u = {s: 0 for s in shifts}, {s: 0 for s in shifts}
    cont = dict(new_irr=0, new_irr_in_max=0, v3max=0, v3max_in_max=0, v3_2017=0, v3_2017_in_max=0, max_px=0)
    v = ds["irrigated"]
    dmax = ds["derived"]["maxIrrExt"]
    for k, r0 in enumerate(range(0, H, STRIP)):
        if k % CHECK_EVERY:
            continue
        r1 = min(r0 + STRIP, H)
        roi = np.asarray(ds["roi_mask"][r0:r1]).astype(bool)
        if not roi.any():
            continue
        acc["n"] += int(roi.sum())
        lab = {}
        for y in years:
            if y >= 2010:
                lab[y] = (np.asarray(v[years.index(y), r0:r1]) == 1) & roi
                irr[y] += int(lab[y].sum())
        for y in NEW_YEARS:
            lab[y] = (nm.read(y, r0, r1) == 1) & roi
            irr[y] += int(lab[y].sum())
        for a, b in pairs:
            inter[(a, b)] = inter.get((a, b), 0) + int((lab[a] & lab[b]).sum())
            union[(a, b)] = union.get((a, b), 0) + int((lab[a] | lab[b]).sum())
        for s in shifts:
            x = (nm.read(2018, r0, r1, *s) == 1) & roi
            sh_i[s] += int((lab[2017] & x).sum())
            sh_u[s] += int((lab[2017] | x).sum())
        mx = (nm.read_max(r0, r1) == 1) & roi
        anynew = lab[2018] | lab[2019] | lab[2020]
        cont["new_irr"] += int(anynew.sum())
        cont["new_irr_in_max"] += int((anynew & mx).sum())
        v3m = (np.asarray(dmax[r0:r1]) == 1) & roi
        cont["v3max"] += int(v3m.sum())
        cont["v3max_in_max"] += int((v3m & mx).sum())
        cont["v3_2017"] += int(lab[2017].sum())
        cont["v3_2017_in_max"] += int((lab[2017] & mx).sum())
        cont["max_px"] += int(mx.sum())
        log(f"   check strip rows {r0}-{r1}")
    rep = dict(
        sampled_cells_in_aquifer=acc["n"],
        irrigated_share={y: irr[y] / acc["n"] for y in sorted(irr) if irr[y] or y in NEW_YEARS},
        jaccard_consecutive_years={f"{a}-{b}": inter[(a, b)] / union[(a, b)] for a, b in pairs},
        shift_test_2017_vs_2018={f"{dr:+d},{dc:+d}": sh_i[(dr, dc)] / sh_u[(dr, dc)] for dr, dc in shifts},
        share_of_2018_2020_irrigated_inside_maxIrrExt_0020=cont["new_irr_in_max"] / max(cont["new_irr"], 1),
        share_of_v3_2017_irrigated_inside_maxIrrExt_0020=cont["v3_2017_in_max"] / max(cont["v3_2017"], 1),
        share_of_v3_maxIrrExt_1997_2017_inside_maxIrrExt_0020=cont["v3max_in_max"] / max(cont["v3max"], 1),
        maxIrrExt_0020_share=cont["max_px"] / acc["n"])
    best = max(rep["shift_test_2017_vs_2018"], key=rep["shift_test_2017_vs_2018"].get)
    rep["best_shift"] = best
    return rep


# ═══════════════════════════════════════════════════════════════════ STEP B  build
def copy_attrs(src, dst, skip=("_FillValue",)):
    for a in src.ncattrs():
        if a not in skip:
            dst.setncattr(a, src.getncattr(a))


def build(ds, nm, report, out_tmp):
    import netCDF4 as nc
    H, W = ds["roi_mask"].shape
    years = [int(v) for v in ds["time"][:]]
    all_years = years + NEW_YEARS
    zkw = dict(zlib=True, complevel=2, shuffle=True)
    with nc.Dataset(out_tmp, "w", format="NETCDF4") as o:
        # STEP B-1  dimensions, coordinates, grid mapping, mask (copied)
        o.createDimension("time", len(all_years))
        o.createDimension("y", H)
        o.createDimension("x", W)
        tv = o.createVariable("time", "i4", ("time",))
        copy_attrs(ds["time"], tv)
        tv[:] = np.asarray(all_years, "i4")
        for name in ("y", "x"):
            v = o.createVariable(name, ds[name].dtype, (name,))
            copy_attrs(ds[name], v)
            v[:] = ds[name][:]
        g = o.createVariable("crs", "i4")
        copy_attrs(ds["crs"], g)
        m = o.createVariable("roi_mask", "i1", ("y", "x"), chunksizes=(512, 512), **zkw)
        copy_attrs(ds["roi_mask"], m)
        rel = o.createVariable("source_release", str, ("time",))
        rel.long_name = "release each annual map comes from"
        for i, y in enumerate(all_years):
            rel[i] = "LANID-US v3.0 (Zenodo 5548555)" if y <= 2017 else f"LANID {y} (pilot test/LANID extended/lanid{y}.tif)"
        # STEP B-2  irrigated(time, y, x): 1997-2017 copied, 2018-2020 by cell centre, aquifer mask applied
        iv = o.createVariable("irrigated", "u1", ("time", "y", "x"), fill_value=FILL, chunksizes=(1, 512, 512), **zkw)
        copy_attrs(ds["irrigated"], iv)
        iv.long_name = "annual irrigation map: 1 irrigated, 0 not irrigated (LANID lanidYYYY)"
        iv.source_files = " ".join(f"lanid{y}.tif" for y in all_years)
        iv.comment = ("1997-2017: exactly the values of LANID_v3_30m_1997-2017_Ogallala.nc. 2018-2020: value of the "
                      "lanidYYYY.tif pixel containing each cell centre; that grid is offset by about -4.99 m (x) and "
                      "-6.06 m (y) from this one, so this is nearest-neighbour resampling. 255: outside the aquifer.")
        # STEP B-3  groups: v3 derived (unchanged) + maxIrrExt_0020; derived_irrAreaGain copied
        dg = o.createGroup("derived")
        copy_attrs(ds["derived"], dg)
        for name, v in ds["derived"].variables.items():
            nv = dg.createVariable(name, v.dtype, v.dimensions, fill_value=getattr(v, "_FillValue", None),
                                   chunksizes=(512, 512), **zkw)
            copy_attrs(v, nv)
        mxv = dg.createVariable("maxIrrExt_0020", "u1", ("y", "x"), fill_value=FILL, chunksizes=(512, 512), **zkw)
        mxv.long_name = "croplands irrigated at least once during 2000-2020 (binary)"
        mxv.units, mxv.grid_mapping = "1", "crs"
        mxv.flag_values, mxv.flag_meanings = np.array([0, 1], "u1"), "never_irrigated irrigated_at_least_once"
        mxv.source_files = "maxIrrExt_0020.tif"
        mxv.comment = "Same lattice as this file: read by a whole-pixel shift, nothing resampled. 255: outside the aquifer."
        ga = ds["derived_irrAreaGain"]
        og = o.createGroup("derived_irrAreaGain")
        copy_attrs(ga, og)
        for dn, d in ga.dimensions.items():
            og.createDimension(dn, len(d))
        for name, v in ga.variables.items():
            nv = og.createVariable(name, v.dtype, v.dimensions, fill_value=getattr(v, "_FillValue", None))
            copy_attrs(v, nv)
            nv[...] = v[...]
        crcs = {}
        for r0 in range(0, H, STRIP):
            r1 = min(r0 + STRIP, H)
            roi = np.asarray(ds["roi_mask"][r0:r1]).astype(bool)
            m[r0:r1] = ds["roi_mask"][r0:r1]
            for t, y in enumerate(years):
                a = np.asarray(ds["irrigated"][t, r0:r1])
                iv[t, r0:r1] = a
                crcs[(t, r0)] = hashlib.md5(a.tobytes()).hexdigest()
            for k, y in enumerate(NEW_YEARS):
                a = nm.read(y, r0, r1)
                a = np.where(roi, np.where(a == 1, 1, 0), FILL).astype("u1")
                iv[len(years) + k, r0:r1] = a
                crcs[(len(years) + k, r0)] = hashlib.md5(a.tobytes()).hexdigest()
            for name, v in ds["derived"].variables.items():
                dg[name][r0:r1] = v[r0:r1]
            a = nm.read_max(r0, r1)
            mxv[r0:r1] = np.where(roi, np.where(a == 1, 1, 0), FILL).astype("u1")
            if (r0 // STRIP) % 8 == 0:
                log(f"   wrote rows {r0}-{r1} of {H}")
        # STEP B-4  global attributes (processing_status written last)
        copy_attrs(ds, o, skip=("processing_status", "history", "years", "title", "summary", "source_md5"))
        o.title = "LANID annual irrigation maps 1997-2020 and derived layers, Ogallala aquifer, 30 m"
        o.summary = ("LANID annual 30 m irrigation maps 1997-2020 on LANID v3's own Albers grid, cut to the aquifer: "
                     "1997-2017 from LANID-US v3.0, 2018-2020 from the separate annual maps (see irrigated:comment), "
                     "the v3 derived layers (1997-2017) and the 2000-2020 maximum irrigated extent.")
        o.years = " ".join(map(str, all_years))
        o.history = (f"{_dt.datetime.now(_dt.timezone.utc):%Y-%m-%dT%H:%M:%SZ} created by "
                     f"local data processing/build_lanid_1997_2020.py from {os.path.basename(SRC_NC)} and "
                     f"lanid2018-2020.tif, maxIrrExt_0020.tif; inputs only read. | v3 history: {ds.history}")
        o.source_md5 = json.dumps(report["inputs_md5"])
        o.checks = json.dumps({k: report[k] for k in report if k != "inputs_md5"}, default=str)
        o.processing_status = "partial"
    return crcs


# ═══════════════════════════════════════════════════════════════════ STEP V  verify by reading back
def verify(out_tmp, ds, nm, crcs):
    import netCDF4 as nc
    H = ds["roi_mask"].shape[0]
    bad = []
    with nc.Dataset(out_tmp) as o:
        o.set_auto_maskandscale(False)
        n = len(o["time"])
        for (t, r0), h in crcs.items():
            r1 = min(r0 + STRIP, H)
            if hashlib.md5(np.asarray(o["irrigated"][t, r0:r1]).tobytes()).hexdigest() != h:
                bad.append((t, r0))
        for r0 in range(0, H, STRIP * 16):
            r1 = min(r0 + STRIP, H)
            if not np.array_equal(np.asarray(o["roi_mask"][r0:r1]), np.asarray(ds["roi_mask"][r0:r1])):
                bad.append(("roi", r0))
            for name in ds["derived"].variables:
                if not np.array_equal(np.asarray(o["derived"][name][r0:r1]), np.asarray(ds["derived"][name][r0:r1])):
                    bad.append((name, r0))
        if n != 24 or [int(v) for v in o["time"][:]] != list(range(1997, 2021)):
            bad.append("time axis")
    return bad


def main(argv):
    import netCDF4 as nc
    os.makedirs(SCRATCH, exist_ok=True)
    if os.path.exists(OUT) and "--check" not in argv:
        raise SystemExit(f"{OUT} exists; remove or rename it first (this script never overwrites)")
    ds = nc.Dataset(SRC_NC)
    ds.set_auto_maskandscale(False)
    x0, y0, dx, dy = nc_grid(ds)
    H, W = ds["roi_mask"].shape
    nm = NewMaps(x0, y0, W)
    log(f"v3 file {H} x {W}; 2018-2020 grid offset dx={((nm.t.c - x0 + 15) % 30) - 15:+.3f} m, "
        f"dy={((nm.t.f - y0 + 15) % 30) - 15:+.3f} m")
    # ══════ STEP C  checks ══════
    rep = check(ds, nm)
    rep["grid_offset_2018_2020_m"] = [((nm.t.c - x0 + 15) % 30) - 15, ((nm.t.f - y0 + 15) % 30) - 15]
    log("checks:\n" + json.dumps(rep, indent=1))
    rep["inputs_md5"] = {os.path.basename(p): md5_full(p) for p in
                         [SRC_NC] + [f"{EXT}/lanid{y}.tif" for y in NEW_YEARS] + [f"{EXT}/maxIrrExt_0020.tif"]}
    with open(f"{SCRATCH}/lanid_checks.json", "w") as f:
        json.dump(rep, f, indent=1, default=str)
    if "--check" in argv:
        return 0
    if rep["best_shift"] != "+0,+0":
        raise SystemExit(f"the 2018-2020 maps align best with an extra shift {rep['best_shift']}: not built")
    # ══════ STEP B  build on scratch ══════
    tmp = f"{SCRATCH}/{os.path.basename(OUT)}.part"
    if os.path.exists(tmp):
        os.remove(tmp)
    crcs = build(ds, nm, rep, tmp)
    # ══════ STEP V  verify, then mark complete and move into place ══════
    bad = verify(tmp, ds, nm, crcs)
    if bad:
        raise SystemExit(f"read-back check failed at {bad[:10]}; the file was not moved into place")
    with nc.Dataset(tmp, "a") as o:
        o.processing_status = "complete"
    os.replace(tmp, OUT) if os.path.splitdrive(tmp)[0] == os.path.splitdrive(OUT)[0] else __import__("shutil").move(tmp, OUT)
    log(f"DONE {OUT} ({os.path.getsize(OUT) / 1e9:.2f} GB); read-back identical for {len(crcs)} strips")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
