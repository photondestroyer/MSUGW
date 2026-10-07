"""Checks on the real source files (skipped when a file is absent).
They read small windows only."""
import numpy as np

from .. import config as C
from .. import climate as K
from .. import drought as D
from .. import sources as S

SKIP = "SKIP"


def _need(*keys):
    return all(C.SRC[k].exists() for k in keys)


def test_lab_drought_file_states_its_windows_and_has_no_gap():
    if not _need("drought"):
        return SKIP
    cube = S.LabCube(C.SRC["drought"], step_days=C.PENTAD_MAX_GAP_DAYS)
    try:
        D.check_long_names(cube)
        year = cube.days.astype("datetime64[Y]").astype(int) + 1970
        for y in range(2015, 2023):
            assert (year == y).sum() == 73, f"{y} has {(year == y).sum()} pentads"
        starts = K.dekad_starts(C.PROFILES["full"].model_t0, C.PROFILES["full"].model_t1)
        j = D.lookup_preceding(cube.days, starts)
        assert ((starts - cube.days[j]).astype(int) >= 1).all()
    finally:
        cube.close()


def test_lab_index_signs_on_the_2012_drought():
    if not _need("drought"):
        return SKIP
    cube = S.LabCube(C.SRC["drought"], step_days=C.PENTAD_MAX_GAP_DAYS)
    try:
        i0, i1 = cube.index("2012-05-01", "2012-10-01")
        ysl, xsl = slice(120, 160), slice(100, 140)
        idx = {n: cube.read(n, i0, i1, ysl, xsl) for n in ("spi30d", "spei30d", "eddi30d", "spi90d")}
    finally:
        cube.close()
    r = lambda a, b: float(np.corrcoef(a.ravel(), b.ravel())[0, 1])
    assert r(idx["spi30d"], idx["spei30d"]) > 0.5
    assert r(idx["eddi30d"], idx["spei30d"]) < -0.5, "EDDI must be positive when it is dry"
    d = D.describe_index(idx["spi90d"])
    assert d["n_levels"] <= 40 and abs(d["minimum"] + 2.09) < 0.01, d


def test_lab_gridmet_decodes_to_physical_units():
    if not _need("gridmet"):
        return SKIP
    cube = S.LabCube(C.SRC["gridmet"], step_days=1)
    try:
        S.assert_covers(cube.days, "1985-01-01", C.PROFILES["full"].model_t1, "gridMET")
        i0, i1 = cube.index("2016-07-01", "2016-07-11")
        ysl, xsl = slice(120, 160), slice(100, 140)
        tmmx = cube.read("tmmx", i0, i1, ysl, xsl)
        pr = cube.read("pr", i0, i1, ysl, xsl)
        assert 290 < np.nanmean(tmmx) < 315, "July maximum temperature is not in Kelvin"
        assert 0 <= np.nanmin(pr) and np.nanmax(pr) < 300
    finally:
        cube.close()


def test_old_drought_axis_is_refused():
    """The Earth Engine drought file had no 2020; its time axis survives in zarr/."""
    p = C.ROOT / "zarr" / "drought" / "time.npy"
    if not p.exists():
        return SKIP
    old = np.load(p).astype("datetime64[D]")
    starts = K.dekad_starts("2019-06-01", "2021-06-01")
    try:
        D.lookup_preceding(old, starts)
    except S.SourceError as e:
        assert "2020" in str(e)
        return
    raise AssertionError("the 2020 hole was not detected")


def test_gfsad_classes_are_cropland_classes():
    if not _need("gfsad"):
        return SKIP
    r = S.Raster(C.SRC["gfsad"], "landcover")
    try:
        cls = r.read_classes(0, slice(None), slice(None))
    finally:
        r.close()
    share = np.bincount(cls[cls >= 0], minlength=6) / (cls >= 0).sum()
    assert cls.max() <= 5
    assert 0.25 < share[3:6].sum() < 0.5, "classes 3-5 (rainfed cropland) should cover about a third of the box"
    assert share[1:3].sum() < 0.15, "irrigated cropland (classes 1-2) is a minority of the box"


def test_ssebop_is_dekadal_and_whole_millimetres():
    if not _need("ssebop"):
        return SKIP
    r = S.Raster(C.SRC["ssebop"], "et")
    try:
        S.assert_cadence(r.days, 11, "SSEBop")
        a = r.read(400, 402, slice(600, 640), slice(300, 340))
    finally:
        r.close()
    v = a[np.isfinite(a)]
    assert np.allclose(v, np.round(v)) and v.min() >= 0
