"""Test helpers: temporary config overrides, the synthetic end-to-end setup, stage folders."""
import contextlib
import copy
import os
import shutil
import tempfile

from .. import config as C

SYNTH_ROOT = os.path.join(tempfile.gettempdir(), "paper1_synthetic")


@contextlib.contextmanager
def override(**kw):
    old = {k: copy.deepcopy(getattr(C, k)) for k in kw}
    for k, v in kw.items():
        setattr(C, k, v)
    try:
        yield
    finally:
        for k, v in old.items():
            setattr(C, k, v)


@contextlib.contextmanager
def patched(obj, name, value):
    old = getattr(obj, name)
    setattr(obj, name, value)
    try:
        yield
    finally:
        setattr(obj, name, old)


def synthetic_settings(root=SYNTH_ROOT):
    from . import synth
    S = synth.build(root)
    ds = copy.deepcopy(C.DATASETS)
    for k in ("gridmet", "drought", "dem", "lanid", "gssurgo", "gde", "aquifer"):
        ds[k]["local"] = S[k]
    for k, tag in (("landsat", synth.LANDSAT_TAG), ("openet", synth.OPENET_TAG), ("cdl", synth.CDL_TAG),
                   ("nwi", synth.NWI_TAG)):
        ds[k]["tag"] = tag                         # the synthetic root keeps its own folder names
    return S, dict(RUN_ENV="local", RUN_MODE="test", LOCAL_DATA_ROOT=root, LOCAL_SCRATCH=f"{root}/_scratch",
                   TEST_TILES=[synth.TILE], TEST_YEARS=list(synth.YEARS), TEST_WINDOW_PX=160, TEST_WINDOW_CENTER=None,
                   ARIDITY_BASE_YEARS=(2005, 2011), BUCKET_SPINUP_FROM="2004-01-01", DRIVER_MONTH_FIRST="2005-01",
                   TEST_SPACE_GROUP_DEG=0.01, TEST_SPACE_FOLDS=2, TEST_YEAR_GROUPS=[list(synth.YEARS)],
                   DRYLAND_MIN_PIXELS=50, MIN_REF_BLOCKS_PER_STRATUM=5, LOCATE_BY_NAME=False, DATASETS=ds,
                   SOURCE_SUBFOLDER="", NWI_MISSING_POLICY="fail", TEST_NWI_MISSING_POLICY="fail", LANDSAT_TAG=synth.LANDSAT_TAG,
                   OPENET_TAG=synth.OPENET_TAG, CDL_TAG=synth.CDL_TAG)


@contextlib.contextmanager
def synthetic(root=SYNTH_ROOT):
    """Config + patches that make every stage run on the synthetic data root without network access."""
    import numpy as np
    from ..core import net
    from ..step0 import hand
    S, kw = synthetic_settings(root)

    def fake_download(url, dest, tries=5, timeout=120, log=None):
        if url == C.COUNTIES_URL:
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.copyfile(S["counties_zip"], dest)
            return dest, dict(url=url, size=os.path.getsize(dest), md5="synthetic")
        raise FileNotFoundError(url)

    def fake_dem_window(top, bottom, left, right, log=None):
        import rasterio
        out = np.full((top - bottom + 1, right - left + 1), np.nan)
        with rasterio.open(S["dem"]) as r:
            k0 = int(round((r.transform.f + r.transform.e * 0.5) * 3600))
            m0 = int(round((r.transform.c + r.transform.a * 0.5) * 3600))
            a = r.read(1).astype("float64")
        for i in range(out.shape[0]):
            rr = k0 - (top - i)
            if 0 <= rr < a.shape[0]:
                c0 = left - m0
                cs = np.arange(out.shape[1]) + c0
                ok = (cs >= 0) & (cs < a.shape[1])
                out[i, ok] = a[rr, cs[ok]]
        return out, ["synthetic_dem"]

    def simple_hand(dem, channel, lat_c):
        """Stand-in when pysheds is missing: height above the lowest channel cell of the same row; no D8."""
        h = np.full(dem.shape, np.nan)
        for i in range(dem.shape[0]):
            c = channel[i] & np.isfinite(dem[i])
            if c.any():
                h[i] = dem[i] - dem[i][c].min()
        return h, np.zeros(dem.shape, "uint8")

    try:
        import pysheds  # noqa: F401
        compute = hand.compute_hand
    except ImportError:
        compute = simple_hand
    with override(**kw), patched(net, "download", fake_download), \
            patched(hand, "hu4_codes", lambda bbox, log=None: ["0000"]), \
            patched(hand, "nhd_subset", lambda hu4, log=None: (S["nhd_gpkg"], dict(source=dict(url="synthetic"),
                                                                                   n_channels=1))), \
            patched(hand, "dem_window", fake_dem_window), patched(hand, "compute_hand", compute):
        yield S


def stage_dir(name):
    """Output folder of a stage for the current settings, or None when the stage is not complete."""
    from ..pipeline import runner
    r = runner()
    d = f"{C.out_root()}/{name}/{r.sig(name)}"
    return d if os.path.exists(f"{d}/_complete.json") else None
