"""Unit tests of the Step 0 pieces (synthetic, seconds)."""
import os
import tempfile

import numpy as np

from .. import config as C
from .. import cdl_groups as CG
from ..step0 import tiles as TL, blocks as BK, drivers as DR, static as ST
from .util import override


def test_frame_offsets_and_block_window():
    zone, left, top = TL.tile_corner("z14e02n29")
    assert (zone, left, top) == (14, 307215.0, 4608015.0)
    F = dict(crs="EPSG:32614", x0=left + 1601 * 30, y0=top - 1603 * 30, w=100, h=50)
    g = TL.geom_from_frame("z14e02n29", F)
    assert (g.r0, g.c0) == (1603, 1601)
    assert (g.br0, g.bc0) == (100, 100) and g.nby == -(-(1603 + 50) // 16) - 100 and g.nbx == -(-(1601 + 100) // 16) - 100
    r, c = g.frame_rows(), g.frame_cols()
    assert r[0] == 100 * 16 - 1603 and c[0] == 100 * 16 - 1601     # block window starts before the frame
    xs, ys = g.pixel_centres()
    assert xs[0] == left + (1600 + 0.5) * 30 and ys[0] == top - (1600 + 0.5) * 30


def test_frame_off_lattice_is_refused():
    _, left, top = TL.tile_corner("z14e02n29")
    try:
        TL.geom_from_frame("z14e02n29", dict(crs="EPSG:32614", x0=left + 10.0, y0=top, w=10, h=10))
        raise AssertionError("accepted a frame off the 30 m lattice")
    except Exception as e:
        assert "lattice" in str(e)


def test_block_stats_equal_brute_force():
    rng = np.random.default_rng(0)
    a = rng.normal(size=(32, 48)).astype("float32")
    a[rng.random(a.shape) < 0.3] = np.nan
    a[:16, :16] = np.nan
    m, s, n = BK.block_stats(a, 2, 3)
    for i in range(2):
        for j in range(3):
            b = a[i * 16:(i + 1) * 16, j * 16:(j + 1) * 16]
            v = b[np.isfinite(b)]
            assert n[i, j] == v.size
            if v.size:
                assert abs(m[i, j] - v.mean()) < 1e-5 and abs(s[i, j] - v.std()) < 1e-5
            else:
                assert np.isnan(m[i, j])


def test_block_ids_are_unique_and_ordered():
    BR, BC = np.meshgrid(np.arange(319, -1, -1), np.arange(320), indexing="ij")
    ids = TL.block_id(14, 2, 29, BR, BC)
    assert np.unique(ids).size == ids.size
    assert ids.dtype == np.int64


def test_primary_copy_is_unique_in_the_zone_overlap():
    from pyproj import Transformer
    lon0, lat0 = -102.0, 40.0                     # the 13/14 zone boundary
    T = {k: [] for k in ("zone", "lon", "lat", "valid", "band_zone", "tile_idx", "br", "bc")}
    tiles = []
    for z in (13, 14):
        x, y = Transformer.from_crs(4326, 32600 + z, always_xy=True).transform(lon0, lat0)
        t = TL.tile_key(z, int((x - 15) // TL.TILE_M), int((y - 15) // TL.TILE_M))
        tiles.append(t)
        _, left, top = TL.tile_corner(t)
        bc0, br0 = int((x - left) // 480), int((top - y) // 480)
        for br in range(br0 - 4, br0 + 5):
            for bc in range(bc0 - 4, bc0 + 5):
                xc, yc = left + bc * 480 + 240, top - br * 480 - 240
                lo, la = Transformer.from_crs(32600 + z, 4326, always_xy=True).transform(xc, yc)
                for k, v in dict(zone=z, lon=lo, lat=la, valid=True, band_zone=int(TL.band_zone(lo)),
                                 tile_idx=len(tiles) - 1, br=br, bc=bc).items():
                    T[k].append(v)
    T = {k: np.array(v) for k, v in T.items()}
    prim = TL.resolve_primary(T, tiles, None)
    own = T["band_zone"] == T["zone"]
    assert prim[own].all()                         # a copy in its own band is always primary
    near = np.hypot((T["lon"] - lon0) * 85.0, (T["lat"] - lat0) * 111.0) < 1.0     # km; other zone covers it
    assert not prim[~own & near].any()             # ... and the other copy is not counted
    assert (~own & near).any()


def _lab_file(d, days, pr):
    import netCDF4 as nc
    p = os.path.join(d, "lab.nc")
    with nc.Dataset(p, "w") as ds:
        ds.createDimension("time", days.size)
        ds.createDimension("lat", 2)
        ds.createDimension("lon", 2)
        t = ds.createVariable("time", "f8", ("time",))
        t.units = "days since 1900-01-01"
        t[:] = (days - np.datetime64("1900-01-01")).astype(int)
        ds.createVariable("lat", "f8", ("lat",))[:] = [41.0, 40.958333]
        ds.createVariable("lon", "f8", ("lon",))[:] = [-101.0, -100.958333]
        for name, (s, o) in dict(pr=(0.1, 0), pet=(0.1, 0), vpd=(0.01, 0), srad=(0.1, 0), tmmn=(0.1, 210),
                                 tmmx=(0.1, 220)).items():
            v = ds.createVariable(name, "i2", ("time", "lat", "lon"), fill_value=32767)
            v.scale_factor, v.add_offset = s, o
            v.set_auto_maskandscale(False)
            val = pr if name == "pr" else np.full(pr.shape, {"pet": 3.0, "vpd": 1.0, "srad": 200.0,
                                                             "tmmn": 270.0, "tmmx": 290.0}[name])
            v[:] = np.round((val - o) / s).astype("i2")
            ds.createVariable(name + "_qc", "u1", ("time",))[:] = np.zeros(days.size, "u1")
        ds.processing_status = "complete"
    return p


def test_rain_windows_end_the_day_before_the_month():
    days = np.arange(np.datetime64("2004-01-01"), np.datetime64("2007-01-01"))
    rng = np.random.default_rng(2)
    pr = np.round(rng.gamma(0.5, 4, (days.size, 2, 2)), 1)
    with tempfile.TemporaryDirectory() as d:
        p = _lab_file(d, days, pr)
        with override(DRIVER_MONTH_FIRST="2005-01", ARIDITY_BASE_YEARS=(2004, 2006), BUCKET_SPINUP_FROM="2004-01-01"):
            R = DR.compute(p, np.array([0, 1]), np.array([0, 1]))
    k = list(R["months"]).index(np.datetime64("2006-03-01"))
    i0 = int((np.datetime64("2006-03-01") - days[0]).astype(int))
    for c, (iy, ix) in enumerate([(0, 0), (1, 1)]):
        assert abs(R["P30d"][k, c] - pr[i0 - 30:i0, iy, ix].sum()) < 1e-3
        assert abs(R["P"][k, c] - pr[i0:i0 + 31, iy, ix].sum()) < 1e-3
        assert abs(R["TMIN"][k, c] - (270 - 273.15)) < 1e-3


def test_bucket_stays_in_bounds_and_balances_mass():
    rng = np.random.default_rng(3)
    cap = np.array([[50.0], [200.0]])
    S = 0.5 * cap
    S0 = S.copy()
    P_tot = E_tot = 0
    aet_tot = np.zeros_like(S)
    drain = np.zeros_like(S)
    for _ in range(400):
        P, E = rng.gamma(0.4, 5), rng.uniform(0, 7)
        before = S + P
        S, aet = DR.bucket_step(S, np.full((1, 1), P), np.full((1, 1), E), cap)
        drain += np.maximum(before - cap, 0)
        aet_tot += aet
        P_tot += P
        assert (S >= -1e-9).all() and (S <= cap + 1e-9).all()
    assert np.allclose(S - S0, P_tot - aet_tot - drain)


def test_cdl_groups_cover_the_legend_once():
    assert CG.LUT[176] == CG.NAMES.index("grass_pasture")
    assert CG.LUT[24] == CG.NAMES.index("winter_wheat")
    assert CG.LUT[1] == CG.NAMES.index("corn")
    assert CG.LUT[111] == CG.NAMES.index("water")
    assert CG.group_of(np.array([255]))[0] == 255
    assert CG.LUT[0] == CG.OTHER


def test_regular_grid_lookup_takes_the_containing_cell():
    g = ST.RegularGrid(np.array([10.0, 9.0, 8.0]), np.array([0.0, 1.0, 2.0]))
    r, c, ok = g.index(np.array([9.6, 9.4, 7.6, 10.6]), np.array([0.4, 0.6, 2.4, 0.0]))
    assert list(r[ok]) == [0, 1, 2] and list(c[ok]) == [0, 1, 2] and not ok[3]


def test_reduce_blocks_statistics():
    a = np.arange(32 * 16, dtype="float64").reshape(32, 16)
    roi = np.ones(a.shape, bool)
    roi[0, 0] = False
    v = ST.reduce_blocks(a, roi, 1, "min")
    assert v[0, 0] == 1 and v[1, 0] == 256
    assert ST.reduce_blocks(a, roi, 1, "n")[0, 0] == 255


def test_hand_follows_d8_to_the_channel():
    try:
        import pysheds  # noqa: F401
    except ImportError:
        return "SKIP pysheds not installed (Colab installs it; locally use the venv named in docs/TESTS.md)"
    from ..step0 import hand as HD
    y, x = np.mgrid[0:60, 0:80]
    dem = 500 + 0.5 * np.abs(x - 40) + 0.05 * y
    channel = x == 40
    h, fd = HD.compute_hand(dem.astype("float64"), channel, 40.0)
    assert np.nanmax(np.abs(h[channel])) < 1e-6
    for (i, j) in [(10, 5), (30, 70), (50, 20)]:
        r, c = i, j
        for _ in range(500):
            if channel[r, c]:
                break
            dr, dc = HD.D8_STEP[int(fd[r, c])]
            r, c = r + dr, c + dc
        assert channel[r, c] and abs(h[i, j] - (dem[i, j] - dem[r, c])) < 1e-6


# ---------------------------------------------------------------- plan v3 (P1.0 and the new Step 0 layers)
def test_block_histograms_count_every_valid_pixel_in_the_right_bin():
    rng = np.random.default_rng(4)
    a = rng.uniform(-0.2, 0.7, (32, 48)).astype("float32")             # beyond both NIRv edges on purpose
    a[rng.random(a.shape) < 0.2] = np.nan
    h = BK.block_hist(a, 2, 3, "NIRv")
    _, _, n = BK.block_stats(a, 2, 3)
    assert h.shape == (C.HIST_BINS, 2, 3) and (h.sum(0) == n).all()
    e = BK.hist_edges("NIRv")
    b = a[16:32, 16:32]
    v = b[np.isfinite(b)]
    k = np.clip(np.searchsorted(e, v, side="right") - 1, 0, C.HIST_BINS - 1)
    assert (np.bincount(k, minlength=C.HIST_BINS) == h[:, 1, 1]).all()
    assert h[0, 1, 1] >= (v < e[0]).sum() and h[-1, 1, 1] >= (v >= e[-1]).sum()


def test_histograms_recover_the_share_of_a_greener_population():
    """P1.0 check in miniature: a block mixing rain-only and subsidised pixels; the share above the rain-only
    95 % value, read from the histogram, minus the 5 % expected by chance, recovers the subsidised share."""
    rng = np.random.default_rng(5)
    shares = [0.0, 0.25, 0.5, 0.75]
    a = np.full((16, 16 * len(shares)), np.nan, "float32")
    for j, sh in enumerate(shares):
        x = rng.normal(0.18, 0.03, 256)
        k = int(round(sh * 256))
        x[:k] = rng.normal(0.40, 0.03, k)
        a[:, j * 16:(j + 1) * 16] = x.reshape(16, 16)
    h = BK.block_hist(a, 1, len(shares), "NIRv")[:, 0, :]
    e = BK.hist_edges("NIRv")
    thr = 0.18 + 1.645 * 0.03                                         # rain-only 95 % quantile
    kb = np.searchsorted(e, thr, side="right") - 1
    frac_in_bin = (e[kb + 1] - thr) / (e[kb + 1] - e[kb])             # uniform within the threshold's bin
    above = (h[kb + 1:].sum(0) + frac_in_bin * h[kb]) / h.sum(0)
    est = (above - 0.05) / 0.95
    assert np.abs(est - np.array(shares)).max() < 0.06, est


def test_observation_support_is_a_pixel_weighted_ratio():
    num = np.zeros((16, 32), "float32")
    den = np.full((16, 32), 2.0, "float32")
    num[:, :16] = 1.0
    den[0, 0] = np.nan
    r = BK.block_ratio(num, den, 1, 2)
    assert abs(r[0, 0] - 0.5) < 1e-6 and r[0, 1] == 0.0
    assert np.isnan(BK.block_ratio(num, np.full((16, 32), np.nan, "float32"), 1, 2)).all()


def test_d8_accumulation_counts_upstream_cells():
    from ..step0 import hand as HD
    # a 3 x 3 grid draining to the bottom-centre cell, which drains out of the grid
    E, S, W, SE, SW = 1, 4, 16, 2, 8
    f = np.array([[SE, S, SW],
                  [E, S, W],
                  [E, S, W]], "uint8")
    acc = HD.d8_accumulation(f, np.ones((3, 3), bool))
    assert acc[2, 1] == 9 and acc[1, 1] == 6 and acc[0, 1] == 1 and acc[2, 0] == 1
    acc2 = HD.d8_accumulation(f, np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], bool))
    assert acc2[2, 1] == 8                                             # an invalid cell passes flow but adds none
    assert (HD.d8_accumulation(np.zeros((2, 2), "uint8"), np.ones((2, 2), bool)) == 1).all()


def test_spi_month_value_is_the_last_record_of_the_month_and_classes():
    import netCDF4 as nc
    d = tempfile.mkdtemp()
    p = os.path.join(d, "drought.nc")
    days = np.arange(np.datetime64("2012-01-01"), np.datetime64("2012-04-01"), 5)
    with nc.Dataset(p, "w") as ds:
        ds.createDimension("time", days.size)
        ds.createDimension("lat", 2)
        ds.createDimension("lon", 3)
        t = ds.createVariable("time", "f8", ("time",))
        t.units = "days since 1900-01-01"
        t[:] = (days - np.datetime64("1900-01-01")).astype(int)
        ds.createVariable("lat", "f8", ("lat",))[:] = [40.0, 40.0 - 1 / 24]
        ds.createVariable("lon", "f8", ("lon",))[:] = [-100.0, -100.0 + 1 / 24, -100.0 + 2 / 24]
        v = ds.createVariable("spi90d", "i2", ("time", "lat", "lon"), fill_value=32767)
        v.set_auto_maskandscale(False)
        v.scale_factor = np.float32(0.01)
        v[:] = (np.arange(days.size)[:, None, None] * 10 + np.zeros((1, 2, 3))).astype("i2")
        q = ds.createVariable("spi90d_qc", "u1", ("time",))
        qq = np.zeros(days.size, "u1")
        qq[-1] = C.QC_NODATA                                            # last record of March is flagged
        q[:] = qq
        ds.processing_status = "complete"
    months = np.arange(np.datetime64("2011-12"), np.datetime64("2012-04"))
    out = DR.spi_monthly(p, np.array([40.0 - 1 / 24]), np.array([-100.0 + 2 / 24]), months)
    last_jan = np.nonzero(days.astype("datetime64[M]") == np.datetime64("2012-01"))[0][-1]
    assert np.isnan(out[0, 0])                                          # December 2011: no record
    assert abs(out[1, 0] - last_jan * 0.1) < 1e-5
    assert np.isnan(out[3, 0])                                          # March: its last record is flagged
    try:
        DR.spi_monthly(p, np.array([41.0]), np.array([-100.0]), months)
        raise AssertionError("a cell off the drought grid was accepted")
    except Exception as e:
        assert "grid" in str(e)
    assert list(DR.drought_class([-1.0, -0.8, 0.0, 0.8, np.nan])) == [0, 0, 1, 2, 1]
