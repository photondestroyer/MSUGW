"""Blinded ET extraction: cell-month OpenET values for cohort and rainfed cells, through the blinding gate.

Physical idea: ET is the water the field actually used. Before the freeze only two kinds of it may be seen:
rainfed cells (the baseline) and cohort cells in non-dry months (to estimate noise). blinding.gate() enforces
that on every tile-year before anything is written, so dry cohort ET exists only in memory.
Decisions: D20 (dry), D21 (model coverage), D29 (ET value).
Deliverable: et_cell_month.parquet (blinded), parts/<tile>_<year>.parquet (resumable units, already gated).

iter_tile_months() is shared with Task 1.4 (coverage): the same files, frames and cell-id rasters.
"""
import os

import numpy as np

FILL = -32768
NO_COUNT = 255


def scoped_cells(ctx):
    """Cells table restricted to the run's scope (smoke county / stage2 sample / all), and the full table."""
    from ..common import cells as CL, report
    from ..phase0 import t0_2_cells as T02
    cells = T02.load(ctx.dep("t0_2"))
    sample = report.read(ctx.dep("sample"), "sample_cells").cell_id if "sample" in ctx.stage["deps"] else None
    return cells, CL.scope(cells, sample)


def iter_tile_months(ctx, cells, in_scope, years, months):
    """Yield (tile, year, month, ids, ds, (t, r0, r1, c0, c1)) for every OpenET tile-month holding a cell in scope.
    ids: cell-id raster of the WINDOW (the bounding box of in-scope cells in the tile), out-of-scope cells -1.
    ds is open; the caller reads ds[var][t, r0:r1, c0:c1] (read_at)."""
    import geopandas as gpd
    from .. import config as C
    from ..common import io, cells as CL
    from ..phase0 import t0_2_cells as T02
    geo = T02.load_geometry(ctx.dep("t0_2"))
    geo = geo.set_index("cell_id").loc[cells.cell_id].reset_index()
    geo["x5070"] = cells["x5070"].to_numpy()
    frames = io.tile_frames("openet")
    sig = CL.cells_signature(cells)
    cache = f"{C.shared_root()}/cellids"
    keep = np.where(in_scope)[0]
    for tile in C.KS_TILES:
        F = frames.get(tile)
        if F is None:
            continue
        ids = CL.tile_ids(geo, F, cache, sig)
        m = np.isin(ids, keep)
        if not m.any():
            continue
        rr, cc = np.where(m)
        r0, r1, c0, c1 = rr.min(), rr.max() + 1, cc.min(), cc.max() + 1
        ids = np.where(m, ids, -1)[r0:r1, c0:c1]
        for y in years:
            p = io.openet_file(tile, y)
            if not os.path.exists(p):
                ctx.log.warning("OpenET %s %d missing: %s", tile, y, p)
                continue
            ds = io.open_tile(p)
            io.tile_frame_of(ds, frames, tile)
            yy, mm = io.openet_months(ds)
            for mo in months:
                k = np.where((yy == y) & (mm == mo))[0]
                if len(k):
                    yield tile, y, mo, ids, ds, (int(k[0]), r0, r1, c0, c1)
            ds.close()


def read_at(ds, var, at, tries=4):
    """ds[var][t, r0:r1, c0:c1]. Reads through the Google Drive mount sometimes fail with a transient HDF error;
    they are retried with a pause, and a read that keeps failing stops the task (it resumes from its parts)."""
    import time
    from ..core.errors import SourceError
    t, r0, r1, c0, c1 = at
    for k in range(tries):
        try:
            return np.asarray(ds[var][t, r0:r1, c0:c1])
        except RuntimeError as e:
            last = e
            time.sleep(3 * (k + 1))
    raise SourceError(f"{ds.filepath()}: {var}[{t}] unreadable after {tries} tries ({last})")


def valid_pixels(ds, model, at, et=None):
    """Pixels where a model holds a value backed by at least one cloud-free observation (count > 0).
    et: the model's raw ET window when the caller has already read it (it is not read a second time)."""
    if et is None:
        et = read_at(ds, f"et_{model}", at)
    cnt = read_at(ds, f"count_{model}", at)
    return (et != FILL) & (cnt != NO_COUNT) & (cnt > 0)


def save_part(ctx, df, parts, key):
    """One resumable unit (a tile-year) -> <parts>/<key>.parquet. Written in scratch, then delivered: on Colab the
    scratch disk and Google Drive are different file systems, where a plain rename fails."""
    from ..core.drive import deliver
    tmp = os.path.join(ctx.scratch, f"{key}.parquet")
    df.to_parquet(tmp, index=False)
    deliver(tmp, f"{parts}/{key}.parquet")


def run(ctx):
    import pandas as pd
    from .. import config as C
    from .. import blinding
    from ..common import report, zonal

    # ══════ STEP E-1  cohort and rainfed cells in scope; dry flags; per-model coverage (Task 1.4) ══════
    cells, scope = scoped_cells(ctx)
    scope = scope & cells.role.isin(["cohort", "rainfed"]).to_numpy()
    dry = report.read(ctx.dep("dry"), "dry_flags")
    dry["county"] = dry["county"].astype(str)
    cov = report.read(ctx.dep("t1_4"), "et_coverage")
    ok_cols = [f"ok_{m}" for m in C.ET_MODELS]
    cov = cov.set_index(["cell", "year", "month"])[ok_cols + ["supported"]]
    n = len(cells)
    parts = f"{ctx.out}/parts"
    os.makedirs(parts, exist_ok=True)
    ctx.log.info("STEP E-1  %d cohort/rainfed cells in scope", int(scope.sum()))

    # ══════ STEP E-2  per tile-year: cell means of the ensemble and each model, gated before writing ══════
    totals = dict(rows_in=0, rows_kept=0, cohort_dry_dropped=0, other_dropped=0)
    done = set(os.listdir(parts))
    buf, cur = [], None

    def flush(key):
        if not buf:
            return
        et = pd.concat(buf, ignore_index=True)
        out, dropped = blinding.gate(et, cells[["cell_id", "county", "role"]], dry)
        for k in totals:
            totals[k] += dropped[k]
        save_part(ctx, out, parts, key)
        buf.clear()

    for tile, y, mo, ids, ds, t in iter_tile_months(ctx, cells, scope, C.PANEL_YEARS, C.SEASON_MONTHS):
        key = f"{tile}_{y}"
        if f"{key}.parquet" in done:
            continue
        if cur is not None and key != cur:
            flush(cur)
        cur = key
        present = np.unique(ids[ids >= 0])
        rec = dict(cell=cells.cell_id.to_numpy()[present], year=y, month=mo)
        # sums and counts, not means: a cell that crosses a tile edge is added up from its pieces in STEP E-3
        ens = read_at(ds, C.ET_ENSEMBLE, t).astype("float64")
        ens[ens == FILL] = np.nan
        s, c = zonal.total(ids, ens, n)
        rec[f"{C.ET_ENSEMBLE}_sum"], rec[f"{C.ET_ENSEMBLE}_n"] = s[present], c[present]
        for m in C.ET_MODELS:
            raw = read_at(ds, f"et_{m}", t)
            v = raw.astype("float64")
            v[~valid_pixels(ds, m, t, et=raw)] = np.nan
            s, c = zonal.total(ids, v, n)
            rec[f"et_{m}_sum"], rec[f"et_{m}_n"] = s[present], c[present]
        buf.append(pd.DataFrame(rec))
        ctx.check_budget()
    flush(cur)

    # ══════ STEP E-3  add up tile pieces per cell-month, apply D21 (D29), gate again, one blinded table ══════
    files = sorted(os.listdir(parts))
    allp = pd.concat([pd.read_parquet(os.path.join(parts, f)) for f in files], ignore_index=True)
    allp = allp.drop(columns=[c for c in ("county", "role", "dry", "blinding_gate") if c in allp.columns])
    allp = allp.groupby(["cell", "year", "month"], as_index=False).sum()
    vars_ = [C.ET_ENSEMBLE] + [f"et_{m}" for m in C.ET_MODELS]
    for v in vars_:
        with np.errstate(invalid="ignore", divide="ignore"):
            allp[v] = np.where(allp[f"{v}_n"] > 0, allp[f"{v}_sum"] / allp[f"{v}_n"], np.nan)
    allp = allp.drop(columns=[f"{v}_{s}" for v in vars_ for s in ("sum", "n")])
    k = cov.reindex(pd.MultiIndex.from_arrays([allp.cell, allp.year, allp.month]))
    for m in C.ET_MODELS:                                         # D29: a model's value only where D21 holds
        allp.loc[~k[f"ok_{m}"].astype("boolean").fillna(False).to_numpy(), f"et_{m}"] = np.nan
    supported = k["supported"].astype("boolean").fillna(False).to_numpy()
    allp.loc[~supported, C.ET_ENSEMBLE] = np.nan
    allp["et_supported"] = supported
    out, dropped = blinding.write_et(ctx, allp, cells[["cell_id", "county", "role"]], dry, "et_cell_month",
                                     units={"et_*": "mm per month"},
                                     notes="cell mean of valid pixels where D21 holds (D29)")
    ok, problems, nchk = blinding.audit(ctx.out)
    if not ok:
        raise blinding.BlindingError("audit failed on this task's outputs: " + "; ".join(problems[:5]))
    ctx.complete(dict(rows=int(len(out)), parts=len(files), gate_last=dropped, gate_parts=totals,
                      cohort_rows=int((out.role == "cohort").sum()), rainfed_rows=int((out.role == "rainfed").sum()),
                      audit_tables=nchk))
