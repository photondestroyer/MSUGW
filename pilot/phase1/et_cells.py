"""Blinded ET extraction: cell-month OpenET values for cohort and rainfed cells, through the blinding gate.

Physical idea: ET is the water the field actually used. Before the freeze only two kinds of it may be seen:
rainfed cells (the baseline) and cohort cells in non-dry months (to estimate noise). blinding.gate() enforces
that on every tile-year before anything is written, so dry cohort ET exists only in memory.
Decisions: D20 (dry), D21 (model coverage), D29 (ET value).
Deliverable: et_cell_month.parquet (blinded), parts/<tile>_<year>.parquet (resumable units, already gated).

tile_jobs(), open_job() and run_jobs() are shared with Task 1.4 (coverage): the same files, frames and cell-id
rasters. Tile-years are read by worker processes (config.N_WORKERS); dry cohort ET passes from a worker to
run() in memory only and is dropped by the gate there, so it still never exists on disk.
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


def tile_jobs(ctx, cells, in_scope, years, months, done=(), **extra):
    """One job (a dict a worker process can take) per OpenET tile-year that holds a cell in scope and has no part
    file yet. The cell-id raster of the tile is built here once and read by the workers from its cache file."""
    from .. import config as C
    from ..common import io, cells as CL, grids
    from ..phase0 import t0_2_cells as T02
    geo = T02.load_geometry(ctx.dep("t0_2"))
    geo = geo.set_index("cell_id").loc[cells.cell_id].reset_index()
    geo["x5070"] = cells["x5070"].to_numpy()
    frames = io.tile_frames("openet")
    sig = CL.cells_signature(cells)
    cache = f"{C.shared_root()}/cellids"
    keep = np.where(in_scope)[0]
    jobs = []
    for tile in C.KS_TILES:
        F = frames.get(tile)
        if F is None:
            continue
        m = np.isin(CL.tile_ids(geo, F, cache, sig), keep)
        if not m.any():
            continue
        rr, cc = np.where(m)
        window = (int(rr.min()), int(rr.max()) + 1, int(cc.min()), int(cc.max()) + 1)
        for y in years:
            p = io.openet_file(tile, y)
            if f"{tile}_{y}.parquet" in done:
                continue
            if not os.path.exists(p):
                ctx.log.warning("OpenET %s %d missing: %s", tile, y, p)
                continue
            jobs.append(dict(tile=tile, year=y, path=p, frame=F, window=window, keep=keep, months=list(months),
                             ids_path=os.path.join(cache, f"ids_{sig}_{grids.frame_key(F)}.npz"), **extra))
    return jobs


def open_job(job):
    """In a worker: (ids, ds, [(month, time index)]) of one tile-year. ids is the cell-id raster of the WINDOW
    (the bounding box of in-scope cells in the tile), out-of-scope cells -1; the caller reads
    ds[var][t, r0:r1, c0:c1] (read_at) and closes ds."""
    from ..common import io
    r0, r1, c0, c1 = job["window"]
    ids = np.load(job["ids_path"])["ids"][r0:r1, c0:c1]
    ids = np.where(np.isin(ids, job["keep"]), ids, -1)
    ds = io.open_tile(job["path"])
    io.tile_frame_of(ds, {job["tile"]: job["frame"]}, job["tile"])
    yy, mm = io.openet_months(ds)
    ts = [(mo, int(k[0])) for mo in job["months"] for k in [np.where((yy == job["year"]) & (mm == mo))[0]] if len(k)]
    return ids, ds, ts


def run_jobs(ctx, worker, jobs):
    """Yield (key, table) per tile-year. Reading the tiles is the slow step, so config.N_WORKERS processes read
    them side by side; the caller alone writes (for ET, after the blinding gate)."""
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from .. import config as C
    if C.N_WORKERS <= 1 or len(jobs) <= 1:
        for j in jobs:
            yield f"{j['tile']}_{j['year']}", worker(j)
            ctx.check_budget()
        return
    with ProcessPoolExecutor(min(C.N_WORKERS, len(jobs))) as ex:
        futs = {ex.submit(worker, j): f"{j['tile']}_{j['year']}" for j in jobs}
        for f in as_completed(futs):
            yield futs[f], f.result()
            ctx.check_budget()


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


def et_part(job):
    """Worker: per cell and month of one tile-year, sums and pixel counts of the ensemble and of each model.
    Not gated: the rows go back to run() in memory, which gates them before anything is written."""
    import pandas as pd
    from .. import config as C
    from ..common import zonal
    ids, ds, ts = open_job(job)
    n = job["n"]
    present = np.unique(ids[ids >= 0])
    out = []
    for mo, k in ts:
        t = (k, *job["window"])
        rec = dict(cell=job["cell_id"][present], year=job["year"], month=mo)
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
        out.append(pd.DataFrame(rec))
    ds.close()
    return pd.concat(out, ignore_index=True) if out else None


def run(ctx):
    import pandas as pd
    from .. import config as C
    from .. import blinding
    from ..common import report

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
    unblind = C.RUN_MODE == "main"          # after the freeze only: gate() itself refuses without a valid record
    jobs =tile_jobs(ctx, cells, scope, C.PANEL_YEARS, C.SEASON_MONTHS, done=set(os.listdir(parts)), n=n,
                     cell_id=cells.cell_id.to_numpy())
    for key, et in run_jobs(ctx, et_part, jobs):
        if et is None:
            continue
        out, dropped = blinding.gate(et, cells[["cell_id", "county", "role"]], dry, unblind)
        for k in totals:
            totals[k] += dropped[k]
        save_part(ctx, out, parts, key)

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
                                     notes="cell mean of valid pixels where D21 holds (D29)", unblind=unblind)
    ok, problems, nchk = blinding.audit(ctx.out)
    if not ok:
        raise blinding.BlindingError("audit failed on this task's outputs: " + "; ".join(problems[:5]))
    ctx.complete(dict(rows=int(len(out)), parts=len(files), gate_last=dropped, gate_parts=totals,
                      cohort_rows=int((out.role == "cohort").sum()), rainfed_rows=int((out.role == "rainfed").sum()),
                      audit_tables=nchk))
