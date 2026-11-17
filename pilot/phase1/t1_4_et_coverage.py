"""Task 1.4  ET coverage.

"For every month of May-September and each of the six OpenET models, compute the share of each cell covered by
valid ET pixels with a positive observation count. Fix the share a model needs to count for a cell-month. Export
these coverage summaries before any ET value. OpenET documents a coverage gap near the Colorado-Kansas border."
Physical idea: a monthly ET total is only as good as the cloud-free Landsat views behind it; a month filled
entirely from neighbouring dates says little about that month's water use. Agreement of at least four models
guards against one model's quirks being read as a field's response.
Decisions: D21.   Deliverables: et_coverage (parquet; shares, no ET values), coverage_summary.csv,
coverage_by_month.png.
"""
import os

import numpy as np


def coverage_part(job):
    """Worker: per cell, month and model of one tile-year, the share of the cell's pixels valid with count > 0."""
    import pandas as pd
    from .. import config as C
    from ..common import zonal
    from .et_cells import open_job, valid_pixels, read_at, FILL
    ids, ds, ts = open_job(job)
    n = job["n"]
    present = np.unique(ids[ids >= 0])
    npx = zonal.count(ids, n)[present]
    expected = (job["area_ha"][present] * 1e4 / C.PIXEL_M ** 2).round().astype(int)
    out = []
    for mo, k in ts:
        t = (k, *job["window"])
        rec = dict(cell=job["cell_id"][present], year=job["year"], month=mo, n_px=npx, expected_px=expected)
        rec["share_ensemble"] = zonal.count(ids, n, read_at(ds, C.ET_ENSEMBLE, t) != FILL)[present] / expected
        for m in C.ET_MODELS:
            rec[f"share_{m}"] = zonal.count(ids, n, valid_pixels(ds, m, t))[present] / expected
        out.append(pd.DataFrame(rec))
    ds.close()
    return pd.concat(out, ignore_index=True) if out else None


def run(ctx):
    import pandas as pd
    from .. import config as C
    from .. import definitions as DEF
    from ..common import report
    from .et_cells import scoped_cells, tile_jobs, run_jobs, save_part

    # ══════ STEP 1.4-1  cells in scope ══════
    cells, scope = scoped_cells(ctx)
    n = len(cells)
    parts = f"{ctx.out}/parts"
    os.makedirs(parts, exist_ok=True)
    ctx.log.info("STEP 1.4-1  %d cells in scope", int(scope.sum()))

    # ══════ STEP 1.4-2  per tile-month and model: share of the cell's pixels valid with count > 0 ══════
    jobs = tile_jobs(ctx, cells, scope, C.PANEL_YEARS, C.SEASON_MONTHS, done=set(os.listdir(parts)), n=n,
                     cell_id=cells.cell_id.to_numpy(), area_ha=cells.area_ha.to_numpy())
    for key, part in run_jobs(ctx, coverage_part, jobs):
        if part is not None:
            save_part(ctx, part, parts, key)
    files = sorted(os.listdir(parts))
    cov = pd.concat([pd.read_parquet(os.path.join(parts, f)) for f in files], ignore_index=True)
    # a cell split over two tiles: add its pieces
    agg = {c: "sum" for c in ["n_px"]}
    share_cols = ["share_ensemble"] + [f"share_{m}" for m in C.ET_MODELS]
    for c in share_cols:
        cov[c] = cov[c] * cov["expected_px"]
    cov = cov.groupby(["cell", "year", "month"]).agg({**agg, "expected_px": "first", **{c: "sum" for c in share_cols}})
    cov = cov.reset_index()
    for c in share_cols:
        cov[c] = (cov[c] / cov["expected_px"]).clip(upper=1.0)

    # ══════ STEP 1.4-3  which models count (D21) and supported cell-months; cells with no OpenET pixel ══════
    for m in C.ET_MODELS:
        cov[f"ok_{m}"] = DEF.et_model_ok(cov[f"share_{m}"])
    cov["n_models"] = cov[[f"ok_{m}" for m in C.ET_MODELS]].sum(axis=1)
    cov["supported"] = DEF.et_supported(cov["n_models"])
    sc = cells.loc[scope, ["cell_id", "county", "role"]]
    seen = set(cov.cell)
    no_px = sc[~sc.cell_id.isin(seen)]
    cov = cov.merge(cells[["cell_id", "role", "county"]], left_on="cell", right_on="cell_id", how="left").drop(
        columns="cell_id")
    summ = (cov[cov.role == "cohort"].groupby(["year", "month"])
            .agg(cohort_cells=("cell", "size"), share_supported=("supported", "mean"),
                 **{f"share_ok_{m}": (f"ok_{m}", "mean") for m in C.ET_MODELS}).reset_index())
    ctx.log.info("STEP 1.4-3  cohort cell-months supported: %.1f%%; %d cells without any OpenET pixel",
                 100 * cov[cov.role == "cohort"].supported.mean(), len(no_px))

    # ══════ STEP 1.4-4  write (coverage only: no ET value leaves this task) ══════
    report.table(ctx, cov, "et_coverage", parquet=True, csv=False, units={"share_*": "fraction of cell pixels"},
                 notes="valid = value present and observation count > 0 (count 0 = filled from other dates)")
    report.table(ctx, summ, "coverage_summary", units={"share_*": "fraction of cohort cells"}, mask="cohort cells")
    report.table(ctx, no_px, "cells_without_openet", notes="e.g. the strip west of about -102.027 deg (zone-14 edge)")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 3.5))
    for m in C.ET_MODELS:
        s = summ.groupby("month")[f"share_ok_{m}"].mean()
        ax.plot(s.index, s.values, marker="o", label=m)
    s = summ.groupby("month").share_supported.mean()
    ax.plot(s.index, s.values, color="k", lw=2.5, label=">= 4 models")
    ax.axhline(C.TARGETS["et_cell_share"], ls="--", color="grey")
    ax.set_xlabel("month")
    ax.set_ylabel("share of cohort cell-months")
    ax.legend(fontsize=7, ncol=2)
    report.figure(ctx, fig, "coverage_by_month", units={"share": "fraction"}, mask="cohort cells, 2016-2024")
    ctx.complete(dict(cell_months=int(len(cov)), cohort_supported=float(cov[cov.role == "cohort"].supported.mean())
                      if (cov.role == "cohort").any() else None, cells_without_openet=int(len(no_px)),
                      parts=len(files)))
