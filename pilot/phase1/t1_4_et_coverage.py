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


def run(ctx):
    import pandas as pd
    from .. import config as C
    from .. import definitions as DEF
    from ..common import report, zonal
    from .et_cells import scoped_cells, iter_tile_months, valid_pixels, read_at, save_part, FILL

    # ══════ STEP 1.4-1  cells in scope ══════
    cells, scope = scoped_cells(ctx)
    n = len(cells)
    parts = f"{ctx.out}/parts"
    os.makedirs(parts, exist_ok=True)
    done = set(os.listdir(parts))
    ctx.log.info("STEP 1.4-1  %d cells in scope", int(scope.sum()))

    # ══════ STEP 1.4-2  per tile-month and model: share of the cell's pixels valid with count > 0 ══════
    buf, cur = [], None

    def flush(key):
        if buf:
            save_part(ctx, pd.concat(buf, ignore_index=True), parts, key)
            buf.clear()

    for tile, y, mo, ids, ds, t in iter_tile_months(ctx, cells, scope, C.PANEL_YEARS, C.SEASON_MONTHS):
        key = f"{tile}_{y}"
        if f"{key}.parquet" in done:
            continue
        if cur is not None and key != cur:
            flush(cur)
        cur = key
        present = np.unique(ids[ids >= 0])
        npx = zonal.count(ids, n)[present]
        rec = dict(cell=cells.cell_id.to_numpy()[present], year=y, month=mo, n_px=npx,
                   expected_px=(cells.area_ha.to_numpy()[present] * 1e4 / C.PIXEL_M ** 2).round().astype(int))
        ens_valid = read_at(ds, C.ET_ENSEMBLE, t) != FILL
        rec["share_ensemble"] = zonal.count(ids, n, ens_valid)[present] / rec["expected_px"]
        for m in C.ET_MODELS:
            rec[f"share_{m}"] = zonal.count(ids, n, valid_pixels(ds, m, t))[present] / rec["expected_px"]
        buf.append(pd.DataFrame(rec))
        ctx.check_budget()
    flush(cur)
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
