"""Draft paper items the pilot can already produce (pilot Sections 9-10), regenerated from stored tables only.

Paper 1: Table 1 (data sources, versions, periods), Table 2 (sample per mask), Figure 2 (study area and data
support), Figure 1 rainfed curve (copied from Check 2.5).
Paper 2 (Kansas): Table 1 (well-years by source and year), Table 2 (validation statistics), Figure 1 (well
network), Figure 2 (thickness and SE, from Check 2.4), Figure 3 (validation, from Check 2.4).
Colours and labels are config.COLORS everywhere (thin, thick, rainfed, cohort).
"""
import os
import shutil


def run(ctx):
    import numpy as np
    import pandas as pd
    from .. import config as C
    from ..common import report
    from ..phase0 import counties as CTY, t0_2_cells as T02
    # ══════ STEP P-1  Paper 1 Table 1: data sources, versions, periods (Task 0.1 manifest) ══════
    man = report.read(ctx.dep("t0_1"), "manifest")
    t1 = man[man.exists.astype(bool)][[c for c in ("key", "title", "period", "status", "path") if c in man.columns]]
    report.table(ctx, t1, "P1_Table1_data_sources")
    # ══════ STEP P-2  Paper 1 Table 2: sample in each mask (Task 1.7) ══════
    report.table(ctx, report.read(ctx.dep("t1_7"), "mask_summary"), "P1_Table2_sample",
                 units=dict(area_ha="ha"), mask="the four sample masks")
    # ══════ STEP P-3  Paper 1 Figure 2: cohort cells, thickness, wells, footprints ══════
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cty = CTY.load(ctx.dep("counties"))
    geo = T02.load_geometry(ctx.dep("t0_2"))
    th = report.read(ctx.dep("c2_4"), "thickness_cell_year")
    am = th[th.year.isin(C.PANEL_YEARS)].groupby("cell_id").A_m.mean()
    wells = report.read(ctx.dep("t1_5"), "well_year")
    fp = report.read(ctx.dep("t1_2"), "footprints")
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    cty.boundary.plot(ax=ax[0], color="#777777", linewidth=0.4)
    g = geo[geo.role == "cohort"].assign(A=lambda d: d.cell_id.map(am))
    g.plot(ax=ax[0], column="A", cmap="viridis", linewidth=0, legend=True, legend_kwds=dict(label="mean A 2016-2024 (m)"))
    w = wells[wells.use.astype(bool) & (wells.year == 2020)]
    ax[0].scatter(w.x5070, w.y5070, s=1, color="k")
    ax[0].set_title("Cohort cells by aquifer thickness; wells measured in the 2020 pre-season", fontsize=9)
    cty.boundary.plot(ax=ax[1], color="#777777", linewidth=0.4)
    ref = fp[fp.is_reference_any.astype(bool)]
    ax[1].scatter(ref.x5070, ref.y5070, s=60, marker="s", facecolor="none", edgecolor=C.COLORS["footprint"])
    geo[geo.role == "rainfed"].plot(ax=ax[1], color=C.COLORS["rainfed"], linewidth=0)
    geo[geo.role == "cohort"].plot(ax=ax[1], color=C.COLORS["cohort"], linewidth=0)
    ax[1].set_title("Reference footprints (squares), cohort (purple), rainfed (green)", fontsize=9)
    for a in ax:
        a.set_axis_off()
    report.figure(ctx, fig, "P1_Figure2_study_area", units={"A": "m"}, mask="study cells",
                  caption="Paper 1 Figure 2 draft")
    # ══════ STEP P-4  Paper 2 Tables 1-2 and Figure 1; copies of the drafts made by the checks ══════
    report.table(ctx, report.read(ctx.dep("t1_5"), "wells_by_source_year"), "P2_Table1_well_years",
                 units=dict(well_years="count"))
    report.table(ctx, report.read(ctx.dep("c2_4"), "validation_stats"), "P2_Table2_validation", units=dict(rmse_m="m"))
    copies = [("t1_5", "wells_support.png", "P2_Figure1_wells.png"),
              ("c2_4", "thickness_map.png", "P2_Figure2_thickness.png"),
              ("c2_4", "validation.png", "P2_Figure3_validation.png"),
              ("c2_5", "rainfed_curve.png", "P1_Figure1_rainfed_curve.png")]
    for task, name, target in copies:
        src = f"{ctx.dep(task)}/{name}"
        if os.path.exists(src):
            shutil.copy2(src, f"{ctx.out}/{target}")
    ctx.complete(dict(items=sorted(f for f in os.listdir(ctx.out) if f.startswith(("P1_", "P2_")))))
