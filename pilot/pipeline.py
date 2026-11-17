"""The pilot's task graph (GW_ET_High_Plains_pilot.docx Sections 5-7, 9, 10).

Each task: name, title (pilot ID), deps (tasks whose outputs it reads), keys (config settings + decision IDs
that change its result: part of its signature), decisions (must be FIXED before real data are read),
modules (code it runs: part of its signature), func ("module:function"), sources (config.DATASETS keys it reads,
plus "smap" for the SMAP file in use: the resolved file is part of its signature, so a task runs again when the
file behind a key changes).
The runner (core/stages.py) turns this into signatures, output folders, resume and the decision check.
"""
DEF = ["definitions.py"]
IO = ["common/io.py", "common/report.py", "common/grids.py"]
CELLS = ["common/cells.py", "common/zonal.py"]


def T(name, title, func, deps=(), keys=(), decisions=(), modules=(), sources=()):
    return dict(name=name, title=title, func=func, deps=list(deps),
                keys=list(keys) + list(decisions) + [f"SRC:{s}" for s in sources],
                decisions=list(decisions), modules=list(modules), sources=list(sources))


STAGES = [
    # ------------------------------------------------------------------ Phase 0: foundations
    T("t0_1", "0.1 Manifest", "phase0.t0_1_manifest:run", keys=["DATASETS", "EARLIER_OUTPUTS"],
      modules=["phase0/t0_1_manifest.py"] + IO),
    T("t0_4", "0.4 Blinding audit", "phase0.t0_4_blinding_audit:run", keys=["EARLIER_OUTPUTS"],
      modules=["phase0/t0_4_blinding_audit.py", "common/report.py"]),
    T("counties", "Study counties (D02)", "phase0.counties:run", keys=["STATE_FIPS", "REP_REGIONS"],
      decisions=["D02"], modules=["phase0/counties.py"] + IO, sources=["counties_zip", "gmd", "aquifer_shp"]),
    T("t0_2", "0.2 Cells and masks", "phase0.t0_2_cells:run", deps=["counties"],
      keys=["LABEL_YEARS", "CORN_YEARS", "LATER_IRR_YEARS", "CDL_CLASSES", "CROPLAND_CODES", "KS_TILES",
            "SMOKE_COUNTY"],
      decisions=["D03", "D04", "D05", "D06", "D07", "D08", "D09"],
      modules=["phase0/t0_2_cells.py"] + DEF + IO + CELLS,
      sources=["plss_sections", "gmd", "aquifer_shp", "gssurgo", "wimas_pods", "kgs_sections", "usgs_bottom", "dem",
               "kgs_extent", "aim_dir", "lanid", "cdl_tiles"]),
    T("t0_3", "0.3 Benchmark reproduction", "phase0.t0_3_benchmark:run",
      modules=["phase0/t0_3_benchmark.py", "phase0/run_mieno.R"]),
    # ------------------------------------------------------------------ Phase 1: data support
    T("t1_1", "1.1 Weather", "phase1.t1_1_weather:run", deps=["counties"],
      keys=["PANEL_YEARS", "SEASON_MONTHS", "PEAK_CANDIDATES"], decisions=["D01", "D10", "D20"],
      modules=["phase1/t1_1_weather.py"] + DEF + IO, sources=["gridmet"]),
    T("t1_2", "1.2 Reference footprints", "phase1.t1_2_footprints:run", deps=["counties"],
      keys=["LABEL_YEARS", "FOOTPRINT_CHECK_YEARS", "SMAP_ORDER", "EASE_M36", "EASE_M09"],
      decisions=["D11", "D12", "D13", "D14", "D15"], modules=["phase1/t1_2_footprints.py"] + DEF + IO + CELLS,
      sources=["smap", "aim_dir", "lanid_1820_dir", "gridmet", "gssurgo"]),
    T("t1_3", "1.3 SMAP dates", "phase1.t1_3_smap_dates:run", deps=["counties", "t1_2"],
      keys=["PANEL_YEARS", "SEASON_MONTHS", "SMAP_ORDER", "SMAP_GAPS_EXPECTED"],
      decisions=["D11", "D15", "D16", "D17", "D18", "D19"], modules=["phase1/t1_3_smap_dates.py"] + DEF + IO,
      sources=["smap"]),
    T("dry", "Dry flags (blinding)", "blinding:run", deps=["t1_1", "t1_3"], keys=["PANEL_YEARS", "SEASON_MONTHS"],
      decisions=["D20"], modules=["blinding.py"] + DEF),
    T("t1_5", "1.5 Well inventory", "phase1.t1_5_wells:run", deps=["counties", "t0_2"],
      keys=["WELL_YEARS", "PRESEASON_MONTHS", "HPA_CODE_PREFIXES"], decisions=["D22"],
      modules=["phase1/t1_5_wells.py"] + DEF + IO + ["common/stats.py"],
      sources=["kgs_winter", "wells_ts", "kgs_extent", "dem"]),
    T("sample", "1.6 Test counties and cells (D27)", "phase1.t1_6_processing_test:select",
      deps=["counties", "t0_2", "t1_1"], decisions=["D27"], modules=["phase1/t1_6_processing_test.py"] + IO,
      sources=["kgs_satthick_1921"]),
    T("t1_4", "1.4 ET coverage", "phase1.t1_4_et_coverage:run", deps=["t0_2", "sample"],
      keys=["KS_TILES", "ET_MODELS", "ET_ENSEMBLE", "PANEL_YEARS", "SEASON_MONTHS"], decisions=["D21"],
      modules=["phase1/t1_4_et_coverage.py", "phase1/et_cells.py"] + DEF + IO + CELLS, sources=["openet_tiles"]),
    T("et_cells", "Blinded ET extraction", "phase1.et_cells:run", deps=["t0_2", "sample", "t1_4", "dry"],
      keys=["KS_TILES", "ET_MODELS", "ET_ENSEMBLE", "PANEL_YEARS", "SEASON_MONTHS"],
      decisions=["D20", "D21", "D29"], modules=["phase1/et_cells.py", "blinding.py"] + DEF + IO + CELLS,
      sources=["openet_tiles"]),
    T("t1_6", "1.6 Processing test log", "phase1.t1_6_processing_test:run",
      deps=["sample", "t0_2", "t1_4", "et_cells", "dry"], modules=["phase1/t1_6_processing_test.py", "blinding.py"]),
    T("t1_inv", "Preliminary county-year inventory", "phase1.t1_7_support_table:inventory",
      deps=["counties", "t0_2", "t1_1", "t1_5"], keys=["PANEL_YEARS"], decisions=["D25"],
      modules=["phase1/t1_7_support_table.py"] + DEF + IO),
    T("t1_7", "1.7 Complete support table", "phase1.t1_7_support_table:run",
      deps=["counties", "t0_2", "t1_1", "t1_2", "t1_3", "dry", "t1_4", "t1_5", "sample"],
      keys=["PANEL_YEARS", "SEASON_MONTHS"], decisions=["D21", "D25"],
      modules=["phase1/t1_7_support_table.py"] + DEF + IO),
    T("rep", "Replication region screen (metadata)", "phase1.rep_screen:run", deps=["counties"],
      keys=["REP_REGIONS", "PANEL_YEARS"], decisions=["D01", "D12", "D20"],
      modules=["phase1/rep_screen.py"] + DEF + IO + CELLS,
      sources=["gridmet", "aim_dir", "wells_ts", "haacker_dir"]),
    # ------------------------------------------------------------------ Phase 2: feasibility checks
    T("c2_1", "2.1 SMAP coverage in dry conditions", "phase2.c2_1_smap_dry:run",
      deps=["t1_1", "t1_3", "dry", "t1_7"], keys=["TARGETS"], modules=["phase2/c2_1_smap_dry.py"] + DEF),
    T("c2_2", "2.2 Reference footprints", "phase2.c2_2_footprints:run", deps=["counties", "t1_2", "t1_3", "t1_1"],
      keys=["TARGETS"], modules=["phase2/c2_2_footprints.py"]),
    T("c2_3", "2.3 ET coverage", "phase2.c2_3_et_coverage:run", deps=["t0_2", "t1_1", "t1_4"],
      keys=["TARGETS"], modules=["phase2/c2_3_et_coverage.py"]),
    T("c2_4", "2.4 Aquifer thickness", "phase2.c2_4_thickness:run", deps=["counties", "t0_2", "t1_5", "sample"],
      keys=["WELL_YEARS", "PANEL_YEARS", "TARGETS"], decisions=["D22", "D23", "D24", "D25", "D26"],
      modules=["phase2/c2_4_thickness.py", "common/stats.py"] + DEF + IO,
      sources=["haacker_dir", "kgs_satthick_1921", "usgs_change"]),
    T("c2_5", "2.5 Soil-moisture axis", "phase2.c2_5_sm_axis:run", deps=["t0_2", "t1_1", "t1_3", "et_cells"],
      decisions=["D08"], modules=["phase2/c2_5_sm_axis.py", "common/stats.py"] + DEF),
    T("c2_6", "2.6 Precision", "phase2.c2_6_precision:run",
      deps=["t0_2", "t1_1", "t1_3", "dry", "et_cells", "t1_7", "c2_4"], keys=["TARGETS"], decisions=["D28"],
      modules=["phase2/c2_6_precision.py", "common/stats.py"] + DEF),
    T("c2_7", "2.7 Well capacity", "phase2.c2_7_capacity:run", deps=["t0_2", "c2_4"],
      keys=["CAPACITY_N_START", "CAPACITY_N_EXTEND", "CAPACITY_SEED"], modules=["phase2/c2_7_capacity.py"] + IO,
      sources=["wimas_pods", "master_wells", "wwc5"]),
    T("c2_8", "2.8 Irrigation status after 2020", "phase2.c2_8_irr_status:run", deps=["t0_2"],
      modules=["phase2/c2_8_irr_status.py"] + IO, sources=["wimas_pods"]),
    T("c2_9", "2.9 Model support", "phase2.c2_9_model_support:run", deps=["t0_2", "t1_1", "t1_3", "t1_7", "c2_4"],
      modules=["phase2/c2_9_model_support.py", "common/stats.py"] + DEF),
    # ------------------------------------------------------------------ deliverables
    T("memo", "Feasibility memo", "phase2.memo:run",
      deps=["c2_1", "c2_2", "c2_3", "c2_4", "c2_5", "c2_6", "c2_7", "c2_8", "c2_9", "t1_7"], keys=["TARGETS"],
      modules=["phase2/memo.py"]),
    T("paper", "Paper items (pilot drafts)", "phase2.paper_items:run",
      deps=["t0_1", "counties", "t0_2", "t1_2", "t1_5", "t1_7", "c2_4", "c2_5"], keys=["COLORS"],
      modules=["phase2/paper_items.py"] + IO),
    T("fig_thick", "ET by aquifer thickness, simulated (Mieno Figs 3-5)", "phase2.fig_thickness:run",
      deps=["t0_2", "t1_1", "t1_7", "c2_4", "c2_6"], keys=["COLORS"], decisions=["D28"],
      modules=["phase2/fig_thickness.py", "phase2/c2_6_precision.py", "common/stats.py", "common/report.py"]),
    # ------------------------------------------------------------------ Main analysis: after the freeze only
    T("m_models", "Main: models C1, P1, J1, T1 and the two contrasts", "phase3.main_models:run",
      deps=["t0_2", "t1_1", "t1_7", "c2_4", "c2_8", "et_cells"], decisions=["D28", "D29", "D30", "D31", "D32"],
      modules=["phase3/main_models.py", "freeze.py", "common/stats.py"] + DEF),
    T("m_figs", "Main: Figure 1, Figure 3, ET by aquifer thickness (Mieno Figs 3-5)", "phase3.main_figs:run",
      deps=["m_models", "t0_2", "t1_1", "t1_7", "c2_4", "c2_8", "et_cells"], keys=["COLORS"],
      decisions=["D29", "D31", "D32"],
      modules=["phase3/main_figs.py", "phase3/main_models.py", "phase2/fig_thickness.py", "freeze.py"] + DEF),
]
MAIN = ["m_models", "m_figs"]

ALIASES = {"0.1": "t0_1", "0.2": "t0_2", "0.3": "t0_3", "0.4": "t0_4", "1.1": "t1_1", "1.2": "t1_2",
           "1.3": "t1_3", "1.4": "t1_4", "1.5": "t1_5", "1.6": "t1_6", "1.7": "t1_7", "2.1": "c2_1", "2.2": "c2_2",
           "2.3": "c2_3", "2.4": "c2_4", "2.5": "c2_5", "2.6": "c2_6", "2.7": "c2_7", "2.8": "c2_8", "2.9": "c2_9"}
GROUPS = {
    "phase0": ["t0_1", "t0_4", "counties", "t0_2"],
    "phase1": ["t1_1", "t1_2", "t1_3", "dry", "t1_5", "sample", "t1_4", "et_cells", "t1_6", "t1_inv", "t1_7", "rep"],
    "phase2": ["c2_1", "c2_2", "c2_3", "c2_4", "c2_5", "c2_6", "c2_7", "c2_8", "c2_9"],
    "no_et": ["t0_1", "t0_4", "counties", "t0_2", "t1_1", "t1_2", "t1_3", "dry", "t1_5", "sample", "t1_inv", "rep",
              "c2_2", "c2_4", "c2_7", "c2_8"],
    "all": [s["name"] for s in STAGES if s["name"] != "t0_3" and s["name"] not in MAIN],
    "main": MAIN,                    # the tasks refuse outside main mode and without a valid freeze record
}


def resolve(names):
    out = []
    for n in names:
        n = ALIASES.get(n, n)
        out += GROUPS.get(n, [n])
    known = {s["name"] for s in STAGES}
    bad = [n for n in out if n not in known]
    if bad:
        raise KeyError(f"unknown task(s) {bad}; tasks: {sorted(known)} or aliases {sorted(ALIASES)} or {sorted(GROUPS)}")
    return out


def runner():
    from .core.stages import Runner
    return Runner(STAGES)
