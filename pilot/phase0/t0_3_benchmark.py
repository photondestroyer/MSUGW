"""Task 0.3  Benchmark reproduction.

"Reproduce the corn drought-response and irrigated-share figures of Mieno et al. (2024) from their public code and
processed data. Record sample counts, coefficients and any difference from the published figures."
Physical idea: Mieno et al. found that drought losses on irrigated land grow as the aquifer thins - the same
mechanism at county scale. Reproducing it first shows the inputs and code path behave before the pilot adds
satellite ET.
Decisions: none.   Deliverables: comparison_note.md, the reproduced figures (PDF), the exported numbers (CSV).

Safety (the user's condition): R runs from an isolated conda environment (R_ENV_PREFIX, created by
`python -m pilot r-env`), on a working copy of their repository in scratch, with Rscript --vanilla. Their
repository and every pilot input are only read.
"""
import os
import shutil
import subprocess

R_ENV_PREFIX = "H:/MSU_GWB/r_env_mieno"
N_BOOT = 1000          # their value; lower it only with a recorded reason
PUBLISHED = [          # statements of the paper compared against the reproduction
    dict(item="corn: yield loss when the aquifer thins from ~70 m to ~40 m (230 to 130 ft) at a 950 mm deficit",
         value="about 15 bu/acre", sat_from=70, sat_to=40, balance_mm=950, expected_bu=15.0),
    dict(item="corn: yield loss when the aquifer thins from ~100 m to ~70 m (330 to 230 ft) at a 950 mm deficit",
         value="negligible", sat_from=100, sat_to=70, balance_mm=950, expected_bu=0.0),
]


def rscript():
    for p in (f"{R_ENV_PREFIX}/Scripts/Rscript.exe", f"{R_ENV_PREFIX}/Lib/R/bin/Rscript.exe",
              f"{R_ENV_PREFIX}/bin/Rscript"):
        if os.path.exists(p):
            return p
    return None


def run(ctx):
    import pandas as pd
    import numpy as np
    from .. import config as C
    from ..common import io, report
    from ..core.errors import PipelineStop

    # ══════ STEP 0.3-1  the isolated R environment ══════
    rs = rscript()
    if rs is None:
        raise PipelineStop(f"no R environment at {R_ENV_PREFIX}. Create it after checking the package list: "
                           f"python -m pilot r-env --dry-run, then python -m pilot r-env --create")

    # ══════ STEP 0.3-2  working copy of their repository (theirs is only read) ══════
    src = io.path("mieno_repo")
    work = os.path.join(ctx.scratch, "mieno_work")
    if os.path.exists(work):
        shutil.rmtree(work)
    shutil.copytree(src, work, ignore=shutil.ignore_patterns(".git", "renv", "Results"))
    raw = os.path.join(work, "Data", "data-raw")
    os.makedirs(raw, exist_ok=True)
    aq = os.path.dirname(C.src("aquifer_shp"))
    for f in os.listdir(aq):
        if f.startswith("hp_bound2010."):
            shutil.copy2(os.path.join(aq, f), raw)
    ctx.log.info("STEP 0.3-2  working copy at %s", work)

    # ══════ STEP 0.3-3  run their analysis and figures (Rscript --vanilla) ══════
    out = os.path.join(ctx.scratch, "mieno_out")
    os.makedirs(out, exist_ok=True)
    log_path = os.path.join(ctx.scratch, "mieno_R.log")
    runner = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_mieno.R")
    env = dict(os.environ, R_LIBS_USER="", R_PROFILE_USER="", R_ENVIRON_USER="")
    with open(log_path, "w", encoding="utf-8") as lg:
        p = subprocess.run([rs, "--vanilla", runner, work, str(N_BOOT), out], stdout=lg, stderr=subprocess.STDOUT,
                           env=env, timeout=36 * 3600)
    report.text(ctx, "mieno_R.log", open(log_path, encoding="utf-8", errors="ignore").read())
    if p.returncode != 0:
        raise PipelineStop(f"R exited with {p.returncode}; see {ctx.out}/mieno_R.log")

    # ══════ STEP 0.3-4  numbers against the published statements ══════
    counts = pd.read_csv(f"{out}/mieno_sample_counts.csv")
    avg = pd.read_csv(f"{out}/mieno_avg_yield_sum.csv")
    info = pd.read_csv(f"{out}/mieno_run_info.csv").iloc[0].to_dict()
    rows = []
    corn = avg[avg["crop"] == "corn"]
    for s in PUBLISHED:
        b = corn["balance"].to_numpy()
        bal = b[np.argmin(np.abs(b - s["balance_mm"]))]
        at = corn[corn["balance"] == bal]

        def y(sat):
            r = at.iloc[np.argmin(np.abs(at["sat"].to_numpy() - sat))]
            return float(r["avg_yield"]), float(r["sat"])
        y0, s0 = y(s["sat_from"])
        y1, s1 = y(s["sat_to"])
        rows.append(dict(statement=s["item"], published=s["value"], balance_used_mm=float(bal), sat_from_m=s0,
                         sat_to_m=s1, reproduced_loss_bu_acre=round(y0 - y1, 2), expected_bu_acre=s["expected_bu"]))
    comp = pd.DataFrame(rows)
    report.table(ctx, counts, "mieno_sample_counts")
    report.table(ctx, comp, "mieno_comparison", units={"loss": "bu/acre", "balance": "mm", "sat": "m"})
    report.table(ctx, avg, "mieno_avg_yield_sum", units={"avg_yield": "bu/acre", "balance": "mm", "sat": "m"})
    figs = []
    fdir = os.path.join(work, "Results", "Figures")
    if os.path.isdir(fdir):
        for f in sorted(os.listdir(fdir)):
            shutil.copy2(os.path.join(fdir, f), os.path.join(ctx.out, f))
            figs.append(f)

    # ══════ STEP 0.3-5  comparison note ══════
    L = ["# Benchmark reproduction: Mieno et al. (2024)", "",
         "Code: github.com/tmieno2/Drought-Production-Risk-Aquifer (their processed data in Data/data-processed).",
         f"R: {info.get('r_version')}; fixest {info.get('fixest')}; mgcv {info.get('mgcv')}; bootstrap draws "
         f"{int(info.get('n_boot'))} (theirs: 1000), run serially (theirs: parallel), so bootstrap intervals "
         "agree in distribution, not draw for draw. Point estimates do not depend on the bootstrap.", "",
         "## Sample", "", report.md(counts), "", "## Against the paper", "",
         report.md(comp), "", f"## Figures reproduced ({len(figs)})", ""] + [f"- {f}" for f in figs] + \
        ["", "Differences to look at: any statement whose reproduced loss differs from the published one by more "
         "than the rounding of the paper (the paper rounds to whole bushels)."]
    report.text(ctx, "comparison_note.md", "\n".join(L))
    ctx.complete(dict(samples=counts.to_dict("records"), comparison=comp.to_dict("records"), figures=figs,
                      figures_ok=bool(info.get("figures_ok"))))
