"""The isolated R environment for Task 0.3 (python -m pilot r-env --dry-run | --create).

Safety: a conda environment created BY PREFIX at R_ENV_PREFIX from conda-forge binaries. No admin rights, no
registry or PATH change, nothing outside that folder; deleting the folder removes it completely.
--dry-run prints the exact package list conda would install and changes nothing.
"""
import shutil
import subprocess

from .phase0.t0_3_benchmark import R_ENV_PREFIX

PACKAGES = ["r-base=4.4", "r-here", "r-dplyr", "r-purrr", "r-ggplot2", "r-stringr", "r-tidyr", "r-data.table",
            "r-sf", "r-mgcv", "r-fixest", "r-cowplot", "r-broom", "r-lubridate", "r-scales", "r-patchwork",
            "r-tidycensus", "r-knitr"]
CRAN_ONLY = ["nombre"]       # not on conda-forge: a small pure-R package, installed from CRAN into the env's library


def conda():
    c = shutil.which("conda") or shutil.which("conda.bat")
    if not c:
        raise RuntimeError("conda not found on PATH")
    return c


def command(dry_run):
    cmd = [conda(), "create", "--yes", "--prefix", R_ENV_PREFIX, "--override-channels", "-c", "conda-forge"] + PACKAGES
    return cmd + (["--dry-run"] if dry_run else [])


def cran_step():
    from .phase0.t0_3_benchmark import rscript
    rs = rscript()
    lib = f"{R_ENV_PREFIX}/Lib/R/library"
    expr = (f"install.packages(c({', '.join(repr(p) for p in CRAN_ONLY)}), lib='{lib}', "
            f"repos='https://cloud.r-project.org', type='source')")
    return [rs or f"{R_ENV_PREFIX}/Scripts/Rscript.exe", "--vanilla", "-e", expr]


def run(dry_run=True):
    cmd = command(dry_run)
    print(" ".join(cmd))
    p = subprocess.run(cmd, capture_output=True, text=True)
    print(p.stdout[-20000:])
    if p.returncode != 0:
        print(p.stderr[-5000:])
        return p.returncode
    step2 = cran_step()
    if dry_run:
        print("then, from CRAN into the environment's own library:\n  " + " ".join(step2))
        return 0
    q = subprocess.run(step2, capture_output=True, text=True)
    print(q.stdout[-4000:], q.stderr[-4000:])
    return q.returncode
