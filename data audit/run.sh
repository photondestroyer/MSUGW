#!/bin/bash
# usage: run.sh script.py [args]   (conda env MSUGW; pyexpat must be imported before GDAL, hence site_fix)
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$HERE/site_fix"
[ -f "$HERE/site_fix/sitecustomize.py" ] || echo "import pyexpat" > "$HERE/site_fix/sitecustomize.py"
export PATH="C:/Users/AlienX/anaconda3/envs/MSUGW/Library/bin:$PATH"
export PYTHONPATH="$HERE/site_fix:$HERE"
export PYTHONIOENCODING=utf-8
# one BLAS/OpenMP thread per process: with 20 threads, importing numpy alone reserves ~0.6 GB of commit
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 GDAL_NUM_THREADS=1
cd "$HERE" && exec C:/Users/AlienX/anaconda3/envs/MSUGW/python.exe "$@"
