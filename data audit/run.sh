#!/bin/bash
# usage: run.sh script.py [args]   (conda env MSUGW with pyexpat fix)
SP="C:/Users/AlienX/AppData/Local/Temp/claude/G--MSU-GWB-datasets/d2155122-1714-4320-b37e-67df7844f446/scratchpad"
export PATH="C:/Users/AlienX/anaconda3/envs/MSUGW/Library/bin:$PATH"
export PYTHONPATH="$SP/test/site_fix"
export PYTHONIOENCODING=utf-8
exec C:/Users/AlienX/anaconda3/envs/MSUGW/python.exe "$@"
