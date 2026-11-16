"""End to end: the whole task chain in smoke mode (real data, one county), then the blinding audit.

This replaces the synthetic data root planned at first: the pilot reads 20 sources in 12 formats, and a synthetic
copy of all of them would mostly test the generator. The smoke county exercises the real formats, frames, flags
and joins; it is SKIPPED where the data are absent (e.g. Colab without the synced inputs).
Smoke outputs are code tests only; they are never pilot results.
"""
import os
import subprocess
import sys


def test_smoke_chain():
    from pilot import config as C
    if not os.path.exists(C.DATASETS["gridmet"][0]):
        return "SKIP: real inputs not on this machine"
    env = dict(os.environ, PILOT_RUN_MODE="smoke")
    if not os.environ.get("PILOT_SMOKE_OPENET") and not os.path.isdir(C.DATASETS["openet_tiles"][0]):
        env["PILOT_SMOKE_OPENET"] = f"{C.LOCAL_DRIVE}/OPENET_v2_0_monthly_30m_tiles"   # code path only, when the
        #                                                                                v2.1 Kansas folder is absent
    p = subprocess.run([sys.executable, "-m", "pilot", "run", "all"], cwd=C.REPO_ROOT, env=env, capture_output=True,
                       text=True, timeout=6 * 3600)
    assert p.returncode == 0, p.stdout[-3000:] + p.stderr[-3000:]
    a = subprocess.run([sys.executable, "-m", "pilot", "audit"], cwd=C.REPO_ROOT, env=env, capture_output=True,
                       text=True)
    assert a.returncode == 0, a.stdout
