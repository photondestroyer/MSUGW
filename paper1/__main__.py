"""Command line:  python -m paper1 <command>

    preflight            check config, inputs (identity, presence, sizes), scratch, libraries, GPU, Drive speed
    status               stage table: signature, complete?, units done, last error, summary
    run <what> [--force] run stages and everything they need; <what> = stage names, "p11", "p12" or "all"
    test [unit|provenance|torch|all] [-k word]
    docs                 regenerate docs/CONFIG.md from config.py

RUN_ENV / RUN_MODE come from config.py or the environment variables P1_RUN_ENV / P1_RUN_MODE.
"""
import json
import os
import shutil
import sys
import time

from . import config as C


def preflight(write=True):
    from .core import paths, provenance as P, drive
    from .core.errors import SourceError
    from . import checks
    C.check()
    rep = dict(run_env=C.RUN_ENV, run_mode=C.RUN_MODE, data_root=C.data_root(), out_root=C.out_root(),
               problems=[], inputs={}, libraries={})
    if not os.path.isdir(C.data_root()):
        rep["problems"].append(f"data root {C.data_root()} does not exist")
    for key, spec in C.DATASETS.items():
        if spec["kind"] != "file":
            continue
        try:
            p = paths.dataset_path(key)
            info = dict(path=p, size=os.path.getsize(p))
            if key in checks.FILE_CHECKS:
                probs = checks.FILE_CHECKS[key](p)
                rep["problems"] += probs
                info["problems"] = probs
            if C.RUN_ENV == "colab":
                info["drive_md5"] = drive.drive_md5(p)
            rep["inputs"][key] = info
        except SourceError as e:
            rep["problems"].append(str(e))
            rep["inputs"][key] = dict(error=str(e))
    try:
        from .step0 import tiles as TL
        frames = TL.load_frames()
        tl = TL.selected_tiles(frames)
        rep["tiles"] = tl
        try:
            rep["frames_agree"] = TL.check_frames_agree(frames, tl)
        except Exception as e:
            rep["problems"].append(str(e))
        inv = checks.inventory(tl, C.agg_years())
        rep["tile_files"] = dict(present=sum(r["exists"] for r in inv), missing=[r["path"] for r in inv if not r["exists"]],
                                 gb=round(sum(r["size"] or 0 for r in inv) / 1e9, 2))
        checked = set()
        for r in inv:
            if not r["exists"] or (C.RUN_MODE == "full" and r["dataset"] in checked):
                continue                                   # full mode: one file per dataset is checked
            rep["problems"] += checks.check_tile_file(r["dataset"], r["path"], frames.get(r["tile"]))
            checked.add(r["dataset"])
        miss = [r for r in inv if not r["exists"]]
        tolerated = [r for r in miss if r["dataset"] == "nwi" and C.nwi_missing_policy() == "pass"]
        if tolerated:
            rep["notes"] = rep.get("notes", []) + [f"{len(tolerated)} NWI tile(s) missing: the NWI rule passes there "
                                                   f"(NWI missing policy 'pass' in {C.RUN_MODE} mode)"]
        if len(miss) > len(tolerated):
            rep["problems"].append(f"{len(miss) - len(tolerated)} tile files missing (see tile_files.missing)")
        # Drive / disk read speed on one Landsat file
        first = next((r for r in checks.inventory(tl[:1], C.agg_years()[:1]) if r["exists"]), None)
        if first:
            t0 = time.time()
            with open(first["path"], "rb") as f:
                n = len(f.read(64 * 1024 * 1024))
            rep["read_mb_per_s"] = round(n / 1e6 / max(time.time() - t0, 1e-6), 1)
    except SourceError as e:
        rep["problems"].append(str(e))
    try:
        os.makedirs(C.scratch_root(), exist_ok=True)
        rep["scratch_free_gb"] = round(shutil.disk_usage(C.scratch_root()).free / 1e9, 1)
    except OSError as e:
        rep["problems"].append(f"scratch {C.scratch_root()}: {e}")
    for lib in ("numpy", "netCDF4", "pyproj", "rasterio", "geopandas", "pyogrio", "shapely", "scipy", "pysheds",
                "requests", "matplotlib", "torch"):
        try:
            m = __import__(lib)
            rep["libraries"][lib] = getattr(m, "__version__", "ok")
        except ImportError:
            rep["libraries"][lib] = None
    if rep["libraries"].get("torch"):
        import torch
        rep["gpu"] = torch.cuda.get_device_name(0) if torch.cuda.is_available() else None
    if write:
        try:
            drive.write_text(f"{C.out_root()}/preflight.json", json.dumps(rep, indent=1, default=str))
        except OSError:
            pass
    return rep


def print_status():
    from .pipeline import runner
    rows = runner().status()
    print(f"RUN_ENV={C.RUN_ENV} RUN_MODE={C.RUN_MODE} outputs: {C.out_root()}")
    print(f"{'stage':12} {'signature':9} {'complete':8} {'units':>6}  last error / summary")
    for r in rows:
        tail = r["last_error"] or (json.dumps(r["summary"], default=str)[:110] if r["summary"] else "")
        print(f"{r['stage']:12} {r['signature']:9} {str(r['complete']):8} {r['units_done']:>6}  {tail}")
        if r["other_signatures"]:
            print(f"{'':12} older outputs with other settings/code: {', '.join(r['other_signatures'])}")


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == "preflight":
        rep = preflight()
        print(json.dumps({k: v for k, v in rep.items() if k != "inputs"}, indent=1, default=str))
        print("\nPROBLEMS:" if rep["problems"] else "\nno problems found")
        for p in rep["problems"]:
            print("  -", p)
        return 1 if rep["problems"] else 0
    if cmd == "status":
        print_status()
        return 0
    if cmd == "run":
        from .pipeline import runner, P11, P12
        force = "--force" in rest
        names = [x for x in rest if not x.startswith("--")]
        expand = dict(p11=P11, p12=P12, all=P11 + P12)
        targets = [s for n in names for s in expand.get(n, [n])]
        if not targets:
            print("run what? stage names, p11, p12 or all")
            return 2
        runner().run(targets, force=force)
        print_status()
        return 0
    if cmd == "test":
        from .tests.run import main as tmain
        return tmain(rest)
    if cmd == "docs":
        from .docs_gen import write_config_md
        print(write_config_md())
        return 0
    print(f"unknown command {cmd!r}\n{__doc__}")
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
