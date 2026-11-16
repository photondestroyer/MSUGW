"""Command line:  python -m pilot <command>

    status                         every task: signature, complete?, blocked by a proposed decision, last error
    run <tasks> [--force]          run tasks and everything they need. <tasks>: names (t1_1), pilot IDs (1.1),
                                   or groups phase0 / phase1 / phase2 / no_et / all
    decisions                      the decision log (id, status, value)
    decisions show D04             one decision in full
    decisions fix D04 D05 --by "name"            proposed -> fixed
    decisions amend D04 --value '{"min_share": 0.6}' --reason "..." --by "name"
    audit                          the blinding audit of every table under the output root (STEP B-3)
    preflight                      which inputs of this mode exist, and which tasks can run now
    test [unit|definitions|blinding|e2e|all] [-k word]
    docs                           regenerate docs/pilot/DEFINITIONS.md and docs/pilot/DECISIONS.md
    sync-inputs [--dry-run]        copy the local-only inputs and the code to Google Drive for Colab

Mode and environment: PILOT_RUN_MODE (test | smoke | stage2 | stage3), PILOT_RUN_ENV (local | colab).
"""
import json
import os
import sys

from . import config as C


def _arg(argv, flag, default=None):
    if flag in argv:
        i = argv.index(flag)
        return argv[i + 1]
    return default


def cmd_status():
    from . import pipeline, decisions as D
    r = pipeline.runner()
    titles = {s["name"]: s["title"] for s in pipeline.STAGES}
    needs = {s["name"]: s["decisions"] for s in pipeline.STAGES}
    print(f"mode={C.RUN_MODE} env={C.RUN_ENV} out={C.out_root()}")
    for row in r.status():
        n = row["stage"]
        proposed = [d for d in needs[n] if D.status(d) != "fixed"]
        state = "done" if row["complete"] else ("blocked " + ",".join(proposed) if proposed and C.RUN_MODE in
                                                ("stage2", "stage3") else "to run")
        err = f"  last error: {row['last_error'][:90]}" if row["last_error"] and not row["complete"] else ""
        print(f"  {n:<9} {row['signature']}  {state:<22} {titles[n]}{err}")


def cmd_decisions(argv):
    from . import decisions as D
    if not argv:
        print(D.table())
        n = sum(e["status"] != "fixed" for e in D.entries().values())
        print(f"\n{n} of {len(D.entries())} decisions still proposed. Fix with: python -m pilot decisions fix <ids> --by <name>")
        return 0
    sub = argv[0]
    if sub == "show":
        print(json.dumps(D.entry(argv[1]), indent=1, default=str))
    elif sub == "fix":
        by = _arg(argv, "--by")
        ids = [a for a in argv[1:] if a.startswith("D") and a != "--by" and a != by]
        if not by:
            print("give your name: --by <name>")
            return 2
        if ids == ["Dall"] or "all" in argv[1:]:
            ids = list(D.entries())
        print("fixed:", ", ".join(D.fix(ids, by)) or "nothing (already fixed)")
    elif sub == "amend":
        did, val, why, by = argv[1], _arg(argv, "--value"), _arg(argv, "--reason"), _arg(argv, "--by")
        if not (val and why and by):
            print("amend needs --value '<json>' --reason \"...\" --by <name>")
            return 2
        D.amend(did, json.loads(val), why, by)
        print(f"{did} amended; old value kept in its history")
    else:
        print(__doc__)
        return 2
    return 0


def cmd_run(argv):
    from . import pipeline
    names = pipeline.resolve([a for a in argv if not a.startswith("--")])
    pipeline.runner().run(names, force="--force" in argv)
    return 0


def cmd_audit():
    from . import blinding
    ok, problems, n = blinding.audit()
    print(f"blinding audit of {C.out_root()}: {n} ET tables checked -> {'PASS' if ok else 'FAIL'}")
    for p in problems:
        print("  ", p)
    return 0 if ok else 1


def cmd_preflight():
    from .common import io
    print(f"mode={C.RUN_MODE} env={C.RUN_ENV}")
    for key in C.DATASETS:
        p = C.src(key)
        print(f"  {'ok ' if os.path.exists(p) else '-- '} {key:<20} {p}")
    try:
        p, label = io.smap_source()
        print(f"SMAP used: {label} ({p})")
    except Exception as e:
        print(f"SMAP: {e}")
    return 0


def cmd_docs():
    from . import docs_gen
    for p in docs_gen.write_all():
        print("wrote", p)
    return 0


def cmd_sync(argv):
    from . import sync
    sync.run(dry_run="--dry-run" in argv)
    return 0


def quiet():
    """Library warnings that say nothing about the data (fragmented frames, rasterio window arithmetic)."""
    import warnings
    warnings.filterwarnings("ignore", message="DataFrame is highly fragmented")
    warnings.filterwarnings("ignore", message="divide by zero encountered in scalar divide")
    warnings.filterwarnings("ignore", category=RuntimeWarning, module="rasterio")


def main(argv):
    C.check()
    quiet()
    if not argv:
        print(__doc__)
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == "status":
        cmd_status()
        return 0
    if cmd == "decisions":
        return cmd_decisions(rest)
    if cmd == "run":
        return cmd_run(rest)
    if cmd == "audit":
        return cmd_audit()
    if cmd == "preflight":
        return cmd_preflight()
    if cmd == "test":
        from .tests import run as T
        return T.main(rest)
    if cmd == "docs":
        return cmd_docs()
    if cmd == "sync-inputs":
        return cmd_sync(rest)
    if cmd == "r-env":
        from . import renv_setup
        if "--create" in rest:
            return renv_setup.run(dry_run=False)
        return renv_setup.run(dry_run=True)
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
