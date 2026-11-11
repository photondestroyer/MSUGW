"""Test runner (pytest is not installed locally; Colab has it, and every test is also a valid pytest test).

    python -m paper1 test                 unit + e2e + provenance (+ torch where torch exists)
    python -m paper1 test unit            synthetic unit tests (seconds)
    python -m paper1 test e2e             every stage up to P1.1 on the synthetic data root (minutes)
    python -m paper1 test provenance      real files in the configured data root (SKIP when absent)
    python -m paper1 test torch           M1, and train/predict/verify on the synthetic root (needs torch)
    python -m paper1 test -k word         tests whose name contains `word`

A test returns "SKIP" (with a printed reason) when what it needs is not there. Exit code 1 on any failure.
"""
import importlib
import sys
import time
import traceback

GROUPS = dict(unit=["unit_core", "unit_step0", "unit_p11", "unit_p12"], e2e=["e2e"], provenance=["provenance"],
              torch=["torch_m1"])
ORDER = ["unit", "e2e", "provenance", "torch"]


def main(argv):
    pattern = None
    if "-k" in argv:
        i = argv.index("-k")
        pattern = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]
    groups = argv or ORDER
    if "all" in groups:
        groups = ORDER
    passed = failed = skipped = 0
    for grp in groups:
        for name in GROUPS[grp]:
            try:
                mod = importlib.import_module(f"paper1.tests.test_{name}")
            except ImportError as e:
                print(f"-- test_{name}: not run in this environment ({e.name} missing)")
                continue
            print(f"-- test_{name}")
            fns = [n for n in mod.__dict__ if n.startswith("test_") and callable(mod.__dict__[n])]
            for fn_name in fns:
                if pattern and pattern not in fn_name:
                    continue
                t0 = time.time()
                try:
                    out = getattr(mod, fn_name)()
                except Exception:
                    failed += 1
                    print(f"   FAIL  {fn_name}")
                    traceback.print_exc()
                    continue
                if isinstance(out, str) and out.startswith("SKIP"):
                    skipped += 1
                    print(f"   skip  {fn_name}  {out[4:].strip(': ')}")
                else:
                    passed += 1
                    print(f"   ok    {fn_name}  ({time.time() - t0:.1f}s)")
    print(f"{passed} passed, {failed} failed, {skipped} skipped")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
