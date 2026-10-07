"""Minimal test runner (no pytest in either environment).

    python -m residual_et.tests.run            # every module whose imports are available
    python -m residual_et.tests.run core data  # selected modules
    python -m residual_et.tests.run -k resume  # tests whose name contains a word

Exit code 1 if any test fails. Functions are also valid pytest tests.
"""
import importlib
import sys
import time
import traceback

MODULES = ("core", "data", "torch")


def main(argv):
    pattern = None
    if "-k" in argv:
        i = argv.index("-k")
        pattern = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]
    names = argv or MODULES
    passed = failed = skipped = 0
    for name in names:
        try:
            mod = importlib.import_module(f"residual_et.tests.test_{name}")
        except ImportError as e:
            print(f"-- test_{name}: not run in this environment ({e.name} missing)")
            continue
        print(f"-- test_{name}")
        for fn_name in [n for n in dir(mod) if n.startswith("test_")]:
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
            if out == "SKIP":
                skipped += 1
                print(f"   skip  {fn_name}")
            else:
                passed += 1
                print(f"   ok    {fn_name}  ({time.time() - t0:.1f}s)")
    print(f"{passed} passed, {failed} failed, {skipped} skipped")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
