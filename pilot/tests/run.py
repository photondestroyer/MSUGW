"""Test runner (pytest is not installed locally; every test is also a valid pytest test).

    python -m pilot test                  definitions + unit + blinding
    python -m pilot test e2e              the smoke chain on real data (one county; SKIP when data are absent)
    python -m pilot test all              everything
    python -m pilot test -k word          tests whose name contains `word`

A test returns "SKIP: reason" when what it needs is not there. Exit code 1 on any failure.
"""
import importlib
import sys
import time
import traceback

GROUPS = dict(definitions=["definitions"], unit=["unit"], blinding=["blinding"], e2e=["e2e"])
DEFAULT = ["definitions", "unit", "blinding"]


def main(argv):
    pattern = None
    if "-k" in argv:
        i = argv.index("-k")
        pattern = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]
    groups = argv or DEFAULT
    if "all" in groups:
        groups = list(GROUPS)
    passed = failed = skipped = 0
    for grp in groups:
        for name in GROUPS[grp]:
            mod = importlib.import_module(f"pilot.tests.test_{name}")
            print(f"-- test_{name}")
            for fn_name in [n for n in mod.__dict__ if n.startswith("test_") and callable(mod.__dict__[n])]:
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
