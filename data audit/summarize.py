"""One line per scanned variable with its statistics and every flag the scan raised."""
import json
import sys

from audit_paths import OUT

name = sys.argv[1] if len(sys.argv) > 1 else "flagged"
r = json.load(open(f"{OUT}/scan_{name}.json"))


def f(x):
    return "-" if x is None else (f"{x:.4g}" if isinstance(x, (int, float)) else str(x))


for x in sorted(r, key=lambda z: (z["file"], z["var"])):
    if "fatal" in x:
        print(x["file"], x["var"], "FATAL", x["fatal"])
        continue
    if "valid_frac" not in x:
        print(f"{x['file'][:34]:34s} {x['var'][:34]:34s} (text) empty={x.get('n_empty')} unique={x.get('n_unique_read')}")
        continue
    flags = []
    if x["errors"]:
        flags.append(f"READERR={len(x['errors'])}")
    if x["n_inf"]:
        flags.append(f"INF={x['n_inf']}")
    if x["sentinels"]:
        flags.append(f"sent={x['sentinels']}")
    for k in ("n_abs_gt_1e20", "n_near_9999_not_fill", "n_tiny_negative", "n_missing_value"):
        if x.get(k):
            flags.append(f"{k[2:]}={x[k]}")
    if x["n_outside_phys"]:
        flags.append(f"outPhys{x['phys_bounds']}={x['n_outside_phys']}")
    if x["n_outside_valid_range_attr"]:
        flags.append(f"outValidRange={x['n_outside_valid_range_attr']}")
    if x["n_outside_valid_minmax_attr"]:
        flags.append(f"outValidMinMax={x['n_outside_valid_minmax_attr']}")
    if x["fill"] is not None and x["fill"] == x["fill"] and x["n_nan"]:
        flags.append(f"NaN-with-nonNaN-fill={x['n_nan']}")
    p = x.get("pcts", {})
    print(f"{x['file'][:34]:34s} {x['var'][:34]:34s} {x['dtype']:7s} valid={x['valid_frac']:.4f} min={f(x['min'])} "
          f"p1={f(p.get('1'))} med={f(p.get('50'))} p99={f(p.get('99'))} max={f(x['max'])} {' '.join(flags)}")
