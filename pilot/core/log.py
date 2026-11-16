"""Logging to the console (Colab cell output) and to <out_root>/logs/run_<UTC date>.log."""
import datetime as _dt
import logging
import os

from .. import config as C

_FMT = "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s"
_configured = set()


def get(name="pilot"):
    log = logging.getLogger(f"pilot.{name}" if name != "pilot" else "pilot")
    root = logging.getLogger("pilot")
    key = (C.RUN_ENV, C.RUN_MODE, C.out_root())
    if key not in _configured:
        for h in list(root.handlers):
            h.close()
        root.handlers.clear()
        _configured.clear()
        root.setLevel(logging.INFO)
        root.propagate = False
        sh = logging.StreamHandler()
        sh.setFormatter(logging.Formatter(_FMT, "%H:%M:%S"))
        root.addHandler(sh)
        try:
            d = os.path.join(C.out_root(), "logs")
            os.makedirs(d, exist_ok=True)
            fh = logging.FileHandler(os.path.join(d, f"run_{_dt.datetime.utcnow():%Y%m%d}.log"), encoding="utf-8")
            fh.setFormatter(logging.Formatter(_FMT, "%Y-%m-%d %H:%M:%S"))
            root.addHandler(fh)
        except OSError:
            pass                                   # no writable output root yet (preflight reports it)
        for noisy in ("rasterio", "pyogrio", "fiona", "urllib3", "matplotlib", "numba"):
            logging.getLogger(noisy).setLevel(logging.ERROR)
        _configured.add(key)
    return log
