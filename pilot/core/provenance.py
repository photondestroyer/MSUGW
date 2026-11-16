"""Provenance: what every output was made from, and by which code and settings.

fingerprint(path)   size + mtime + sampled MD5 (first 4 MB and 64 x 64 KB at fixed offsets, ~8 MB read)
                    + header digest (NetCDF: dims/variables/attributes; GeoTIFF: profile/tags);
                    full MD5 on request (cached per path+size+mtime, so a multi-GB file is hashed once).
                    mtime changes when a file is downloaded, so file IDENTITY is size + sampled MD5 +
                    header digest (+ full MD5 when available); mtime is informational only.
code_hash(modules)  MD5 over the source of the modules a stage declares.
signature(...)      8 hex characters naming a stage's output folder: code hash + the config keys the
                    stage declares + the signatures of the stages it reads.
Ledger              append-only JSON lines: one line per unit of work (start/done/error).
"""
import datetime as _dt
import hashlib
import json
import os
import platform
import socket

from .. import config as C

SAMPLE_HEAD = 4 * 1024 * 1024
SAMPLE_CHUNK = 64 * 1024
SAMPLE_N = 64
PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def utcnow():
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (str, int, float, bool)) or x is None:
        return x
    if hasattr(x, "tolist"):
        return x.tolist()
    return str(x)


def dumps(x):
    return json.dumps(_jsonable(x), sort_keys=True, separators=(",", ":"))


def sampled_md5(path):
    size = os.path.getsize(path)
    h = hashlib.md5()
    with open(path, "rb") as f:
        h.update(f.read(SAMPLE_HEAD))
        if size > SAMPLE_HEAD:
            for k in range(SAMPLE_N):
                off = SAMPLE_HEAD + (size - SAMPLE_HEAD - SAMPLE_CHUNK) * k // max(SAMPLE_N - 1, 1)
                f.seek(max(off, 0))
                h.update(f.read(SAMPLE_CHUNK))
    return h.hexdigest()


def full_md5(path, cache=True):
    st = os.stat(path)
    key = hashlib.md5(f"{os.path.abspath(path)}|{st.st_size}|{int(st.st_mtime)}".encode()).hexdigest()
    cdir = os.path.join(C.scratch_root(), "md5_cache")
    cfile = os.path.join(cdir, key + ".txt")
    if cache and os.path.exists(cfile):
        return open(cfile).read().strip()
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(16 * 1024 * 1024), b""):
            h.update(chunk)
    md5 = h.hexdigest()
    if cache:
        os.makedirs(cdir, exist_ok=True)
        with open(cfile, "w") as f:
            f.write(md5)
    return md5


def header_digest(path):
    """MD5 of the file's self-description; changes when variables, types, shapes or attributes change."""
    ext = os.path.splitext(path)[1].lower()
    desc = {}
    try:
        if ext == ".nc":
            import netCDF4 as nc
            with nc.Dataset(path) as ds:
                def grp(g):
                    return dict(dims={k: len(v) for k, v in g.dimensions.items()},
                                attrs={a: str(g.getncattr(a)) for a in g.ncattrs() if a not in ("history",)},
                                vars={n: dict(dtype=str(v.dtype), dims=list(v.dimensions),
                                              attrs={a: str(v.getncattr(a)) for a in v.ncattrs()})
                                      for n, v in g.variables.items()},
                                groups={n: grp(s) for n, s in g.groups.items()})
                desc = grp(ds)
        elif ext in (".tif", ".tiff"):
            import rasterio
            with rasterio.open(path) as r:
                desc = dict(profile={k: str(v) for k, v in r.profile.items()}, tags=r.tags(),
                            descriptions=list(r.descriptions))
        else:
            return None
    except Exception as e:                      # a digest must never break a run; record why it is missing
        return f"unavailable: {type(e).__name__}: {e}"
    return hashlib.md5(dumps(desc).encode()).hexdigest()


def fingerprint(path, full=False):
    st = os.stat(path)
    fp = dict(path=str(path).replace("\\", "/"), size=st.st_size,
              mtime=_dt.datetime.fromtimestamp(st.st_mtime, _dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
              sampled_md5=sampled_md5(path), header=header_digest(path))
    if full:
        fp["md5"] = full_md5(path)
    return fp


def same_file(fp_a, fp_b):
    """Identity of two fingerprints (mtime ignored; full MD5 compared when both have it)."""
    keys = ["size", "sampled_md5", "header"] + (["md5"] if "md5" in fp_a and "md5" in fp_b else [])
    return all(fp_a.get(k) == fp_b.get(k) for k in keys)


def code_hash(modules):
    """modules: paths relative to the package directory, e.g. 'step0/blocks.py'."""
    h = hashlib.md5()
    for m in sorted(modules):
        with open(os.path.join(PKG_DIR, m), "rb") as f:
            h.update(m.encode() + b"\0" + f.read().replace(b"\r\n", b"\n"))
    return h.hexdigest()


def config_subset(keys):
    snap = C.snapshot()
    out = {}
    for k in keys:
        if k.startswith("SRC:"):                # a source the stage reads: which file, not only which setting
            out[k] = C.source_id(k[4:])
            continue
        if k not in snap:
            raise KeyError(f"config has no setting {k!r}")
        out[k] = snap[k]
    return out


def signature(stage, keys, modules, upstream):
    parts = dict(stage=stage, config=config_subset(keys), upstream=dict(sorted(upstream.items())),
                 mode=C.RUN_MODE)
    if C.CODE_IN_SIGNATURE:
        parts["code"] = code_hash(modules)
    return hashlib.md5(dumps(parts).encode()).hexdigest()[:8]


def record(stage, sig, inputs, keys, modules, upstream, extra=None):
    """The provenance dictionary written into every output (as a JSON global attribute)."""
    return dict(stage=stage, signature=sig, created=utcnow(), run_env=C.RUN_ENV, run_mode=C.RUN_MODE,
                host=socket.gethostname(), python=platform.python_version(),
                code_hash=code_hash(modules), modules=sorted(modules), config=config_subset(keys),
                upstream=upstream, inputs={k: v for k, v in sorted(inputs.items())}, **(extra or {}))


class Ledger:
    """<out_root>/ledger.jsonl - one JSON line per event, flushed and fsynced."""

    def __init__(self, path):
        self.path = path
        os.makedirs(os.path.dirname(path), exist_ok=True)

    def write(self, **kw):
        kw = dict(time=utcnow(), **kw)
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(dumps(kw) + "\n")
            f.flush()
            os.fsync(f.fileno())

    def read(self):
        if not os.path.exists(self.path):
            return []
        out = []
        with open(self.path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        out.append(json.loads(line))
                    except json.JSONDecodeError:
                        pass                      # a line cut by a crash is ignored, not fatal
        return out
