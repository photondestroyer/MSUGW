"""The decision log: load, read values, block tasks on proposed decisions, fix and amend.

    value("D04")              the value dict of a decision
    values()                  {id: value} of every decision (part of every task signature, see config.snapshot)
    require(ids, task)        raises PipelineStop in stage2/stage3 when any of ids is still "proposed"
    fix(ids, by)              status proposed -> fixed, with name and date
    amend(id, value, reason)  a fixed value changes; the old value is appended to its history
    log_hash()                MD5 of the log file (written into every deliverable's sidecar)

The file is rewritten by fix/amend only (yaml.safe_dump, field order kept). Comments at the top are kept.
"""
import datetime as _dt
import hashlib
import json
import os

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "decisions.yaml")
FIELDS = ["id", "title", "value", "reason", "tasks", "status", "proposed_by", "proposed_on", "fixed_by",
          "fixed_on", "history"]
_cache = {}


def _load(path=PATH):
    import yaml
    st = os.stat(path)
    key = (path, st.st_size, st.st_mtime_ns)
    if key not in _cache:
        with open(path, encoding="utf-8") as f:
            entries = yaml.safe_load(f) or []
        ids = [e["id"] for e in entries]
        if len(ids) != len(set(ids)):
            raise ValueError(f"{path}: duplicate decision ids")
        for e in entries:
            missing = [k for k in FIELDS if k not in e]
            if missing:
                raise ValueError(f"{path}: {e.get('id')} lacks {missing}")
            if e["status"] not in ("proposed", "fixed"):
                raise ValueError(f"{path}: {e['id']} status must be proposed or fixed, not {e['status']!r}")
        _cache.clear()
        _cache[key] = {e["id"]: e for e in entries}
    return _cache[key]


def entries(path=PATH):
    return _load(path)


def entry(did, path=PATH):
    d = _load(path)
    if did not in d:
        raise KeyError(f"no decision {did} in {path}")
    return d[did]


def value(did, path=PATH):
    return entry(did, path)["value"]


def values(path=PATH):
    return {k: v["value"] for k, v in _load(path).items()}


def status(did, path=PATH):
    return entry(did, path)["status"]


def log_hash(path=PATH):
    with open(path, "rb") as f:
        return hashlib.md5(f.read().replace(b"\r\n", b"\n")).hexdigest()


def require(ids, task, path=PATH):
    """Stop a task on real pilot data (stage2/stage3) while any decision it uses is still proposed.
    Test (synthetic) and smoke (code test) modes run with proposed values and say so."""
    from . import config as C
    from .core.errors import PipelineStop
    proposed = [i for i in ids if status(i, path) != "fixed"]
    if proposed and C.RUN_MODE in ("stage2", "stage3", "main"):
        raise PipelineStop(f"task {task} uses decisions that are still proposed: {', '.join(proposed)}. "
                           f"Review them (python -m pilot decisions) and fix them "
                           f"(python -m pilot decisions fix {' '.join(proposed)} --by <name>) before real data are read.")
    return proposed


def _write(entries_by_id, path=PATH):
    import yaml
    with open(path, encoding="utf-8") as f:
        head = []
        for line in f:
            if line.startswith("#"):
                head.append(line)
            else:
                break
    body = yaml.safe_dump([{k: e[k] for k in FIELDS} for e in entries_by_id.values()], sort_keys=False,
                          allow_unicode=True, width=110)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        f.write("".join(head) + body)
    os.replace(tmp, path)
    _cache.clear()


def _today():
    return _dt.date.today().isoformat()


def fix(ids, by, path=PATH):
    d = dict(_load(path))
    done = []
    for i in ids:
        e = dict(entry(i, path))
        if e["status"] == "fixed":
            continue
        e.update(status="fixed", fixed_by=by, fixed_on=_today())
        d[i] = e
        done.append(i)
    if done:
        _write(d, path)
    return done


def amend(did, new_value, reason, by, path=PATH):
    """Change a FIXED value. The old value, its fixing and the reason for the change go into history."""
    d = dict(_load(path))
    e = dict(entry(did, path))
    if e["status"] != "fixed":
        raise ValueError(f"{did} is still proposed: edit its value in {path} directly, then fix it")
    if not reason or not reason.strip():
        raise ValueError("an amendment needs a reason")
    hist = list(e["history"] or [])
    hist.append(dict(value=e["value"], fixed_by=e["fixed_by"], fixed_on=e["fixed_on"], replaced_on=_today(),
                     replaced_by=by, reason_for_change=reason))
    e.update(value=new_value, history=hist, fixed_by=by, fixed_on=_today())
    d[did] = e
    _write(d, path)
    return e


def table(path=PATH):
    rows = []
    for e in _load(path).values():
        v = json.dumps(e["value"], separators=(",", ":"))
        rows.append(f"{e['id']}  {e['status']:<8}  {e['title'][:42]:<42}  {v[:70]}")
    return "\n".join(rows)
