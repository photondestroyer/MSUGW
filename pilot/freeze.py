"""The freeze (pilot Section 11: "With the supervisor, freeze the sample, exclusions, model formulas and thresholds.
The main analysis then begins").

    python -m pilot freeze --by "<name>" --supervisor "<name>" [--deviation "<reason>"]

refuses unless every decision is fixed and pilot/ has no uncommitted change; then writes pilot/freeze.json (who,
when, the supervisor, the SHA-256 of decisions.yaml, the commit of the frozen code), commits it and tags the commit
pilot-freeze-<date>. check() is what the main analysis calls: no record, or a decision log that differs from the
frozen one, stops it. A change after the freeze is a deviation: freeze again with --deviation, which keeps the
earlier record in "history" and is reported with the results.
"""
import datetime as _dt
import hashlib
import json
import os
import subprocess

from .core.errors import PipelineStop

PKG = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(PKG, "freeze.json")
REPO = os.path.dirname(PKG)


def decisions_hash(path=None):
    """SHA-256 of the decision log, line endings normalised (git may check it out with CRLF)."""
    from . import decisions as D
    with open(path or D.PATH, "rb") as f:
        return hashlib.sha256(f.read().replace(b"\r\n", b"\n")).hexdigest()


def record(path=PATH):
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def check(path=PATH, decisions_path=None):
    """The freeze record, if it exists and the decision log is the frozen one; otherwise PipelineStop."""
    rec = record(path)
    if rec is None:
        raise PipelineStop("no freeze record (pilot/freeze.json): the main analysis runs only after "
                           "`python -m pilot freeze --by ... --supervisor ...`")
    if decisions_hash(decisions_path) != rec["decisions_sha256"]:
        raise PipelineStop("decisions.yaml differs from the frozen one: a change after the freeze is a deviation; "
                           "record it with `python -m pilot freeze --deviation \"<reason>\" ...`")
    return rec


def make_record(by, supervisor, commit, previous=None, deviation=None, today=None, decisions_path=None):
    """The freeze record (pure: no git, no file)."""
    from . import decisions as D
    entries = D._load(decisions_path or D.PATH)
    proposed = sorted(k for k, e in entries.items() if e["status"] != "fixed")
    if proposed:
        raise PipelineStop(f"cannot freeze while decisions are proposed: {proposed}")
    if not by or not supervisor:
        raise PipelineStop("a freeze needs --by and --supervisor")
    if previous is not None and not deviation:
        raise PipelineStop("already frozen: a new freeze is a deviation and needs --deviation \"<reason>\"")
    day = today or _dt.date.today().isoformat()
    n = len((previous or {}).get("history", [])) + (1 if previous else 0)
    rec = dict(frozen_on=day, by=by, supervisor=supervisor, decisions=sorted(entries),
               decisions_sha256=decisions_hash(decisions_path), code_commit=commit,
               tag=f"pilot-freeze-{day}" + (f"-dev{n}" if n else ""), deviation=deviation,
               history=((previous or {}).get("history", []) + [{k: v for k, v in previous.items() if k != "history"}])
               if previous else [])
    return rec


def _git(*args):
    p = subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True)
    if p.returncode != 0:
        raise PipelineStop(f"git {' '.join(args)} failed: {p.stderr.strip()}")
    return p.stdout.strip()


def freeze(by, supervisor, deviation=None):
    """Write the record, commit it and tag the commit. Run by the user from the command line."""
    dirty = [l for l in _git("status", "--porcelain", "--", "pilot").splitlines()
             if not l.endswith("freeze.json") and "/results/" not in l and "__pycache__" not in l]
    if dirty:
        raise PipelineStop("commit the pilot code first; uncommitted:\n  " + "\n  ".join(dirty[:20]))
    # a record git does not track is an attempt that stopped before its commit, not an earlier freeze
    tracked = subprocess.run(["git", "ls-files", "--error-unmatch", PATH], cwd=REPO, capture_output=True).returncode == 0
    rec = make_record(by, supervisor, _git("rev-parse", "HEAD"), previous=record() if tracked else None,
                      deviation=deviation)
    with open(PATH, "w", encoding="utf-8") as f:
        json.dump(rec, f, indent=1)
    _git("add", "-f", PATH)                           # -f: the repository ignores **/*.json; this one must be tracked
    _git("commit", "-m", f"Pilot freeze {rec['frozen_on']} (by {by}, supervisor {supervisor})"
         + (f"\n\nDeviation: {deviation}" if deviation else ""))
    _git("tag", "-a", rec["tag"], "-m", f"Pilot freeze: decisions {rec['decisions_sha256'][:12]}")
    return rec
