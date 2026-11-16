"""Stage runner: dependency order, signatures, resume, lock, time budget, status.

A stage is declared in pilot/pipeline.py with
    deps     stages whose outputs it reads
    keys     config settings it depends on (part of its signature)
    modules  package files whose code it runs (part of its signature when CODE_IN_SIGNATURE)
    func     "module:function", called with a Context
    gate     (optional) name of the gate the stage decides, written by it as summary["gate"] =
             dict(name, passes, reasons) in _complete.json (plan v3: G1a by verify_p11, G1b by verify_p12)
    gated_by (optional) stages whose gate must pass before this stage runs. In `focus` and `full` mode a failed
             gate halts the stage (PipelineStop) unless the gate's name is in config.GATE_OVERRIDE; in `test`
             mode it is only logged. Gate states, overrides included, are written into the stage's provenance.
Its output folder is <out_root>/<stage>/<signature>/. A stage is complete when that folder holds
_complete.json. Units of work inside a stage (a tile, a tile-year, a fold) are resumed by the stage
itself: it skips a unit whose output file is complete for the same signature (ncio.is_complete).
"""
import importlib
import json
import os
import socket
import time
import traceback

from .. import config as C
from . import provenance as P
from . import log as L
from .drive import write_text
from .errors import PipelineStop


class Context:
    def __init__(self, runner, stage):
        self.runner, self.stage = runner, stage
        self.name = stage["name"]
        self.sig = runner.sig(self.name)
        self.out = f"{C.out_root()}/{self.name}/{self.sig}"
        os.makedirs(self.out, exist_ok=True)
        self.log = L.get(self.name)
        self.ledger = runner.ledger
        self.t0 = runner.t0
        self.upstream = {d: runner.sig(d) for d in stage["deps"]}
        self.scratch = os.path.join(C.scratch_root(), "work", self.name)
        os.makedirs(self.scratch, exist_ok=True)
        self.gates = runner.gate_states(stage)

    def dep(self, name):
        """Output folder of an upstream stage (raises if it is not complete)."""
        if name not in self.stage["deps"]:
            raise KeyError(f"stage {self.name} does not declare {name} as a dependency")
        d = f"{C.out_root()}/{name}/{self.runner.sig(name)}"
        if not os.path.exists(f"{d}/_complete.json"):
            raise PipelineStop(f"{self.name} needs stage {name}, which is not complete ({d})")
        return d

    def check_budget(self):
        if (time.time() - self.t0) / 3600 > C.SESSION_BUDGET_H:
            raise PipelineStop(f"time budget of {C.SESSION_BUDGET_H} h used up; rerun to continue")

    def record(self, inputs, extra=None):
        extra = dict(extra or {})
        if self.gates:
            extra["gates"] = self.gates
        return P.record(self.name, self.sig, inputs, self.stage["keys"], self.stage["modules"], self.upstream, extra)

    def unit(self, unit, status, **kw):
        self.ledger.write(stage=self.name, sig=self.sig, unit=unit, status=status, **kw)

    def complete(self, summary):
        write_text(f"{self.out}/_complete.json",
                   json.dumps(dict(provenance=self.record({}), summary=summary), indent=1, default=str))


class Runner:
    def __init__(self, stages):
        C.check()
        self.stages = {s["name"]: s for s in stages}
        self.order = [s["name"] for s in stages]
        self._sig = {}
        self.t0 = time.time()
        self.ledger = P.Ledger(f"{C.out_root()}/ledger.jsonl")
        self.log = L.get()

    # ------------------------------------------------------------ signatures
    def sig(self, name):
        if name not in self._sig:
            s = self.stages[name]
            up = {d: self.sig(d) for d in s["deps"]}
            self._sig[name] = P.signature(name, s["keys"], s["modules"], up)
        return self._sig[name]

    def is_complete(self, name):
        return os.path.exists(f"{C.out_root()}/{name}/{self.sig(name)}/_complete.json")

    # ------------------------------------------------------------ gates
    def gate_states(self, stage):
        """{gate name: dict(stage, passes, reasons, overridden)} of the stages that gate `stage`."""
        out = {}
        for up in stage.get("gated_by", []):
            p = f"{C.out_root()}/{up}/{self.sig(up)}/_complete.json"
            gate = None
            if os.path.exists(p):
                try:
                    gate = json.load(open(p)).get("summary", {}).get("gate")
                except (OSError, json.JSONDecodeError):
                    gate = None
            name = self.stages[up].get("gate", up)
            if gate is None:
                out[name] = dict(stage=up, passes=None, reasons=["gate not written (stage incomplete or older code)"],
                                 overridden=name in C.GATE_OVERRIDE)
            else:
                out[name] = dict(stage=up, passes=bool(gate.get("passes")), reasons=gate.get("reasons", []),
                                 overridden=name in C.GATE_OVERRIDE)
        return out

    def check_gates(self, stage):
        for name, g in self.gate_states(stage).items():
            if g["passes"]:
                continue
            why = "; ".join(g["reasons"]) or "failed"
            if C.RUN_MODE == "test":
                self.log.warning("gate %s (%s) does not pass - test mode, %s runs anyway: %s", name, g["stage"],
                                 stage["name"], why)
            elif g["overridden"]:
                self.log.warning("gate %s (%s) does not pass - overridden by GATE_OVERRIDE, %s runs: %s", name,
                                 g["stage"], stage["name"], why)
            else:
                raise PipelineStop(f"gate {name} ({g['stage']}) does not pass, so {stage['name']} does not run: {why}. "
                                   f"Read the {g['stage']} report; to continue anyway add {name!r} to GATE_OVERRIDE "
                                   f"(recorded in the provenance of every later output).")

    def closure(self, names):
        """Stages needed for `names`, in pipeline order."""
        need, todo = set(), list(names)
        while todo:
            n = todo.pop()
            if n not in self.stages:
                raise KeyError(f"unknown stage {n!r}; stages: {self.order}")
            if n not in need:
                need.add(n)
                todo.extend(self.stages[n]["deps"])
        return [n for n in self.order if n in need]

    # ------------------------------------------------------------ lock
    def _lock_path(self):
        return f"{C.out_root()}/.lock"

    def acquire(self):
        p = self._lock_path()
        me = dict(host=socket.gethostname(), pid=os.getpid(), time=time.time())
        if os.path.exists(p):
            try:
                other = json.load(open(p))
            except (OSError, json.JSONDecodeError):
                other = None
            if other and (other.get("host"), other.get("pid")) != (me["host"], me["pid"]) \
                    and (time.time() - other.get("time", 0)) / 3600 < C.LOCK_STALE_H:
                raise PipelineStop(f"another session holds the lock {p} ({other}); delete it if that session is gone")
        write_text(p, json.dumps(me))

    def release(self):
        try:
            os.remove(self._lock_path())
        except OSError:
            pass

    # ------------------------------------------------------------ run
    def run(self, names, force=False):
        todo = self.closure(names)
        self.acquire()
        try:
            for n in todo:
                if self.is_complete(n) and not (force and n in names):
                    self.log.info("stage %-12s %s  complete - skipped", n, self.sig(n))
                    continue
                s = self.stages[n]
                self.check_gates(s)
                from .. import decisions as D           # pilot: a task on real data needs its decisions fixed
                proposed = D.require(s.get("decisions", []), n)
                if proposed:
                    self.log.warning("task %s runs with PROPOSED decisions %s (%s mode: a code test, not a result)",
                                     n, ",".join(proposed), C.RUN_MODE)
                ctx = Context(self, s)
                mod, fn = s["func"].split(":")
                func = getattr(importlib.import_module(f"{__package__.rsplit(chr(46), 1)[0]}.{mod}"), fn)
                self.log.info("stage %-12s %s  running -> %s", n, ctx.sig, ctx.out)
                ctx.unit("_stage", "start")
                t = time.time()
                try:
                    func(ctx)
                except PipelineStop as e:
                    ctx.unit("_stage", "stopped", note=str(e))
                    self.log.warning("stage %s stopped: %s", n, e)
                    raise
                except Exception as e:
                    ctx.unit("_stage", "error", note=f"{type(e).__name__}: {e}",
                             trace=traceback.format_exc()[-2000:])
                    raise
                ctx.unit("_stage", "done", seconds=round(time.time() - t, 1))
                if not self.is_complete(n):
                    raise RuntimeError(f"stage {n} returned without calling ctx.complete()")
        finally:
            self.release()

    # ------------------------------------------------------------ status
    def status(self):
        events = self.ledger.read()
        rows = []
        for n in self.order:
            sig = self.sig(n)
            ev = [e for e in events if e.get("stage") == n and e.get("sig") == sig]
            units_done = {e["unit"] for e in ev if e.get("status") == "done" and e.get("unit") != "_stage"}
            errors = [e for e in ev if e.get("status") in ("error", "stopped")]
            summary = None
            p = f"{C.out_root()}/{n}/{sig}/_complete.json"
            if os.path.exists(p):
                try:
                    summary = json.load(open(p)).get("summary")
                except (OSError, json.JSONDecodeError):
                    pass
            rows.append(dict(stage=n, signature=sig, complete=os.path.exists(p), units_done=len(units_done),
                             last_error=(errors[-1].get("note") if errors else None), summary=summary,
                             other_signatures=sorted(d for d in _listdir(f"{C.out_root()}/{n}") if d != sig)))
        return rows


def _listdir(d):
    try:
        return [x for x in os.listdir(d) if os.path.isdir(os.path.join(d, x))]
    except OSError:
        return []
