"""docs/pilot/DEFINITIONS.md and docs/pilot/DECISIONS.md, generated from definitions.TERMS and decisions.yaml.
A test fails when the files on disk differ from what this module writes (python -m pilot docs)."""
import json
import os

from . import config as C
from . import decisions as D
from . import definitions as DEF

DOCS = os.path.join(C.REPO_ROOT, "docs", "pilot")


def definitions_md():
    L = ["# Pilot definitions (written label rules)", "",
         "Generated from `pilot/definitions.py` and `pilot/decisions.yaml` by `python -m pilot docs`. Each term is "
         "one function; every task calls it, none re-implements it.", ""]
    for term, wording, funcs, dids, task in DEF.TERMS:
        L += [f"## {term}", "", f"> {wording}", "", f"- Built in: Task {task}",
              f"- Functions: {', '.join(f'`definitions.{f}()`' for f in funcs)}"]
        for did in dids:
            e = D.entry(did)
            L.append(f"- {did} ({e['status']}) {e['title']}: `{json.dumps(e['value'])}`")
        for f in funcs:
            doc = (getattr(DEF, f).__doc__ or "").strip().split("\n")
            L.append(f"  - `{f}`: " + " ".join(x.strip() for x in doc))
        L.append("")
    return "\n".join(L)


def decisions_md():
    L = ["# Pilot decision log", "", "Generated from `pilot/decisions.yaml` (`python -m pilot docs`).", "",
         "| ID | status | decision | value | tasks |", "|---|---|---|---|---|"]
    for e in D.entries().values():
        L.append(f"| {e['id']} | {e['status']} | {e['title']} | `{json.dumps(e['value'])}` | "
                 f"{', '.join(map(str, e['tasks']))} |")
    L.append("")
    for e in D.entries().values():
        L += [f"## {e['id']} {e['title']}", "", e["reason"], "",
              f"Status: {e['status']} (proposed by {e['proposed_by']} on {e['proposed_on']}"
              + (f"; fixed by {e['fixed_by']} on {e['fixed_on']}" if e["status"] == "fixed" else "") + ")", ""]
        for h in e["history"] or []:
            L.append(f"- replaced {h['replaced_on']}: was `{json.dumps(h['value'])}` ({h['reason_for_change']})")
    return "\n".join(L)


def write_all():
    os.makedirs(DOCS, exist_ok=True)
    out = []
    for name, text in (("DEFINITIONS.md", definitions_md()), ("DECISIONS.md", decisions_md())):
        p = os.path.join(DOCS, name)
        with open(p, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
        out.append(p)
    return out


def current():
    """True when the docs on disk equal what write_all() would write."""
    for name, text in (("DEFINITIONS.md", definitions_md()), ("DECISIONS.md", decisions_md())):
        p = os.path.join(DOCS, name)
        if not os.path.exists(p) or open(p, encoding="utf-8").read() != text:
            return False
    return True
