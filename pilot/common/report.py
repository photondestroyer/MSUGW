"""Writing outputs: every table as CSV (and parquet for internal tables), every figure as PNG, each with a
sidecar JSON that states units, sample mask, git commit, decision-log hash and the task signature.

Output standards of the pilot (Section 10): "Produce every figure and table with a script that reads the
stored tables ... Save each item as a data file (CSV) and a figure file. State the units and the sample mask
on each. Use the same colours and labels for thin, thick and rainfed cells in every figure."
"""
import json
import os
import subprocess

from .. import config as C
from .. import decisions as D
from ..core.drive import deliver, write_text
from ..core import provenance as P

_git = None


def git_tag():
    """Commit of the repository and whether the working tree has uncommitted changes."""
    global _git
    if _git is None:
        try:
            sha = subprocess.run(["git", "-C", C.REPO_ROOT, "rev-parse", "--short", "HEAD"], capture_output=True,
                                 text=True, timeout=20).stdout.strip()
            dirty = subprocess.run(["git", "-C", C.REPO_ROOT, "status", "--porcelain", "--", "pilot"],
                                   capture_output=True, text=True, timeout=60).stdout.strip()
            _git = dict(commit=sha or "unknown", pilot_code_uncommitted=bool(dirty))
        except (OSError, subprocess.SubprocessError):
            _git = dict(commit="unknown", pilot_code_uncommitted=None)
    return _git


def sidecar(ctx, name, units=None, mask=None, notes=None, extra=None):
    return dict(item=name, task=ctx.name, signature=ctx.sig, run_mode=C.RUN_MODE, created=P.utcnow(),
                units=units or {}, sample_mask=mask or "none", notes=notes or "", git=git_tag(),
                decision_log_md5=D.log_hash(), decisions_proposed=[k for k, e in D.entries().items()
                                                                    if e["status"] != "fixed"],
                **(extra or {}))


def _put(ctx, name, writer):
    tmp = os.path.join(ctx.scratch, name)
    os.makedirs(os.path.dirname(tmp), exist_ok=True)
    writer(tmp)
    return deliver(tmp, f"{ctx.out}/{name}")


def table(ctx, df, name, units=None, mask=None, notes=None, parquet=False, csv=True):
    """Write a table (CSV deliverable and/or parquet for internal use) and its sidecar."""
    paths = {}
    if csv:
        paths["csv"] = _put(ctx, f"{name}.csv", lambda p: df.to_csv(p, index=False))
    if parquet:
        paths["parquet"] = _put(ctx, f"{name}.parquet", lambda p: df.to_parquet(p, index=False))
    meta = sidecar(ctx, name, units, mask, notes, dict(rows=int(len(df)), columns=list(map(str, df.columns))))
    write_text(f"{ctx.out}/{name}.meta.json", json.dumps(meta, indent=1, default=str))
    return paths


def text(ctx, name, body):
    write_text(f"{ctx.out}/{name}", body)
    return f"{ctx.out}/{name}"


def figure(ctx, fig, name, units=None, mask=None, caption=""):
    """Save a matplotlib figure as PNG with units and sample mask written under it, and a sidecar."""
    import matplotlib
    matplotlib.use("Agg")
    foot = " | ".join(x for x in [caption, f"units: {', '.join(f'{k} {v}' for k, v in (units or {}).items())}"
                                  if units else "", f"sample: {mask}" if mask else "",
                                  f"{ctx.name} {ctx.sig} {git_tag()['commit']}"] if x)
    fig.text(0.01, 0.005, foot, fontsize=6, color="#555555", ha="left", va="bottom", wrap=True)
    p = _put(ctx, f"{name}.png", lambda q: fig.savefig(q, dpi=160, bbox_inches="tight"))
    import matplotlib.pyplot as plt
    plt.close(fig)
    write_text(f"{ctx.out}/{name}.meta.json", json.dumps(sidecar(ctx, name, units, mask, caption), indent=1,
                                                         default=str))
    return p


def md(df, digits=3):
    """A DataFrame as a markdown table (no `tabulate` dependency). Floats rounded to `digits` significant places."""
    import numpy as np

    def cell(v):
        if isinstance(v, (float, np.floating)):
            return "" if np.isnan(v) else f"{v:.{digits}g}"
        return str(v).replace("|", "/")
    cols = [str(c) for c in df.columns]
    rows = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    rows += ["| " + " | ".join(cell(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(rows)


def read(stage_dir, name):
    """Read a table written by table(): parquet if present, else CSV."""
    import pandas as pd
    p = f"{stage_dir}/{name}.parquet"
    return pd.read_parquet(p) if os.path.exists(p) else pd.read_csv(f"{stage_dir}/{name}.csv", low_memory=False)
