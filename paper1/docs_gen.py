"""Generate docs/CONFIG.md from config.py (every setting, its default, its comment, the stages it changes)."""
import ast
import os

from . import config as C

HERE = os.path.dirname(os.path.abspath(__file__))


def _settings():
    src = open(os.path.join(HERE, "config.py"), encoding="utf-8").read()
    lines = src.splitlines()
    tree = ast.parse(src)
    out, section = [], "General"
    sec_lines = {i + 1: l.strip("# =").strip() for i, l in enumerate(lines) if l.startswith("# ====")}
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name) and node.targets[0].id.isupper():
            for ln in sorted(sec_lines):
                if ln < node.lineno:
                    section = sec_lines[ln]
            name = node.targets[0].id
            seg = ast.get_source_segment(src, node.value) or ""
            comment = ""
            first = lines[node.lineno - 1]
            if "#" in first.split("=", 1)[1]:
                comment = first.split("#", 1)[1].strip() if first.count("#") else ""
            default = " ".join(seg.split())
            if len(default) > 70:
                default = default[:67] + "..."
            out.append((section, name, default, comment))
    return out


def render_config_md():
    from .pipeline import STAGES
    used = {}
    for s in STAGES:
        for k in s["keys"]:
            used.setdefault(k, []).append(s["name"])
    L = ["# paper1 configuration reference", "",
         "Generated from `config.py` by `python -m paper1 docs` - do not edit by hand (a test checks it is current).", "",
         "Changing a setting listed under *stages* gives those stages (and every stage after them) a new signature,",
         "so they recompute into a new folder; outputs made with the old value stay where they are.", ""]
    cur = None
    for section, name, default, comment in _settings():
        if section != cur:
            L += ["", f"## {section}", "", "| setting | default | meaning | stages |", "|---|---|---|---|"]
            cur = section
        d = default.replace("|", "\\|")
        c = comment.replace("|", "\\|")
        L.append(f"| `{name}` | `{d}` | {c} | {', '.join(used.get(name, [])) or '-'} |")
    return "\n".join(L) + "\n"


def write_config_md():
    p = os.path.join(HERE, "docs", "CONFIG.md")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write(render_config_md())
    return p
