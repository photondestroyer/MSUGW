"""Build the Papers 1 and 2 detailed plan PDF.
Run:  python figs.py ; python fetch_refs.py (needs network, once) ; python build_pdf.py
Numbers come from checks/*.json, which checks/check_*.py write from the files in merged_datasets."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine as E                                                     # noqa: E402
from content_a import part_a, paper1                                   # noqa: E402
from content_b import paper2, back                                     # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "Groundwater_Buffering_Papers_1_2_Detailed_Plan_v3_Ogallala.pdf")
# rain by aridity decile (mm per year, 1985-2014), printed by checks/check_support.py
P_DEC = [411, 439, 443, 449, 472, 511, 519, 568, 635, 706]
# revision-2 references cited in this document (first author, year as printed there)
V2_CITED = [("Abatzoglou", "2013"), ("Butler", "2016"), ("Chernozhukov", "2018"), ("Deines", "2019a"), ("Deines", "2019b"),
            ("Drysdale", "2018"), ("Fan", "2017"), ("Jasechko", "2021"), ("Lipsitch", "2010"), ("Ma", "2026"),
            ("Melton", "2022"), ("Meyer", "2021"), ("Rohde", "2024b"), ("Scanlon", "2012"), ("VanderWeele", "2017"),
            ("Volk", "2024"), ("Whittemore", "2016"), ("Xie", "2021")]


def pivot_geometry(radius=400.0, n=16, px=30.0, steps=161):
    """Share of random block offsets for which a 480 m block is wholly (or at least 80 %) inside a pivot circle."""
    c = (np.arange(n) - (n - 1) / 2) * px
    off = np.linspace(-n * px / 2, n * px / 2, steps, endpoint=False) + n * px / steps / 2
    inside_all = inside80 = 0
    for dx in off:
        for dy in off:
            frac = (((c[None, :] + dx) ** 2 + (c[:, None] + dy) ** 2) <= radius ** 2).mean()
            inside_all += frac == 1.0
            inside80 += frac >= 0.8
    return dict(p_inside=100 * inside_all / steps ** 2, p80=100 * inside80 / steps ** 2)


def load():
    J = lambda n: json.load(open(f"{HERE}/checks/{n}.json", encoding="utf-8"))
    sup, mi, ww = J("support"), J("misc"), J("wtd_wells")
    pan = mi["two_drought_panel"]
    return dict(S=sup["rows"], SS=sup["summary"], MI=mi, WW=ww, geom=pivot_geometry(), figs=f"{HERE}/fig", P_dec=P_DEC,
                n_both=sum(v["n"] for v in pan.values()))


def references():
    new = json.load(open(f"{HERE}/refs_new.json", encoding="utf-8"))
    old = json.load(open(f"{HERE}/refs_v2.json", encoding="utf-8"))
    items = [(r["first"], str(r["year"]), "* " + r["text"]) for r in new]
    for first, year in V2_CITED:
        hit = [r for r in old if r["first"] == first and r["year"] == year]
        assert len(hit) == 1, (first, year, len(hit))
        items.append((first, year, hit[0]["text"]))
    items.sort(key=lambda t: (t[0].lower(), t[1]))
    fix = lambda t: t.replace("‐", "-").replace("‑", "-").replace("–", "-")
    out = E.H1("References")
    out.append(E.P("Entries marked * are added in this document and were checked against Crossref by DOI on 10 November "
                   "2026 (title, authors and year agree). The others are cited from revision 2, where they were checked on "
                   "20 October 2026; revision 2 holds the full list of 144 works.", "meta"))
    out += [E.P(fix(t), "ref") for _, _, t in items]
    return out, len(new), len(V2_CITED)


def main():
    D = load()
    story = part_a(D) + paper1(D) + paper2(D) + back(D)
    refs, n_new, n_old = references()
    story += refs
    missing = E.check_glyphs()
    if missing:
        raise SystemExit("characters missing from the font: " + " ".join(f"U+{ord(c):04X}({c})" for c in missing))
    pages = E.build(OUT, story, "Groundwater buffering in the High Plains aquifer: Papers 1 and 2, detailed analysis plan",
                    "Groundwater buffering in the High Plains aquifer  ·  Papers 1 and 2, detailed analysis plan  ·  "
                    "10 November 2026")
    print(f"written {OUT}\n  pages {pages}  references {n_new} new + {n_old} from revision 2  "
          f"pivot geometry {D['geom']}  size {os.path.getsize(OUT) / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
