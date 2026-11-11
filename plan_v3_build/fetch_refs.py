"""Fetch the 19 references added in this plan from Crossref by DOI and write them, formatted, to refs_new.json.
Also extract the reference list of the revision-2 PDF to refs_v2.json (those entries were verified on 2026-10-20)."""
import html
import json
import os
import re
import subprocess
import sys
import urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
DOIS = """10.1007/BF00329030 10.1111/j.1745-6584.2007.00302.x 10.1198/016214504000000089 10.1214/23-AOS2276
10.1111/rssb.12212 10.1016/S0304-4076(98)00084-0 10.1016/j.wace.2015.07.003 10.1016/j.agwat.2020.106061
10.3133/sir20125177 10.1038/d41586-018-01023-3 10.2307/2981697 10.1257/aeri.20210236 10.1093/restud/rdad018
10.1046/j.1365-2745.2002.00682.x 10.1016/j.rse.2019.02.015 10.1214/16-STS558 10.1177/0049124115585360
10.1016/S0022-1694(01)00504-2 10.3133/sir20175040""".split()


def initials(given):
    parts = re.split(r"[\s\-]+", (given or "").replace(".", " ").strip())
    return " ".join(p[0] + "." for p in parts if p)


def fetch(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/")
    r = subprocess.run(["curl", "-s", "-m", "30", "-A", "reference-check/1.0", url], capture_output=True)
    m = json.loads(r.stdout.decode("utf-8", "replace"))["message"]
    au = m.get("author", [])
    names = [f"{a.get('family', a.get('name', '?'))}, {initials(a.get('given'))}".strip().rstrip(",") for a in au]
    if len(names) > 4:
        who = ", ".join(names[:3]) + ", et al."
    elif len(names) > 1:
        who = ", ".join(names[:-1]) + ", & " + names[-1]
    else:
        who = names[0]
    year = (m.get("published-print") or m.get("published") or m.get("issued"))["date-parts"][0][0]
    title = html.unescape(re.sub(r"\s+", " ", (m.get("title") or [""])[0])).replace("Consquences", "Consequences")
    journal = html.unescape((m.get("container-title") or [""])[0])
    vol, page = m.get("volume"), m.get("page") or m.get("article-number")
    tail = journal + (f", {vol}" if vol else "") + (f", {page}" if page and m.get("type") == "journal-article" else "")
    if m.get("type") != "journal-article":
        tail = "U.S. Geological Survey " + journal + " " + doi.split("sir")[-1][:4] + "-" + doi.split("sir")[-1][4:]
    return dict(doi=doi, first=au[0].get("family"), year=year, n_auth=len(au), type=m.get("type"),
                text=f"{who} ({year}). {title}. {tail}. doi:{doi}")


def v2_refs():
    pdf = os.path.join(os.path.dirname(HERE), "Groundwater_Buffering_Research_Plan_v2_Ogallala.pdf")
    txt = subprocess.run(["pdftotext", "-enc", "UTF-8", pdf, "-"], capture_output=True).stdout.decode("utf-8", "replace")
    txt = txt.replace("\r", "").replace("\x0c", "\n")
    body = txt[re.search(r"\n\s*References\s*\n", txt).end():]
    lines = [l.strip() for l in body.splitlines()]
    lines = [l for l in lines if l and not l.startswith("Groundwater buffering in the High Plains aquifer") and not re.fullmatch(r"Page \d+", l)]
    entries, cur = [], ""
    start = re.compile(r"^[A-Z][^()]{2,140}\((19|20)\d\d[a-z]?\)\.")
    for l in lines:
        if start.match(l) and cur:
            entries.append(cur); cur = l
        else:
            cur = (cur + " " + l).strip()
    entries.append(cur)
    out = []
    for e in entries:
        m = re.match(r"^([^,]+),.*?\(((?:19|20)\d\d[a-z]?)\)\.", e)
        if m:
            e = re.sub(r"\s+", " ", e).replace("https://doi.org/", "doi:").replace("- ", "-") if "doi.org" in e else e
            out.append(dict(first=m.group(1).strip(), year=m.group(2), text=e))
    return out


if __name__ == "__main__":
    new = [fetch(d) for d in DOIS]
    json.dump(new, open(f"{HERE}/refs_new.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    old = v2_refs()
    json.dump(old, open(f"{HERE}/refs_v2.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    sys.stdout.reconfigure(encoding="utf-8")
    print(len(new), "new references;", len(old), "entries parsed from revision 2")
    for r in new:
        print(" ", r["text"][:170])
