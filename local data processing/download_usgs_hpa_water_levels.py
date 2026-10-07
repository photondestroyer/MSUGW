"""
Download every depth-to-water field measurement USGS holds for wells of the High Plains aquifer.

Source: USGS Water Data API (OGC API), https://api.waterdata.usgs.gov/ogcapi/v0
  monitoring-locations  national_aquifer_code = N100HGHPLN, site_type_code = GW   -> wells
  field-measurements    national_aquifer_code = N100HGHPLN, parameter_code = 72019 -> depth to water below land
                                                                                     surface, ft
The older service waterservices.usgs.gov/nwis/gwlevels was retired in late 2025.

Raw pages are stored unchanged (gzip JSON) in RAW so the run can be resumed and re-read; nothing is filtered here.
"""
import gzip
import json
import sys
import time
from pathlib import Path

import requests

API = "https://api.waterdata.usgs.gov/ogcapi/v0/collections"
RAW = Path(r"G:/MSU_GWB/_usgs_wl_raw")
RAW.mkdir(parents=True, exist_ok=True)
AQUIFER = "N100HGHPLN"
PARAM = "72019"
LIMIT = 10000
FIRST_YEAR, LAST_YEAR = 1900, 2026
S = requests.Session()
S.headers["User-Agent"] = "msu-gwb-research/1.0 (High Plains groundwater study)"


def log(m):
    print(time.strftime("%H:%M:%S"), m, flush=True)


def get(url, params=None):
    last = None
    for attempt in range(1, 9):
        try:
            r = S.get(url, params=params, timeout=(30, 300))
            if r.status_code == 429 or r.status_code >= 500:
                wait = int(r.headers.get("Retry-After", 30 * attempt))
                log(f"  HTTP {r.status_code}, waiting {wait} s")
                time.sleep(min(wait, 900))
                continue
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last = e
            time.sleep(15 * attempt)
    raise RuntimeError(f"gave up on {url}: {last}")


def pages(collection, params, stem):
    """Fetch a query page by page (following the server's 'next' link). Finished queries are skipped."""
    done = RAW / f"{stem}.done"
    if done.exists():
        return int(done.read_text())
    for old in RAW.glob(f"{stem}_p*.json.gz"):
        old.unlink()
    url, p, n, k = f"{API}/{collection}/items", dict(params, f="json", limit=LIMIT), 0, 0
    while url:
        j = get(url, p)
        feats = j.get("features", [])
        with gzip.open(RAW / f"{stem}_p{k:03d}.json.gz", "wt", encoding="utf-8") as f:
            json.dump(feats, f)
        n += len(feats)
        k += 1
        url = next((l["href"] for l in j.get("links", []) if l.get("rel") == "next"), None)
        p = None
        time.sleep(0.4)
    done.write_text(str(n))
    return n


if __name__ == "__main__":
    n = pages("monitoring-locations", dict(national_aquifer_code=AQUIFER, site_type_code="GW"), "wells")
    log(f"wells: {n}")
    total = pages("field-measurements", dict(national_aquifer_code=AQUIFER, parameter_code=PARAM,
                                             time=f"../{FIRST_YEAR - 1}-12-31T23:59:59Z", skipGeometry="true"),
                  "fm_pre1900")
    log(f"before {FIRST_YEAR}: {total}")
    for y in range(FIRST_YEAR, LAST_YEAR + 1):
        k = pages("field-measurements", dict(national_aquifer_code=AQUIFER, parameter_code=PARAM,
                                             time=f"{y}-01-01T00:00:00Z/{y}-12-31T23:59:59Z", skipGeometry="true"),
                  f"fm_{y}")
        total += k
        log(f"{y}: {k} measurements (total {total})")
    log(f"DONE: {total} measurements")
    sys.exit(0)
