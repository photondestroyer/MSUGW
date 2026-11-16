"""Resumable downloads of public files (NHDPlus HR, Copernicus DEM tiles, Census counties).

Each download goes to <dest>.part with HTTP Range resume, is checked against the server's
Content-Length, then renamed. A JSON sidecar <dest>.json records URL, size, Last-Modified/ETag and MD5,
and is what provenance cites.
"""
import hashlib
import json
import os
import time

USER_AGENT = "pilot-groundwater-et/1.0 (research pipeline)"


def head(url, timeout=60):
    import requests
    r = requests.head(url, allow_redirects=True, timeout=timeout, headers={"User-Agent": USER_AGENT})
    return r.status_code, r.headers


def download(url, dest, tries=5, timeout=120, log=None):
    import requests
    meta_path = dest + ".json"
    if os.path.exists(dest) and os.path.exists(meta_path):
        meta = json.load(open(meta_path))
        if meta.get("size") == os.path.getsize(dest):
            return dest, meta
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    part = dest + ".part"
    status, hdr = head(url)
    if status == 404:
        raise FileNotFoundError(url)
    total = int(hdr.get("Content-Length", 0)) or None
    for k in range(tries):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if total and have >= total:
            break
        h = {"User-Agent": USER_AGENT}
        if have:
            h["Range"] = f"bytes={have}-"
        try:
            with requests.get(url, headers=h, stream=True, timeout=timeout) as r:
                if r.status_code == 200 and have:
                    have = 0                                  # server ignored Range: start again
                    open(part, "wb").close()
                r.raise_for_status()
                with open(part, "ab") as f:
                    for chunk in r.iter_content(8 * 1024 * 1024):
                        f.write(chunk)
        except Exception as e:                                # network hiccup: retry from where we are
            if log:
                log.warning("download %s: %s (try %d)", url, e, k + 1)
            time.sleep(5 * (k + 1))
            continue
        if not total or os.path.getsize(part) >= total:
            break
    size = os.path.getsize(part)
    if total and size != total:
        raise IOError(f"download of {url} incomplete: {size} of {total} bytes")
    h = hashlib.md5()
    with open(part, "rb") as f:
        for chunk in iter(lambda: f.read(16 * 1024 * 1024), b""):
            h.update(chunk)
    os.replace(part, dest)
    meta = dict(url=url, size=size, md5=h.hexdigest(), last_modified=hdr.get("Last-Modified"),
                etag=hdr.get("ETag"), downloaded=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=1)
    if log:
        log.info("downloaded %s (%.1f MB)", url, size / 1e6)
    return dest, meta
