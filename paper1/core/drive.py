"""Google Drive helpers (Colab) that degrade to plain file operations locally.

deliver(tmp, dst)   move a finished file into place. On Colab: copy to dst + '.tmp' on Drive, fsync,
                    wait, check size and MD5 (up to 4 tries), rename - the robust_drive_copy pattern of
                    the dataset notebooks. Locally: os.replace.
drive_md5(path)     the MD5 Google keeps for a Drive file (Drive API, no file read). Colab only; None
                    elsewhere or when the API is not reachable.
"""
import hashlib
import os
import shutil
import time

from .. import config as C
from .errors import DataIntegrityError

DRIVE_SETTLE_S = 10
TRIES = 4


def in_colab():
    try:
        import google.colab  # noqa: F401
        return True
    except ImportError:
        return False


def mount():
    """Mount Drive in Colab (no-op when already mounted or not in Colab)."""
    if C.RUN_ENV == "colab" and not os.path.ismount("/content/drive") and in_colab():
        from google.colab import drive
        drive.mount("/content/drive")


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(16 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def deliver(tmp, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if not (C.RUN_ENV == "colab" and dst.startswith("/content/drive")):
        os.replace(tmp, dst)
        return dst
    want_size, want_md5 = os.path.getsize(tmp), _md5(tmp)
    part = dst + ".tmp"
    for k in range(TRIES):
        shutil.copyfile(tmp, part)
        with open(part, "rb+") as f:
            os.fsync(f.fileno())
        time.sleep(DRIVE_SETTLE_S * (k + 1) if k else 2)
        if os.path.getsize(part) == want_size and _md5(part) == want_md5:
            os.replace(part, dst)
            os.remove(tmp)
            return dst
    raise DataIntegrityError(f"copy to Drive did not verify after {TRIES} tries: {dst}")


def write_text(dst, text):
    """Small text/JSON outputs (reports, manifests): written in place with fsync."""
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    tmp = dst + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, dst)


_service = None


def drive_md5(path):
    """MD5 that Google Drive stores for a file under /content/drive/MyDrive (None if unavailable)."""
    global _service
    if not (C.RUN_ENV == "colab" and in_colab() and path.startswith(C.COLAB_MYDRIVE)):
        return None
    try:
        if _service is None:
            from google.colab import auth
            from googleapiclient.discovery import build
            auth.authenticate_user()
            _service = build("drive", "v3", cache_discovery=False)
        parent = "root"
        parts = path[len(C.COLAB_MYDRIVE):].strip("/").split("/")
        for i, name in enumerate(parts):
            q = f"name = '{name.replace(chr(39), chr(92) + chr(39))}' and '{parent}' in parents and trashed = false"
            res = _service.files().list(q=q, fields="files(id, md5Checksum, mimeType)", pageSize=10).execute()
            files = res.get("files", [])
            if len(files) != 1:
                return None
            parent = files[0]["id"]
            if i == len(parts) - 1:
                return files[0].get("md5Checksum")
    except Exception:
        return None
    return None
