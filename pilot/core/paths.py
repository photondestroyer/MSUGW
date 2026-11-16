"""Folders of a run: task outputs, shared caches, scratch (pilot version of paper1/core/paths.py)."""
import os

from .. import config as C


def stage_dir(stage, sig):
    return f"{C.out_root()}/{stage}/{sig}"


def shared_dir(*parts):
    d = "/".join([C.shared_root(), *parts])
    os.makedirs(d, exist_ok=True)
    return d


def scratch_dir(*parts):
    d = "/".join([C.scratch_root(), *parts])
    os.makedirs(d, exist_ok=True)
    return d
