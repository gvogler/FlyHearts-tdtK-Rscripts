"""Fallback for all other formats through Bio-Formats command line tools (bftools, needs Java).

Only used when 'showinf' and 'bfconvert' are on the PATH. The movie is converted
once to OME-TIFF in '<output>/import_cache' and then read like any OME-TIFF."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys

import numpy as np

from .base import ImportOptions, MovieInfo, output_name
from .tiff_reader import parse_ome

CACHE_DIR = {"path": None}   # set by the scanner (output folder)


def _tool(name: str):
    for cand in (name, name + ".bat") if sys.platform == "win32" else (name,):
        p = shutil.which(cand)
        if p:
            return p
    return None


def available() -> bool:
    return bool(_tool("showinf") and _tool("bfconvert"))


def probe(path: str, rel: str) -> list[MovieInfo]:
    showinf = _tool("showinf")
    if not showinf:
        raise RuntimeError("unsupported format (install Bio-Formats 'bftools' to read it)")
    txt = subprocess.run([showinf, "-nopix", "-omexml-only", "-no-upgrade", path],
                         capture_output=True, text=True, timeout=600).stdout
    start = txt.find("<?xml") if "<?xml" in txt else txt.find("<OME")
    if start < 0:
        raise RuntimeError("Bio-Formats could not read the file")
    images = parse_ome(txt[start:])
    out = []
    for i, im in enumerate(images):
        s = im["sizes"]
        t_ax = "T" if s["T"] > 1 else ("Z" if s["Z"] > 1 else "T")
        info = MovieInfo(path=path, rel=rel, format="Bio-Formats: " + os.path.splitext(path)[1].lstrip("."),
                         reader="bioformats", series=i, n_series=len(images), series_name=im["name"],
                         size_t=s[t_ax], size_y=s["Y"], size_x=s["X"], size_c=s["C"],
                         size_z=s["Z"] if t_ax == "T" else 1, dtype=im["type"] or "uint16",
                         time_interval=im["time_interval"], pixel_size=im["pixel_size"],
                         created_unix=im["created"], channel_names=im["channels"], extra={"t_axis": t_ax})
        info.name = output_name(os.path.basename(path), im["name"], len(images), i)
        out.append(info)
    return out


def read(info: MovieInfo, opts: ImportOptions, max_frames: int | None = None) -> np.ndarray:
    from . import tiff_reader

    cache = CACHE_DIR["path"] or os.path.join(os.path.dirname(info.path), ".tdtk_import_cache")
    os.makedirs(cache, exist_ok=True)
    key = hashlib.sha1(f"{info.path}|{os.path.getmtime(info.path)}|{info.series}".encode()).hexdigest()[:12]
    target = os.path.join(cache, f"{key}.ome.tif")
    if not os.path.exists(target):
        subprocess.run([_tool("bfconvert"), "-overwrite", "-series", str(info.series), "-bigtiff",
                        info.path, target], check=True, capture_output=True, timeout=24 * 3600)
    tinfo = tiff_reader.probe(target, info.rel)[0]
    tinfo.channel, tinfo.channel_names = info.channel, info.channel_names
    return tiff_reader.read(tinfo, opts, max_frames)
