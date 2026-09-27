"""Hamamatsu HCImage / SimplePCI .cxd (see ../cxd.py)."""

from __future__ import annotations

import os

import numpy as np

from ..cxd import CxdError, read_cxd_frames, read_cxd_info
from .base import ImportOptions, MovieInfo


def probe(path: str, rel: str) -> list[MovieInfo]:
    ci = read_cxd_info(path)
    info = MovieInfo(path=path, rel=rel, format="Hamamatsu CXD", reader="cxd", name=os.path.basename(path),
                     size_t=ci.size_t, size_y=ci.size_y, size_x=ci.size_x, size_c=ci.size_c, size_z=ci.size_z,
                     dtype=ci.pixel_type, bits=ci.bits_per_pixel, time_interval=ci.time_interval,
                     created_unix=ci.created_unix)
    try:
        info.pixel_size = ci.resolution          # includes the R script's calibration rules
    except CxdError:
        info.notes.append("no calibration ('factor') in the file")
    info.extra = {"cxd": ci}
    return [info]


def read(info: MovieInfo, opts: ImportOptions, max_frames: int | None = None) -> np.ndarray:
    ci = info.extra.get("cxd") or read_cxd_info(info.path)
    return read_cxd_frames(ci, n_frames=max_frames)
