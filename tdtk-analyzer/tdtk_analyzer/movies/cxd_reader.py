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
                     created_unix=ci.created_unix if ci.last_field_date is not None else None)
    info.sources = {"time_interval": "file time stamps"}
    if ci.last_field_date is not None:
        info.sources["created"] = "file"
    try:
        info.pixel_size = ci.resolution          # includes the R script's calibration rules
        info.sources["pixel_size"] = ("assumed: file not calibrated (factor = 1), 0.65 µm × magnification "
                                      "used as in the R script") if ci.scale_factor == 1 else "file"
    except CxdError:
        pass                                     # reported as missing pixel size
    info.set_timing([ci.time_from_start[k] for k in sorted(ci.time_from_start)])
    info.extra = {"cxd": ci}
    return [info]


def read(info: MovieInfo, opts: ImportOptions, max_frames: int | None = None) -> np.ndarray:
    ci = info.extra.get("cxd") or read_cxd_info(info.path)
    return read_cxd_frames(ci, n_frames=max_frames)
