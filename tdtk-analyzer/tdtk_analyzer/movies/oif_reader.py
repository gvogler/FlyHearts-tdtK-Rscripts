"""Olympus FluoView .oib / .oif (package 'oiffile')."""

from __future__ import annotations

import os

import numpy as np

from .base import ImportOptions, MovieInfo, pick_time_axis, select_planes
from .tiff_reader import to_s, to_um


def _axis_params(main: dict):
    """Frame interval (s) and pixel size (um) from the FluoView 'Axis n Parameters Common' sections."""
    interval = pixel = None
    for key, sec in main.items():
        if not key.startswith("Axis") or not isinstance(sec, dict):
            continue
        code = str(sec.get("AxisCode", "")).strip('"')
        try:
            n = int(sec.get("MaxSize", 0))
            start, end = float(sec.get("StartPosition", 0)), float(sec.get("EndPosition", 0))
        except (TypeError, ValueError):
            continue
        unit = str(sec.get("UnitName", "")).strip('"')
        if code == "T" and n > 1 and end > start:
            interval = to_s((end - start) / (n - 1), unit or "ms")
        if code == "X" and n > 1 and end > start:
            pixel = to_um((end - start) / (n - 1), unit or "um")
    return interval, pixel


def probe(path: str, rel: str) -> list[MovieInfo]:
    import oiffile

    with oiffile.OifFile(path) as oif:
        axes, shape = oif.axes.upper(), tuple(oif.shape)
        interval, pixel = _axis_params(oif.mainfile)
        dtype = str(oif.series[0].dtype) if oif.series else "uint16"
    t_ax, note = pick_time_axis(axes, shape)
    size = lambda a: shape[axes.index(a)] if a in axes else 1  # noqa: E731
    info = MovieInfo(path=path, rel=rel, format="Olympus OIB/OIF", reader="oif", name=os.path.basename(path),
                     size_t=size(t_ax), size_y=size("Y"), size_x=size("X"), size_c=size("C"),
                     size_z=size("Z") if t_ax != "Z" else 1, dtype=dtype, time_interval=interval,
                     pixel_size=pixel, created_unix=os.path.getmtime(path),
                     extra={"axes": axes, "t_axis": t_ax})
    info.sources = {"created": "file modification time"}
    if interval:
        info.sources["time_interval"] = "file"
    if pixel:
        info.sources["pixel_size"] = "file"
    if note:
        info.notes.append(note)
    return [info]


def read(info: MovieInfo, opts: ImportOptions, max_frames: int | None = None) -> np.ndarray:
    import oiffile

    with oiffile.OifFile(info.path) as oif:
        arr = oif.asarray()
    movie = select_planes(arr, info.extra["axes"], info.extra["t_axis"], info.channel_index(opts), opts.z_plane)
    return movie[:max_frames] if max_frames else movie
