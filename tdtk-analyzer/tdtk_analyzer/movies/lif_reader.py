"""Leica .lif (package 'readlif'). Every image in the project is a separate movie."""

from __future__ import annotations

import os

import numpy as np

from .base import ImportOptions, MovieInfo, output_name


def probe(path: str, rel: str) -> list[MovieInfo]:
    from readlif.reader import LifFile

    lif = LifFile(path)
    imgs = list(lif.get_iter_image())
    out = []
    for i, img in enumerate(imgs):
        x, y, z, t, _m = img.dims
        t_ax = "T" if t > 1 else ("Z" if z > 1 else "T")
        sc = img.scale  # (px/um x, px/um y, px/um z, frames/s)
        info = MovieInfo(path=path, rel=rel, format="Leica LIF", reader="lif", series=i, n_series=len(imgs),
                         series_name=img.name, size_t=t if t_ax == "T" else z, size_y=y, size_x=x,
                         size_c=img.channels, size_z=z if t_ax == "T" else 1,
                         dtype="uint16" if max(img.bit_depth) > 8 else "uint8", bits=max(img.bit_depth),
                         time_interval=(1.0 / sc[3]) if sc[3] else None,
                         pixel_size=(1.0 / sc[0]) if sc[0] else None,
                         created_unix=os.path.getmtime(path), extra={"t_axis": t_ax})
        info.sources = {"created": "file modification time"}
        if info.time_interval:
            info.sources["time_interval"] = "file"
        if info.pixel_size:
            info.sources["pixel_size"] = "file"
        if t_ax == "Z":
            info.notes.append("no time axis - the Z axis is used as time")
        if y == 1:
            info.notes.append("line scan (x-t) - not a movie")
        info.name = output_name(os.path.basename(path), img.name, len(imgs), i)
        out.append(info)
    return out


def read(info: MovieInfo, opts: ImportOptions, max_frames: int | None = None) -> np.ndarray:
    from readlif.reader import LifFile

    img = LifFile(info.path).get_image(info.series)
    ch = info.channel_index(opts)
    n = info.size_t if not max_frames else min(info.size_t, max_frames)
    zsel = str(opts.z_plane).lower()
    frames = []
    for t in range(n):
        if info.extra["t_axis"] == "Z":
            frames.append(np.asarray(img.get_frame(z=t, t=0, c=ch)))
        elif zsel == "max" and info.size_z > 1:
            frames.append(np.max([np.asarray(img.get_frame(z=z, t=t, c=ch)) for z in range(info.size_z)], axis=0))
        else:
            frames.append(np.asarray(img.get_frame(z=int(zsel or 0), t=t, c=ch)))
    return np.stack(frames)
