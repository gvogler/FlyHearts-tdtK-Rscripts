"""Zeiss .czi (package 'pylibCZIrw'). Every scene is a separate movie."""

from __future__ import annotations

import os

import numpy as np

from .base import ImportOptions, MovieInfo, output_name


def _get(d, *keys):
    for k in keys:
        if isinstance(d, dict):
            d = d.get(k)
        elif isinstance(d, list) and isinstance(k, int) and k < len(d):
            d = d[k]
        else:
            return None
    return d


def _as_list(v):
    return v if isinstance(v, list) else ([] if v is None else [v])


def probe(path: str, rel: str) -> list[MovieInfo]:
    from pylibCZIrw import czi as pyczi

    out = []
    with pyczi.open_czi(path) as d:
        bb = d.total_bounding_box
        md = d.metadata or {}
        scenes = sorted(d.scenes_bounding_rectangle) or [None]
        size = {k: (v[1] - v[0]) for k, v in bb.items()}
        t_ax = "T" if size.get("T", 1) > 1 else ("Z" if size.get("Z", 1) > 1 else "T")
        pixel = None
        for item in _as_list(_get(md, "ImageDocument", "Metadata", "Scaling", "Items", "Distance")):
            if isinstance(item, dict) and item.get("@Id") == "X" and item.get("Value"):
                pixel = float(item["Value"]) * 1e6
        interval = None
        tinc = _get(md, "ImageDocument", "Metadata", "Information", "Image", "Dimensions", "T", "Positions",
                    "Interval", "Increment")
        if tinc:
            try:
                interval = float(tinc)
            except (TypeError, ValueError):
                interval = None
        names = []
        for c in _as_list(_get(md, "ImageDocument", "Metadata", "Information", "Image", "Dimensions",
                               "Channels", "Channel")):
            if isinstance(c, dict):
                names.append(c.get("@Name") or c.get("Fluor") or c.get("@Id", ""))
        created = None
        acq = _get(md, "ImageDocument", "Metadata", "Information", "Image", "AcquisitionDateAndTime")
        if acq:
            from .tiff_reader import parse_iso
            created = parse_iso(acq)
        for si, sc in enumerate(scenes):
            rect = d.scenes_bounding_rectangle[sc] if sc is not None else d.total_bounding_rectangle
            info = MovieInfo(path=path, rel=rel, format="Zeiss CZI", reader="czi", series=si,
                             n_series=len(scenes), series_name=f"Scene{sc + 1}" if (sc is not None and len(scenes) > 1) else "",
                             size_t=size.get(t_ax, 1), size_y=rect.h, size_x=rect.w, size_c=size.get("C", 1),
                             size_z=size.get("Z", 1) if t_ax != "Z" else 1,
                             dtype=str(d.get_channel_pixel_type(0)), time_interval=interval, pixel_size=pixel,
                             created_unix=created or os.path.getmtime(path), channel_names=names,
                             extra={"t_axis": t_ax, "scene": sc, "roi": (rect.x, rect.y, rect.w, rect.h)})
            info.sources = {"created": "file" if created else "file modification time"}
            if interval:
                info.sources["time_interval"] = "file"
            if pixel:
                info.sources["pixel_size"] = "file"
            if t_ax == "Z":
                info.notes.append("no time axis - the Z axis is used as time")
            info.name = output_name(os.path.basename(path), info.series_name, len(scenes), si)
            out.append(info)
    return out


def read(info: MovieInfo, opts: ImportOptions, max_frames: int | None = None) -> np.ndarray:
    from pylibCZIrw import czi as pyczi

    ch = info.channel_index(opts)
    t_ax = info.extra["t_axis"]
    n = info.size_t if not max_frames else min(info.size_t, max_frames)
    zsel = str(opts.z_plane).lower()
    frames = []
    with pyczi.open_czi(info.path) as d:
        def plane(t, z):
            p = {"C": ch, "T": t if t_ax == "T" else 0, "Z": z}
            a = d.read(roi=info.extra["roi"], plane=p, scene=info.extra["scene"])
            return a[..., 0] if a.ndim == 3 else a

        for t in range(n):
            if t_ax == "Z":
                frames.append(plane(0, t))
            elif zsel == "max" and info.size_z > 1:
                frames.append(np.max([plane(t, z) for z in range(info.size_z)], axis=0))
            else:
                frames.append(plane(t, int(zsel or 0)))
    return np.stack(frames)
