"""Nikon NIS-Elements .nd2 (package 'nd2'). Every XY position is a separate movie."""

from __future__ import annotations

import os

import numpy as np

from .base import ImportOptions, MovieInfo, output_name


def probe(path: str, rel: str) -> list[MovieInfo]:
    import nd2

    out = []
    with nd2.ND2File(path) as f:
        sizes = dict(f.sizes)
        n_pos = sizes.get("P", 1)
        t_ax = "T" if sizes.get("T", 1) > 1 else ("Z" if sizes.get("Z", 1) > 1 else "T")
        try:
            vx = f.voxel_size().x
        except Exception:
            vx = None
        try:
            names = [c.channel.name for c in (f.metadata.channels or [])]
        except Exception:
            names = []
        interval = _interval(f, sizes, t_ax)
        created = None
        try:
            created = os.path.getmtime(path)
        except OSError:
            pass
        for p in range(n_pos):
            info = MovieInfo(path=path, rel=rel, format="Nikon ND2", reader="nd2", series=p, n_series=n_pos,
                             series_name=f"P{p + 1}" if n_pos > 1 else "",
                             size_t=sizes.get(t_ax, 1), size_y=sizes.get("Y", 0), size_x=sizes.get("X", 0),
                             size_c=sizes.get("C", 1) if not f.is_rgb else 3,
                             size_z=sizes.get("Z", 1) if t_ax != "Z" else 1,
                             dtype=str(f.dtype), bits=np.dtype(f.dtype).itemsize * 8,
                             time_interval=interval, pixel_size=vx if vx and vx > 0 else None,
                             created_unix=created, channel_names=names or (["red", "green", "blue"] if f.is_rgb else []),
                             extra={"sizes": sizes, "t_axis": t_ax})
            if t_ax == "Z":
                info.notes.append("no time loop - the Z axis is used as time")
            info.notes.append("date taken from the file's modification time")
            info.name = output_name(os.path.basename(path), info.series_name, n_pos, p)
            out.append(info)
    return out


def _interval(f, sizes, t_ax):
    """Mean frame interval from the per-frame time stamps (fallback: the time loop period)."""
    if t_ax != "T" or sizes.get("T", 1) < 2:
        return None
    try:
        seqs = [i for i, idx in enumerate(f.loop_indices)
                if idx.get("P", 0) == 0 and idx.get("Z", 0) == 0]
        first = f.frame_metadata(seqs[0]).channels[0].time.relativeTimeMs
        last = f.frame_metadata(seqs[-1]).channels[0].time.relativeTimeMs
        if last > first:
            return (last - first) / (len(seqs) - 1) / 1000.0
    except Exception:
        pass
    try:
        for loop in f.experiment:
            if loop.type in ("TimeLoop", "NETimeLoop"):
                period = getattr(loop.parameters, "periodMs", None)
                if period is None and getattr(loop.parameters, "periods", None):
                    period = loop.parameters.periods[0].periodMs
                if period:
                    return float(period) / 1000.0
    except Exception:
        pass
    return None


def read(info: MovieInfo, opts: ImportOptions, max_frames: int | None = None) -> np.ndarray:
    import nd2

    ch = info.channel_index(opts)
    t_ax = info.extra["t_axis"]
    with nd2.ND2File(info.path) as f:
        frames = {}
        for seq, idx in enumerate(f.loop_indices):
            if idx.get("P", 0) != info.series:
                continue
            t = idx.get(t_ax, 0)
            if max_frames and t >= max_frames:
                continue
            z = idx.get("Z", 0) if t_ax != "Z" else 0
            zsel = str(opts.z_plane).lower()
            if t_ax != "Z" and zsel != "max" and z != int(zsel or 0):
                continue
            img = np.asarray(f.read_frame(seq))        # (C, Y, X), (Y, X) or (Y, X, S)
            if img.ndim == 3 and f.is_rgb:
                img = img[..., min(ch, img.shape[-1] - 1)]
            elif img.ndim == 3:
                img = img[min(ch, img.shape[0] - 1)]
            frames[t] = np.maximum(frames[t], img) if t in frames else img
    return np.stack([frames[t] for t in sorted(frames)])
