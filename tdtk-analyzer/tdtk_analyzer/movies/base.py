"""Format-independent movie description and helpers shared by all readers.

Every reader turns a file (or one series of a multi-series file) into a
``MovieInfo`` and can load it as a numpy array [T, Y, X] of the selected
channel / z-plane. Everything downstream only sees this.
"""

from __future__ import annotations

import os
import re
from dataclasses import asdict, dataclass, field

import numpy as np

# channel names that most likely hold tdTomato (the heart signal)
RED_CHANNEL_HINTS = ("tdt", "tom", "mcherry", "rfp", "red", "dsred", "561", "568", "594", "tritc", "cy3", "texas")


@dataclass
class ImportOptions:
    """How to turn a multi-dimensional recording into the [T, Y, X] movie."""
    channel: str = "auto"          # "auto", a 0-based index ("1") or part of a channel name ("mCherry")
    z_plane: str = "0"             # 0-based index or "max" (maximum projection)
    rotate: str = "0"              # "0", "90", "180", "270" (counter-clockwise) or "auto"
    default_interval_ms: float = 0.0   # used when a file has no frame interval (0 = none)
    default_pixel_um: float = 0.0      # used when a file has no pixel size (0 = none)


@dataclass
class MovieInfo:
    path: str                       # absolute path of the file
    rel: str                        # path relative to the movie folder
    format: str                     # e.g. "Nikon ND2"
    reader: str                     # reader key
    name: str = ""                  # output base name ('<name>_peak_1_at Xpos_...tiff')
    series: int = 0
    series_name: str = ""
    n_series: int = 1
    size_t: int = 0
    size_y: int = 0
    size_x: int = 0
    size_c: int = 1
    size_z: int = 1
    dtype: str = "uint16"
    bits: int = 0
    time_interval: float | None = None   # seconds
    pixel_size: float | None = None      # micrometre per pixel
    created_unix: float | None = None
    channel_names: list = field(default_factory=list)
    notes: list = field(default_factory=list)   # e.g. "time axis taken from Z"
    extra: dict = field(default_factory=dict)   # reader specific (axis layout, ...)
    # per-movie settings (from the import table); None = use the global ImportOptions
    use: bool = True
    channel: str | None = None
    rotate: str | None = None
    error: str = ""

    def to_row(self) -> dict:
        d = asdict(self)
        d.pop("extra")
        return d

    # ------------------------------------------------------------------ helpers
    def channel_index(self, opts: ImportOptions) -> int:
        spec = (self.channel if self.channel not in (None, "") else opts.channel) or "auto"
        return resolve_channel(spec, self.channel_names, self.size_c)

    def rotation(self, opts: ImportOptions) -> str:
        return str(self.rotate if self.rotate not in (None, "") else opts.rotate or "0")

    def effective_interval(self, opts: ImportOptions) -> float | None:
        if self.time_interval and self.time_interval > 0:
            return self.time_interval
        return opts.default_interval_ms / 1000.0 if opts.default_interval_ms > 0 else None

    def effective_pixel(self, opts: ImportOptions) -> float | None:
        if self.pixel_size and self.pixel_size > 0:
            return self.pixel_size
        return opts.default_pixel_um if opts.default_pixel_um > 0 else None


def resolve_channel(spec: str, names: list, size_c: int) -> int:
    spec = str(spec).strip()
    if size_c <= 1:
        return 0
    if spec.lower() in ("", "auto"):
        for i, n in enumerate(names):
            if any(h in str(n).lower() for h in RED_CHANNEL_HINTS):
                return i
        return 0
    if re.fullmatch(r"\d+", spec):
        i = int(spec)
        if i >= size_c:
            raise ValueError(f"channel {i} requested but the file has {size_c} channels (0-{size_c - 1})")
        return i
    for i, n in enumerate(names):
        if spec.lower() in str(n).lower():
            return i
    raise ValueError(f"no channel named like '{spec}' (channels: {', '.join(map(str, names))})")


def select_planes(arr: np.ndarray, axes: str, t_axis: str, channel: int, z_plane: str) -> np.ndarray:
    """Reduce an array with named axes to [T, Y, X].

    ``t_axis`` is the letter used as time; 'C'/'S' is the channel; 'Z' is
    reduced by ``z_plane``; any other axis is reduced to its first element."""
    axes = axes.upper()
    for ax in list(axes):
        if ax in (t_axis, "Y", "X"):
            continue
        i = axes.index(ax)
        if ax in ("C", "S"):
            arr = np.take(arr, min(channel, arr.shape[i] - 1), axis=i)
        elif ax == "Z" and str(z_plane).lower() == "max":
            arr = arr.max(axis=i)
        elif ax == "Z":
            arr = np.take(arr, min(int(z_plane or 0), arr.shape[i] - 1), axis=i)
        else:
            arr = np.take(arr, 0, axis=i)
        axes = axes.replace(ax, "", 1)
    if t_axis not in axes:              # a single image
        arr = arr[None]
        axes = t_axis + axes
    order = [axes.index(t_axis), axes.index("Y"), axes.index("X")]
    return np.transpose(arr, order)


def pick_time_axis(axes: str, shape: tuple) -> tuple[str, str]:
    """Return (time axis letter, note). Stacks without a T axis often store time as I, Q or Z."""
    axes = axes.upper()
    if "T" in axes and shape[axes.index("T")] > 1:
        return "T", ""
    for cand in ("I", "Q", "Z"):
        if cand in axes and shape[axes.index(cand)] > 1:
            return cand, f"no time axis in the file - the {cand} axis ({shape[axes.index(cand)]} planes) is used as time"
    return ("T" if "T" in axes else "T"), ""


def rotate_movie(movie: np.ndarray, rotation: str) -> np.ndarray:
    k = {"0": 0, "90": 1, "180": 2, "270": 3}.get(str(rotation), 0)
    return np.rot90(movie, k=k, axes=(1, 2)) if k else movie


def auto_rotation(movie: np.ndarray) -> str:
    """'0' if the heart tube lies along X (as the analysis expects), '90' if along Y.

    The moving heart walls have a high temporal SD. For a horizontal heart the
    SD is concentrated in a narrow band of rows but spread over many columns."""
    t = movie.shape[0]
    step = max(1, t // 200)
    sd = movie[::step].astype(np.float32).std(axis=0)
    rows, cols = sd.sum(axis=1), sd.sum(axis=0)

    def spread(p):
        p = np.sort(p)[::-1]
        c = np.cumsum(p) / max(p.sum(), 1e-12)
        return (np.searchsorted(c, 0.8) + 1) / p.size    # fraction of lines holding 80 % of the SD

    return "0" if spread(rows) < spread(cols) else "90"


def output_name(file_name: str, series_name: str, n_series: int, series: int) -> str:
    """Name used for all output files of a movie; keeps the file extension so that
    '<name>_peak_...' and '<name>_directionmarks.csv' can be traced back."""
    if n_series <= 1:
        return file_name
    root, ext = split_movie_ext(file_name)
    s = re.sub(r"[^\w\-]+", "-", series_name).strip("-") if series_name else f"S{series + 1:02d}"
    # series named like the flies (e.g. 'MAYO0001_1_1wf') are used as they are
    if series_name and series_name.count("_") >= 2:
        return f"{s}{ext}"
    return f"{root}_{s}{ext}"


def split_movie_ext(file_name: str) -> tuple[str, str]:
    low = file_name.lower()
    for ext in (".ome.tiff", ".ome.tif"):
        if low.endswith(ext):
            return file_name[: -len(ext)], file_name[-len(ext):]
    root, ext = os.path.splitext(file_name)
    return root, ext
