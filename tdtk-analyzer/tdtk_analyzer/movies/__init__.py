"""Movie import: reads scientific movie formats into one common model.

    infos = scan_movies(folder)                 # one MovieInfo per movie (series)
    movie = load_movie(info, ImportOptions())   # numpy [T, Y, X]

Supported without extra software:
    Hamamatsu .cxd, OME-TIFF, ImageJ / Micro-Manager / plain TIFF stacks, Zeiss .lsm and
    .czi, MetaMorph .stk, Nikon .nd2, Leica .lif, Olympus .oib/.oif, video (.avi .mp4 .mov .mkv).
With Bio-Formats' command line tools (bftools) installed, every other Bio-Formats
format (Imaris .ims, DeltaVision .dv, Olympus .vsi, .ics, ...) as well.

The import table ('movie_import.csv' in the output folder) stores per-movie
choices (use, output name, channel, rotation, frame interval, pixel size).
"""

from __future__ import annotations

import importlib
import os
import re
from typing import Callable

import numpy as np
import pandas as pd

from .base import ImportOptions, MovieInfo, auto_rotation, rotate_movie  # noqa: F401
from .metadata import Issue, check_metadata, format_date, parse_date  # noqa: F401

READER_BY_EXT = {
    ".cxd": "cxd",
    ".ome.tif": "tiff", ".ome.tiff": "tiff", ".tif": "tiff", ".tiff": "tiff", ".btf": "tiff", ".tf8": "tiff",
    ".lsm": "tiff", ".stk": "tiff",
    ".nd2": "nd2", ".czi": "czi", ".lif": "lif", ".oib": "oif", ".oif": "oif",
    ".avi": "video", ".mp4": "video", ".mov": "video", ".mkv": "video",
}
BIOFORMATS_EXT = {".ims", ".dv", ".r3d", ".vsi", ".ics", ".ids", ".zvi", ".lei", ".sld", ".nd", ".ser",
                  ".dcimg", ".mvd2", ".obf", ".msr", ".lof", ".xlef", ".scn", ".ipl", ".liff", ".dm3", ".dm4",
                  ".apl", ".mrc", ".sif", ".fli", ".pic", ".mea", ".oir", ".vws"}
FORMAT_NAMES = {
    "cxd": "Hamamatsu HCImage (.cxd)", "tiff": "TIFF / OME-TIFF / ImageJ / Micro-Manager / LSM / STK",
    "nd2": "Nikon NIS-Elements (.nd2)", "czi": "Zeiss ZEN (.czi)", "lif": "Leica LAS X (.lif)",
    "oif": "Olympus FluoView (.oib, .oif)", "video": "Video (.avi, .mp4, .mov, .mkv)",
    "bioformats": "other formats via Bio-Formats (needs bftools)",
}
# files written by this program next to the movies - never treated as movies
GENERATED = re.compile(r"(_peak_\d+_at Xpos_\d+\.tiff|_SD and peaklines\.tiff|_traced\.jpg|\.tiff\.csv)$")
SKIP_DIRS = {"import_cache", ".tdtk_import_cache", "TIFFs", "balled", "__MACOSX"}
MANIFEST = "movie_import.csv"


def _ext(name: str) -> str:
    low = name.lower()
    for e in (".ome.tiff", ".ome.tif"):
        if low.endswith(e):
            return e
    return os.path.splitext(low)[1]


def reader_for(path: str) -> str | None:
    e = _ext(path)
    if e in READER_BY_EXT:
        return READER_BY_EXT[e]
    if e in BIOFORMATS_EXT:
        from .bioformats_reader import available
        return "bioformats" if available() else None
    return None


def _module(key: str):
    return importlib.import_module(f"{__name__}.{key}_reader")


def probe_file(path: str, rel: str | None = None) -> list[MovieInfo]:
    """Describe the movie(s) in a file. Problems are reported in MovieInfo.error."""
    rel = rel or os.path.basename(path)
    key = reader_for(path)
    if key is None:
        return [MovieInfo(path=path, rel=rel, format=_ext(path), reader="", name=os.path.basename(path),
                          error="unsupported file format")]
    try:
        infos = _module(key).probe(path, rel)
        if not infos:
            raise ValueError("no image data found")
        for m in infos:
            for attr, k in (("time_interval", "time_interval"), ("pixel_size", "pixel_size"), ("created_unix", "created")):
                if getattr(m, attr) is not None and k not in m.sources:
                    m.sources[k] = "file"
            m.remember_original()
        return infos
    except ImportError as e:
        err = f"reader library missing ({e.name}); install it with pip"
    except Exception as e:
        err = f"{type(e).__name__}: {e}"
    # 'use' stays the user's choice; the error itself keeps the movie out of the analysis
    return [MovieInfo(path=path, rel=rel, format=FORMAT_NAMES.get(key, key), reader=key,
                      name=os.path.basename(path), error=err)]


def find_movie_files(folder: str, skip: list[str] | None = None) -> list[str]:
    """Relative paths of all movie files below ``folder`` (sorted)."""
    skip_abs = {os.path.abspath(s) for s in (skip or []) if s}
    exts = set(READER_BY_EXT) | BIOFORMATS_EXT
    out = []
    for root, dirs, files in os.walk(folder):
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS and not d.startswith(".")
                         and os.path.abspath(os.path.join(root, d)) not in skip_abs)
        for f in files:
            if f.startswith(".") or GENERATED.search(f):
                continue
            if _ext(f) in exts:
                out.append(os.path.relpath(os.path.join(root, f), folder).replace(os.sep, "/"))
    return sorted(out)


def scan_movies(folder: str, output_dir: str | None = None,
                progress: Callable[[int, int, str], None] | None = None) -> list[MovieInfo]:
    files = find_movie_files(folder, skip=[output_dir] if output_dir else None)
    if output_dir:
        from . import bioformats_reader
        bioformats_reader.CACHE_DIR["path"] = os.path.join(output_dir, "import_cache")
    infos: list[MovieInfo] = []
    for i, rel in enumerate(files):
        if progress:
            progress(i, len(files), rel)
        infos.extend(probe_file(os.path.join(folder, rel), rel))
    if progress:
        progress(len(files), len(files), "")
    _unique_names(infos)
    return infos


def _unique_names(infos: list[MovieInfo]) -> None:
    """Output names must be unique within a folder (outputs are written next to the movie)."""
    seen = {}
    for m in infos:
        key = (os.path.dirname(m.rel), m.name)
        if key in seen:
            root, ext = os.path.splitext(m.name)
            m.name = f"{root}_{seen[key] + 1}{ext}"
        seen[key] = seen.get(key, 0) + 1


def load_movie(info: MovieInfo, opts: ImportOptions, max_frames: int | None = None) -> tuple[np.ndarray, str]:
    """Load the movie as [T, Y, X]; returns (movie, rotation applied)."""
    movie = _module(info.reader).read(info, opts, max_frames)
    if movie.ndim != 3:
        raise ValueError(f"expected a [T, Y, X] movie, got shape {movie.shape}")
    rot = info.rotation(opts)
    if rot == "auto":
        rot = auto_rotation(movie)
    return rotate_movie(movie, rot), rot


# --------------------------------------------------------------------------- status

def movie_issues(info: MovieInfo, opts: ImportOptions) -> list[Issue]:
    if info.error:
        return [Issue("file", "error", info.error)]
    issues = check_metadata(info, opts)
    try:
        info.channel_index(opts)
    except ValueError as e:
        issues.append(Issue("channel", "error", str(e)))
    return issues


def movie_status(info: MovieInfo, opts: ImportOptions, min_frames: int = 200,
                 max_interval_s: float | None = 0.010) -> tuple[str, str]:
    """(level, message) with level 'error', 'skip', 'warn' or 'ok'."""
    issues = movie_issues(info, opts)
    errors = [i for i in issues if i.level == "error"]
    warns = [i for i in issues if i.level == "warn"]
    if errors:                                   # errors first, but do not hide the warnings
        return "error", "; ".join(str(i) for i in errors + warns)
    if not info.use:
        return "skip", "not selected"
    if info.size_t < min_frames:
        return "skip", f"only {info.size_t} frames (minimum {min_frames})"
    ti = info.effective_interval(opts)
    if max_interval_s and ti and ti > max_interval_s:
        return "skip", f"frame interval {ti * 1000:.2f} ms is longer than {max_interval_s * 1000:.2f} ms"
    if warns:
        return "warn", "; ".join(str(i) for i in warns)
    return "ok", "ready"


# --------------------------------------------------------------------------- import table

TABLE_COLUMNS = ["use", "file", "series", "series_name", "name", "format", "frames", "height", "width",
                 "channels", "channel_names", "z_planes", "frame_interval_ms", "pixel_size_um", "recorded",
                 "channel", "rotate", "frame_interval_source", "pixel_size_source", "recorded_source", "issues"]


def to_table(infos: list[MovieInfo]) -> pd.DataFrame:
    rows = []
    for m in infos:
        rows.append({
            "use": "yes" if m.use else "no", "file": m.rel, "series": m.series, "series_name": m.series_name,
            "name": m.name, "format": m.format, "frames": m.size_t, "height": m.size_y, "width": m.size_x,
            "channels": m.size_c, "channel_names": "|".join(map(str, m.channel_names)), "z_planes": m.size_z,
            "frame_interval_ms": round(m.time_interval * 1000, 6) if m.time_interval else None,
            "pixel_size_um": round(m.pixel_size, 6) if m.pixel_size else None,
            "recorded": format_date(m.created_unix),
            "channel": m.channel or "", "rotate": m.rotate or "",
            "frame_interval_source": m.sources.get("time_interval", "missing"),
            "pixel_size_source": m.sources.get("pixel_size", "missing"),
            "recorded_source": m.sources.get("created", "missing"),
            "issues": m.error or "; ".join(str(i) for i in check_metadata(m, ImportOptions())),
        })
    return pd.DataFrame(rows, columns=TABLE_COLUMNS)


def save_table(infos: list[MovieInfo], output_dir: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, MANIFEST)
    to_table(infos).to_csv(path, index=False)
    return path


def apply_table(infos: list[MovieInfo], table: pd.DataFrame) -> None:
    """Apply the user's choices from an import table to freshly probed movies."""
    table = table.copy()
    table["series"] = pd.to_numeric(table["series"], errors="coerce").fillna(0).astype(int)
    rows = {(str(r["file"]), int(r["series"])): r for _, r in table.iterrows()}
    for m in infos:
        r = rows.get((m.rel, m.series))
        if r is None:
            continue
        m.use = str(r.get("use", "yes")).strip().lower() not in ("no", "n", "false", "0", "")
        if isinstance(r.get("name"), str) and r["name"].strip():
            from ..rio import has_movie_extension
            from .base import split_movie_ext
            new = r["name"].strip()
            if not has_movie_extension(new):      # outputs are traced back through the extension
                new += split_movie_ext(os.path.basename(m.path))[1]
            m.name = new
        # a value is taken from the table when the user entered it (source 'entered')
        # or when it differs from what the file says (e.g. edited in Excel)
        for col, attr, scale, key in (("frame_interval_ms", "time_interval", 1e-3, "frame_interval_source"),
                                      ("pixel_size_um", "pixel_size", 1.0, "pixel_size_source")):
            v = pd.to_numeric(r.get(col), errors="coerce")
            entered = str(r.get(key, "")).strip() == "entered"
            if pd.notna(v) and v > 0:
                v = float(v) * scale
                cur = getattr(m, attr)
                if entered or cur is None or abs(v - cur) > 1e-6 * max(abs(cur), 1e-9) + 1e-12:
                    m.set_value(attr, v)
        rec = str(r.get("recorded", "") or "").strip()
        if rec and (str(r.get("recorded_source", "")).strip() == "entered" or rec != format_date(m.created_unix)):
            try:
                m.set_value("created_unix", parse_date(rec))
            except ValueError:
                pass
        for col in ("channel", "rotate"):
            v = r.get(col)
            if isinstance(v, (int, float)) and pd.notna(v):
                v = str(int(v))
            if isinstance(v, str) and v.strip():
                setattr(m, col, v.strip())


def load_table(output_dir: str) -> pd.DataFrame | None:
    path = os.path.join(output_dir, MANIFEST)
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, dtype=str, keep_default_na=False)
