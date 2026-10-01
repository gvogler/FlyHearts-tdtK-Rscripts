"""File helpers that write files the way the R pipeline did.

* CSV files like ``write.csv``/``write.table(sep = ",")``: quoted header and
  strings, ``NA`` for missing values, up to 15 significant digits.
* 8-bit TIFFs like ``EBImage::writeImage`` (values in [0, 1], clipped).
* R matrices are indexed [x, y]; image files are [row, col] = [y, x], so a
  matrix written by EBImage corresponds to ``m.T`` here.
"""

from __future__ import annotations

import csv
import math
import os
import re

import numpy as np
import pandas as pd


def fmt_num(v) -> str:
    if v is None:
        return "NA"
    if isinstance(v, (bool, np.bool_)):
        return "TRUE" if v else "FALSE"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if math.isnan(f):
        return "NA"
    if math.isinf(f):
        return "Inf" if f > 0 else "-Inf"
    if f.is_integer() and abs(f) < 1e15:
        return str(int(f))
    return "%.15g" % f


def _cell(v) -> tuple[str, bool]:
    """Return (text, quote)."""
    if isinstance(v, str):
        return v, True
    if v is None or (isinstance(v, float) and math.isnan(v)) or v is pd.NA or v is pd.NaT:
        return "NA", False
    if isinstance(v, (np.floating, float, int, np.integer, bool, np.bool_)):
        return fmt_num(v), False
    return str(v), True


def write_r_csv(df: pd.DataFrame, path: str) -> None:
    """write.csv(df, path, row.names = FALSE)."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        fh.write(",".join('"%s"' % str(c).replace('"', '""') for c in df.columns) + "\n")
        cols = [df[c].tolist() for c in df.columns]
        for row in zip(*cols) if cols else []:
            out = []
            for v in row:
                text, quote = _cell(v)
                out.append('"%s"' % text.replace('"', '""') if quote else text)
            fh.write(",".join(out) + "\n")


def read_r_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, na_values=["NA", "NaN", "Inf", "-Inf"], keep_default_na=True)


def write_tiff8(values: np.ndarray, path: str) -> None:
    """EBImage::writeImage(x, file) for a grey image; ``values`` is [row, col]."""
    import tifffile

    v = np.clip(np.nan_to_num(np.asarray(values, dtype=float), nan=0.0), 0.0, 1.0)
    tifffile.imwrite(path, np.floor(v * 255.0 + 0.5).astype(np.uint8), photometric="minisblack")


def read_tiff_gray(path: str) -> np.ndarray | None:
    """Read a TIFF like RBioFormats::read.image (normalized to [0, 1]).

    Returns None for RGB/multi-plane images (which the R script skipped)."""
    import tifffile

    a = tifffile.imread(path)
    if a.ndim != 2:
        return None
    if np.issubdtype(a.dtype, np.integer):
        return a.astype(float) / float(np.iinfo(a.dtype).max)
    return a.astype(float)


def write_rgb_jpeg(red: np.ndarray, green: np.ndarray, blue: np.ndarray, path: str) -> None:
    from PIL import Image

    rgb = np.stack([red, green, blue], axis=-1)
    rgb = np.floor(np.clip(np.nan_to_num(rgb), 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
    Image.fromarray(rgb).save(path, quality=100)       # uint8 H x W x 3 -> RGB


MOVIE_EXTENSIONS = sorted(
    ["ome.tiff", "ome.tif", "cxd", "tiff", "tif", "btf", "tf8", "lsm", "stk", "nd2", "czi", "lif", "oib", "oif",
     "avi", "mp4", "mov", "mkv", "ims", "dv", "r3d", "vsi", "ics", "ids", "zvi", "lei", "sld", "nd", "ser",
     "dcimg", "mvd2", "obf", "msr", "lof", "xlef", "scn", "ipl", "liff", "dm3", "dm4", "apl", "mrc", "sif",
     "fli", "pic", "mea", "oir", "vws"], key=len, reverse=True)
# first '.<movie extension>' that is followed by '_' or the end: 'fly.nd2_peak_1_at Xpos_5.tiff'
_MOVIE_RE = re.compile(r"\.(?:" + "|".join(re.escape(e) for e in MOVIE_EXTENSIONS) + r")(?=_|$)", re.IGNORECASE)


def _movie_match(name: str):
    m = _MOVIE_RE.search(name)
    if m:
        return m.start(), m.end()
    m = re.search(r".cxd", name)          # the R script's rule
    return (m.start(), m.start() + 4) if m else None


def movie_prefix(name: str) -> str:
    """Text up to the '.' of the movie extension (R: substr(x, 1, gregexpr('.cxd', x)[[1]][1]))."""
    m = _movie_match(name)
    return name[: m[0] + 1] if m else name


def movie_stem(name: str) -> str:
    """Movie file name at the start of an output name (R: substr(x, 1, gregexpr('.cxd', x) + 3))."""
    m = _movie_match(name)
    return name[: m[1]] if m else name


def has_movie_extension(name: str) -> bool:
    return bool(_MOVIE_RE.search(name))


# names used by the R script
cxd_prefix = movie_prefix
cxd_stem = movie_stem


def meta_csv_for(name: str) -> str:
    return movie_prefix(name) + "_new_meta_data.csv"


def write_meta_csv(pairs: list[tuple[str, object]], path: str) -> None:
    """write.csv(melt(cxd_meta.data)) - two columns 'value' and 'L1', all quoted."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, quoting=csv.QUOTE_ALL, lineterminator="\n")
        w.writerow(["value", "L1"])
        for k, v in pairs:
            w.writerow([v if isinstance(v, str) else fmt_num(v), k])


def read_meta_csv(path: str) -> dict:
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    return dict(zip(df["L1"], df["value"]))


def meta_float(meta: dict, key: str) -> float:
    try:
        return float(meta.get(key, "nan"))
    except ValueError:
        return float("nan")


def list_files(folder: str, pattern: str, recursive: bool = False, full_names: bool = False) -> list[str]:
    """list.files(folder, pattern, recursive, full.names) - sorted like R (C locale)."""
    rx = re.compile(pattern)
    out = []
    if recursive:
        for root, _dirs, files in os.walk(folder):
            for f in files:
                if rx.search(f):
                    rel = os.path.relpath(os.path.join(root, f), folder)
                    out.append(rel.replace(os.sep, "/"))
    else:
        for f in os.listdir(folder):
            if os.path.isfile(os.path.join(folder, f)) and rx.search(f):
                out.append(f)
    out.sort()
    if full_names:
        return [os.path.join(folder, f) for f in out]
    return out
