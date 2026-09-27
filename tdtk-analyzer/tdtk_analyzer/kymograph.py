"""Step 1: movie -> kymograph TIFFs + beat direction table.

Port of "MAYO Screen Script No.1" of tdtK_Full_Analysis_script (v0.7.4/0.7.5).
Array convention: the movie is [T, Y, X]; R's image_[x, y, t] is movie[t, y, x].
All positions reported (peak/X positions) are 1-based like in R.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import rio
from .movies import ImportOptions, MovieInfo, load_movie, movie_status
from .rstats import (
    find_peaks_m, lm_slope_adjr2, quantile7, r_sd, rollmean, rollmedian, trapz_unit, which_max_n,
)


@dataclass
class KymographSettings:
    min_file_size: int = 150_000_000   # bytes; smaller .cxd movies are skipped (as in R)
    min_frames: int = 200               # shorter movies are skipped
    max_frame_interval: float = 0.010   # s; slower movies are skipped (high-FPS only)
    max_loop_iterations: int = 1000     # safety cap for the window-size loops


@dataclass
class StepResult:
    name: str
    status: str          # "done", "skipped", "error"
    message: str = ""


def _stripe_borders(v: np.ndarray):
    """Find the bright/dark stripe borders along X (R lines ~350-410).

    Returns the list of 1-based borders or None when the movie must be skipped."""
    mv = np.mean(v)
    bright = v >= mv
    edge = np.where(bright, 100.0, 1.0)
    typ = np.where(bright, "bright", "dark")
    # duplicate the last row 29 times, rolling median over 30
    edge_ext = np.concatenate([edge, np.repeat(edge[-1], 29)])
    typ = np.concatenate([typ, np.repeat(typ[-1], 29)])
    borders = rollmedian(edge_ext, 30)
    edges = np.nonzero(borders == 50.5)[0] + 1
    if edges.size == 0:
        return None
    e1 = edges[0]
    if np.unique(typ[:e1]).size != 1:
        typ[:e1] = typ[e1 - 1]
    borders_ = list(edges) if typ[0] == "dark" else [1] + list(edges)
    if len(borders_) < 2:
        return None
    # while (diff(borders_[1:2]) < 60) borders_ <- borders_[3:length(borders_)]
    while borders_[1] - borders_[0] < 60:
        rest = borders_[2:]
        if len(rest) >= 2:
            borders_ = rest
            continue
        return None  # R ends up with < 2 usable borders (or NA) -> skipped / crashed
    return borders_


def _peak_positions(x: np.ndarray, borders_: list):
    """Kymograph X positions (1-based), R lines ~416-487."""
    ncol = x.shape[0]
    bright_area = np.array([trapz_unit(row) for row in x])
    b = borders_
    ba1 = bright_area[b[0] - 1:b[1]].copy()
    if len(b) <= 3:
        ba2 = np.array([0.0])
    else:
        ba2 = bright_area[b[2] - 1:b[3]].copy()
    for ba in (ba1, ba2):
        q4 = quantile7(ba, 0.75)
        ba[ba < q4] = q4
    peaks1 = pd.unique(find_peaks_m(ba1, 3))
    peaks2 = pd.unique(find_peaks_m(ba2, 3))
    b = list(b)
    if peaks1.size == 0:
        peaks1, peaks2 = peaks2, peaks2[:0]
        ba1 = ba2
        if len(b) > 3:
            b[0], b[1] = b[2], b[3]
    if peaks1.size == 0:
        return None
    peaklist = [peaks1[np.argmax(ba1[peaks1 - 1])] + b[0]]
    if peaks1.size > 1:
        peaklist.append(peaks1[which_max_n(ba1[peaks1 - 1], 2) - 1] + b[0])
    if peaks1.size > 2:
        peaklist.append(peaks1[which_max_n(ba1[peaks1 - 1], 3) - 1] + b[0])
    if peaks2.size != 0:
        peaklist.append(peaks2[np.argmax(ba2[peaks2 - 1])] + b[2])
        if peaks2.size > 1:
            peaklist.append(peaks2[which_max_n(ba2[peaks2 - 1], 2) - 1] + b[2])
        if peaks2.size > 2:
            peaklist.append(peaks2[which_max_n(ba2[peaks2 - 1], 3) - 1] + b[2])
    peaklist = np.array(peaklist, dtype=int)
    peaklist[peaklist > ncol] = ncol
    return peaklist


def count_peaks_in_windows(peaks_positions: np.ndarray, n_frames: int, window_size: int) -> np.ndarray:
    """sum(peaks_positions %in% i:(i + window_size)) for i in 1:(n_frames - window_size)."""
    i = np.array(rcolon_list(1, n_frames - window_size), dtype=np.int64)
    lo = np.minimum(i, i + window_size)
    hi = np.maximum(i, i + window_size)
    pp = peaks_positions[peaks_positions >= 1].astype(np.int64)
    nb = int(max(hi.max(initial=1), pp.max(initial=1), 1))
    cs = np.concatenate([[0], np.cumsum(np.bincount(pp, minlength=nb + 1)[1:nb + 1])])
    return cs[np.maximum(hi, 0)] - cs[np.maximum(lo - 1, 0)]


def rcolon_list(a: int, b: int) -> list:
    return list(range(a, b + 1)) if b >= a else list(range(a, b - 1, -1))


def binned_by_zero(x: np.ndarray):
    """Mid positions of zero-count tracks (R's binned_by_zero); None if none."""
    n = x.size
    if n < 3:
        return None
    i = np.arange(2, n)  # 2..n-1 (1-based)
    z = x == 0
    starts = i[~z[i - 2] & z[i - 1]]
    ends = i[z[i - 1] & ~z[i]]
    if ends.size == 0 or starts.size == 0:
        return None
    ends = ends.astype(float)
    if ends.size < starts.size:
        ends = np.concatenate([ends, np.full(starts.size - ends.size, np.nan)])
    else:
        ends = ends[: starts.size]
    miss = np.isnan(ends)
    ends[miss] = starts[miss]
    return np.round((starts + ends) / 2.0)  # R round(): half to even, same as numpy


def _direction_table(movie: np.ndarray, peaklist: np.ndarray, time_interval: float,
                     resolution_x: float, settings: KymographSettings) -> pd.DataFrame:
    """Beat direction per contraction ('_directionmarks.csv'), R lines ~530-835."""
    T, Y, X = movie.shape
    step = X // 24
    k = list(peaklist)
    if step > 0:
        k += list(range(step, X // 2 + 1, step))
    k = np.sort(np.array(k, dtype=int))
    k[k < 21] = 21
    k[k > X - 21] = X - 21
    k = pd.unique(k)

    extent = int(math.floor(Y * 0.35))
    ycols = np.concatenate([np.arange(0, extent), np.arange(Y - extent - 1, Y)])
    series = []
    for xpos in k:
        sub = movie[:, :, xpos - 21:xpos + 20]
        s = sub[:, ycols, :].mean(axis=(1, 2), dtype=np.float64)
        s = rollmean(s, 10)
        s = s - s.min()
        s = s / s.max()
        series.append(1.0 - s)
    peaks_position_list = [find_peaks_m(s, 10) for s in series]
    names = [str(int(v)) for v in k]
    peaks_positions = np.sort(np.concatenate(peaks_position_list)) if peaks_position_list else np.empty(0)
    L = series[0].size
    window_size = int(math.floor(r_sd(np.diff(peaks_positions))))

    def bins_for(ws):
        counted = count_peaks_in_windows(peaks_positions, L, ws)
        bz = binned_by_zero(counted)
        return np.concatenate([[0.0], bz if bz is not None else [], [float(T)]])

    bins_final = bins_for(window_size)
    target_bin = math.floor(np.mean([p.size for p in peaks_position_list]))
    if bins_final.size > 1.2 * target_bin:
        n = 1
        guard = 0
        while bins_final.size > 1.2 * target_bin or n < 11:
            n += 1
            window_size += 1
            bins_final = bins_for(window_size)
            guard += 1
            if guard >= settings.max_loop_iterations:
                break
    if bins_final.size < 0.8 * target_bin:
        n = 1
        guard = 0
        while bins_final.size < 1.2 * target_bin or n < 11 or window_size <= 0:
            n += 1
            window_size -= 1
            bins_final = bins_for(window_size)
            guard += 1
            if guard >= settings.max_loop_iterations:  # R would loop forever here
                break

    nk, nb = len(k), bins_final.size
    reg = np.full((nk, nb), np.nan)
    for i, pk in enumerate(peaks_position_list):
        for p in pk:
            j = np.searchsorted(bins_final, p, side="right")  # .bincode(right = FALSE)
            if 1 <= j <= nb - 1:
                reg[i, j - 1] = p
    keep = ~np.all(np.isnan(reg), axis=0)
    reg = reg[:, keep]
    slopes, rsq = [], []
    rows = np.arange(1, nk + 1, dtype=float)
    for c in range(reg.shape[1]):
        s, r = lm_slope_adjr2(reg[:, c], rows)
        slopes.append(s)
        rsq.append(r)
    slopes, rsq = np.array(slopes), np.array(rsq)
    good = np.sum(~np.isnan(reg), axis=0)
    thr = np.mean(np.unique(good)) if good.size else np.nan
    sel = good >= thr
    fb = reg[:, sel].T                      # rows = bins, cols = X positions
    slopes, rsq = slopes[sel], rsq[sel]
    xpos_nas = np.sum(np.isnan(fb), axis=0)
    if not np.any(xpos_nas == 0):
        n_na = np.sum(np.isnan(fb), axis=1)
        drop = n_na > 4
        if drop.any():
            fb, slopes, rsq = fb[~drop], slopes[~drop], rsq[~drop]
        else:  # R's x[-which(...)] with no match selects nothing
            fb, slopes, rsq = fb[:0], slopes[:0], rsq[:0]
        xpos_nas = np.sum(np.isnan(fb), axis=0)

    if np.any(xpos_nas == 0):
        ok = np.nonzero(xpos_nas < 2)[0]
        c1, c2 = ok.min(), ok.max()
        branch2 = False
    else:
        c1, c2 = 0, fb.shape[1] - 1
        branch2 = True
    n1, n2 = names[c1], names[c2]
    if c1 == c2:
        n2 = n2 + ".1"  # data.frame makes duplicated names unique
    a, b = fb[:, c1], fb[:, c2]
    delta = a - b
    with np.errstate(invalid="ignore"):
        direction = np.where(np.isnan(a) | np.isnan(b), np.nan, (a > b).astype(float))
    if branch2:
        nad = np.isnan(direction)
        direction[nad & (slopes < 0)] = 1
        direction[nad & (slopes > 0)] = 0
    distance = float(n2) - float(n1) * resolution_x   # (sic) operator precedence as in R
    with np.errstate(divide="ignore", invalid="ignore"):
        speed = distance / (delta * time_interval)
        neg = (delta == 0) & (slopes < 0) & (rsq > 0.15)
        pos = (delta == 0) & (slopes > 0) & (rsq > 0.15)
    direction[pos] = 1
    direction[neg] = 0
    return pd.DataFrame({n1: a, n2: b, "delta": delta, "direction": direction,
                         "slope": slopes, "speed": speed, "rsq": rsq})


def meta_pairs(info: MovieInfo, opts: ImportOptions, shape: tuple, rotation: str) -> list[tuple[str, object]]:
    """The '<movie>._new_meta_data.csv' rows: R's 21 values + the movie file and its format."""
    T, Y, X = shape
    interval, pixel = info.effective_interval(opts), info.effective_pixel(opts)
    if info.reader == "cxd" and "cxd" in info.extra:
        ci = info.extra["cxd"]
        pairs = [(k, v) for k, v in ci.meta_table()[:18]]
    else:
        pairs = [("sizeX", X), ("sizeY", Y), ("sizeZ", info.size_z), ("sizeC", 1), ("sizeT", T),
                 ("pixelType", info.dtype), ("bitsPerPixel", info.bits), ("imageCount", T),
                 ("dimensionOrder", "XYCZT"), ("orderCertain", "XYCZT"), ("rgb", "false"),
                 ("littleEndian", "true"), ("interleaved", "false"), ("falseColor", 0),
                 ("metadataComplete", 1), ("thumbnail", 0), ("series", info.series + 1), ("resolutionLevel", 1)]
    pairs[0], pairs[1] = ("sizeX", X), ("sizeY", Y)      # size of the movie as analyzed (after rotation)
    created = info.created_unix if info.created_unix is not None else float("nan")
    def src(key, own):
        return info.sources.get(key, "file") if own else "default setting"

    return pairs + [("time_interval", interval), ("resolution", pixel), ("created_unix_from_file", created),
                    ("movie_file", info.name), ("source_format", info.format),
                    ("source_series", info.series_name or str(info.series)), ("rotation", rotation),
                    ("time_interval_source", src("time_interval", bool(info.time_interval))),
                    ("resolution_source", src("pixel_size", bool(info.pixel_size))),
                    ("created_source", info.sources.get("created", "missing"))]


def process_movie(info: MovieInfo, opts: ImportOptions | None = None,
                  settings: KymographSettings | None = None) -> StepResult:
    """Process one movie (any supported format); outputs are written next to the movie file."""
    settings = settings or KymographSettings()
    opts = opts or ImportOptions()
    folder = os.path.dirname(info.path)
    name = info.name
    label = info.rel + (f" [{info.series_name or info.series + 1}]" if info.n_series > 1 else "")
    base = os.path.join(folder, name)
    try:
        tiffs = [f for f in os.listdir(folder) if (name + "_peak_") in f]
        if tiffs and os.path.exists(base + "_directionmarks.csv"):
            return StepResult(label, "skipped", "already processed")
        if info.reader == "cxd" and os.path.getsize(info.path) < settings.min_file_size:
            return StepResult(label, "skipped", "file smaller than the minimum size")
        level, why = movie_status(info, opts, settings.min_frames, settings.max_frame_interval)
        if level in ("error", "skip"):
            return StepResult(label, "skipped" if level == "skip" else "error", why)
        time_interval = info.effective_interval(opts)
        resolution = info.effective_pixel(opts)

        movie, rotation = load_movie(info, opts)
        T, Y, X = movie.shape
        rio.write_meta_csv(meta_pairs(info, opts, movie.shape, rotation),
                           os.path.join(folder, rio.meta_csv_for(name)))
        vmin, vmax = float(movie.min()), float(movie.max())
        if vmax <= vmin:
            return StepResult(label, "error", "the selected channel is empty (constant)")

        # SD of every pixel over time -> map [X, Y]
        sd = np.empty((Y, X))
        rows = max(1, int(256e6 // (T * X * 8)))
        for y0 in range(0, Y, rows):
            chunk = movie[:, y0:y0 + rows, :].astype(np.float64)
            sd[y0:y0 + rows] = chunk.std(axis=0, ddof=1)
        x = sd.T / sd.max()
        x[x < quantile7(x.ravel(), 0.75)] = 0
        v = x.sum(axis=1)

        borders_ = _stripe_borders(v)
        if borders_ is None:
            return StepResult(label, "skipped", "no heart stripes found (is the heart horizontal? see 'rotate')")
        peaklist = _peak_positions(x, borders_)
        if peaklist is None:
            return StepResult(label, "skipped", "no kymograph positions found")

        q = x.copy()
        q[peaklist - 1, :] = 1
        rio.write_tiff8(q.T, os.path.join(folder, rio.movie_prefix(name) + "_SD and peaklines.tiff"))
        peaklist = pd.unique(peaklist)
        rng = vmax - vmin
        for i, p in enumerate(peaklist, start=1):
            kymo = (movie[:, :, p - 1].astype(np.float64) - vmin) / rng   # [T, Y]
            rio.write_tiff8(kymo.T, f"{base}_peak_{i}_at Xpos_{int(p)}.tiff")

        if os.path.exists(base + "_directionmarks.csv"):
            return StepResult(label, "done", f"{len(peaklist)} kymographs")
        table = _direction_table(movie, np.asarray(peaklist), time_interval, resolution, settings)
        rio.write_r_csv(table, base + "_directionmarks.csv")
        return StepResult(label, "done", f"{len(peaklist)} kymographs")
    except Exception as e:  # keep the batch running; report the file
        return StepResult(label, "error", f"{type(e).__name__}: {e}")
