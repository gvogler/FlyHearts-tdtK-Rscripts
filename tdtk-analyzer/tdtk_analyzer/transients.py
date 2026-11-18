"""Step 3a: beat (transient) analysis of one traced kymograph.

Port of the per-fly loop of "MAYO Screen Script No.4" and of the spline
metrics ("Transient analysis pt.1") of tdtK_Full_Analysis_script.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import rio
from .rstats import (
    NA, SmoothSpline, baseline_rolling_ball, find_peaks_simple, r_mean_narm, r_median, trapz_unit,
)

METRIC_NAMES = ["peaked", "tt10p", "tt25p", "tt50p", "tt75p", "tt90p",
                "tt90r", "tt75r", "tt50r", "tt25r", "tt10r"]
AVGT = 400  # rolling-mean window for the peak/valley threshold


@dataclass
class TraceResult:
    csv: str
    ok: bool = True
    message: str = ""
    timescale: float = NA
    total_time: float = NA
    reversals: float = NA
    anterior_beat_percent: float = NA
    speed_anterograde: float = NA
    speed_retrograde: float = NA
    edd: float = NA            # median end-diastolic "diameter" (inverted distance, R's EDD_this_fly)
    esd: float = NA
    fs: float = NA
    intervals: pd.DataFrame | None = None
    transients: pd.DataFrame | None = None
    n_identicals: int = 0
    time_used: float = NA
    max_velocity: list = field(default_factory=list)
    max_neg_velocity: list = field(default_factory=list)
    metrics: list = field(default_factory=list)   # one list of 11 values per transient


def spline_metrics(t: np.ndarray, d: np.ndarray) -> list:
    """smooth.spline + predict on a 0.2 ms grid + metrics() of the R script."""
    s = SmoothSpline(t, d)
    m = float(np.max(s.ux))
    n = int(math.floor(m / 0.0002 + 1e-10)) + 1
    gx = np.arange(n) * 0.0002
    gy = s.predict(gx)
    mv = int(np.argmax(gy)) + 1            # which.max (1-based)
    contraction = gy[:mv]
    relaxation = gy[mv:] if mv < gy.size else np.array([np.nan, gy[-1]])

    def first(mask, offset=0):
        w = np.nonzero(mask)[0]
        if w.size == 0:
            return NA
        j = w[0] + 1 + offset               # 1-based position in gx
        return float(gx[j - 1]) if j <= gx.size else NA

    out = [float(gx[mv - 1])]
    for thr in (0.1, 0.25, 0.5, 0.75, 0.9):
        out.append(first(contraction > thr))
    with np.errstate(invalid="ignore"):
        for thr in (0.9, 0.75, 0.5, 0.25, 0.10):
            out.append(first(relaxation < thr, mv))
    return out


def _directions(dfile: str | None, res: TraceResult) -> None:
    if not dfile or not os.path.exists(dfile):
        return
    d = rio.read_r_csv(dfile)
    if d.shape[1] < 6:
        return
    direction = pd.to_numeric(d["direction"], errors="coerce").to_numpy(float, copy=True)
    speed = pd.to_numeric(d.iloc[:, 5], errors="coerce").to_numpy(float, copy=True)
    speed[np.isinf(speed)] = np.nan
    valid = direction[~np.isnan(direction)]
    # reversals: as.numeric(as.factor(direction)) -> table(abs(diff(.)))
    levels = np.unique(valid)
    codes = np.array([np.nan if np.isnan(v) else float(np.searchsorted(levels, v) + 1) for v in direction])
    dd = np.abs(np.diff(codes))
    dd = dd[~np.isnan(dd)]
    tab = np.unique(dd, return_counts=True)
    if tab[0].size == 0:
        res.reversals = NA
    elif tab[0].size == 1:
        res.reversals = 0.0
    else:
        res.reversals = float(tab[1][1])
    total = valid.size
    ante = float(np.sum(valid == 1))
    res.anterior_beat_percent = ante / total if total else NA
    if np.any(direction == 1):
        res.speed_anterograde = r_mean_narm(speed[direction == 1])
    if np.any(direction == 0):
        res.speed_retrograde = r_mean_narm(speed[direction == 0])


def analyze_trace(folder: str, csv_name: str, direction_file: str | None) -> TraceResult:
    """Analyze one traced kymograph CSV ('<...>.tiff.csv')."""
    res = TraceResult(csv=csv_name)
    try:
        raw = rio.read_r_csv(os.path.join(folder, csv_name))
        time = raw.iloc[:, 0].to_numpy(float, copy=True)
        dist = raw.iloc[:, 1].to_numpy(float, copy=True)
        n = time.size
        res.timescale = float(np.median(np.diff(time)))
        res.total_time = n * res.timescale
        _directions(direction_file, res)

        x2 = np.round(-dist, 2)                           # inverted transient
        corr = x2 - baseline_rolling_ball(x2, 20, 20)     # baseline correction

        valleys = find_peaks_simple(-corr)
        hills = find_peaks_simple(corr)
        if hills.size == 0 or valleys.size == 0:
            raise ValueError("no beats found")
        peak = np.array([None] * n, dtype=object)
        peak[hills - 1] = "Peak"
        peak[valleys - 1] = "Valley"

        avgt = AVGT if n >= AVGT else int(n / 3)
        rmv = np.convolve(corr, np.ones(avgt) / avgt, mode="valid") if avgt > 0 else corr.copy()
        # right-aligned rolling mean; the tail (incl. the last value) is the median of the last 6
        thr = np.empty(n)
        m = rmv.size
        thr[:m] = rmv
        thr[m - 1:] = np.median(rmv[-6:])

        vid = np.nonzero(peak == "Valley")[0]
        peak[vid[corr[vid] >= thr[vid]]] = None
        pid = np.nonzero(peak == "Peak")[0]
        peak[pid[corr[pid] <= thr[pid]]] = None

        all_pos = np.nonzero(peak != None)[0]              # noqa: E711 (object array)
        kinds = peak[all_pos]
        nk = kinds.size
        dbl = np.nonzero((kinds[:-1] == "Peak") & (kinds[1:] == "Peak"))[0] if nk > 1 else np.empty(0, int)
        cleaned = np.delete(all_pos, dbl + 1) if dbl.size else all_pos
        kinds = peak[cleaned]
        m3 = kinds.size
        if m3 >= 3:
            vpv = np.nonzero((kinds[:-2] == "Valley") & (kinds[1:-1] == "Peak") & (kinds[2:] == "Valley"))[0]
        else:
            vpv = np.empty(0, int)
        if vpv.size == 0:
            raise ValueError("no complete beats (valley-peak-valley) found")
        starts, peaks_, ends = cleaned[vpv], cleaned[vpv + 1], cleaned[vpv + 2]   # 0-based

        # diameters at start (end-diastolic) and peak (end-systolic) - from the uncorrected transient
        DD = x2[starts]
        SD = x2[peaks_]
        with np.errstate(divide="ignore", invalid="ignore"):
            FS = (DD - SD) / DD
        res.edd, res.esd, res.fs = r_median(DD), r_median(SD), r_median(FS)
        mean_edd = -res.edd

        ts, tp, te = time[starts], time[peaks_], time[ends]
        SI = te - ts
        dSI = np.append(np.diff(SI), NA)
        dSIp = np.full(SI.size, NA)
        dSIp[: SI.size - 1] = dSI[1:]
        HP = np.append(np.diff(ts), NA)
        res.intervals = pd.DataFrame({"starts": ts, "peaks": tp, "ends": te, "SI": SI, "deltaSI": dSI,
                                      "deltaSIplus": dSIp, "HP": HP, "DI": HP - SI})

        beats = []
        for b, (s0, e0) in enumerate(zip(starts, ends), start=1):
            t = time[s0:e0 + 1]
            d = x2[s0:e0 + 1]
            diameter = -d
            timestamp = t.copy()
            t = t - t[0]
            d = d - d[0]
            if np.max(d) != 0:
                d = d / np.max(d)
            delta = (diameter - diameter[0]) * -1
            from_edd = mean_edd - delta
            with np.errstate(divide="ignore", invalid="ignore"):
                vel = np.concatenate([[0.0], np.diff(delta) / np.diff(t)])
            nz = np.nonzero(vel != 0)[0]
            if nz.size:
                k0 = max(nz[0] - 1, 0)          # start one sample before the first movement
                sl = slice(k0, None)
                t, d, diameter, timestamp, delta, from_edd, vel = (
                    a[sl] for a in (t, d, diameter, timestamp, delta, from_edd, vel))
            t = t - t[0]
            beats.append(pd.DataFrame({"time": t, "distance": d, "diameter": diameter,
                                       "timestamp": timestamp, "delta": delta, "from_EDD": from_edd,
                                       "velocity": vel, "beat": b}))
        res.transients = pd.concat(beats, ignore_index=True)

        auc = np.array([trapz_unit(bt["distance"].to_numpy()) for bt in beats])
        lo, hi = np.nanmin(auc), np.nanmax(auc)
        identicals = [bt for bt, a in zip(beats, auc) if lo <= a <= hi and len(bt) > 10]
        if not identicals:
            raise ValueError("no usable beats (all shorter than 11 samples)")
        res.n_identicals = len(identicals)
        res.time_used = float(sum(bt["time"].max() for bt in identicals))
        for bt in identicals:
            v = bt["velocity"].to_numpy(float, copy=True)
            res.max_velocity.append(float(v[np.argmax(v)]))
            res.max_neg_velocity.append(float(v[np.argmin(v)]))
            try:
                res.metrics.append(spline_metrics(bt["time"].to_numpy(float, copy=True), bt["distance"].to_numpy(float, copy=True)))
            except Exception:
                res.metrics.append([NA] * len(METRIC_NAMES))
        return res
    except Exception as e:
        res.ok = False
        res.message = f"{type(e).__name__}: {e}"
        return res
