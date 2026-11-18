"""Step 2: background subtraction, edge tracing of kymographs and quality control.

Ports of the Fiji macro (see background.py), "MAYO Screen Script No.2"
(tracing) and "No.3" (quality control) of tdtK_Full_Analysis_script.

A kymograph image is [Y, T] (rows = position across the heart, columns = time),
which is what R's ``image_`` holds after ``t()``.
"""

from __future__ import annotations

import math
import os
import shutil

import numpy as np
import pandas as pd
from scipy import ndimage

from . import rio
from .background import subtract_background_8bit
from .kymograph import StepResult
from .rstats import NA, quantile7, r_cor, r_sd, r_median_narm, rollmean, trapz_unit_cols


# ----------------------------------------------------------------------- Fiji step

def subtract_background_file(src: str, dst: str, radius: float = 50.0) -> StepResult:
    import tifffile

    name = os.path.basename(src)
    try:
        img = tifffile.imread(src)
        if img.ndim != 2:
            shutil.copyfile(src, dst)
            return StepResult(name, "skipped", "not a grey image; copied unchanged")
        if img.dtype != np.uint8:
            # the kymographs are written as 8-bit; scale anything else to 8-bit first
            img = np.clip(img.astype(float) / np.iinfo(img.dtype).max * 255 + 0.5, 0, 255).astype(np.uint8)
        tifffile.imwrite(dst, subtract_background_8bit(img, radius), photometric="minisblack")
        return StepResult(name, "done")
    except Exception as e:
        return StepResult(name, "error", f"{type(e).__name__}: {e}")


# ----------------------------------------------------------------------- tracing

def _gblur(img: np.ndarray, sigma: float = 1.5) -> np.ndarray:
    """EBImage::gblur(x, sigma): gaussian brush of size 2*ceiling(3*sigma)+1, circular boundary."""
    size = 2 * math.ceil(3 * sigma) + 1
    half = size // 2
    ax = np.arange(-half, half + 1)
    k = np.exp(-(ax[:, None] ** 2 + ax[None, :] ** 2) / (2 * sigma * sigma))
    k /= k.sum()
    return ndimage.convolve(img, k, mode="wrap")


def _first_above(col_mask: np.ndarray) -> np.ndarray:
    """which(col)[1] (1-based) for every column; NaN when none."""
    any_ = col_mask.any(axis=0)
    idx = np.argmax(col_mask, axis=0).astype(float) + 1
    idx[~any_] = np.nan
    return idx


def _nth_above(col: np.ndarray, thr: float, g: int) -> float:
    w = np.nonzero(col > thr)[0]
    return float(w[g - 1] + 1) if g <= w.size else np.nan


def trace_kymograph(tiff_name: str, folder: str) -> StepResult:
    """Trace one background-subtracted kymograph; writes '<tiff>.csv' and '<tiff>_traced.jpg'."""
    path = os.path.join(folder, tiff_name)
    try:
        if os.path.exists(path + "_traced.jpg"):
            return StepResult(tiff_name, "skipped", "already traced")
        img = rio.read_tiff_gray(path)
        if img is None:
            return StepResult(tiff_name, "skipped", "not a grey image")
        image_ = img / img.max()
        image_ = _gblur(image_, 1.5)
        ny = image_.shape[0]
        p = np.linspace(0.0, math.pi, ny) if ny > 1 else np.zeros(1)
        image_ = image_ * np.sin(p)[:, None]

        rm = rollmean(image_, 4, axis=1)                  # [Y, T-3]
        rm = rm / rm.max(axis=0, keepdims=True)           # normalize every time point
        background = quantile7(rm.ravel(), 0.0)
        bg = rm - background                              # image_mean_bg [Y, T-3]

        col_mean = bg.mean(axis=0)
        above = bg > col_mean[None, :]
        up = _first_above(above)
        down = ny - _first_above(above[::-1, :])
        positions = np.arange(1, bg.shape[1] + 1)

        # sanity check: re-pick boundaries that are far off
        limits = np.nanmean(up) - r_sd(up) if not np.isnan(up).any() else np.nan
        tbc = np.nonzero(up <= limits)[0]
        g = 1
        while tbc.size > 1:
            tbc = np.nonzero(up <= limits)[0]
            g += 1
            for i in tbc:
                up[i] = _nth_above(bg[:, i], np.mean(bg[:, i]), g)
        limits2 = np.nanmean(down) + r_sd(down) if not np.isnan(down).any() else np.nan
        tbc = np.nonzero(down >= limits2)[0]
        g = 1
        while tbc.size > 1:
            tbc = np.nonzero(down >= limits2)[0]
            g += 1
            for i in tbc:
                col = bg[::-1, i]
                down[i] = ny - _nth_above(col, np.mean(bg[:, i]), g)

        if np.isnan(up).any() or np.isnan(down).any():
            return StepResult(tiff_name, "error", "could not trace the heart edges")

        meta = rio.read_meta_csv(os.path.join(folder, rio.meta_csv_for(tiff_name)))
        time_interval = rio.meta_float(meta, "time_interval")
        resolution = rio.meta_float(meta, "resolution")

        n = up.size
        distance = (down - up) * resolution
        time = np.concatenate([[0.0], time_interval * np.arange(1, n)])

        up_mid = int(np.max(up) - r_sd(up))        # as.integer truncates
        down_mid = int(np.max(down) - r_sd(down))
        rows_up = np.arange(1, up_mid + 1) if up_mid >= 1 else np.array([1])
        rows_down = np.arange(down_mid, ny + 1) if down_mid <= ny else np.arange(ny, down_mid + 1)
        rows_down = rows_down[(rows_down >= 1) & (rows_down <= ny)]
        up_bright = trapz_unit_cols(bg[rows_up - 1, :])
        down_bright = trapz_unit_cols(bg[rows_down - 1, :])

        with np.errstate(invalid="ignore", divide="ignore"):
            down_norm = (down - down.min()) / (down.max() - down.min())
            up_norm = (up - up.min()) / (up.max() - up.min())
        out = pd.DataFrame({
            "time": time, "distance": distance, "up": up, "down": down,
            "up_bright": up_bright, "down_bright": down_bright,
            "down_norm": down_norm, "up_norm": up_norm, "distance_norm": down_norm - up_norm,
        })
        rio.write_r_csv(out, path + ".csv")

        # overlay image: red/blue = smoothed kymograph, green = traced edges (flipped like R)
        edges = np.zeros_like(bg)
        for arr in (up, down):
            ok = (arr >= 1) & (arr <= ny)
            edges[arr[ok].astype(int) - 1, positions[ok] - 1] = 1.0
        red = rm[::-1, :]
        rio.write_rgb_jpeg(red, edges[::-1, :], red, path + "_traced.jpg")
        return StepResult(tiff_name, "done")
    except Exception as e:
        return StepResult(tiff_name, "error", f"{type(e).__name__}: {e}")


# ----------------------------------------------------------------------- quality control

def _autocorr(v: np.ndarray) -> float:
    d = np.diff(v)
    return r_cor(d[:-1], d[1:])


def quality_control(folder: str) -> pd.DataFrame:
    """Sort traced kymographs into 'excellent', 'good' and 'bad traces' (R Script No.3)."""
    csvs = rio.list_files(folder, r"\.tiff.csv$")
    rows = []
    for idx, name in enumerate(csvs, start=1):
        r = {"csv": name, "smoothie_up": float(idx), "smoothie_up_norm": float(idx),
             "smoothie_down": float(idx), "smoothie_down_norm": float(idx), "sd": float(idx),
             "dist_cor": float(idx), "norm_dist_cor": float(idx),
             "autocorr_up": NA, "autocorr_down": NA, "autocorr_sd": NA}
        try:
            x = rio.read_r_csv(os.path.join(folder, name))
        except Exception:
            x = pd.DataFrame()
        if x.shape[1] >= 4:
            up, down = x["up"].to_numpy(float, copy=True), x["down"].to_numpy(float, copy=True)
            sdu, sdd = r_sd(np.diff(up)), r_sd(np.diff(down))
            if sdd != 0 and sdu != 0 and not (math.isnan(sdd) or math.isnan(sdu)):
                r["smoothie_up"], r["smoothie_down"] = sdu, sdd
                r["autocorr_up"], r["autocorr_down"] = _autocorr(up), _autocorr(down)
                r["autocorr_sd"] = r_sd([r["autocorr_up"], r["autocorr_down"]])
                r["sd"] = r_sd([sdd, sdu])
                r["dist_cor"] = _autocorr(x["distance"].to_numpy(float, copy=True))
                r["norm_dist_cor"] = _autocorr(x["distance_norm"].to_numpy(float, copy=True))
            else:
                r.update(smoothie_up=10.0, smoothie_down=10.0, autocorr_up=0.0, autocorr_down=0.0,
                         autocorr_sd=10.0, sd=100.0, dist_cor=0.0, norm_dist_cor=0.0)
        rows.append(r)
    qc = pd.DataFrame(rows, columns=["csv", "smoothie_up", "smoothie_up_norm", "smoothie_down",
                                     "smoothie_down_norm", "sd", "dist_cor", "norm_dist_cor",
                                     "autocorr_up", "autocorr_down", "autocorr_sd"])
    if qc.empty:
        return qc
    names = qc["csv"]
    qc["CODE"] = names.str.split("_").str[0]
    qc["ID"] = names.str.split("_").str[1]
    qc["Xpos"] = pd.to_numeric(names.str.extract(r"Xpos_(.*?)\.tiff", expand=False), errors="coerce")
    au, ad, sd = qc["autocorr_up"], qc["autocorr_down"], qc["sd"]
    quality = np.array(["n.t."] * len(qc), dtype=object)
    with np.errstate(invalid="ignore"):
        exc = ((au > r_median_narm(au) * 0.9) & (ad > r_median_narm(ad) * 0.9)
               & (sd < r_median_narm(sd) * 1.5)).to_numpy()
    quality[exc] = "excellent"
    good = (quality == "n.t.") & (qc["dist_cor"] > 0.6).to_numpy() & (sd < 4).to_numpy() \
        & (qc["norm_dist_cor"] > 0.6).to_numpy()
    quality[good] = "good"
    qc["quality"] = quality
    qc["tiff"] = names.str.replace(r"(\.tiff).*$", r"\1", regex=True)
    qc["jpeg"] = qc["tiff"] + "_traced.jpg"

    dirs = {"n.t.": "bad traces", "good": "good traces", "excellent": "excellent traces"}
    for q, d in dirs.items():
        target = os.path.join(folder, d)
        os.makedirs(target, exist_ok=True)
        for _, row in qc[qc["quality"] == q].iterrows():
            for f in (row["csv"], row["jpeg"]):
                src = os.path.join(folder, f)
                if os.path.exists(src) and not os.path.exists(os.path.join(target, f)):
                    shutil.copy(src, target)
    rio.write_r_csv(qc, os.path.join(folder, "Quality_control.csv"))

    # movies without any excellent trace: rescue their 'good' traces
    qc["cxd"] = qc["jpeg"].map(rio.movie_stem)
    counts = qc.groupby(["cxd", "quality"]).size().unstack(fill_value=0)
    if "excellent" in counts and "good" in counts:
        need = counts.index[(counts["excellent"] == 0) & (counts["good"] > 0)]
        rescue = qc[(qc["quality"] == "good") & qc["cxd"].isin(need)]
        target = os.path.join(folder, dirs["excellent"])
        for _, row in rescue.iterrows():
            for f in (row["csv"], row["jpeg"]):
                src = os.path.join(folder, f)
                if os.path.exists(src) and not os.path.exists(os.path.join(target, f)):
                    shutil.copy(src, target)

    # traces the user added to / removed from 'excellent traces' on the Review tab
    from .curation import apply_decisions
    apply_decisions(folder)
    return qc
