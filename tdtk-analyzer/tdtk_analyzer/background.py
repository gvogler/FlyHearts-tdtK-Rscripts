"""ImageJ "Process > Subtract Background..." (rolling ball), replacing the Fiji macro.

The Fiji macro ran ``run("Subtract Background...", "rolling=50")`` on every
8-bit kymograph. In macro mode that means: rolling ball (not paraboloid),
dark background, 3x3 mean pre-smoothing, subtract. This is a port of
ij/plugin/filter/BackgroundSubtracter.java for exactly that case; the
ball-rolling itself is a grey-scale opening with the ball as a non-flat
structuring element, which gives the same result as ImageJ's loop.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import ndimage


class RollingBall:
    def __init__(self, radius: float):
        if radius <= 10:
            self.shrink, trim = 1, 24
        elif radius <= 30:
            self.shrink, trim = 2, 24
        elif radius <= 100:
            self.shrink, trim = 4, 32
        else:
            self.shrink, trim = 8, 40
        small = radius / self.shrink
        if small < 1:
            small = 1
        rsq = small * small
        xtrim = int(int(trim * small) / 100)          # (int)(arcTrimPer*r)/100 in Java
        half = int(math.floor(small - xtrim + 0.5))    # Math.round
        self.width = 2 * half + 1
        yy, xx = np.mgrid[-half:half + 1, -half:half + 1]
        temp = rsq - xx * xx - yy * yy
        self.data = np.where(temp > 0, np.sqrt(np.maximum(temp, 0)), 0.0).astype(np.float32)


def _filter3_mean(img: np.ndarray) -> np.ndarray:
    """filter3x3(MEAN): 1-D 3-point means along rows, then columns (edges replicated)."""
    out = img.astype(np.float32).copy()
    for axis in (1, 0):
        p = np.pad(out, [(1, 1) if a == axis else (0, 0) for a in range(2)], mode="edge")
        a = np.take(p, range(0, p.shape[axis] - 2), axis=axis)
        b = np.take(p, range(1, p.shape[axis] - 1), axis=axis)
        c = np.take(p, range(2, p.shape[axis]), axis=axis)
        out = ((a + b + c) * np.float32(0.333333333)).astype(np.float32)
    return out


def _shrink(img: np.ndarray, f: int) -> np.ndarray:
    h, w = img.shape
    sh, sw = (h + f - 1) // f, (w + f - 1) // f
    pad = np.full((sh * f, sw * f), np.inf, dtype=np.float32)
    pad[:h, :w] = img
    return pad.reshape(sh, f, sw, f).min(axis=(1, 3))


def _roll_ball(img: np.ndarray, ball: RollingBall) -> np.ndarray:
    r = ball.width // 2
    se = ball.data
    padded = np.pad(img, r, mode="constant", constant_values=np.inf)
    # z(center) = min over patch (pixel - ball); ball centres may lie outside the image
    z = ndimage.grey_erosion(padded, structure=se, mode="constant", cval=np.inf)
    # background(pixel) = max over centres (z + ball)
    zp = np.pad(z, r, mode="constant", constant_values=-np.inf)
    bg = ndimage.grey_dilation(zp, structure=se, mode="constant", cval=-np.inf)
    return bg[2 * r:-2 * r, 2 * r:-2 * r].astype(np.float32)


def _interp_arrays(length: int, small_length: int, f: int):
    idx = np.empty(length, dtype=int)
    wts = np.empty(length, dtype=np.float32)
    for i in range(length):
        si = int((i - f // 2) / f)  # Java int division truncates toward zero
        if si >= small_length - 1:
            si = small_length - 2
        idx[i] = si
        dist = (i + 0.5) / f - (si + 0.5)
        wts[i] = 1.0 - dist
    return idx, wts


def _enlarge(small: np.ndarray, shape: tuple, f: int) -> np.ndarray:
    h, w = shape
    sh, sw = small.shape
    xi, xw = _interp_arrays(w, sw, f)
    yi, yw = _interp_arrays(h, sh, f)
    lines = small[:, xi] * xw[None, :] + small[:, xi + 1] * (1 - xw[None, :])
    return (lines[yi, :] * yw[:, None] + lines[yi + 1, :] * (1 - yw[:, None])).astype(np.float32)


def rolling_ball_background(img: np.ndarray, radius: float = 50.0, presmooth: bool = True) -> np.ndarray:
    """Background estimated by ImageJ's rolling ball (dark background)."""
    ball = RollingBall(radius)
    fp = img.astype(np.float32)
    if presmooth:
        fp = _filter3_mean(fp)
    shrink = ball.shrink > 1
    small = _shrink(fp, ball.shrink) if shrink else fp
    if shrink and min(small.shape) < 2:
        shrink = False
        small = fp
    bg = _roll_ball(small, ball)
    if shrink:
        bg = _enlarge(bg, fp.shape, ball.shrink)
    return bg


def subtract_background_8bit(img: np.ndarray, radius: float = 50.0) -> np.ndarray:
    """Subtract Background (rolling=radius) on an 8-bit image, returns uint8."""
    bg = rolling_ball_background(img, radius)
    value = img.astype(np.float32) - bg + np.float32(0.5)
    return np.clip(value, 0, 255).astype(np.uint8)
