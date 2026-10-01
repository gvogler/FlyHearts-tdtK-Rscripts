"""The vectorized rolling ball must equal ImageJ's loop (BackgroundSubtracter.rollBall)."""

import numpy as np

from tdtk_analyzer.background import RollingBall, _roll_ball, subtract_background_8bit


def java_roll_ball(pixels: np.ndarray, ball: RollingBall) -> np.ndarray:
    """Literal port of ImageJ's rollBall() (slow, for testing only)."""
    height, width = pixels.shape
    src = pixels.astype(np.float32).ravel()
    out = np.full(src.size, -np.finfo(np.float32).max, dtype=np.float32)
    z_ball = ball.data.ravel()
    bw = ball.width
    radius = bw // 2
    for y in range(-radius, height + radius):
        y0 = max(y - radius, 0)
        y_ball0 = y0 - y + radius
        yend = min(y + radius, height - 1)
        for x in range(-radius, width + radius):
            z = np.finfo(np.float32).max
            x0 = max(x - radius, 0)
            x_ball0 = x0 - x + radius
            xend = min(x + radius, width - 1)
            for yp, yb in zip(range(y0, yend + 1), range(y_ball0, bw)):
                for xp, bp in zip(range(x0, xend + 1), range(x_ball0 + yb * bw, bw * bw)):
                    zr = src[xp + yp * width] - z_ball[bp]
                    if z > zr:
                        z = zr
            for yp, yb in zip(range(y0, yend + 1), range(y_ball0, bw)):
                for xp, bp in zip(range(x0, xend + 1), range(x_ball0 + yb * bw, bw * bw)):
                    zmin = z + z_ball[bp]
                    p = xp + yp * width
                    if out[p] < zmin:
                        out[p] = zmin
    return out.reshape(height, width)


def test_roll_ball_matches_imagej_loop():
    rng = np.random.default_rng(1)
    img = (rng.random((14, 23)) * 200).astype(np.float32)
    img[5:9, 8:15] += 60
    ball = RollingBall(8)          # no shrinking for radius <= 10
    np.testing.assert_allclose(_roll_ball(img, ball), java_roll_ball(img, ball), rtol=0, atol=1e-4)


def test_ball_geometry_radius_50():
    b = RollingBall(50)
    assert b.shrink == 4
    assert b.width == 19           # halfWidth = round(12.5 - 4) = 9 (Java Math.round)


def test_subtract_background_keeps_bright_band():
    img = np.full((40, 300), 30, np.uint8)
    img[15:25, :] = 200
    out = subtract_background_8bit(img, 50)
    assert out.dtype == np.uint8
    assert out[20].mean() > 120 and out[2].mean() < 10


def java_shrink(pixels, f):
    h, w = pixels.shape
    sh, sw = (h + f - 1) // f, (w + f - 1) // f
    out = np.empty((sh, sw), np.float32)
    for ys in range(sh):
        for xs in range(sw):
            out[ys, xs] = pixels[f * ys:min(f * ys + f, h), f * xs:min(f * xs + f, w)].min()
    return out


def java_enlarge(small, shape, f):
    h, w = shape
    sh, sw = small.shape

    def arrays(length, small_length):
        idx, wts = [], []
        for i in range(length):
            si = int((i - f // 2) / f)
            if si >= small_length - 1:
                si = small_length - 2
            idx.append(si)
            wts.append(np.float32(1.0) - np.float32((i + 0.5) / f - (si + 0.5)))
        return idx, wts

    xi, xw = arrays(w, sw)
    yi, yw = arrays(h, sh)
    out = np.empty((h, w), np.float32)
    for y in range(h):
        for x in range(w):
            l0 = small[yi[y], xi[x]] * xw[x] + small[yi[y], xi[x] + 1] * (1 - xw[x])
            l1 = small[yi[y] + 1, xi[x]] * xw[x] + small[yi[y] + 1, xi[x] + 1] * (1 - xw[x])
            out[y, x] = l0 * yw[y] + l1 * (1 - yw[y])
    return out


def test_radius_50_path_matches_imagej():
    from tdtk_analyzer.background import _filter3_mean, rolling_ball_background

    rng = np.random.default_rng(2)
    img = (rng.random((40, 90)) * 60 + 20).astype(np.float32)
    img[14:26, :] += 150
    ball = RollingBall(50)
    sm = _filter3_mean(img)
    expected = java_enlarge(java_roll_ball(java_shrink(sm, 4), ball), img.shape, 4)
    np.testing.assert_allclose(rolling_ball_background(img, 50), expected, rtol=0, atol=1e-3)
