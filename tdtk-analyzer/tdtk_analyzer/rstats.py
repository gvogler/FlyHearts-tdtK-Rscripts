"""R-compatible numerical helpers.

The original pipeline was written in R. To reproduce its numbers, the helpers
here follow R's definitions exactly (1-based index conventions are converted
at the call sites):

* ``quantile`` type 7, ``sd``/``var`` with n-1, ``mad`` with constant 1.4826
* ``zoo::rollmean`` / ``rollapply`` ("valid" windows)
* ``pracma::trapz`` with unit spacing
* ``stats::smooth.spline`` (GCV, R's knot selection and its Brent spar search,
  ported from R's sbart.c)
* ``baseline::rollingBall`` (ported line by line from the baseline package)
* ``summary(lm(y ~ x))`` slope / adjusted R^2 as used by the script
"""

from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import cho_factor, cho_solve, LinAlgError

NA = float("nan")


# --------------------------------------------------------------------------- basics

def quantile7(x, probs) -> np.ndarray:
    """R quantile(x, probs) (type 7). NaNs are an error in R; here they are dropped."""
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if x.size == 0:
        return np.full(np.shape(probs), NA)
    return np.quantile(x, probs)


def r_sd(x) -> float:
    x = np.asarray(x, dtype=float)
    if x.size < 2 or np.isnan(x).any():
        return NA
    return float(np.std(x, ddof=1))


def r_sd_narm(x) -> float:
    x = np.asarray(x, dtype=float)
    return r_sd(x[~np.isnan(x)])


def r_mean(x) -> float:
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return NA
    return float(np.mean(x))


def r_mean_narm(x) -> float:
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    return float(np.mean(x)) if x.size else NA


def r_median(x) -> float:
    x = np.asarray(x, dtype=float)
    if x.size == 0 or np.isnan(x).any():
        return NA
    return float(np.median(x))


def r_median_narm(x) -> float:
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    return float(np.median(x)) if x.size else NA


def r_mad_narm(x) -> float:
    """R mad(x, na.rm = TRUE) (constant 1.4826)."""
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if x.size == 0:
        return NA
    m = np.median(x)
    return float(1.4826 * np.median(np.abs(x - m)))


def r_cor(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size < 2 or np.isnan(a).any() or np.isnan(b).any():
        return NA
    sa, sb = np.std(a), np.std(b)
    if sa == 0 or sb == 0:
        return NA
    return float(np.corrcoef(a, b)[0, 1])


def trapz_unit(y) -> float:
    """pracma::trapz(y): trapezoid rule with x = 1..n."""
    y = np.asarray(y, dtype=float)
    if y.size < 2:
        return 0.0
    return float(np.sum((y[:-1] + y[1:]) * 0.5))


def trapz_unit_cols(m) -> np.ndarray:
    """apply(m, 2, trapz) for a 2-D array (columns)."""
    m = np.asarray(m, dtype=float)
    if m.shape[0] < 2:
        return np.zeros(m.shape[1])
    return np.sum((m[:-1] + m[1:]) * 0.5, axis=0)


def rollmean(x, k: int, axis: int = -1) -> np.ndarray:
    """zoo::rollmean(x, k) / rollapply(x, k, mean): valid windows only."""
    x = np.asarray(x, dtype=float)
    if x.shape[axis] < k:
        shape = list(x.shape)
        shape[axis] = 0
        return np.empty(shape)
    xv = np.moveaxis(x, axis, -1)
    win = np.lib.stride_tricks.sliding_window_view(xv, k, axis=-1)
    return np.moveaxis(win.mean(axis=-1), -1, axis)


def rollmedian(x, k: int) -> np.ndarray:
    """rollapply(x, k, median) (valid windows)."""
    x = np.asarray(x, dtype=float)
    if x.size < k:
        return np.empty(0)
    win = np.lib.stride_tricks.sliding_window_view(x, k)
    return np.median(win, axis=-1)


def rcolon(a: int, b: int) -> np.ndarray:
    """R's a:b (may be descending)."""
    return np.arange(a, b + 1) if b >= a else np.arange(a, b - 1, -1)


def r_index(vec, idx_1based):
    """vec[idx] with R semantics for a vector of 1-based indices.

    zero indices are dropped, indices beyond the end give NaN, negative indices
    are not supported (raise)."""
    vec = np.asarray(vec, dtype=float)
    idx = np.asarray(idx_1based, dtype=float)
    idx = idx[idx != 0]
    if (idx < 0).any():
        raise IndexError("negative R index")
    out = np.full(idx.shape, NA)
    ok = ~np.isnan(idx) & (idx <= vec.size)
    out[ok] = vec[idx[ok].astype(int) - 1]
    return out


# --------------------------------------------------------------------------- peaks

def find_peaks_m(x, m: int) -> np.ndarray:
    """find_peaks()/findPeaks() of the script (stats.stackexchange 22974).

    Returns 1-based positions."""
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 3:
        return np.empty(0, dtype=int)
    shape = np.diff(np.sign(np.diff(x)))
    out = []
    for i in (np.nonzero(shape < 0)[0] + 1):  # 1-based i
        z = i - m + 1
        z = z if z > 0 else 1
        w = i + m + 1
        w = w if w < n else n
        idx = np.concatenate([rcolon(z, i), rcolon(i + 2, w)])
        if np.all(x[idx - 1] <= x[i]):  # x[i + 1] in R
            out.append(i + 1)
    return np.asarray(out, dtype=int)


def find_peaks_simple(x) -> np.ndarray:
    """findPeaksP / findPeaksV of the script: which(diff(sign(diff(x))) < 0) + 1 (1-based)."""
    x = np.asarray(x, dtype=float)
    if x.size < 3:
        return np.empty(0, dtype=int)
    d = np.diff(np.sign(np.diff(x)))
    return np.nonzero(d < 0)[0] + 2  # which() is 1-based -> +1, then +1


def which_max_n(x, n: int) -> int:
    """First index (1-based) of which.maxN(x, N) as used in the script."""
    x = np.asarray(x, dtype=float)
    length = x.size
    if n > length:
        n = length
    val = np.sort(x)[length - n]
    return int(np.nonzero(x == val)[0][0]) + 1


# --------------------------------------------------------------------------- lm

def lm_slope_adjr2(x, y):
    """z <- summary(lm(y ~ x)); return (z$coefficients[2], z$adj.r.squared).

    Replicates R including the degenerate cases used implicitly by the script:
    when the slope is not estimable, coefficients[2] is the standard error of
    the intercept (R's coefficient matrix only has the intercept row)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    n = x.size
    if n == 0:
        return NA, NA
    xm, ym = x.mean(), y.mean()
    sxx = float(np.sum((x - xm) ** 2))
    tss = float(np.sum((y - ym) ** 2))
    # aliased slope (all x equal or single point): intercept-only model
    if n == 1 or sxx <= 1e-7 * max(1.0, float(np.sum(x ** 2))):
        if n < 2:
            return NA, NA
        se_int = math.sqrt(tss / (n - 1)) / math.sqrt(n)
        return se_int, 0.0
    slope = float(np.sum((x - xm) * (y - ym)) / sxx)
    rss = float(np.sum((y - (ym + slope * (x - xm))) ** 2))
    r2 = 1.0 - rss / tss if tss > 0 else NA
    rdf = n - 2
    if rdf == 0:  # exact fit through 2 points: R gives NaN/-Inf, both fail "> 0.15"
        return slope, NA
    adj = 1.0 - (1.0 - r2) * ((n - 1) / rdf) if not math.isnan(r2) else NA
    return slope, adj


# --------------------------------------------------------------------------- baseline::rollingBall

def baseline_rolling_ball(y, wm: int, ws: int) -> np.ndarray:
    """baseline:::baseline.rollingBall for one spectrum; returns the baseline.

    Literal port of the R code (1-based indices kept via helper y1)."""
    y = np.asarray(y, dtype=float)
    n = y.size

    def s(v, a, b):  # v[a:b] in R (1-based, inclusive, may be descending)
        return v[rcolon(a, b) - 1]

    T1 = np.zeros(n)
    T2 = np.zeros(n)
    basel = np.zeros(n)

    # Minimize
    u1 = math.ceil((wm + 1) / 2) + 1
    T1[0] = np.min(y[:u1])
    for i in range(2, wm + 1):
        u2 = u1 + 1 + (i % 2)
        T1[i - 1] = min(np.min(s(y, u1 + 1, u2)), T1[i - 2])
        u1 = u2
    for i in range(wm + 1, n - wm + 1):
        if (y[u1] <= T1[i - 2]) and (y[u1 - wm - 1] != T1[i - 2]):
            T1[i - 1] = y[u1]
        else:
            T1[i - 1] = np.min(s(y, i - wm, i + wm))
        u1 += 1
    u1 = n - 2 * wm - 1
    for i in range(n - wm + 1, n + 1):
        u2 = u1 + 1 + ((i + 1) % 2)
        if np.min(s(y, u1, u2 - 1)) > T1[i - 2]:
            T1[i - 1] = T1[i - 2]
        else:
            T1[i - 1] = np.min(s(y, u2, n))
        u1 = u2

    # Maximize
    u1 = math.ceil((wm + 1) / 2) + 1
    T2[0] = np.max(T1[:u1])
    for i in range(2, wm + 1):
        u2 = u1 + 1 + (i % 2)
        T2[i - 1] = max(np.max(s(T1, u1 + 1, u2)), T2[i - 2])
        u1 = u2
    for i in range(wm + 1, n - wm + 1):
        if (T1[u1] >= T2[i - 2]) and (T1[u1 - wm - 1] != T2[i - 2]):
            T2[i - 1] = T1[u1]
        else:
            T2[i - 1] = np.max(s(T1, i - wm, i + wm))
        u1 += 1
    u1 = n - 2 * wm - 1
    for i in range(n - wm + 1, n + 1):
        u2 = u1 + 1 + ((i + 1) % 2)
        if np.max(s(T1, u1, u2 - 1)) < T2[i - 2]:
            T2[i - 1] = T2[i - 2]
        else:
            T2[i - 1] = np.max(s(T1, u2, n))
        u1 = u2

    # Smooth
    u1 = math.ceil(ws / 2)
    v = float(np.sum(T2[:u1]))
    for i in range(1, ws + 1):
        u2 = u1 + 1 + (i % 2)
        v = v + float(np.sum(s(T2, u1 + 1, u2)))
        basel[i - 1] = v / u2
        u1 = u2
    v = float(np.sum(T2[: 2 * ws + 1]))
    basel[ws] = v / (2 * ws + 1)
    for i in range(ws + 2, n - ws + 1):
        v = v - T2[i - ws - 2] + T2[i + ws - 1]
        basel[i - 1] = v / (2 * ws + 1)
    u1 = n - 2 * ws + 1
    v = v - T2[u1 - 1]
    basel[n - ws] = v / (2 * ws)
    for i in range(n - ws + 2, n + 1):
        u2 = u1 + 1 + (i + 1) % 2
        v = v - float(np.sum(s(T2, u1, u2 - 1)))
        basel[i - 1] = v / (n - u2 + 1)
        u1 = u2
    return basel


# --------------------------------------------------------------------------- smooth.spline

def _nknots_smspl(n: int) -> int:
    if n < 50:
        return n
    a1, a2, a3, a4 = math.log2(50), math.log2(100), math.log2(140), math.log2(200)
    if n < 200:
        v = 2 ** (a1 + (a2 - a1) * (n - 50) / 150)
    elif n < 800:
        v = 2 ** (a2 + (a3 - a2) * (n - 200) / 600)
    elif n < 3200:
        v = 2 ** (a3 + (a4 - a3) * (n - 800) / 2400)
    else:
        v = 200 + (n - 3200) ** 0.2
    return int(math.trunc(v))


def _r_seq_length(frm: float, to: float, length_out: int) -> np.ndarray:
    """seq.int(from, to, length.out = n)."""
    if length_out == 1:
        return np.array([frm])
    by = (to - frm) / (length_out - 1)
    out = frm + np.arange(length_out) * by
    out[-1] = to
    return out


class SmoothSpline:
    """Minimal port of R's smooth.spline() with default arguments (GCV).

    Only what the tdtK pipeline needs: fit(x, y) and predict(xnew)."""

    def __init__(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        n = x.size
        tol = 1e-6 * (np.quantile(x, 0.75) - np.quantile(x, 0.25))
        if tol > 0:
            xx = np.round((x - x.mean()) / tol)
        else:
            xx = x.copy()
        order = np.argsort(xx, kind="stable")
        xxs = xx[order]
        nd = np.concatenate([[True], xxs[:-1] < xxs[1:]])
        ux = x[order][nd]
        nx = ux.size
        if nx <= 3:
            raise ValueError("need at least four unique 'x' values")
        if nx == n:
            wbar = np.ones(nx)
            ybar = y[order]
            yssw = 0.0
        else:
            grp = np.cumsum(nd) - 1
            ys = y[order]
            wbar = np.bincount(grp).astype(float)
            sy = np.bincount(grp, weights=ys)
            sy2 = np.bincount(grp, weights=ys ** 2)
            ybar = sy / wbar
            yssw = float(np.sum(sy2 - wbar * ybar ** 2))
        self.ux = ux
        r_ux = ux[-1] - ux[0]
        xbar = (ux - ux[0]) / r_ux
        nknots = _nknots_smspl(nx)
        sel = np.trunc(_r_seq_length(1, nx, nknots)).astype(int) - 1
        knot = np.concatenate([np.repeat(xbar[0], 3), xbar[sel], np.repeat(xbar[-1], 3)])
        nk = nknots + 2
        self._x0 = ux[0]
        self._r = r_ux
        self._knot = knot

        eye = np.eye(nk)
        basis = BSpline(knot, eye, 3, extrapolate=False)
        X = basis(xbar)
        X = np.nan_to_num(X)
        # the last point sits on the right boundary: make sure the basis sums to 1 there
        X[xbar >= knot[-1]] = 0.0
        X[xbar >= knot[-1], nk - 1] = 1.0
        w = np.sqrt(wbar)
        Xw = X * (w ** 2)[:, None]
        XtX = X.T @ Xw
        Xty = Xw.T @ ybar

        # SIGMA = int B_i'' B_j'' dt. B'' is linear on each knot interval; R's
        # sgram.f integrates it as w*(a_i a_j + (a_i b_j + a_j b_i)/2 + b_i b_j * .3330)
        # (note .3330, not 1/3) - replicated here to get R's lambda scaling.
        d2 = BSpline(knot, eye, 3, extrapolate=True).derivative(2)
        inner = np.unique(knot)
        left, right = inner[:-1], inner[1:]
        wpt = right - left
        ya = d2(left)
        yb = d2(right) - ya
        sigma = (np.einsum("k,ki,kj->ij", wpt, ya, ya)
                 + 0.5 * (np.einsum("k,ki,kj->ij", wpt, yb, ya) + np.einsum("k,ki,kj->ij", wpt, ya, yb))
                 + 0.3330 * np.einsum("k,ki,kj->ij", wpt, yb, yb))

        t1 = float(np.sum(np.diag(XtX)[2:nk - 3]))
        t2 = float(np.sum(np.diag(sigma)[2:nk - 3]))
        ratio = t1 / t2

        self._fit_args = (X, XtX, Xty, sigma, ybar, w, yssw, ratio)
        self.coef = self._search_spar()

    # sslvrg for a given spar: returns (coef, crit)
    def _sslvrg(self, spar: float):
        X, XtX, Xty, sigma, ybar, w, yssw, ratio = self._fit_args
        lam = ratio * 16.0 ** (spar * 6.0 - 2.0)
        A = XtX + lam * sigma
        try:
            cf = cho_factor(A, lower=False, check_finite=False)
            coef = cho_solve(cf, Xty, check_finite=False)
            Ainv_Xt = cho_solve(cf, (X * (w ** 2)[:, None]).T, check_finite=False)
        except (LinAlgError, ValueError):
            return None, 2e100
        sz = X @ coef
        lev = np.einsum("ij,ji->i", X, Ainv_Xt)
        r = (ybar - sz) * w
        rss = yssw + float(np.sum(r * r))
        df = float(np.sum(lev))
        sumw = float(np.sum(w * w))
        crit = (rss / sumw) / (1.0 - df / sumw) ** 2
        if not math.isfinite(crit):
            crit = 2e100
        return coef, crit

    def _search_spar(self):
        # Forsythe-Malcolm-Moler / Brent fmin exactly as in R's sbart.c
        c_gold = 0.381966011250105151795413165634
        big = 1e100
        tol, eps, maxit = 1e-4, 2e-8, 500
        a, b = -1.5, 1.5
        v = a + c_gold * (b - a)
        w = x = v
        d = e = 0.0
        coef, fx = self._sslvrg(x)
        last_coef = coef
        fv = fw = fx
        it = 0
        while True:
            xm = (a + b) * 0.5
            tol1 = eps * abs(x) + tol / 3.0
            tol2 = tol1 * 2.0
            it += 1
            if abs(x - xm) <= tol2 - (b - a) * 0.5 or it > maxit:
                break
            golden = True
            if not (abs(e) <= tol1 or fx >= big or fv >= big or fw >= big):
                r = (x - w) * (fx - fv)
                q = (x - v) * (fx - fw)
                p = (x - v) * q - (x - w) * r
                q = (q - r) * 2.0
                if q > 0.0:
                    p = -p
                q = abs(q)
                r = e
                e = d
                if not (abs(p) >= abs(0.5 * q * r) or q == 0.0) and not (p <= q * (a - x) or p >= q * (b - x)):
                    d = p / q
                    u = x + d
                    if u - a < tol2 or b - u < tol2:
                        d = math.copysign(tol1, xm - x)
                    golden = False
            if golden:
                e = (a - x) if x >= xm else (b - x)
                d = c_gold * e
            u = x + (d if abs(d) >= tol1 else math.copysign(tol1, d))
            coef, fu = self._sslvrg(u)
            if coef is not None:
                last_coef = coef
            if not math.isfinite(fu):
                fu = 2.0 * big
            if fu <= fx:
                if u >= x:
                    a = x
                else:
                    b = x
                v, fv = w, fw
                w, fw = x, fx
                x, fx = u, fu
            else:
                if u < x:
                    a = u
                else:
                    b = u
                if fu <= fw or w == x:
                    v, fv = w, fw
                    w, fw = u, fu
                elif fu <= fv or v == x or v == w:
                    v, fv = u, fu
        self.spar = x
        self.crit = fx
        # like R, the coefficients are those of the last evaluation
        return last_coef

    def predict(self, xnew) -> np.ndarray:
        xnew = np.asarray(xnew, dtype=float)
        xs = (xnew - self._x0) / self._r
        spl = BSpline(self._knot, self.coef, 3, extrapolate=False)
        out = spl(np.clip(xs, 0.0, 1.0))
        # linear extrapolation outside [0, 1] like predict.smooth.spline
        lo, hi = xs < 0, xs > 1
        if lo.any() or hi.any():
            d1 = spl.derivative(1)
            if lo.any():
                out[lo] = spl(0.0) + d1(0.0) * xs[lo]
            if hi.any():
                out[hi] = spl(1.0) + d1(1.0) * (xs[hi] - 1.0)
        # right boundary: BSpline with extrapolate=False returns nan at exactly 1.0
        at1 = xs == 1.0
        if at1.any():
            out[at1] = self.coef[-1]
        return out
