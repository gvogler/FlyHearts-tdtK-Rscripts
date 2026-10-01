"""R-compatibility of the numerical helpers (reference values produced by real R, see data/)."""

import json
import os

import numpy as np
import pytest

from tdtk_analyzer.rstats import (
    SmoothSpline, baseline_rolling_ball, find_peaks_m, find_peaks_simple, lm_slope_adjr2, quantile7,
    r_mad_narm, rollmean, which_max_n,
)

REF = json.load(open(os.path.join(os.path.dirname(__file__), "data", "r_reference.json")))


@pytest.mark.parametrize("case", REF[:3])
def test_baseline_rolling_ball_matches_r(case):
    b = baseline_rolling_ball(np.array(case["x"]), 20, 20)
    np.testing.assert_allclose(b, case["b"], rtol=0, atol=1e-12)


@pytest.mark.parametrize("case", REF[3:10])
def test_smooth_spline_matches_r(case):
    s = SmoothSpline(case["t"], case["y"])
    assert s.spar == pytest.approx(case["spar"], abs=1e-6)
    n = int(np.floor(max(s.ux) / 0.0002 + 1e-10)) + 1
    pred = s.predict(np.arange(n) * 0.0002)
    idx = np.array(case["idx"], dtype=int) - 1
    np.testing.assert_allclose(pred[idx], case["pred"], rtol=0, atol=1e-8)


def test_lm_matches_r():
    assert lm_slope_adjr2([3, 5, 7, 9, np.nan, 11], [1, 2, 3, 4, 5, 6]) == pytest.approx(
        (REF[10]["slope"], REF[10]["adj"]))
    slope, adj = lm_slope_adjr2([4, 4, 4], [1, 2, 3])        # aliased slope -> SE of intercept
    assert slope == pytest.approx(REF[11]["slope"]) and adj == 0


def test_basic_r_semantics():
    assert quantile7([1, 2, 3, 4, 10], 0.75) == 4                     # type 7
    assert r_mad_narm([1, 2, 3, 4, 100, np.nan]) == pytest.approx(1.4826 * 1)
    np.testing.assert_allclose(rollmean([1, 2, 3, 4, 5], 2), [1.5, 2.5, 3.5, 4.5])
    assert list(find_peaks_simple([0, 1, 0, 2, 0])) == [2, 4]         # 1-based
    assert which_max_n([5, 9, 7, 9], 2) == 2                          # 2nd largest value is 9 (tie)
    x = np.array([0, 1, 2, 5, 2, 1, 0, 0, 3, 0], float)
    assert list(find_peaks_m(x, 3)) == [4, 9]                      # same as R
