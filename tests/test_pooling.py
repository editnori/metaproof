"""Pooling: hand-worked small examples plus a published benchmark.

Benchmark: the BCG vaccine trials (Colditz GA et al. JAMA 1994;271(9):698-702), log risk
ratios, as analysed in Viechtbauer W. Conducting meta-analyses in R with the metafor package.
J Stat Softw 2010;36(3):1-48. This is a published methods benchmark, not NR data."""

import math

import numpy as np
import pytest
from scipy import optimize, stats

from metaproof import pooling as pl

# --- hand-worked: equal within-study variances -------------------------------------
# y = (1, 2, 4), v = 1: ybar = 7/3, sum (y - ybar)^2 = 16/9 + 1/9 + 25/9 = 42/9.
# With equal v every estimator has a closed form:
#   DL = PM = REML: tau2 = (42/9)/2 - 1 = 4/3;   ML: tau2 = (42/9)/3 - 1 = 5/9.
Y3 = np.array([1.0, 2.0, 4.0])
V3 = np.ones(3)


def test_fixed_effect_hand():
    f = pl.fit(Y3, V3, method="FE", hksj=True)
    assert f.estimate == pytest.approx(7 / 3)
    assert f.se[0] == pytest.approx(math.sqrt(1 / 3))
    assert f.QE == pytest.approx(42 / 9)
    assert f.I2 == pytest.approx(100 * (42 / 9 - 2) / (42 / 9))  # 57.142857
    assert f.H2 == pytest.approx(21 / 9)
    assert not f.hksj  # never applied to FE
    assert f.ci_lb[0] == pytest.approx(7 / 3 - stats.norm.ppf(0.975) * math.sqrt(1 / 3))


@pytest.mark.parametrize("method,tau2", [("DL", 4 / 3), ("PM", 4 / 3), ("REML", 4 / 3),
                                         ("ML", 5 / 9)])
def test_tau2_equal_variances(method, tau2):
    f = pl.fit(Y3, V3, method=method, hksj=False)
    assert f.tau2 == pytest.approx(tau2, abs=1e-9)
    assert f.estimate == pytest.approx(7 / 3)
    assert f.se[0] == pytest.approx(math.sqrt((1 + tau2) / 3))


def test_i2_h2_random_effects_hand():
    # typical within variance s2 = (k-1) sum w / ((sum w)^2 - sum w^2) = 2*3/(9-3) = 1
    f = pl.fit(Y3, V3, method="REML", hksj=False)
    assert f.I2 == pytest.approx(100 * (4 / 3) / (4 / 3 + 1))
    assert f.H2 == pytest.approx(7 / 3)


def test_hksj_hand():
    # w* = 1/(1 + 4/3) = 3/7; sum w* r^2 = (3/7)(42/9) = 2; s2 = 2/(k-1) = 1
    # Var_HK = s2 / sum w* = 7/9; CI uses t_2
    f = pl.fit(Y3, V3, method="REML", hksj=True)
    assert f.s2_hk == pytest.approx(1.0)
    assert f.se[0] == pytest.approx(math.sqrt(7 / 9))
    assert f.df == 2
    half = 4.302652729749462 * math.sqrt(7 / 9)
    assert (f.ci_lb[0], f.ci_ub[0]) == pytest.approx((7 / 3 - half, 7 / 3 + half))
    assert f.pval[0] == pytest.approx(2 * stats.t.sf((7 / 3) / math.sqrt(7 / 9), 2))


def test_prediction_interval_hand():
    # mu +/- t_{1, .975} sqrt(tau2 + SE^2) = 7/3 +/- 12.7062047 sqrt(4/3 + 7/9)
    f = pl.fit(Y3, V3, method="REML", hksj=True)
    lb, ub = pl.prediction_interval(f)
    half = 12.706204736174694 * math.sqrt(4 / 3 + 7 / 9)
    assert half == pytest.approx(18.46168746697446)
    assert (lb, ub) == pytest.approx((7 / 3 - half, 7 / 3 + half))
    # Undefined for k < 3 and for FE
    assert all(map(math.isnan, pl.prediction_interval(pl.fit(Y3[:2], V3[:2], method="REML"))))
    assert all(map(math.isnan, pl.prediction_interval(pl.fit(Y3, V3, method="FE"))))


# --- hand-worked: unequal variances, DL ----------------------------------------------
# y = (0, 1, 3), v = (1/4, 1/2, 1): w = (4, 2, 1), mu_FE = 5/7,
# Q = 4(5/7)^2 + 2(2/7)^2 + (16/7)^2 = 364/49, C = 7 - 21/7 = 4, tau2_DL = (364/49 - 2)/4 = 19/14

def test_dl_unequal_hand():
    y, v = np.array([0.0, 1.0, 3.0]), np.array([0.25, 0.5, 1.0])
    fe = pl.fit(y, v, method="FE")
    assert fe.estimate == pytest.approx(5 / 7)
    assert fe.QE == pytest.approx(364 / 49)
    f = pl.fit(y, v, method="DL", hksj=False)
    assert f.tau2 == pytest.approx(19 / 14)
    w = 1 / (v + 19 / 14)
    assert f.estimate == pytest.approx(np.sum(w * y) / np.sum(w))
    assert f.estimate == pytest.approx(1.1427590940288264)
    assert f.se[0] == pytest.approx(0.79431996715221)


def test_pm_solves_estimating_equation():
    y, v = np.array([0.0, 1.0, 3.0, 0.5]), np.array([0.25, 0.5, 1.0, 0.3])
    f = pl.fit(y, v, method="PM")
    w = 1 / (v + f.tau2)
    mu = np.sum(w * y) / np.sum(w)
    assert np.sum(w * (y - mu) ** 2) == pytest.approx(len(y) - 1)


@pytest.mark.parametrize("reml", [True, False])
def test_ml_reml_maximise_likelihood(reml):
    y, v = np.array([0.0, 1.0, 3.0, 0.5, -0.4]), np.array([0.25, 0.5, 1.0, 0.3, 0.2])
    f = pl.fit(y, v, method="REML" if reml else "ML")
    assert f.converged
    brute = optimize.minimize_scalar(lambda t: -pl.loglik(y, v, np.ones((5, 1)), t, reml),
                                     bounds=(0, 50), method="bounded",
                                     options={"xatol": 1e-12})
    assert f.tau2 == pytest.approx(brute.x, abs=1e-6)
    # Fixed-point forms (Viechtbauer 2005):
    w = 1 / (v + f.tau2)
    mu = np.sum(w * y) / np.sum(w)
    fp = np.sum(w**2 * ((y - mu) ** 2 - v)) / np.sum(w**2) + (1 / np.sum(w) if reml else 0)
    assert f.tau2 == pytest.approx(fp, abs=1e-8)


def test_tau2_zero_boundary():
    y, v = np.array([1.0, 1.1, 0.9, 1.0]), np.array([1.0, 1.0, 1.0, 1.0])
    for m in ("DL", "PM", "ML", "REML"):
        assert pl.fit(y, v, method=m).tau2 == 0.0
    # Q(0) = 0.02 is below chi2_{3, 0.025} = 0.216, so both Q-profile bounds are 0
    assert pl.tau2_qprofile_ci(y, v) == (0.0, 0.0)


def test_qprofile_hits_chi2_targets():
    y, v = np.array([0.0, 1.0, 3.0, 0.5, -0.4]), np.array([0.25, 0.5, 1.0, 0.3, 0.2])
    lb, ub = pl.tau2_qprofile_ci(y, v)
    X = np.ones((5, 1))
    if lb > 0:
        assert pl.q_gen(y, v, X, lb) == pytest.approx(stats.chi2.ppf(0.975, 4))
    assert pl.q_gen(y, v, X, ub) == pytest.approx(stats.chi2.ppf(0.025, 4))


def test_k1_pooling():
    f = pl.fit([0.3], [0.04], method="REML", hksj=True)
    assert f.estimate == pytest.approx(0.3)
    assert f.se[0] == pytest.approx(0.2)
    assert not f.hksj and f.tau2 == 0


# --- published benchmark ------------------------------------------------------------

BCG = np.array([[4, 119, 11, 128], [6, 300, 29, 274], [3, 228, 11, 209],
                [62, 13536, 248, 12619], [33, 5036, 47, 5761], [180, 1361, 372, 1079],
                [8, 2537, 10, 619], [505, 87886, 499, 87892], [29, 7470, 45, 7232],
                [17, 1699, 65, 1600], [186, 50448, 141, 27197], [5, 2493, 3, 2338],
                [27, 16886, 29, 17825]], float)
ABLAT = np.array([44, 55, 42, 52, 13, 44, 19, 13, 27, 42, 18, 33, 33], float)


def bcg():
    a, b, c, d = BCG.T
    return np.log((a / (a + b)) / (c / (c + d))), 1 / a - 1 / (a + b) + 1 / c - 1 / (c + d)


def test_bcg_reml_matches_metafor():
    y, v = bcg()
    f = pl.fit(y, v, method="REML", hksj=False)
    assert f.estimate == pytest.approx(-0.7145, abs=5e-5)
    assert f.se[0] == pytest.approx(0.1798, abs=5e-5)
    assert f.tau2 == pytest.approx(0.3132, abs=5e-5)
    assert f.I2 == pytest.approx(92.22, abs=5e-3)
    assert f.H2 == pytest.approx(12.8558, abs=5e-5)
    assert f.QE == pytest.approx(152.2330, abs=5e-5)
    assert pl.tau2_qprofile_ci(y, v) == pytest.approx((0.1197, 1.1115), abs=5e-5)


def test_bcg_other_estimators_match_metafor():
    y, v = bcg()
    assert pl.fit(y, v, method="DL").tau2 == pytest.approx(0.3088, abs=5e-5)
    assert pl.fit(y, v, method="ML").tau2 == pytest.approx(0.2800, abs=5e-5)
    assert pl.fit(y, v, method="FE").estimate == pytest.approx(-0.4303, abs=5e-5)


def test_bcg_metaregression_matches_metafor():
    y, v = bcg()
    X = np.column_stack([np.ones(13), ABLAT])
    f = pl.fit(y, v, X, "REML", hksj=False)
    assert f.beta == pytest.approx([0.2515, -0.0291], abs=5e-5)
    assert f.se == pytest.approx([0.2491, 0.0072], abs=5e-5)
    # metafor prints tau2 0.0764, QM 16.3571 and R^2 75.62% at its default threshold (Fisher
    # scoring stops when the tau2 change falls below 1e-5, at tau2 0.0763547); the exact REML
    # optimum is 0.076348 (QM 16.3582). Hence the 1e-4 and 2e-3 tolerances here.
    assert f.tau2 == pytest.approx(0.0764, abs=1e-4)
    assert f.tau2 == pytest.approx(0.0763480, abs=1e-6)
    assert f.QE == pytest.approx(30.7331, abs=5e-5)
    assert (f.beta[1] / f.se[1]) ** 2 == pytest.approx(16.3571, abs=2e-3)


def test_pooled_row_columns():
    y, v = bcg()
    row, test, _ = pl.pooled_row(y, v, "REML", hksj=True)
    assert list(row) == pl.POOLED_COLUMNS + pl.POOLED_EXTRA_COLUMNS
    assert test["test"] == "t" and test["df"] == 12
    fe_row, _, _ = pl.pooled_row(y, v, "FE", hksj=True)
    assert math.isnan(fe_row["pi_lb"]) and math.isnan(fe_row["tau2_ci_lb"])


# --- numerical stability ------------------------------------------------------------

def test_q_is_unchanged_by_a_common_offset():
    # Astra 19: Q = y'Py cancelled catastrophically, so y = (0, 1, 2) + 1e8 gave Q = -1.24
    # instead of 1.0909. Heterogeneity, tau2 and the Q-profile limits do not depend on a
    # common offset in the effects; the estimate moves by exactly the offset.
    y, v = np.array([0.0, 1.0, 2.0]), np.array([1.0, 2.0, 3.0])
    shift = 1e8
    assert pl.fit(y + shift, v, method="FE").QE == pytest.approx(12 / 11, rel=1e-9)
    yb, vb = bcg()
    X = np.column_stack([np.ones(13), ABLAT])
    for method in ("FE", "DL", "PM", "ML", "REML"):
        for XX in (None, X):
            a = pl.fit(yb, vb, XX, method, hksj=True)
            b = pl.fit(yb + shift, vb, XX, method, hksj=True)
            assert b.QE == pytest.approx(a.QE, rel=1e-6)
            assert b.tau2 == pytest.approx(a.tau2, rel=1e-6, abs=1e-10)
            assert b.se == pytest.approx(a.se, rel=1e-6)
            assert b.beta[0] - shift == pytest.approx(a.beta[0], abs=1e-6)
            assert b.beta[1:] == pytest.approx(a.beta[1:], rel=1e-6)
    assert pl.q_gen(yb + shift, vb, X, 0.1) == pytest.approx(pl.q_gen(yb, vb, X, 0.1), rel=1e-6)
    assert (pl.tau2_qprofile_ci(yb + shift, vb)
            == pytest.approx(pl.tau2_qprofile_ci(yb, vb), rel=1e-6))
