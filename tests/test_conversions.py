"""Hand-worked checks for summary-statistic conversions.

Reference normal quantiles below come from Python's statistics.NormalDist, not scipy, so the
tests do not reuse the code path under test."""

import math
from statistics import NormalDist

import numpy as np
import pytest

from metaproof import conversions as cv

Z = NormalDist()


def test_sd_from_se():
    # SD = SE sqrt(n) = 2 * 5
    assert cv.sd_from_se(2.0, 25) == pytest.approx(10.0)


def test_change_sd_imputation_handbook_6_5_2_8():
    # sqrt(10^2 + 12^2 - 2 * 0.5 * 10 * 12) = sqrt(100 + 144 - 120) = sqrt(124)
    assert cv.sd_change_imputed(10, 12, 0.5) == pytest.approx(math.sqrt(124))
    assert cv.sd_change_imputed(10, 12, 0.5) == pytest.approx(11.135528725660043)
    # r = 0 gives sqrt(SDb^2 + SDf^2); r = 1 with equal SDs gives 0
    assert cv.sd_change_imputed(3, 4, 0.0) == pytest.approx(5.0)
    assert cv.sd_change_imputed(4, 4, 1.0) == pytest.approx(0.0)
    # Inverse formula recovers r
    assert cv.corr_from_sds(10, 12, math.sqrt(124)) == pytest.approx(0.5)


def test_combine_two_groups_table_6_5_a():
    # N = 30; M = (10*5 + 20*8)/30 = 7
    # SS = 9*4 + 19*9 + (10*20/30)*(25 + 64 - 80) = 36 + 171 + 60 = 267; SD = sqrt(267/29)
    n, m, sd = cv.combine_two_groups(10, 5, 2, 20, 8, 3)
    assert n == 30
    assert m == pytest.approx(7.0)
    assert sd == pytest.approx(math.sqrt(267 / 29))
    assert sd == pytest.approx(3.0342868275303405)


def test_combine_groups_matches_raw_data():
    rng = np.random.default_rng(1)
    parts = [rng.normal(loc, 2, size) for loc, size in ((0, 7), (3, 11), (-1, 5))]
    stats_ = [(len(p), p.mean(), p.std(ddof=1)) for p in parts]
    n, m, sd = cv.combine_groups(stats_)
    allx = np.concatenate(parts)
    assert n == len(allx)
    assert m == pytest.approx(allx.mean())
    assert sd == pytest.approx(allx.std(ddof=1))


def test_wan_divisors():
    # eta(50) = 2 Phi^-1((0.75*50 - 0.125)/(50 + 0.25)) = 2 Phi^-1(37.375/50.25)
    assert cv.wan_eta(50) == pytest.approx(2 * Z.inv_cdf(37.375 / 50.25))
    assert cv.wan_eta(50) == pytest.approx(1.3100930154533832)
    # xi(16) = 2 Phi^-1((16 - 0.375)/(16 + 0.25)) = 2 Phi^-1(15.625/16.25)
    assert cv.wan_xi(16) == pytest.approx(3.5376500770374117)
    # eta(n) -> 2 * 0.67449 (the normal IQR in SD units) as n grows
    assert cv.wan_eta(1e7) == pytest.approx(2 * Z.inv_cdf(0.75), rel=1e-6)


def test_wan_scenarios():
    # S2: mean (10 + 15 + 22)/3; SD = (22 - 10)/eta(50)
    mean, sd, _ = cv.median_to_mean_sd(50, 15, q1=10, q3=22, method="wan")
    assert mean == pytest.approx(47 / 3)
    assert sd == pytest.approx(12 / 1.3100930154533832)
    assert sd == pytest.approx(9.159654969877971)
    # S1: mean (2 + 2*8 + 20)/4 = 9.5; SD = 18/xi(16)
    mean, sd, _ = cv.median_to_mean_sd(16, 8, mn=2, mx=20, method="wan")
    assert mean == pytest.approx(9.5)
    assert sd == pytest.approx(5.088123360995052)
    # S3: mean (2 + 10 + 16 + 24 + 20)/8 = 9; SD = 18/(2 xi) + 7/(2 eta)
    mean, sd, _ = cv.median_to_mean_sd(16, 8, q1=5, q3=12, mn=2, mx=20, method="wan")
    assert mean == pytest.approx(9.0)
    assert sd == pytest.approx(5.3889374121386515)


def test_luo_scenarios():
    # S1, n = 16: n^0.75 = 8, w = 4/12 = 1/3; mean = (1/3)(2+20)/2 + (2/3) 8 = 9
    mean, sd, note = cv.median_to_mean_sd(16, 8, mn=2, mx=20, method="luo")
    assert mean == pytest.approx(9.0)
    assert sd == pytest.approx(5.088123360995052)  # Wan SD
    assert "Luo 2018" in note
    # S2, n = 50: w = 0.7 + 0.39/50 = 0.7078; mean = 0.7078 * 16 + 0.2922 * 15 = 15.7078
    mean, _, _ = cv.median_to_mean_sd(50, 15, q1=10, q3=22, method="luo")
    assert mean == pytest.approx(15.7078)
    # S3, n = 16: w1 = 2.2/10.2, w2 = 0.7 - 0.72/16^0.55 (16^0.55 = 2^2.2)
    w1 = 2.2 / 10.2
    w2 = 0.7 - 0.72 / 2**2.2
    mean, _, _ = cv.median_to_mean_sd(16, 8, q1=5, q3=12, mn=2, mx=20, method="luo")
    assert mean == pytest.approx(w1 * 11 + w2 * 8.5 + (1 - w1 - w2) * 8)
    assert mean == pytest.approx(8.91870927283276)


def test_median_needs_spread():
    with pytest.raises(ValueError):
        cv.median_to_mean_sd(20, 5, method="wan")
    with pytest.raises(ValueError):
        cv.median_to_mean_sd(20, 5, q1=1, q3=9, method="hozo")


def test_combine_two_groups_is_translation_invariant():
    # Astra 19: M1^2 + M2^2 - 2 M1 M2 cancels catastrophically for large means; (M1 - M2)^2
    # does not. Two N = 20, SD = 1 groups with means 0 and 1: SS = 19 + 19 + 10 * 1 = 48,
    # SD = sqrt(48 / 39). Shifting both means by 1e8 must leave the SD unchanged.
    expected = math.sqrt(48 / 39)
    assert cv.combine_two_groups(20, 0.0, 1.0, 20, 1.0, 1.0)[2] == pytest.approx(expected)
    n, m, sd = cv.combine_two_groups(20, 1e8, 1.0, 20, 1e8 + 1.0, 1.0)
    assert sd == pytest.approx(expected, rel=1e-12)
    assert m == pytest.approx(1e8 + 0.5, rel=1e-15)


def test_median_mean_reports_scenario_and_method():
    # the mean-only helper used when the SD was reported (CR1-6)
    mean, scen, label = cv.median_mean(50, 15, q1=10, q3=22, method="wan")
    assert (mean, scen, label) == (pytest.approx(47 / 3), "S2", "Wan 2014")
    mean, scen, label = cv.median_mean(50, 15, q1=10, q3=22, method="luo")
    assert (mean, scen, label) == (pytest.approx(15.7078), "S2", "Luo 2018")
