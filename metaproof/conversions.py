"""Summary-statistic conversions applied before effect sizes are computed.

Sources
-------
* SD from SE of a group mean: Cochrane Handbook v6.5, section 6.5.2.2.
* SD of change from baseline from baseline SD, final SD and a correlation:
  Cochrane Handbook v6.5, section 6.5.2.8.
* Combining two groups into one: Cochrane Handbook v6.5, Table 6.5.a.
* Median / IQR / range to mean and SD:
  Wan X, Wang W, Liu J, Tong T. BMC Med Res Methodol 2014;14:135.
  Luo D, Wan X, Liu J, Tong T. Stat Methods Med Res 2018;27(6):1785-1805.
"""

from __future__ import annotations

import math

from scipy.stats import norm


def _finite(*xs) -> bool:
    return all(x is not None and isinstance(x, (int, float)) and math.isfinite(x) for x in xs)


def sd_from_se(se: float, n: float) -> float:
    """SD = SE * sqrt(n) (Handbook 6.5.2.2)."""
    return se * math.sqrt(n)


def sd_change_imputed(sd_baseline: float, sd_final: float, corr: float) -> float:
    """SD_change = sqrt(SD_b^2 + SD_f^2 - 2 r SD_b SD_f) (Handbook 6.5.2.8)."""
    val = sd_baseline**2 + sd_final**2 - 2.0 * corr * sd_baseline * sd_final
    if val < 0:
        raise ValueError("negative variance for imputed change SD")
    return math.sqrt(val)


def corr_from_sds(sd_baseline: float, sd_final: float, sd_change: float) -> float:
    """Back-calculate the pre/post correlation from a study that reports all three SDs
    (Handbook 6.5.2.8): r = (SD_b^2 + SD_f^2 - SD_change^2) / (2 SD_b SD_f)."""
    return (sd_baseline**2 + sd_final**2 - sd_change**2) / (2.0 * sd_baseline * sd_final)


def combine_two_groups(n1, m1, sd1, n2, m2, sd2):
    """Handbook Table 6.5.a: one group from two.

    N = N1 + N2
    M = (N1 M1 + N2 M2) / N
    SD = sqrt( ((N1-1) SD1^2 + (N2-1) SD2^2 + N1 N2 / N (M1^2 + M2^2 - 2 M1 M2)) / (N - 1) )

    M1^2 + M2^2 - 2 M1 M2 is evaluated as (M1 - M2)^2, the same quantity without the
    cancellation that the expanded form suffers when the means are large next to their gap.
    """
    n = n1 + n2
    m = m1 + n2 / n * (m2 - m1)
    num = (n1 - 1) * sd1**2 + (n2 - 1) * sd2**2 + n1 * n2 / n * (m1 - m2) ** 2
    return n, m, math.sqrt(num / (n - 1))


def combine_groups(groups):
    """Apply Table 6.5.a repeatedly to a list of (n, mean, sd) tuples."""
    groups = list(groups)
    if not groups:
        raise ValueError("no groups to combine")
    n, m, sd = groups[0]
    for n2, m2, sd2 in groups[1:]:
        n, m, sd = combine_two_groups(n, m, sd, n2, m2, sd2)
    return n, m, sd


# --- Wan 2014 -----------------------------------------------------------------

def wan_xi(n: float) -> float:
    """xi(n) = 2 Phi^-1((n - 0.375) / (n + 0.25)), Wan 2014 range divisor."""
    return 2.0 * norm.ppf((n - 0.375) / (n + 0.25))


def wan_eta(n: float) -> float:
    """eta(n) = 2 Phi^-1((0.75 n - 0.125) / (n + 0.25)), Wan 2014 IQR divisor."""
    return 2.0 * norm.ppf((0.75 * n - 0.125) / (n + 0.25))


def wan_sd(n, q1=None, q3=None, mn=None, mx=None):
    """Wan 2014 SD for scenarios S1 (min, max), S2 (q1, q3) and S3 (all five)."""
    has_iqr = _finite(q1, q3)
    has_rng = _finite(mn, mx)
    if has_iqr and has_rng:
        return (mx - mn) / (2.0 * wan_xi(n)) + (q3 - q1) / (2.0 * wan_eta(n)), "S3"
    if has_iqr:
        return (q3 - q1) / wan_eta(n), "S2"
    if has_rng:
        return (mx - mn) / wan_xi(n), "S1"
    raise ValueError("need q1/q3 or min/max")


def wan_mean(median, q1=None, q3=None, mn=None, mx=None):
    """Wan 2014 means: S1 (a + 2m + b)/4, S2 (q1 + m + q3)/3, S3 (a + 2q1 + 2m + 2q3 + b)/8."""
    has_iqr = _finite(q1, q3)
    has_rng = _finite(mn, mx)
    if has_iqr and has_rng:
        return (mn + 2 * q1 + 2 * median + 2 * q3 + mx) / 8.0, "S3"
    if has_iqr:
        return (q1 + median + q3) / 3.0, "S2"
    if has_rng:
        return (mn + 2 * median + mx) / 4.0, "S1"
    raise ValueError("need q1/q3 or min/max")


def luo_mean(n, median, q1=None, q3=None, mn=None, mx=None):
    """Luo 2018 optimally weighted means.

    S1: w = 4 / (4 + n^0.75);           mean = w (a+b)/2 + (1-w) m
    S2: w = 0.7 + 0.39 / n;             mean = w (q1+q3)/2 + (1-w) m
    S3: w1 = 2.2 / (2.2 + n^0.75), w2 = 0.7 - 0.72 / n^0.55;
        mean = w1 (a+b)/2 + w2 (q1+q3)/2 + (1 - w1 - w2) m
    """
    has_iqr = _finite(q1, q3)
    has_rng = _finite(mn, mx)
    if has_iqr and has_rng:
        w1 = 2.2 / (2.2 + n**0.75)
        w2 = 0.7 - 0.72 / n**0.55
        return w1 * (mn + mx) / 2 + w2 * (q1 + q3) / 2 + (1 - w1 - w2) * median, "S3"
    if has_iqr:
        w = 0.7 + 0.39 / n
        return w * (q1 + q3) / 2 + (1 - w) * median, "S2"
    if has_rng:
        w = 4.0 / (4.0 + n**0.75)
        return w * (mn + mx) / 2 + (1 - w) * median, "S1"
    raise ValueError("need q1/q3 or min/max")


def median_mean(n, median, q1=None, q3=None, mn=None, mx=None, method="luo"):
    """Return (mean, scenario, label): the Wan 2014 ('wan') or Luo 2018 ('luo') mean."""
    if method == "wan":
        return (*wan_mean(median, q1, q3, mn, mx), "Wan 2014")
    if method == "luo":
        return (*luo_mean(n, median, q1, q3, mn, mx), "Luo 2018")
    raise ValueError(f"unknown median method {method!r}")


def median_to_mean_sd(n, median, q1=None, q3=None, mn=None, mx=None, method="luo"):
    """Return (mean, sd, note). method 'wan' = Wan 2014 mean and SD;
    'luo' = Luo 2018 mean with Wan 2014 SD."""
    mean, scen, label = median_mean(n, median, q1, q3, mn, mx, method)
    sd, _ = wan_sd(n, q1, q3, mn, mx)
    return mean, sd, f"mean/SD from median ({scen}): {label} mean, Wan 2014 SD"
