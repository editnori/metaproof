"""Study-level effect sizes and their sampling variances.

Sign convention: every effect is experimental (is_control = 0) minus control
(is_control = 1); ROM is log(mean_experimental / mean_control).

Sources
-------
* MD and its variance: Cochrane Handbook v6.5, section 6.5.1.
* SMD (Hedges g): Hedges LV. J Educ Stat 1981;6(2):107-128, exact correction
  J(m) = Gamma(m/2) / (sqrt(m/2) Gamma((m-1)/2)), m = n1 + n2 - 2.
* Log ratio of means: Hedges LV, Gurevitch J, Curtis PS. Ecology 1999;80(4):1150-1156.
* Crossover trials: Elbourne DR, Altman DG, Higgins JPT, et al. Int J Epidemiol
  2002;31(1):140-149; Cochrane Handbook v6.5, section 23.2.
* Multi-arm trials: combining groups (recommended) or splitting the shared group,
  Cochrane Handbook v6.5, section 23.3.4 and Table 6.5.a; unit-of-analysis errors, 6.2.
  A crossover gives one NR-vs-placebo comparison. Several NR or placebo periods in the same
  people are refused: Table 6.5.a is for separate groups, and splitting N does not supply the
  covariance of contrasts that share participants. Unequal period n is refused unless a paired
  MD and SE are reported, because the number of complete pairs is not guessed.
* SMD standardisers (Handbook 10.5.2, 23.2.7.2): each SMD row records the SD it is divided by
  (smd_standardiser); a crossover SMD (between-person SD of the measurements) is not entered
  in a change-score SMD analysis.
* Geary's check for the ratio of means: Hedges, Gurevitch & Curtis 1999 (above).
* Median/range conversion limits: Cochrane Handbook v6.5, sections 6.5.2.5 and 6.5.2.6.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.special import gammaln

from . import conversions as cv

NUMERIC_COLUMNS = [
    "is_control", "n", "mean", "sd", "se", "median", "q1", "q3", "min", "max",
    "mean_change", "sd_change", "baseline_mean", "baseline_sd", "corr_prepost",
    "corr_within", "paired_md", "paired_se", "dose_mg_day", "duration_weeks",
    "se_change", "timepoint_weeks",      # optional extension columns
]
TEXT_COLUMNS = ["study_id", "outcome", "timepoint", "design", "arm", "population",
                "risk_of_bias", "notes"]
REQUIRED_COLUMNS = TEXT_COLUMNS[:5] + NUMERIC_COLUMNS[:2]
MODERATOR_COLUMNS = ["dose_mg_day", "duration_weeks", "population", "risk_of_bias"]
# Arm-level moderators are taken only from the experimental rows of a comparison, never from
# the control row (a placebo dose of 0 is not the dose of the NR arm).
ARM_LEVEL_MODERATORS = ("dose_mg_day",)
GEARY_MIN = 3.0


# --- formula layer ---------------------------------------------------------------

def hedges_j(m: float) -> float:
    """Exact small-sample correction J(m) = Gamma(m/2) / (sqrt(m/2) Gamma((m-1)/2))."""
    if m <= 1:
        raise ValueError("J(m) needs m > 1")
    return math.exp(gammaln(m / 2.0) - 0.5 * math.log(m / 2.0) - gammaln((m - 1) / 2.0))


def md_parallel(m1, sd1, n1, m2, sd2, n2):
    """MD = M1 - M2; Var = SD1^2/N1 + SD2^2/N2."""
    return m1 - m2, sd1**2 / n1 + sd2**2 / n2


def smd_parallel(m1, sd1, n1, m2, sd2, n2, vtype: str = "LS"):
    """Hedges g with exact J and one of three variance formulas.

    g = J(m) d, d = (M1 - M2) / S_pooled, m = n1 + n2 - 2,
    S_pooled = sqrt(((n1-1) SD1^2 + (n2-1) SD2^2) / m).

    vtype 'LS'  : 1/n1 + 1/n2 + g^2 / (2 (n1 + n2))            (large-sample; metafor default)
    vtype 'LS2' : J^2 (1/n1 + 1/n2 + d^2 / (2 (n1 + n2)))       (Borenstein et al. 2009)
    vtype 'UB'  : 1/n1 + 1/n2 + (1 - (m - 2) / (m J^2)) g^2     (unbiased; Hedges 1983)
    """
    m = n1 + n2 - 2
    sp = math.sqrt(((n1 - 1) * sd1**2 + (n2 - 1) * sd2**2) / m)
    if sp == 0:
        raise ValueError("pooled SD is 0 (every arm SD is 0), so the SMD is not defined; check "
                         "the source for rounded or missing SDs")
    d = (m1 - m2) / sp
    j = hedges_j(m)
    g = j * d
    if vtype == "LS":
        v = 1 / n1 + 1 / n2 + g**2 / (2 * (n1 + n2))
    elif vtype == "LS2":
        v = j**2 * (1 / n1 + 1 / n2 + d**2 / (2 * (n1 + n2)))
    elif vtype == "UB":
        v = 1 / n1 + 1 / n2 + (1 - (m - 2) / (m * j**2)) * g**2
    else:
        raise ValueError(f"unknown smd_variance {vtype!r}")
    return g, v


def rom_parallel(m1, sd1, n1, m2, sd2, n2, bias_correction: bool = False):
    """ln(M1/M2); Var = SD1^2 / (N1 M1^2) + SD2^2 / (N2 M2^2). Means must be > 0.

    bias_correction adds (SD1^2 / (N1 M1^2) - SD2^2 / (N2 M2^2)) / 2 to the estimate
    (Lajeunesse MJ. Ecology 2015;96(8):2056-2063), leaving the variance unchanged, which is
    what metafor escalc(measure = "ROM", correct = TRUE) does with its default vtype."""
    if m1 <= 0 or m2 <= 0:
        raise ValueError("ROM needs positive means")
    c1, c2 = sd1**2 / (n1 * m1**2), sd2**2 / (n2 * m2**2)
    yi = math.log(m1 / m2) + (0.5 * (c1 - c2) if bias_correction else 0.0)
    return yi, c1 + c2


def crossover_md(me, sde, mc, sdc, n, r):
    """Paired MD from period means: SD_diff = sqrt(SDe^2 + SDc^2 - 2 r SDe SDc); Var = SD_diff^2 / N."""
    sdiff2 = sde**2 + sdc**2 - 2 * r * sde * sdc
    return me - mc, sdiff2 / n


def crossover_smd(md, sde, sdc, n, r, j: float = 1.0):
    """Crossover SMD (Handbook 23.2.7.2; Elbourne 2002).

    d = MD / SD_within, SD_within = sqrt((SDe^2 + SDc^2) / 2) (between-participant SD);
    g = J d; Var = (1/N + g^2 / (2N)) * 2 (1 - r), the plug-in form used for parallel g.
    j = 1 gives the uncorrected d; the engine passes the exact J on N - 1 df
    (crossover_smd_j), so crossover and parallel rows carry the same small-sample correction.
    """
    sdw = math.sqrt((sde**2 + sdc**2) / 2.0)
    if sdw == 0:
        raise ValueError("both period SDs are 0, so the crossover SMD is not defined")
    d = j * md / sdw
    se = math.sqrt(1 / n + d**2 / (2 * n)) * math.sqrt(2 * (1 - r))
    return d, se**2


def crossover_rom(me, sde, mc, sdc, n, r, bias_correction: bool = False):
    """ln(Me/Mc); delta-method paired variance
    (SDe^2/Me^2 + SDc^2/Mc^2 - 2 r SDe SDc / (Me Mc)) / N.

    bias_correction adds (SDe^2 / (N Me^2) - SDc^2 / (N Mc^2)) / 2, the Lajeunesse 2015 term.
    It corrects E[ln Xbar] ~ ln mu - sigma^2 / (2 N mu^2) for each period mean on its own,
    and expectation is linear, so the within-person r does not enter it. metafor's ROMC has no
    correction; tests/test_effects.py checks this one by simulation."""
    if me <= 0 or mc <= 0:
        raise ValueError("ROM needs positive means")
    ce, cc = sde**2 / (n * me**2), sdc**2 / (n * mc**2)
    v = ce + cc - 2 * r * sde * sdc / (n * me * mc)
    return math.log(me / mc) + (0.5 * (ce - cc) if bias_correction else 0.0), v


def corr_within_from_paired(sde, sdc, n, paired_se):
    """Back-calculate the within-person correlation from a crossover trial that reports
    the paired SE and both period SDs: SD_diff = SE * sqrt(N),
    r = (SDe^2 + SDc^2 - SD_diff^2) / (2 SDe SDc)."""
    sdiff = paired_se * math.sqrt(n)
    return (sde**2 + sdc**2 - sdiff**2) / (2 * sde * sdc)


# --- arm preparation --------------------------------------------------------------

@dataclass
class Arm:
    label: str
    is_control: bool
    n: float
    mean: float | None = None
    sd: float | None = None
    mean_change: float | None = None
    sd_change: float | None = None
    baseline_mean: float | None = None
    baseline_sd: float | None = None
    corr_within: float | None = None
    paired_md: float | None = None
    paired_se: float | None = None
    notes: list = field(default_factory=list)
    sd_change_imputed: bool = False
    converted_from_median: bool = False
    sd_from_range_only: bool = False        # the SD itself is estimated from a range alone
    range_only: bool = False                 # the mean or the SD is estimated from a range alone
    median_asymmetry: float | None = None   # largest |quartile or range skewness| used
    change_derived: bool = False             # change mean or SD built from final values
    members: list = field(default_factory=list)

    def uses_conversion(self, basis: str) -> bool:
        """True when a median/IQR/range conversion feeds the statistics used for basis."""
        if not self.converted_from_median:
            return False
        return basis != "change" or self.change_derived


def _val(row, col):
    v = row.get(col)
    if v is None:
        return None
    try:
        v = float(v)
    except (TypeError, ValueError):
        raise ValueError(f"{row.get('arm')}: {col}={v!r} is not a number") from None
    if math.isnan(v):
        return None
    if not math.isfinite(v):
        raise ValueError(f"{row.get('arm')}: {col}={v!r} is not finite")
    return v


NONNEG_COLUMNS = ("sd", "se", "sd_change", "se_change", "baseline_sd")
CORR_COLUMNS = ("corr_prepost", "corr_within")


def _check_row(row) -> None:
    """Reject arm values the formulas would silently accept (a squared negative SD, an r
    outside [-1, 1], a blank is_control read as experimental). Raised as ValueError, so the
    study is excluded with the reason and the other studies still run."""
    label = row.get("arm")
    if _val(row, "is_control") not in (0.0, 1.0):
        raise ValueError(f"{label}: is_control must be 0 or 1 (got {row.get('is_control')!r})")
    n = _val(row, "n")
    if n is None or n < 1:
        raise ValueError(f"{label}: n is missing or < 1 (got {row.get('n')!r})")
    for col in NONNEG_COLUMNS:
        v = _val(row, col)
        if v is not None and v < 0:
            raise ValueError(f"{label}: {col} must be >= 0 (got {v:g})")
    v = _val(row, "paired_se")
    if v is not None and v <= 0:
        raise ValueError(f"{label}: paired_se must be > 0 (got {v:g})")
    for col in CORR_COLUMNS:
        v = _val(row, col)
        if v is not None and not -1 <= v <= 1:
            raise ValueError(f"{label}: {col} must be in [-1, 1] (got {v:g})")
    five = [(c, _val(row, c)) for c in ("min", "q1", "median", "q3", "max")]
    five = [(c, v) for c, v in five if v is not None]
    if any(a[1] > b[1] for a, b in zip(five, five[1:])):
        raise ValueError(f"{label}: five-number summary out of order; need min <= q1 <= "
                         f"median <= q3 <= max (got {dict(five)})")


def _asymmetry(med, q1, q3, mn, mx):
    """Bowley quartile skewness (q1 + q3 - 2 m)/(q3 - q1) and the range analogue
    (min + max - 2 m)/(max - min); 0 for a symmetric five-number summary. Reported only, as
    a pointer to skew; it is not the Shi 2023 test."""
    out = []
    if med is not None and q1 is not None and q3 is not None and q3 > q1:
        out.append((q1 + q3 - 2 * med) / (q3 - q1))
    if med is not None and mn is not None and mx is not None and mx > mn:
        out.append((mn + mx - 2 * med) / (mx - mn))
    return max(out, key=abs) if out else None


def prepare_arm(row: dict, median_method: str, corr_prepost: float) -> Arm:
    """Fill missing SD/mean/change statistics for one arm row. A row-level corr_prepost
    wins over the configured value, so sensitivity runs only move assumed correlations."""
    _check_row(row)
    n = _val(row, "n")
    arm = Arm(label=str(row["arm"]), is_control=_val(row, "is_control") == 1.0, n=n,
              members=[str(row["arm"])])
    for col in ("mean", "sd", "mean_change", "sd_change", "baseline_mean", "baseline_sd",
                "corr_within", "paired_md", "paired_se"):
        setattr(arm, col, _val(row, col))

    se = _val(row, "se")
    if arm.sd is None and se is not None:
        arm.sd = cv.sd_from_se(se, n)
        arm.notes.append(f"{arm.label}: SD from SE (Handbook 6.5.2.2)")
    se_change = _val(row, "se_change")
    if arm.sd_change is None and se_change is not None:
        arm.sd_change = cv.sd_from_se(se_change, n)
        arm.notes.append(f"{arm.label}: change SD from SE of change (Handbook 6.5.2.2)")

    if arm.mean is None or arm.sd is None:
        med = _val(row, "median")
        q1, q3, mn, mx = (_val(row, c) for c in ("q1", "q3", "min", "max"))
        has_iqr = q1 is not None and q3 is not None
        has_rng = mn is not None and mx is not None
        spread = has_iqr or has_rng
        mean_conv = arm.mean is None and med is not None and spread
        sd_conv = arm.sd is None and spread
        if mean_conv:
            arm.mean, scen, who = cv.median_mean(n, med, q1, q3, mn, mx, median_method)
            arm.notes.append(f"{arm.label}: mean from median ({scen}): {who}")
        if sd_conv:
            # The Wan 2014 SD uses only the IQR and/or range and n, not the median.
            arm.sd, scen = cv.wan_sd(n, q1, q3, mn, mx)
            what = {"S1": "range", "S2": "IQR", "S3": "IQR and range"}[scen]
            arm.notes.append(f"{arm.label}: SD from {what} ({scen}): Wan 2014")
        if mean_conv or sd_conv:
            arm.converted_from_median = True
            # Whether a range-only estimate is excluded (allow_range_only) is decided in
            # _one_study, once it is known that the converted value feeds the effect.
            arm.range_only = not has_iqr
            arm.sd_from_range_only = sd_conv and not has_iqr
            arm.median_asymmetry = _asymmetry(med, q1, q3, mn, mx)

    if arm.mean_change is None and arm.mean is not None and arm.baseline_mean is not None:
        arm.mean_change = arm.mean - arm.baseline_mean
        arm.change_derived = True
        arm.notes.append(f"{arm.label}: change mean = final - baseline")

    if arm.sd_change is None and arm.sd is not None and arm.baseline_sd is not None \
            and arm.mean_change is not None:
        row_r = _val(row, "corr_prepost")
        if row_r is not None:
            r, src = row_r, "reported"
        else:
            r, src = corr_prepost, "assumed"
        arm.sd_change = cv.sd_change_imputed(arm.baseline_sd, arm.sd, r)
        arm.sd_change_imputed = True
        arm.change_derived = True
        arm.notes.append(f"{arm.label}: change SD imputed, r={r:g} ({src}; Handbook 6.5.2.8)")
    return arm


def combine_arms(arms: list[Arm], label: str) -> Arm:
    """Merge arms with Handbook Table 6.5.a, separately for final, change and baseline stats."""
    out = Arm(label=label, is_control=arms[0].is_control, n=sum(a.n for a in arms))
    for mcol, scol in (("mean", "sd"), ("mean_change", "sd_change"),
                       ("baseline_mean", "baseline_sd")):
        if all(getattr(a, mcol) is not None and getattr(a, scol) is not None for a in arms):
            _, m, s = cv.combine_groups([(a.n, getattr(a, mcol), getattr(a, scol)) for a in arms])
            setattr(out, mcol, m)
            setattr(out, scol, s)
    out.sd_change_imputed = any(a.sd_change_imputed for a in arms)
    out.converted_from_median = any(a.converted_from_median for a in arms)
    out.sd_from_range_only = any(a.sd_from_range_only for a in arms)
    out.range_only = any(a.range_only for a in arms)
    out.change_derived = any(a.change_derived for a in arms)
    asym = [a.median_asymmetry for a in arms if a.median_asymmetry is not None]
    out.median_asymmetry = max(asym, key=abs) if asym else None
    out.members = [m for a in arms for m in a.members]
    for a in arms:
        out.notes.extend(a.notes)
    out.notes.append(f"{label}: combined {', '.join(a.label for a in arms)} (Handbook Table 6.5.a)")
    return out


# --- per-study effects ------------------------------------------------------------

@dataclass
class EffectSettings:
    measure: str = "MD"
    smd_variance: str = "LS"
    data_basis: str = "auto"         # auto | change | final
    corr_prepost: float = 0.5
    corr_within: float = 0.5
    median_method: str = "luo"
    multi_arm: str = "combine"       # combine (main analysis) | split (dose analyses)
    combine: dict = field(default_factory=dict)  # {study_id: {new_label: [arm, ...]}}
    rom_bias_correction: bool = False
    extra_columns: tuple = ()        # configured moderator/subgroup columns beyond the standard ones
    allow_range_only: bool = True    # the config default is False; see config.py
    crossover_smd_j: bool = True
    rom_geary: str = "flag"          # flag | exclude


def choose_basis(arms: list[Arm], s: EffectSettings) -> str:
    """Pick final values or change from baseline for one parallel study (all its arms).

    data_basis 'final' : final values.
    data_basis 'change': change from baseline, with the change SD imputed where it is
                         missing (Handbook 6.5.2.8). Not for SMD: an imputed change SD is
                         the standardiser, so the assumed r would move g itself (Handbook
                         23.2.7.2), and a study without a reported change SD is excluded.
    data_basis 'auto'  : MD - change when every arm reports a change SD; otherwise final
                         values when every arm has them; otherwise change with an imputed SD.
                         This imputes only when the study gives no assumption-free choice.
                         SMD - final values only, because change and final values should not
                         be combined as SMDs (Handbook 10.5.2).
    ROM always uses final values (a ratio of change scores is not defined)."""
    have_final = all(a.mean is not None and a.sd is not None for a in arms)
    have_change = all(a.mean_change is not None and a.sd_change is not None for a in arms)
    have_change_rep = have_change and not any(a.sd_change_imputed for a in arms)
    if s.measure == "ROM" or s.data_basis == "final" or (s.data_basis == "auto"
                                                        and s.measure == "SMD"):
        if have_final:
            return "final"
        raise ValueError(f"no usable final mean and SD for every arm (measure {s.measure}, "
                         f"data_basis {s.data_basis})")
    if s.data_basis == "change":
        if s.measure == "SMD":
            if have_change_rep:
                return "change"
            raise ValueError("SMD on change scores needs a reported SD of change in every arm; "
                             "an imputed change SD would make g depend on the assumed r "
                             "(Handbook 23.2.7.2). Use data_basis: final for this outcome.")
        if have_change:
            return "change"
        raise ValueError("data_basis=change but change statistics unavailable")
    if have_change_rep:
        return "change"
    if have_final:
        return "final"
    if have_change:
        return "change"
    raise ValueError("no usable mean and SD for both arms")


def geary(a: Arm) -> float | None:
    """Geary's statistic mean sqrt(n) / SD; the log ratio of means and its delta-method
    variance are reliable only when it is at least 3 in every arm (Hedges, Gurevitch &
    Curtis 1999)."""
    if a.mean is None or not a.sd or not a.n:
        return None
    return a.mean * math.sqrt(a.n) / a.sd


def _parallel_effect(e: Arm, c: Arm, s: EffectSettings, notes: list, basis: str):
    measure = s.measure
    if measure == "ROM" and s.data_basis == "change":
        notes.append("ROM uses final values only; change requested but ignored")
    if basis == "change":
        me, sde, mc, sdc = e.mean_change, e.sd_change, c.mean_change, c.sd_change
    else:
        me, sde, mc, sdc = e.mean, e.sd, c.mean, c.sd
    if measure == "MD":
        yi, vi = md_parallel(me, sde, e.n, mc, sdc, c.n)
    elif measure == "SMD":
        yi, vi = smd_parallel(me, sde, e.n, mc, sdc, c.n, s.smd_variance)
        notes.append(f"Hedges g, exact J, variance {s.smd_variance}")
    elif measure == "ROM":
        yi, vi = rom_parallel(me, sde, e.n, mc, sdc, c.n, s.rom_bias_correction)
        if s.rom_bias_correction:
            notes.append("ROM bias-corrected (Lajeunesse 2015)")
    else:
        raise ValueError(f"unknown measure {measure!r}")
    return yi, vi, basis


def _crossover_effect(e: Arm, c: Arm, s: EffectSettings, notes: list, flags: dict):
    """One NR period against one placebo period of a two-condition crossover."""
    paired_md = e.paired_md
    paired_se = e.paired_se
    if s.measure == "MD" and paired_md is not None and paired_se is not None:
        notes.append("crossover: paired MD and SE as reported")
        return paired_md, paired_se**2, "paired"
    if e.n != c.n:
        # min(n_e, n_c) is not the number of people with both periods, and the period
        # summaries may describe different people (Handbook 23.2.7.1 needs complete pairs).
        raise ValueError(f"crossover period n differ ({e.n:g} vs {c.n:g}): the number of "
                         "complete pairs is not known and the kit does not guess it. Enter "
                         "the complete-pair N on both rows with the period summaries of those "
                         "participants, or (MD) the reported paired MD and SE")
    n = e.n
    r_row = e.corr_within if e.corr_within is not None else c.corr_within
    if r_row is not None:
        r, src = r_row, "reported"
    else:
        r, src = s.corr_within, "assumed"
    have_sd = e.sd is not None and c.sd is not None
    marginal = ("paired MD from the period-condition means, which assumes no period effect "
                "(balanced sequences); prefer a period-adjusted paired estimate (Handbook 23.2.5)")
    if s.measure == "MD" and (e.mean is None or c.mean is None or not have_sd):
        only_change = (e.mean_change is not None and c.mean_change is not None
                       and (e.mean is None or c.mean is None))
        if only_change:
            raise ValueError("crossover row carries only within-period change; enter "
                             "paired_md/paired_se or the period means and SDs")
        raise ValueError("crossover needs paired_md/paired_se or both period means and SDs")
    if not have_sd:
        raise ValueError(f"crossover {s.measure} needs period SDs")
    if paired_se is not None and r_row is None and (e.sd == 0 or c.sd == 0):
        # r is undefined with a zero period SD, but for MD the reported paired SE is the SE.
        if s.measure == "MD":
            notes.append(f"crossover: {marginal}; variance from the reported paired SE "
                         "(a period SD is 0)")
            return e.mean - c.mean, paired_se ** 2, "paired"
        raise ValueError(f"crossover {s.measure}: a period SD is 0, so the within-person r "
                         "cannot be derived from the paired SE")
    if paired_se is not None and r_row is None:
        # Handbook 23.2.7.1: a reported paired SE with both period SDs fixes r.
        r = corr_within_from_paired(e.sd, c.sd, n, paired_se)
        src = "derived from paired SE"
        if not -1 <= r <= 1:
            raise ValueError(f"within-person r derived from paired SE is {r:.3g}, outside [-1, 1]")
    flags["corr_within_assumed"] = src == "assumed"
    if s.measure == "MD":
        notes.append(f"crossover: {marginal}; r_within={r:g} ({src})")
        yi, vi = crossover_md(e.mean, e.sd, c.mean, c.sd, n, r)
        return yi, vi, "paired"
    if s.measure == "SMD":
        if paired_md is not None:
            md = paired_md
        elif e.mean is None or c.mean is None:
            raise ValueError("crossover SMD needs both period means or the reported paired MD")
        else:
            md = e.mean - c.mean
            notes.append(f"crossover: {marginal}")
        j = hedges_j(n - 1) if s.crossover_smd_j else 1.0
        yi, vi = crossover_smd(md, e.sd, c.sd, n, r, j)
        notes.append(f"crossover SMD, {'J on N - 1 df' if s.crossover_smd_j else 'no J correction'}"
                     f", r_within={r:g} ({src}), standardised by the between-person SD of the "
                     "period measurements (Handbook 23.2.7.2)")
        return yi, vi, "paired"
    if s.measure == "ROM":
        if e.mean is None or c.mean is None:
            raise ValueError("crossover ROM needs period means")
        yi, vi = crossover_rom(e.mean, e.sd, c.mean, c.sd, n, r, s.rom_bias_correction)
        notes.append(f"crossover ROM (delta method), r_within={r:g} ({src})")
        if s.rom_bias_correction:
            notes.append("ROM bias-corrected (Lajeunesse 2015, per period mean; independent "
                         "of r_within)")
        return yi, vi, "paired"
    raise ValueError(f"unknown measure {s.measure!r}")


def _first(rows: pd.DataFrame, col: str):
    if col not in rows:
        return None
    vals = rows[col].dropna()
    return vals.iloc[0] if len(vals) else None


def moderator_values(rows: pd.DataFrame, members: list, carried: list) -> dict:
    """Moderators for one comparison, taken from its experimental rows (members). A value
    that differs between combined arms is left blank so the comparison drops out of models
    on that moderator. Study-level columns fall back to any row of the study; arm-level
    ones (dose) never take the control row's value."""
    src = rows[rows["arm"].astype(str).isin([str(m) for m in members])]
    out = {}
    for m in carried:
        vals = list(dict.fromkeys(src[m].dropna())) if m in src else []
        if len(vals) == 1:
            out[m] = vals[0]
        elif len(vals) > 1:
            out[m] = None
        elif m in ARM_LEVEL_MODERATORS:
            out[m] = None
        else:
            out[m] = _first(rows, m)
    if "dose_mg_day" in src:
        d = pd.to_numeric(src["dose_mg_day"], errors="coerce").dropna()
        out["dose_min"] = float(d.min()) if len(d) else None
        out["dose_max"] = float(d.max()) if len(d) else None
    return out


EFFECT_COLUMNS = ["study_id", "yi", "vi", "sei", "method_notes", "source_study", "arm_e",
                  "arm_c", "n_e", "n_c", "design", "data_basis", "smd_standardiser",
                  "sd_change_imputed",
                  "corr_within_assumed", "converted_from_median", "sd_from_range_only",
                  "median_asymmetry", "digitised", "imputed_any", "geary_min", "geary_ok",
                  "dose_min", "dose_max"]


def study_effects(df: pd.DataFrame, s: EffectSettings):
    """Compute one row per comparison. Returns (effects DataFrame, excluded DataFrame)."""
    out, excluded = [], []
    for sid, rows in df.groupby("study_id", sort=False):
        try:
            out.extend(_one_study(str(sid), rows, s))
        except (ValueError, ZeroDivisionError) as exc:
            excluded.append({"study_id": sid, "reason": str(exc)})
    cols = EFFECT_COLUMNS + _carried(s)
    return (pd.DataFrame(out, columns=cols),
            pd.DataFrame(excluded, columns=["study_id", "reason"]))


def _carried(s: EffectSettings) -> list:
    return MODERATOR_COLUMNS + [c for c in s.extra_columns if c not in MODERATOR_COLUMNS
                                and c not in ("design",)]


DIGITISED = ("digitis", "digitiz", "from figure", "from a figure", "read from fig")


def _digitised(rows: pd.DataFrame) -> bool:
    text = " ".join(str(x) for x in rows.get("notes", pd.Series(dtype=str)).dropna()).lower()
    return any(w in text for w in DIGITISED)


CROSSOVER_MULTI = (
    "multi-period crossover ({what} in the same participants): the comparisons share "
    "participants and the covariance of their contrasts is not available from summary data, "
    "so the periods can be neither combined with Table 6.5.a (which is for separate groups) "
    "nor entered as separate comparisons; splitting N does not supply the covariance "
    "(Handbook 23.2.7 and 23.3.4). Keep one prespecified NR-vs-placebo pair of rows, or one "
    "paired estimate the trial reports")
RANGE_ONLY = (
    "{labels}: only a range is reported (Wan 2014 scenario S1); a mean or SD estimated from a "
    "range alone is unstable (Handbook 6.5.2.5-6.5.2.6), so it is excluded unless "
    "allow_range_only: true (use for a sensitivity analysis)")


def _one_study(sid: str, rows: pd.DataFrame, s: EffectSettings):
    designs = {str(d).strip().lower() for d in rows["design"].dropna()}
    if len(designs) != 1 or not designs <= {"parallel", "crossover"}:
        raise ValueError(f"design must be one of parallel/crossover, got {sorted(designs)}")
    design = designs.pop()
    dup = rows["arm"].astype(str)[rows["arm"].astype(str).duplicated()].unique()
    if len(dup):
        raise ValueError(f"arm label {list(dup)} appears more than once (several timepoints "
                         "or duplicate rows); one row per arm is allowed (Handbook 6.2.4)")
    if s.multi_arm not in ("split", "combine"):
        raise ValueError(f"unknown multi_arm {s.multi_arm!r}")
    arms = [prepare_arm(r, s.median_method, s.corr_prepost) for r in rows.to_dict("records")]
    on_control = [a.label for a in arms if a.is_control
                  and (a.paired_md is not None or a.paired_se is not None)]
    if on_control:
        # The data dictionary puts them on the NR row; a control-row value was reused for every
        # comparison of the study (CR1-3).
        raise ValueError(f"paired_md/paired_se belong on the NR row, not the control row "
                         f"{on_control}")
    if design == "crossover":
        return _crossover_study(sid, rows, arms, s)

    for new_label, members in (s.combine.get(sid) or {}).items():
        group = [a for a in arms if a.label in set(map(str, members))]
        if len(group) != len(members):
            raise ValueError(f"combine_arms: arms {members} not all found")
        if len({a.is_control for a in group}) != 1:
            raise ValueError("combine_arms: cannot merge control with experimental arms")
        arms = [a for a in arms if a not in group] + [combine_arms(group, str(new_label))]

    controls = [a for a in arms if a.is_control]
    exps = [a for a in arms if not a.is_control]
    if not exps:
        raise ValueError("need at least one experimental arm")
    paired_only = s.measure == "MD" and len(exps) > 1 and any(_has_paired(a) for a in exps)
    if paired_only and (not controls or s.multi_arm == "split"
                        or not all(a.mean is not None and a.sd is not None for a in exps)):
        # Handbook 23.3.4: entering several comparisons that share one comparator "must be
        # avoided". Reported between-group MDs against a shared control cannot be split or
        # combined from their SEs alone.
        raise ValueError("more than one reported between-group MD (paired_md/paired_se) "
                         "against one control: the comparisons share the control group "
                         "(Handbook 23.3.4). Keep one prespecified comparison, or supply "
                         "arm-level data so the arms can be combined")
    if not controls:
        # A parallel trial may be extracted as one row holding the reported between-group
        # MD and its SE (for example an ANCOVA-adjusted MD), used directly by the generic
        # inverse-variance method (Handbook 6.3 and 10.5.2). Only MD can be passed through.
        if s.measure == "MD" and all(_has_paired(a) for a in exps):
            return _paired_only(sid, rows, exps, _carried(s))
        raise ValueError("need at least one control and one experimental arm "
                         "(or, for MD, paired_md and paired_se on every experimental row)")
    if len(controls) > 1:
        controls = [combine_arms(controls, "+".join(a.label for a in controls))]
    ctrl = controls[0]
    if len(exps) > 1 and s.multi_arm == "combine":
        exps = [combine_arms(exps, "+".join(a.label for a in exps))]

    k_arms = len(exps)
    # A reported between-group MD needs no arm moments, so the basis (which requires them)
    # is chosen only when an effect is built from arm summaries (Astra 4).
    from_arms = not (s.measure == "MD" and all(_has_paired(e) for e in exps))
    study_basis = choose_basis([ctrl] + exps, s) if from_arms else None
    result = []
    for e in exps:
        notes = list(ctrl.notes) + list(e.notes)
        reported = s.measure == "MD" and _has_paired(e)
        if reported:
            yi, vi, basis = e.paired_md, e.paired_se**2, "paired"
            notes.append("reported between-group MD and SE used directly (paired_md/paired_se)")
            n_c = ctrl.n
        else:
            c = ctrl
            if k_arms > 1:
                c = Arm(**{**ctrl.__dict__, "n": ctrl.n / k_arms, "notes": []})
                notes.append(f"shared control n {ctrl.n:g} split /{k_arms} (Handbook 23.3.4); "
                             "the comparisons stay correlated")
            yi, vi, basis = _parallel_effect(e, c, s, notes, study_basis)
            n_c = c.n
        result.append(_comparison(
            sid, f"{sid}|{e.label}" if k_arms > 1 else sid, "parallel", e, ctrl, yi, vi, basis,
            notes, {"corr_within_assumed": False}, e.n, n_c, reported, s, rows))
    return result


def _crossover_study(sid: str, rows: pd.DataFrame, arms: list[Arm], s: EffectSettings):
    """A crossover trial gives one comparison: one NR period against one placebo period."""
    if s.combine.get(sid):
        raise ValueError("combine_arms cannot merge crossover periods: they are the same "
                         "participants, and Handbook Table 6.5.a is for separate groups")
    controls = [a for a in arms if a.is_control]
    exps = [a for a in arms if not a.is_control]
    if not exps:
        raise ValueError("need at least one experimental arm")
    if not controls:
        raise ValueError("need at least one control and one experimental arm "
                         "(or, for MD, paired_md and paired_se on every experimental row)")
    if len(exps) > 1 or len(controls) > 1:
        what = f"{len(exps)} NR and {len(controls)} control periods"
        raise ValueError(CROSSOVER_MULTI.format(what=what))
    if s.measure == "SMD" and s.data_basis == "change":
        raise ValueError("a crossover SMD is standardised by the between-person SD of the "
                         "period measurements, not by an SD of change, so it cannot be pooled "
                         "with change-score SMDs (Handbook 10.5.2 and 23.2.7.2). Use "
                         "data_basis: final or auto for SMD analyses with crossover trials")
    e, c = exps[0], controls[0]
    notes = list(c.notes) + list(e.notes)
    flags = {"corr_within_assumed": False}
    yi, vi, basis = _crossover_effect(e, c, s, notes, flags)
    reported = s.measure == "MD" and _has_paired(e)
    return [_comparison(sid, sid, "crossover", e, c, yi, vi, basis, notes, flags, e.n, c.n,
                        reported, s, rows)]


def _comparison(sid, study_id, design, e, ctrl, yi, vi, basis, notes, flags, n_e, n_c,
                reported, s, rows) -> dict:
    """Checks shared by every comparison, and its row of study_effects."""
    if not (math.isfinite(yi) and math.isfinite(vi) and vi > 0):
        raise ValueError(f"non-finite effect or non-positive variance (yi={yi}, vi={vi})")
    conv_arms = [] if reported else [a for a in (e, ctrl) if a.uses_conversion(basis)]
    range_arms = [a.label for a in conv_arms if a.range_only]
    if range_arms and not s.allow_range_only:
        raise ValueError(RANGE_ONLY.format(labels=", ".join(range_arms)))
    g_min = None
    if s.measure == "ROM":
        gs = [geary(a) for a in (e, ctrl)]
        g_min = min(gs) if all(g is not None for g in gs) else None
        if g_min is not None and g_min < GEARY_MIN:
            msg = (f"Geary's check failed (min mean sqrt(n)/SD = {g_min:.2f} < 3): the log "
                   "ratio and its variance are unreliable (Hedges 1999)")
            if s.rom_geary == "exclude":
                raise ValueError(msg)
            notes.append(msg)
    asym = [a.median_asymmetry for a in conv_arms if a.median_asymmetry is not None]
    sd_imp = bool(basis == "change" and (e.sd_change_imputed or ctrl.sd_change_imputed))
    return {
        "study_id": study_id,
        "yi": yi, "vi": vi, "sei": math.sqrt(vi),
        "method_notes": "; ".join(dict.fromkeys(notes)),
        "source_study": sid, "arm_e": e.label, "arm_c": ctrl.label,
        "n_e": n_e, "n_c": n_c, "design": design, "data_basis": basis,
        "smd_standardiser": _standardiser(s, basis),
        "sd_change_imputed": sd_imp,
        "corr_within_assumed": bool(flags.get("corr_within_assumed")),
        "converted_from_median": bool(conv_arms),
        "sd_from_range_only": any(a.sd_from_range_only for a in conv_arms),
        "median_asymmetry": max(asym, key=abs) if asym else None,
        "digitised": _digitised(rows),
        "imputed_any": bool(sd_imp or flags.get("corr_within_assumed")),
        "geary_min": g_min, "geary_ok": None if g_min is None else bool(g_min >= GEARY_MIN),
        **moderator_values(rows, e.members, _carried(s)),
    }


def _standardiser(s: EffectSettings, basis: str) -> str | None:
    """SD an SMD is divided by: 'sd_change' for parallel change scores, 'sd_final' for
    parallel final values and for crossover period measurements (the between-person SD,
    Handbook 23.2.7.2). SMDs with different standardisers are not pooled (Handbook 10.5.2)."""
    if s.measure != "SMD":
        return None
    return "sd_change" if basis == "change" else "sd_final"


def _has_paired(a: Arm) -> bool:
    return a.paired_md is not None and a.paired_se is not None


def _paired_only(sid: str, rows: pd.DataFrame, exps: list[Arm], carried: list):
    """A row with no control arm whose reported MD and SE are used as they stand. Only one
    such comparison per study is allowed (several would share the unseen control group)."""
    if len(exps) > 1:
        raise ValueError("more than one reported between-group MD in one study with no "
                         "control row; the comparisons share a control group (Handbook "
                         "23.3.4). Keep one prespecified comparison")
    e = exps[0]
    vi = e.paired_se**2
    if not (math.isfinite(e.paired_md) and vi > 0):
        raise ValueError(f"non-finite paired_md or non-positive paired_se for {e.label}")
    return [{
        "study_id": sid, "yi": e.paired_md, "vi": vi, "sei": math.sqrt(vi),
        "method_notes": "reported between-group MD and SE used directly (paired_md/paired_se)",
        "source_study": sid, "arm_e": e.label, "arm_c": None, "n_e": e.n, "n_c": None,
        "design": "parallel", "data_basis": "paired", "smd_standardiser": None,
        "sd_change_imputed": False,
        "corr_within_assumed": False, "converted_from_median": False,
        "sd_from_range_only": False, "median_asymmetry": None, "digitised": _digitised(rows),
        "imputed_any": False, "geary_min": None, "geary_ok": None,
        **moderator_values(rows, e.members, carried),
    }]


def coerce_extraction(df: pd.DataFrame) -> pd.DataFrame:
    """Validate columns and coerce numeric types. Missing optional columns are added empty.
    A blank cell is missing; any other value that is not a finite number is an error."""
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"extraction CSV missing required columns: {missing}")
    df = df.copy()
    blank_id = df["study_id"].isna() | df["arm"].isna() \
        | df["study_id"].astype(str).str.strip().eq("") | df["arm"].astype(str).str.strip().eq("")
    if blank_id.any():
        raise ValueError(f"extraction CSV has a blank study_id or arm on CSV lines "
                         f"{[i + 2 for i in np.flatnonzero(blank_id)]}")
    bad = []
    for c in NUMERIC_COLUMNS:
        if c not in df.columns:
            df[c] = np.nan
        raw = df[c]
        text = raw.astype(str).str.strip()
        blank = raw.isna() | text.eq("")
        num = pd.to_numeric(raw.where(~blank), errors="coerce")
        for i in np.flatnonzero((num.isna() & ~blank).to_numpy()):
            bad.append(f"{df['study_id'].iloc[i]}/{df['arm'].iloc[i]}: {c}={raw.iloc[i]!r} "
                       "is not a number")
        for i in np.flatnonzero(np.isinf(num.to_numpy(float))):
            bad.append(f"{df['study_id'].iloc[i]}/{df['arm'].iloc[i]}: {c}={raw.iloc[i]!r} "
                       "is not finite")
        df[c] = num.astype(float)
    if bad:
        # A typo read as missing would silently change which data a study contributes.
        raise ValueError("extraction CSV has values that are not numbers (leave a cell blank "
                         "when it is not reported): " + "; ".join(bad))
    for c in TEXT_COLUMNS:
        if c not in df.columns:
            df[c] = None
    for c in ("study_id", "outcome", "arm", "timepoint", "design"):
        df[c] = df[c].astype(str).str.strip()
    return df
