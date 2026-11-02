"""Hand-worked checks for effect sizes and the arm-to-comparison logic."""

import math

import numpy as np
import pandas as pd
import pytest

from metaproof import effects as ef
from metaproof.effects import hedges_j, EffectSettings, study_effects

J10 = 24 / (math.sqrt(5) * 6.5625 * math.sqrt(math.pi))  # Gamma(5)/(sqrt(5) Gamma(4.5))


def test_hedges_j_exact():
    # Gamma(4.5) = 3.5 * 2.5 * 1.5 * 0.5 * sqrt(pi) = 6.5625 sqrt(pi)
    assert ef.hedges_j(10) == pytest.approx(J10)
    assert ef.hedges_j(10) == pytest.approx(0.922745608053087)
    # Close to, but not equal to, the approximation 1 - 3/(4m - 1)
    assert abs(ef.hedges_j(10) - (1 - 3 / 39)) > 1e-4
    assert ef.hedges_j(4000) == pytest.approx(1 - 3 / (4 * 4000 - 1), rel=1e-7)


def test_md_parallel():
    # 10 - 8 = 2; 4/20 + 9/30 = 0.5
    assert ef.md_parallel(10, 2, 20, 8, 3, 30) == pytest.approx((2.0, 0.5))


def test_smd_parallel_variances():
    # n1 = n2 = 6, m = 10, SD 2 in both arms -> S_pooled = 2, d = (12 - 10)/2 = 1, g = J(10)
    g, v = ef.smd_parallel(12, 2, 6, 10, 2, 6, "LS")
    assert g == pytest.approx(J10)
    assert v == pytest.approx(1 / 6 + 1 / 6 + J10**2 / 24)
    assert v == pytest.approx(0.36881081071588584)
    _, v2 = ef.smd_parallel(12, 2, 6, 10, 2, 6, "LS2")
    assert v2 == pytest.approx(J10**2 * (1 / 3 + 1 / 24))
    _, v3 = ef.smd_parallel(12, 2, 6, 10, 2, 6, "UB")
    assert v3 == pytest.approx(1 / 3 + (1 - 8 / (10 * J10**2)) * J10**2)


def test_smd_pooled_sd_unequal():
    # S_pooled = sqrt((9*4 + 19*9)/28) = sqrt(207/28)
    g, _ = ef.smd_parallel(5, 2, 10, 4, 3, 20)
    assert g == pytest.approx(ef.hedges_j(28) * 1 / math.sqrt(207 / 28))


def test_rom_parallel():
    # ln(20/10); 16/(10*400) + 4/(10*100) = 0.008
    y, v = ef.rom_parallel(20, 4, 10, 10, 2, 10)
    assert y == pytest.approx(math.log(2))
    assert v == pytest.approx(0.008)
    with pytest.raises(ValueError):
        ef.rom_parallel(-1, 1, 10, 2, 1, 10)


def test_crossover_md_from_arms():
    # SD_diff^2 = 16 + 9 - 2*0.5*4*3 = 13; Var = 13/16
    y, v = ef.crossover_md(12, 4, 10, 3, 16, 0.5)
    assert y == pytest.approx(2.0)
    assert v == pytest.approx(0.8125)


def test_crossover_smd():
    # SD_within = 4, d = 2/4 = 0.5; SE^2 = (1/16 + 0.25/32) * 2 (1 - 0.5) = 0.0703125
    d, v = ef.crossover_smd(2, 4, 4, 16, 0.5)
    assert d == pytest.approx(0.5)
    assert v == pytest.approx(0.0703125)


def test_crossover_rom():
    # (16/400 + 4/100 - 2*0.5*4*2/(20*10)) / 10 = (0.04 + 0.04 - 0.04)/10 = 0.004
    y, v = ef.crossover_rom(20, 4, 10, 2, 10, 0.5)
    assert y == pytest.approx(math.log(2))
    assert v == pytest.approx(0.004)


def test_corr_within_from_paired_se():
    # SD_diff = 1 * sqrt(16) = 4; r = (16 + 16 - 16)/(2*16) = 0.5
    assert ef.corr_within_from_paired(4, 4, 16, 1.0) == pytest.approx(0.5)


# --- study-level assembly --------------------------------------------------------

def frame(rows):
    df = pd.DataFrame(rows)
    df.insert(1, "outcome", "SYNTHETIC_x")
    df.insert(2, "timepoint", "8 weeks")
    df["notes"] = "SYNTHETIC test row"
    return ef.coerce_extraction(df)


def arm(sid, label, ctl, n, design="parallel", **kw):
    return {"study_id": sid, "design": design, "arm": label, "is_control": ctl, "n": n, **kw}


def test_parallel_final_md():
    df = frame([arm("S1", "pbo", 1, 30, mean=8, sd=3), arm("S1", "nr", 0, 20, mean=10, sd=2)])
    eff, exc = study_effects(df, EffectSettings(measure="MD"))
    assert exc.empty
    r = eff.iloc[0]
    assert (r.yi, r.vi, r.sei) == pytest.approx((2.0, 0.5, math.sqrt(0.5)))
    assert r.data_basis == "final"


def test_change_imputed_and_auto_basis():
    # Each arm: change SD = sqrt(100 + 144 - 120) = sqrt(124) at r = 0.5
    rows = [arm("S1", "pbo", 1, 20, baseline_mean=100, baseline_sd=10, mean=101, sd=12),
            arm("S1", "nr", 0, 20, baseline_mean=100, baseline_sd=10, mean=97, sd=12)]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="MD", corr_prepost=0.5,
                                                       data_basis="change"))
    r = eff.iloc[0]
    assert r.data_basis == "change"
    assert r.yi == pytest.approx(-3 - 1)
    assert r.vi == pytest.approx(2 * 124 / 20)
    assert bool(r.sd_change_imputed)
    # auto: final values are complete, so no change SD is imputed (same rule as run_meta.R)
    effA, _ = study_effects(frame(rows), EffectSettings(measure="MD", corr_prepost=0.5))
    assert effA.iloc[0].data_basis == "final"
    assert effA.iloc[0].vi == pytest.approx(2 * 144 / 20)
    # auto falls back to an imputed change SD only when final values are incomplete
    rows_nf = [dict(r_, sd=12, mean=None, mean_change=c) for r_, c in zip(rows, (1, -3))]
    effN, _ = study_effects(frame(rows_nf), EffectSettings(measure="MD", corr_prepost=0.5))
    assert effN.iloc[0].data_basis == "change"
    assert effN.iloc[0].vi == pytest.approx(2 * 124 / 20)
    # r = 0.8: SD^2 = 100 + 144 - 192 = 52
    eff8, _ = study_effects(frame(rows), EffectSettings(measure="MD", corr_prepost=0.8,
                                                        data_basis="change"))
    assert eff8.iloc[0].vi == pytest.approx(2 * 52 / 20)
    # final basis ignores change data: Var = 2 * 144/20
    effF, _ = study_effects(frame(rows), EffectSettings(measure="MD", data_basis="final"))
    assert effF.iloc[0].vi == pytest.approx(2 * 144 / 20)
    assert effF.iloc[0].yi == pytest.approx(-4)


def test_reported_corr_wins_over_config():
    rows = [arm("S1", "pbo", 1, 20, baseline_mean=0, baseline_sd=10, mean=0, sd=12, corr_prepost=0.8),
            arm("S1", "nr", 0, 20, baseline_mean=0, baseline_sd=10, mean=1, sd=12, corr_prepost=0.8)]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="MD", corr_prepost=0.2,
                                                       data_basis="change"))
    assert eff.iloc[0].vi == pytest.approx(2 * 52 / 20)


def test_reported_change_stats_used_directly():
    rows = [arm("S1", "pbo", 1, 25, mean_change=1, sd_change=5),
            arm("S1", "nr", 0, 25, mean_change=-2, sd_change=5)]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="MD"))
    assert (eff.iloc[0].yi, eff.iloc[0].vi) == pytest.approx((-3, 2.0))
    assert not bool(eff.iloc[0].sd_change_imputed)


def test_rom_forces_final_values():
    rows = [arm("S1", "pbo", 1, 10, mean=10, sd=2, mean_change=1, sd_change=1),
            arm("S1", "nr", 0, 10, mean=20, sd=4, mean_change=5, sd_change=1)]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="ROM"))
    assert eff.iloc[0].yi == pytest.approx(math.log(2))
    assert eff.iloc[0].data_basis == "final"


def test_sd_from_se_in_study():
    rows = [arm("S1", "pbo", 1, 25, mean=8, se=0.6), arm("S1", "nr", 0, 25, mean=10, se=0.4)]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="MD"))
    # SDs 3 and 2 -> Var = 9/25 + 4/25 = 0.52
    assert eff.iloc[0].vi == pytest.approx(0.52)


def test_median_rows_converted():
    rows = [arm("S1", "pbo", 1, 50, median=15, q1=10, q3=22),
            arm("S1", "nr", 0, 50, mean=18, sd=9.159654969877971)]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="MD", median_method="wan"))
    assert eff.iloc[0].yi == pytest.approx(18 - 47 / 3)
    assert eff.iloc[0].vi == pytest.approx(2 * 9.159654969877971**2 / 50)


def test_multi_arm_split_shared_control():
    # Control n 30 split over two comparisons -> 15 each
    rows = [arm("S1", "pbo", 1, 30, mean=10, sd=3),
            arm("S1", "low", 0, 20, mean=11, sd=2), arm("S1", "high", 0, 20, mean=12, sd=2)]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="MD", multi_arm="split"))
    assert list(eff.study_id) == ["S1|low", "S1|high"]
    assert list(eff.yi) == pytest.approx([1, 2])
    assert list(eff.vi) == pytest.approx([4 / 20 + 9 / 15] * 2)
    assert list(eff.n_c) == [15, 15]


def test_multi_arm_combine():
    # Experimental arms (20, 11, 2) and (20, 13, 2): N 40, mean 12,
    # SD^2 = (19*4 + 19*4 + 10*(121 + 169 - 286))/39 = (152 + 40)/39
    rows = [arm("S1", "pbo", 1, 30, mean=10, sd=3),
            arm("S1", "low", 0, 20, mean=11, sd=2), arm("S1", "high", 0, 20, mean=13, sd=2)]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="MD", multi_arm="combine"))
    assert len(eff) == 1
    assert eff.iloc[0].yi == pytest.approx(2)
    assert eff.iloc[0].vi == pytest.approx(192 / 39 / 40 + 9 / 30)


def test_two_controls_are_combined():
    rows = [arm("S1", "pbo_a", 1, 10, mean=5, sd=2), arm("S1", "pbo_b", 1, 20, mean=8, sd=3),
            arm("S1", "nr", 0, 30, mean=9, sd=3)]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="MD"))
    # combined control: N 30, mean 7, SD^2 = 267/29 (see test_conversions)
    assert eff.iloc[0].yi == pytest.approx(2)
    assert eff.iloc[0].vi == pytest.approx(9 / 30 + 267 / 29 / 30)


def test_configured_arm_combination():
    rows = [arm("S1", "pbo", 1, 30, mean=7, sd=3),
            arm("S1", "nr_men", 0, 10, mean=5, sd=2), arm("S1", "nr_women", 0, 20, mean=8, sd=3)]
    s = EffectSettings(measure="MD", combine={"S1": {"nr": ["nr_men", "nr_women"]}})
    eff, _ = study_effects(frame(rows), s)
    assert eff.iloc[0].study_id == "S1"
    assert eff.iloc[0].yi == pytest.approx(0)
    assert eff.iloc[0].vi == pytest.approx(267 / 29 / 30 + 9 / 30)


def test_crossover_paired_and_from_arms():
    paired = [arm("X1", "nr", 0, 16, "crossover", mean=12, sd=4, paired_md=1.5, paired_se=0.5),
              arm("X1", "pbo", 1, 16, "crossover", mean=10, sd=3)]
    eff, _ = study_effects(frame(paired), EffectSettings(measure="MD"))
    assert (eff.iloc[0].yi, eff.iloc[0].vi) == pytest.approx((1.5, 0.25))
    arms_only = [arm("X2", "nr", 0, 16, "crossover", mean=12, sd=4),
                 arm("X2", "pbo", 1, 16, "crossover", mean=10, sd=3)]
    eff, _ = study_effects(frame(arms_only), EffectSettings(measure="MD", corr_within=0.5))
    assert (eff.iloc[0].yi, eff.iloc[0].vi) == pytest.approx((2, 0.8125))
    eff, _ = study_effects(frame(arms_only), EffectSettings(measure="MD", corr_within=0.0))
    assert eff.iloc[0].vi == pytest.approx(25 / 16)


def test_crossover_smd_uses_corr_from_paired_se():
    # paired_se = 1, N = 16, SDs 4 and 4 -> r = 0.5 (see test_corr_within_from_paired_se)
    rows = [arm("X1", "nr", 0, 16, "crossover", mean=12, sd=4, paired_md=2, paired_se=1),
            arm("X1", "pbo", 1, 16, "crossover", mean=10, sd=4)]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="SMD", corr_within=0.9,
                                                       crossover_smd_j=False))
    assert (eff.iloc[0].yi, eff.iloc[0].vi) == pytest.approx((0.5, 0.0703125))
    # default: Hedges J on N - 1 = 15 df, as for parallel g (F08b)
    eff, _ = study_effects(frame(rows), EffectSettings(measure="SMD", corr_within=0.9))
    g = 0.5 * hedges_j(15)
    assert (eff.iloc[0].yi, eff.iloc[0].vi) == pytest.approx((g, (1 / 16 + g**2 / 32) * 2 * 0.5))


def test_crossover_impossible_derived_corr_excluded():
    # paired_se 3 -> SD_diff 12 > SDe + SDc, r = (16 + 16 - 144)/32 = -3.5
    rows = [arm("X1", "nr", 0, 16, "crossover", mean=12, sd=4, paired_md=2, paired_se=3),
            arm("X1", "pbo", 1, 16, "crossover", mean=10, sd=4)]
    eff, exc = study_effects(frame(rows), EffectSettings(measure="SMD"))
    assert eff.empty and "outside [-1, 1]" in exc.reason.iloc[0]


def test_exclusions_are_reported():
    rows = [arm("S1", "nr", 0, 10, mean=1, sd=1), arm("S1", "nr2", 0, 10, mean=1, sd=1),
            arm("S2", "pbo", 1, 10, mean=1), arm("S2", "nr", 0, 10, mean=2)]
    eff, exc = study_effects(frame(rows), EffectSettings(measure="MD"))
    assert eff.empty
    assert set(exc.study_id) == {"S1", "S2"}


def test_missing_required_column():
    with pytest.raises(ValueError, match="missing required"):
        ef.coerce_extraction(pd.DataFrame({"study_id": ["S1"]}))


def test_effects_carry_moderators():
    rows = [arm("S1", "pbo", 1, 10, mean=1, sd=1, dose_mg_day=np.nan, population="SYNTHETIC_a"),
            arm("S1", "nr", 0, 10, mean=2, sd=1, dose_mg_day=1000, population="SYNTHETIC_a")]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="MD"))
    assert eff.iloc[0].dose_mg_day == 1000
    assert eff.iloc[0].population == "SYNTHETIC_a"


# --- review fixes (FX1): crossover structure, standardisers, routing, input checks ---------

XO = "crossover"


def xo_pair(sid, n_e=16, n_c=16, **kw):
    return [arm(sid, "nr", 0, n_e, XO, mean=12, sd=4, **kw),
            arm(sid, "pbo", 1, n_c, XO, mean=10, sd=4)]


@pytest.mark.parametrize("multi_arm", ["combine", "split"])
def test_multi_period_crossover_is_refused(multi_arm):
    # Astra 1: two NR periods and one placebo period in the same people. Splitting N does not
    # supply the covariance of the two contrasts, so they are not pooled as separate studies.
    rows = [arm("X1", "low", 0, 18, XO, mean=44, sd=10), arm("X1", "high", 0, 18, XO, mean=49, sd=11),
            arm("X1", "pbo", 1, 18, XO, mean=40, sd=10.5)] + xo_pair("X2")
    for measure in ("MD", "SMD", "ROM"):
        eff, exc = study_effects(frame(rows), EffectSettings(measure=measure, multi_arm=multi_arm))
        assert list(eff.study_id) == ["X2"]
        assert list(exc.study_id) == ["X1"]
        assert "multi-period crossover" in exc.reason.iloc[0]
        assert "same participants" in exc.reason.iloc[0]


def test_multi_period_crossover_with_reported_paired_mds_is_refused():
    rows = [arm("X1", "low", 0, 18, XO, paired_md=4, paired_se=1),
            arm("X1", "high", 0, 18, XO, paired_md=9, paired_se=1),
            arm("X1", "pbo", 1, 18, XO, mean=40, sd=10.5)]
    eff, exc = study_effects(frame(rows), EffectSettings(measure="MD"))
    assert eff.empty and "multi-period crossover" in exc.reason.iloc[0]


def test_crossover_repeated_placebo_periods_are_not_combined_as_separate_groups():
    # Table 6.5.a is for separate groups; two placebo periods are the same people.
    rows = [arm("X1", "nr", 0, 16, XO, mean=12, sd=4), arm("X1", "pbo1", 1, 16, XO, mean=10, sd=4),
            arm("X1", "pbo2", 1, 16, XO, mean=11, sd=4)]
    eff, exc = study_effects(frame(rows), EffectSettings(measure="MD"))
    assert eff.empty and "multi-period crossover" in exc.reason.iloc[0]


def test_configured_combine_refused_for_crossover_periods():
    rows = [arm("X1", "low", 0, 16, XO, mean=11, sd=4), arm("X1", "high", 0, 16, XO, mean=13, sd=4),
            arm("X1", "pbo", 1, 16, XO, mean=10, sd=4)]
    s = EffectSettings(measure="MD", combine={"X1": {"nr": ["low", "high"]}})
    eff, exc = study_effects(frame(rows), s)
    assert eff.empty and "combine_arms" in exc.reason.iloc[0]
    assert "same participants" in exc.reason.iloc[0]


def test_crossover_unequal_period_n_is_not_guessed():
    # Astra 2: min(n_e, n_c) is not the number of complete pairs.
    for measure in ("MD", "SMD", "ROM"):
        eff, exc = study_effects(frame(xo_pair("X1", n_e=16, n_c=14)),
                                 EffectSettings(measure=measure))
        assert eff.empty, measure
        assert "complete pairs" in exc.reason.iloc[0]
    # a reported paired MD and SE needs no N, so it is still used
    rows = [arm("X1", "nr", 0, 16, XO, paired_md=1.5, paired_se=0.5),
            arm("X1", "pbo", 1, 14, XO, mean=10, sd=3)]
    eff, exc = study_effects(frame(rows), EffectSettings(measure="MD"))
    assert exc.empty and (eff.iloc[0].yi, eff.iloc[0].vi) == pytest.approx((1.5, 0.25))


def test_smd_standardiser_recorded_and_change_basis_refuses_crossover():
    # Astra 3 / Handbook 10.5.2: parallel change-score SMDs are standardised by the SD of
    # change; crossover SMDs by the SD of the period measurements. They are not pooled together.
    par = [arm("P1", "pbo", 1, 25, mean_change=1, sd_change=5),
           arm("P1", "nr", 0, 25, mean_change=-2, sd_change=5)]
    eff, exc = study_effects(frame(par + xo_pair("X1")),
                             EffectSettings(measure="SMD", data_basis="change"))
    assert list(eff.study_id) == ["P1"] and eff.iloc[0].smd_standardiser == "sd_change"
    assert list(exc.study_id) == ["X1"] and "standardised" in exc.reason.iloc[0]
    par_final = [arm("P2", "pbo", 1, 25, mean=10, sd=5), arm("P2", "nr", 0, 25, mean=12, sd=5)]
    eff, exc = study_effects(frame(par_final + xo_pair("X1")), EffectSettings(measure="SMD"))
    assert exc.empty and list(eff.smd_standardiser) == ["sd_final", "sd_final"]
    eff, _ = study_effects(frame(par_final), EffectSettings(measure="MD"))
    assert eff.smd_standardiser.isna().all()


def test_adjusted_md_used_when_control_row_has_arm_data():
    # Astra 4: the NR row carries an ANCOVA-adjusted MD and SE; the control row only its arm
    # summary. The reported MD is used; arm moments are not required for it.
    rows = [arm("S1", "nr", 0, 30, paired_md=2, paired_se=0.5),
            arm("S1", "pbo", 1, 30, mean=10, sd=3)]
    eff, exc = study_effects(frame(rows), EffectSettings(measure="MD"))
    assert exc.empty
    assert (eff.iloc[0].yi, eff.iloc[0].vi, eff.iloc[0].data_basis) == (2, 0.25, "paired")


def test_unused_range_summary_does_not_exclude_a_study():
    # Astra 10: change statistics are reported; a median/range for the (unused) final values
    # must not trigger the range-only exclusion.
    rows = [arm("S1", "pbo", 1, 20, mean_change=1, sd_change=1, median=10, min=2, max=20),
            arm("S1", "nr", 0, 20, mean_change=2, sd_change=1, median=12, min=3, max=22)]
    eff, exc = study_effects(frame(rows), EffectSettings(measure="MD", allow_range_only=False))
    assert exc.empty and eff.iloc[0].data_basis == "change"
    assert (eff.iloc[0].yi, eff.iloc[0].vi) == pytest.approx((1, 0.1))
    assert not bool(eff.iloc[0].converted_from_median)
    # the same range is used, and so excluded, when the analysis needs final values
    eff, exc = study_effects(frame(rows), EffectSettings(measure="MD", data_basis="final",
                                                         allow_range_only=False))
    assert eff.empty and "only a range" in exc.reason.iloc[0]
    # a reported between-group MD needs no arm summary at all
    rows = [arm("S2", "nr", 0, 30, paired_md=2, paired_se=0.5, median=12, min=3, max=22),
            arm("S2", "pbo", 1, 30, median=10, min=2, max=20)]
    eff, exc = study_effects(frame(rows), EffectSettings(measure="MD", allow_range_only=False))
    assert exc.empty and eff.iloc[0].yi == 2


def test_range_conversion_labels_say_what_was_converted():
    # CR1-6: SD reported, mean from median + range
    a = ef.prepare_arm({"arm": "a", "is_control": 1, "n": 40, "sd": 5.0, "median": 10,
                        "min": 2, "max": 20}, "luo", 0.5)
    assert a.sd == 5.0 and not a.sd_from_range_only and a.range_only
    assert any("mean from median (S1): Luo 2018" in x for x in a.notes)
    assert not any("SD" in x for x in a.notes)
    # mean reported, SD from the IQR
    a = ef.prepare_arm({"arm": "b", "is_control": 1, "n": 50, "mean": 16.0, "median": 15,
                        "q1": 10, "q3": 22}, "luo", 0.5)
    assert a.mean == 16.0 and not a.range_only
    assert a.notes == ["b: SD from IQR (S2): Wan 2014"]
    rows = [arm("S1", "pbo", 1, 40, sd=5.0, median=10, min=2, max=20),
            arm("S1", "nr", 0, 40, mean=12, sd=5.0)]
    eff, exc = study_effects(frame(rows), EffectSettings(measure="MD", allow_range_only=False))
    assert eff.empty and "mean or SD estimated from a range alone" in exc.reason.iloc[0]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="MD", allow_range_only=True))
    assert not bool(eff.iloc[0].sd_from_range_only)   # K7 drops range-only SDs, not this row


@pytest.mark.parametrize("bad, why", [
    ({"is_control": np.nan}, "is_control must be 0 or 1"),
    ({"is_control": 2}, "is_control must be 0 or 1"),
    ({"n": np.nan}, "n is missing"),
    ({"n": 0}, "n is missing or < 1"),
    ({"sd": -3}, "sd must be >= 0"),
    ({"se": -1, "sd": np.nan}, "se must be >= 0"),
    ({"paired_se": 0, "paired_md": 1}, "paired_se must be > 0"),
    ({"corr_prepost": 1.1}, "corr_prepost must be in [-1, 1]"),
    ({"q1": 9, "q3": 4, "median": 6}, "q1 <= median <= q3"),
])
def test_invalid_arm_inputs_exclude_the_study_with_a_reason(bad, why):
    # CR1-1, CR1-2, Astra 11: bad rows are excluded with a reason, other studies still run.
    rows = [{**arm("S1", "pbo", 1, 30, mean=8, sd=3), **bad},
            arm("S1", "nr", 0, 20, mean=10, sd=2),
            arm("S2", "pbo", 1, 30, mean=8, sd=3), arm("S2", "nr", 0, 20, mean=10, sd=2)]
    eff, exc = study_effects(frame(rows), EffectSettings(measure="MD"))
    assert list(eff.study_id) == ["S2"]
    assert list(exc.study_id) == ["S1"] and why in exc.reason.iloc[0]


def test_corr_within_out_of_range_is_excluded():
    rows = [arm("X1", "nr", 0, 16, XO, mean=12, sd=10, corr_within=1.2),
            arm("X1", "pbo", 1, 16, XO, mean=10, sd=3)]
    eff, exc = study_effects(frame(rows), EffectSettings(measure="MD"))
    assert eff.empty and "corr_within must be in [-1, 1]" in exc.reason.iloc[0]


def test_coerce_reports_unparseable_numbers_and_blank_ids():
    df = pd.DataFrame([arm("S1", "pbo", 1, 30, mean="8,5", sd=3),
                       arm("S1", "nr", 0, 20, mean=10, sd=" ")])
    df.insert(1, "outcome", "SYNTHETIC_x")
    df.insert(2, "timepoint", "8 weeks")
    with pytest.raises(ValueError, match=r"S1/pbo: mean='8,5'"):
        ef.coerce_extraction(df)
    ok = df.assign(mean=[8.5, 10])
    assert np.isnan(ef.coerce_extraction(ok).sd.iloc[1])       # blank stays missing
    with pytest.raises(ValueError, match="not finite"):
        ef.coerce_extraction(ok.assign(sd=[3, "inf"]))
    with pytest.raises(ValueError, match="blank study_id or arm"):
        ef.coerce_extraction(ok.assign(arm=["pbo", None]))


def test_control_row_paired_data_is_refused():
    # CR1-3: paired_md/paired_se belong on the NR row; a control-row value was reused for
    # every comparison (crossover) or replaced the combined-arm MD (parallel).
    par = [arm("S1", "pbo", 1, 30, mean=10, sd=3, paired_md=5, paired_se=1),
           arm("S1", "low", 0, 20, mean=11, sd=2), arm("S1", "high", 0, 20, mean=13, sd=2)]
    xo = [arm("X1", "nr", 0, 16, XO, mean=12, sd=4),
          arm("X1", "pbo", 1, 16, XO, mean=10, sd=4, paired_md=5, paired_se=1)]
    for rows in (par, xo):
        eff, exc = study_effects(frame(rows), EffectSettings(measure="MD"))
        assert eff.empty and "belong on the NR row" in exc.reason.iloc[0]


def test_crossover_md_uses_r_implied_by_reported_paired_se():
    # CR1-4: n 16, SDs 4 and 4, paired_se 1 -> r = 0.5 and Var = paired_se^2 = 1 (as in R),
    # not (16 + 16) / 16 = 2 from the assumed r = 0.
    rows = [arm("X1", "nr", 0, 16, XO, mean=12, sd=4, paired_se=1.0),
            arm("X1", "pbo", 1, 16, XO, mean=10, sd=4)]
    eff, _ = study_effects(frame(rows), EffectSettings(measure="MD", corr_within=0.0))
    assert (eff.iloc[0].yi, eff.iloc[0].vi) == pytest.approx((2.0, 1.0))
    assert not bool(eff.iloc[0].corr_within_assumed)
    assert "derived from paired SE" in eff.iloc[0].method_notes


def test_crossover_rom_bias_correction():
    # CR1-5 / Astra 15: E[ln Xbar] ~ ln mu - sigma^2 / (2 N mu^2) for each period mean whatever
    # the within-person r, so the Lajeunesse 2015 correction has the parallel form with N per
    # period: (SDe^2/(N Me^2) - SDc^2/(N Mc^2)) / 2 = (16/4000 - 16/1000) / 2 = -0.006.
    y0, v0 = ef.crossover_rom(20, 4, 10, 4, 10, 0.5)
    y1, v1 = ef.crossover_rom(20, 4, 10, 4, 10, 0.5, bias_correction=True)
    assert y1 - y0 == pytest.approx(-0.006) and v1 == v0
    rows = [arm("X1", "nr", 0, 10, XO, mean=20, sd=4), arm("X1", "pbo", 1, 10, XO, mean=10, sd=4)]
    on, _ = study_effects(frame(rows), EffectSettings(measure="ROM", rom_bias_correction=True))
    off, _ = study_effects(frame(rows), EffectSettings(measure="ROM"))
    assert on.iloc[0].yi - off.iloc[0].yi == pytest.approx(-0.006)
    assert "bias-corrected" in on.iloc[0].method_notes


def test_crossover_rom_bias_correction_removes_bias_in_simulation():
    # Independent check: bivariate normal periods, N = 10, r = 0.6. The uncorrected log ratio
    # is biased by about -0.006 + O(1/N^2); the corrected one is unbiased to within MC error.
    rng = np.random.default_rng(20261101)
    reps, n, mu_e, mu_c, sd, r = 200_000, 10, 20.0, 10.0, 4.0, 0.6
    cov = [[sd**2, r * sd**2], [r * sd**2, sd**2]]
    x = rng.multivariate_normal([mu_e, mu_c], cov, size=(reps, n))
    me, mc = x[..., 0].mean(1), x[..., 1].mean(1)
    se_, sc = x[..., 0].std(1, ddof=1), x[..., 1].std(1, ddof=1)
    raw = np.log(me / mc)
    corr = raw + 0.5 * (se_**2 / (n * me**2) - sc**2 / (n * mc**2))
    truth = math.log(mu_e / mu_c)
    mc_se = raw.std() / math.sqrt(reps)
    assert raw.mean() - truth == pytest.approx(0.006, abs=0.001)
    assert abs(corr.mean() - truth) < max(4 * mc_se, 0.0005)
