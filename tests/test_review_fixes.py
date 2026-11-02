"""Regression tests for the methods review of 2026-11-01 (findings F01-F17), plus the
engine-level property tests asked for in F19. Every row built here is SYNTHETIC: invented
numbers, not trial data."""

import json
import math

import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import stats

from metaproof import diagnostics as dg
from metaproof import effects as ef
from metaproof import pooling as pl
from metaproof.config import DEFAULTS, load_config
from metaproof.effects import EffectSettings, study_effects
from metaproof.pipeline import run, select_timepoint, tp_weeks


def arm(sid, label, ctl, n, design="parallel", outcome="SYNTHETIC_x", timepoint="8 weeks",
        **kw):
    return {"study_id": sid, "outcome": outcome, "timepoint": timepoint, "design": design,
            "arm": label, "is_control": ctl, "n": n, "notes": "SYNTHETIC test row", **kw}


def frame(rows):
    return ef.coerce_extraction(pd.DataFrame(rows))


def pair(sid, ne, me, sde, nc, mc, sdc, **kw):
    return [arm(sid, "SYNTHETIC NR", 0, ne, mean=me, sd=sde, **kw),
            arm(sid, "SYNTHETIC placebo", 1, nc, mean=mc, sd=sdc, **kw)]


def three_arm(sid="SYN-3A", **extra):
    return [arm(sid, "NR300", 0, 20, mean=41.0, sd=9.0, dose_mg_day=300, **extra),
            arm(sid, "NR1000", 0, 22, mean=55.0, sd=12.0, dose_mg_day=1000, **extra),
            arm(sid, "PBO", 1, 24, mean=30.0, sd=8.0, dose_mg_day=0, **extra)]


def run_config(tmp_path, rows, outcomes, defaults=None):
    csv = tmp_path / "SYNTHETIC_rows.csv"
    pd.DataFrame(rows).to_csv(csv, index=False)
    cfg = tmp_path / "SYNTHETIC_config.yaml"
    cfg.write_text(yaml.safe_dump({"defaults": defaults or {}, "outcomes": outcomes}))
    out = tmp_path / "out"
    return run(csv, cfg, out), out


# --- F01: HKSJ with few studies and tau2 = 0 ----------------------------------------------

def test_f01_revman_rule_uses_wald_at_k2_and_tau2_zero():
    y, v = np.array([0.50, 0.51]), np.array([0.04, 0.04])
    always = pl.fit(y, v, method="REML", hksj=True, hksj_rule="always")
    revman = pl.fit(y, v, method="REML", hksj=True, hksj_rule="revman")
    fe = pl.fit(y, v, method="FE")
    # the failure the review reproduced: untruncated HKSJ ~4 times narrower than FE
    assert (always.ci_ub[0] - always.ci_lb[0]) < 0.3 * (fe.ci_ub[0] - fe.ci_lb[0])
    assert not revman.hksj
    assert revman.ci_lb[0] == pytest.approx(0.228, abs=1e-3)
    assert revman.ci_ub[0] == pytest.approx(0.782, abs=1e-3)
    assert DEFAULTS["hksj_rule"] == "revman"


def test_f01_truncated_hksj_never_narrower_than_wald():
    y, v = np.array([0.50, 0.51, 0.49]), np.array([0.04, 0.05, 0.03])
    tr = pl.fit(y, v, method="REML", hksj=True, hksj_rule="truncated")
    wald = pl.fit(y, v, method="REML", hksj=False)
    assert tr.hksj and tr.s2_hk >= 1.0
    assert (tr.ci_ub[0] - tr.ci_lb[0]) >= (wald.ci_ub[0] - wald.ci_lb[0])


def test_f01_pooled_row_reports_all_three_intervals():
    y, v = np.array([0.2, 0.6, 0.9, 0.1]), np.array([0.04, 0.05, 0.03, 0.06])
    row, _, f = pl.pooled_row(y, v, "REML", hksj_rule="revman")
    for c in ("wald_ci_lb", "wald_ci_ub", "hksj_ci_lb", "hksj_ci_ub", "hksj_trunc_ci_lb",
              "hksj_trunc_ci_ub", "ci_method"):
        assert c in row
    wald = pl.fit(y, v, method="REML", hksj=False)
    assert row["wald_ci_lb"] == pytest.approx(wald.ci_lb[0])
    assert row["ci_method"].startswith("HKSJ") == f.hksj


# --- F02 / F05: multi-arm trials ---------------------------------------------------------

def test_f02_combine_is_default_and_counts_the_control_once():
    assert DEFAULTS["multi_arm"] == "combine"
    eff, _ = study_effects(frame(three_arm()), EffectSettings())
    assert len(eff) == 1
    r = eff.iloc[0]
    assert r.n_c == 24 and r.n_e == 42 and r.arm_e == "NR300+NR1000"


def test_f02_split_reports_studies_and_leave_one_out_drops_whole_trials(tmp_path):
    rows = three_arm()
    rows += pair("SYN-B", 18, 52.0, 13.0, 17, 31.0, 9.0)
    rows += pair("SYN-C", 30, 46.0, 11.0, 29, 33.0, 10.0)
    _, out = run_config(tmp_path, rows, {"SYNTHETIC_split": {
        "outcome": "SYNTHETIC_x", "multi_arm": "split"}})
    pooled = pd.read_csv(out / "SYNTHETIC_split" / "pooled.csv")
    assert pooled.k.iloc[0] == 4 and pooled.k_studies.iloc[0] == 3
    assert "split multi-arm" in pooled.note.iloc[0]
    loo = pd.read_csv(out / "SYNTHETIC_split" / "leave1out.csv")
    assert set(loo.omitted) == {"SYN-3A", "SYN-B", "SYN-C"}
    assert loo.set_index("omitted").loc["SYN-3A", "k"] == 2


def test_f05_dose_never_taken_from_the_control_row():
    split, _ = study_effects(frame(three_arm()), EffectSettings(multi_arm="split"))
    assert sorted(split.dose_mg_day) == [300, 1000]
    comb, _ = study_effects(frame(three_arm()), EffectSettings(multi_arm="combine"))
    r = comb.iloc[0]
    assert pd.isna(r.dose_mg_day) and r.dose_min == 300 and r.dose_max == 1000
    same = three_arm()
    same[1]["dose_mg_day"] = 300
    comb, _ = study_effects(frame(same), EffectSettings(multi_arm="combine"))
    assert comb.iloc[0].dose_mg_day == 300


# --- F03: reported between-group MDs against a shared control -----------------------------

def test_f03_several_paired_mds_against_one_control_are_refused():
    rows = three_arm()
    rows[0].update(paired_md=11.0, paired_se=2.6)
    rows[1].update(paired_md=25.0, paired_se=2.8)
    for mode in ("split", "combine"):
        rr = [dict(r) for r in rows]
        if mode == "combine":   # no arm-level data left, so the arms cannot be combined
            for r in rr[:2]:
                r.pop("mean"), r.pop("sd")
        eff, excl = study_effects(frame(rr), EffectSettings(multi_arm=mode))
        assert len(eff) == 0
        assert "share the control group" in excl.reason.iloc[0]


def test_f03_combine_falls_back_to_arm_level_data():
    rows = three_arm()
    rows[0].update(paired_md=11.0, paired_se=2.6)
    rows[1].update(paired_md=25.0, paired_se=2.8)
    eff, excl = study_effects(frame(rows), EffectSettings(multi_arm="combine"))
    assert len(eff) == 1 and len(excl) == 0 and eff.iloc[0].n_c == 24


# --- F04 / F12: one timepoint per study ----------------------------------------------------

def test_f04_list_rule_keeps_the_first_listed_timepoint_per_study():
    rows = (pair("SYN-T1", 20, 15.0, 4.0, 20, 11.0, 4.2, timepoint="wk6")
            + pair("SYN-T1", 20, 18.0, 4.0, 20, 11.0, 4.2, timepoint="wk12")
            + pair("SYN-T2", 22, 16.0, 4.5, 22, 12.0, 4.0, timepoint="wk12"))
    kept, _, excl = select_timepoint(frame(rows), ["wk6", "wk12"])
    assert kept.groupby("study_id").timepoint.nunique().max() == 1
    assert set(kept[kept.study_id == "SYN-T1"].timepoint) == {"wk6"}
    assert set(kept[kept.study_id == "SYN-T2"].timepoint) == {"wk12"}
    assert not excl


def test_f04_duplicate_arm_labels_are_refused():
    rows = (pair("SYN-T1", 20, 15.0, 4.0, 20, 11.0, 4.2, timepoint="wk6")
            + pair("SYN-T1", 20, 18.0, 4.0, 20, 11.0, 4.2, timepoint="wk12"))
    eff, excl = study_effects(frame(rows), EffectSettings())
    assert len(eff) == 0 and "appears more than once" in excl.reason.iloc[0]


def test_f12_timepoint_units_unknown_labels_and_missing_labels():
    assert tp_weeks("day 21") == pytest.approx(3.0)
    assert tp_weeks("wk6") == 6 and tp_weeks("12 weeks") == 12 and math.isnan(tp_weeks("end"))
    rows = (pair("SYN-U1", 20, 30.0, 6.0, 20, 27.0, 6.0, timepoint="day 21")
            + pair("SYN-U1", 20, 33.0, 6.0, 20, 27.0, 6.0, timepoint="wk6")
            + pair("SYN-U2", 18, 31.0, 6.5, 18, 28.0, 6.0, timepoint="end")
            + pair("SYN-U2", 18, 32.0, 6.5, 18, 28.0, 6.0, timepoint="wk6"))
    kept, _, excl = select_timepoint(frame(rows), "last")
    assert set(kept.timepoint) == {"wk6"} and set(kept.study_id) == {"SYN-U1"}
    assert [e["study_id"] for e in excl] == ["SYN-U2"]
    kept, _, excl = select_timepoint(frame(rows), "day 21")
    assert set(kept.study_id) == {"SYN-U1"} and [e["study_id"] for e in excl] == ["SYN-U2"]


def test_f12_timepoint_weeks_column_wins_and_exclusions_reach_excluded_csv(tmp_path):
    rows = (pair("SYN-V", 19, 30.5, 6.0, 19, 28.0, 6.0, timepoint="visit 1", timepoint_weeks=4)
            + pair("SYN-V", 19, 33.5, 6.0, 19, 28.0, 6.0, timepoint="visit 2",
                   timepoint_weeks=12)
            + pair("SYN-W", 18, 31.0, 6.5, 18, 28.0, 6.0, timepoint="end")
            + pair("SYN-W", 18, 32.0, 6.5, 18, 28.0, 6.0, timepoint="wk6"))
    _, out = run_config(tmp_path, rows, {"SYNTHETIC_tp": {"outcome": "SYNTHETIC_x"}})
    eff = pd.read_csv(out / "SYNTHETIC_tp" / "study_effects.csv")
    assert eff.yi.tolist() == [pytest.approx(5.5)]
    excl = pd.read_csv(out / "SYNTHETIC_tp" / "excluded.csv")
    assert excl.study_id.tolist() == ["SYN-W"]
    assert "excluded SYN-W" in (out / "SYNTHETIC_tp" / "log.txt").read_text()


# --- F06: minimum k for pooling -------------------------------------------------------------

def test_f06_two_studies_are_shown_but_not_pooled(tmp_path):
    rows = (pair("SYN-A", 20, 10.0, 3.0, 20, 8.0, 3.2)
            + pair("SYN-B", 30, 9.0, 4.0, 31, 8.5, 3.6))
    res, out = run_config(tmp_path, rows, {"SYNTHETIC_k2": {"outcome": "SYNTHETIC_x"}},
                          defaults={"bh": True})
    pooled = pd.read_csv(out / "SYNTHETIC_k2" / "pooled.csv")
    assert pooled.estimate.isna().all() and (pooled.ci_method == "not pooled").all()
    assert "min_k_pool" in pooled.note.iloc[0]
    assert (out / "SYNTHETIC_k2" / "forest.png").exists()
    assert len(pd.read_csv(out / "SYNTHETIC_k2" / "study_effects.csv")) == 2
    # nothing pooled, so the BH family is empty; the analysis is still returned (CR3-09)
    assert len(pd.read_csv(out / "bh_across_outcomes.csv")) == 0
    assert res.analysis.tolist() == ["SYNTHETIC_k2"] and res.p_bh.isna().all()


# --- F07: prediction interval and tau2 CI with few studies ----------------------------------

def test_f07_pi_and_tau2_ci_flagged_below_five_studies(tmp_path):
    rows = []
    for i, (me, n) in enumerate([(14.0, 25), (11.0, 40), (16.5, 18)]):
        rows += pair(f"SYN-{i}", n, me, 4.0, n + 1, 10.0, 4.2)
    _, out = run_config(tmp_path, rows, {"SYNTHETIC_k3": {"outcome": "SYNTHETIC_x"}})
    note = pd.read_csv(out / "SYNTHETIC_k3" / "pooled.csv").note.iloc[0]
    assert "not interpretable" in note and "tau2 CI uninformative" in note


def test_f07_pi_type_revman_uses_t_on_k_minus_1():
    y, v = np.array([0.2, 0.6, 0.9, 0.1, 0.4]), np.array([0.04, 0.05, 0.03, 0.06, 0.05])
    f = pl.fit(y, v, method="REML", hksj=True, hksj_rule="revman")
    lb, ub = pl.prediction_interval(f, 0.95, "revman")
    half = stats.t.ppf(0.975, 4) * math.sqrt(f.tau2 + f.se[0] ** 2)
    assert lb == pytest.approx(f.estimate - half) and ub == pytest.approx(f.estimate + half)
    lb_r, _ = pl.prediction_interval(f, 0.95, "riley")
    assert lb_r < lb       # t(3) is wider than t(4)


# --- F08: SMD and change scores, crossover J, mixing warning -------------------------------

def change_rows(sid, sd_change=None, r=None):
    e = arm(sid, "NR", 0, 25, mean=128.0, sd=13.0, baseline_mean=134.0, baseline_sd=14.0,
            mean_change=-6.0)
    c = arm(sid, "PBO", 1, 24, mean=132.0, sd=12.5, baseline_mean=133.5, baseline_sd=13.0,
            mean_change=-1.5)
    if sd_change is not None:
        e["sd_change"], c["sd_change"] = sd_change
    return [e, c]


def test_f08a_smd_never_standardises_by_an_imputed_change_sd():
    s = EffectSettings(measure="SMD", data_basis="change")
    eff, excl = study_effects(frame(change_rows("SYN-C")), s)
    assert len(eff) == 0 and "imputed change SD" in excl.reason.iloc[0]
    g = []
    for r in (0.2, 0.8):
        e2, _ = study_effects(frame(change_rows("SYN-C", sd_change=(8.0, 8.5))),
                              EffectSettings(measure="SMD", data_basis="change", corr_prepost=r))
        g.append(e2.yi.iloc[0])
    assert g[0] == g[1]       # reported change SDs: the assumed r cannot move g
    auto, _ = study_effects(frame(change_rows("SYN-C")), EffectSettings(measure="SMD"))
    assert auto.data_basis.iloc[0] == "final"


def test_f08b_crossover_smd_carries_hedges_j():
    rows = [arm("SYN-X", "NR", 0, 12, "crossover", mean=45.0, sd=9.5, corr_within=0.5),
            arm("SYN-X", "PBO", 1, 12, "crossover", mean=41.5, sd=9.0, corr_within=0.5)]
    with_j, _ = study_effects(frame(rows), EffectSettings(measure="SMD"))
    no_j, _ = study_effects(frame(rows), EffectSettings(measure="SMD", crossover_smd_j=False))
    assert with_j.yi.iloc[0] / no_j.yi.iloc[0] == pytest.approx(ef.hedges_j(11))
    assert ef.hedges_j(11) == pytest.approx(0.93, abs=0.01)


def test_f08c_smd_rows_carry_one_standardiser_or_are_refused(tmp_path, monkeypatch):
    """Handbook 10.5.2: SMDs divided by an SD of change and by an SD of final values (or of
    crossover period measurements) are on different scales and are not pooled together. The
    old test asserted only that a warning was absent, and that warning could never fire."""
    change_only = [arm("SYN-C", "NR", 0, 20, mean_change=4.0, sd_change=3.0),
                   arm("SYN-C", "PBO", 1, 20, mean_change=1.0, sd_change=3.1)]
    xo = [arm("SYN-X", "NR", 0, 12, "crossover", mean=45.0, sd=9.5, corr_within=0.5),
          arm("SYN-X", "PBO", 1, 12, "crossover", mean=41.5, sd=9.0, corr_within=0.5)]
    rows = pair("SYN-P1", 20, 10.0, 3.0, 20, 8.0, 3.2) + change_only + xo
    # auto: every SMD is divided by a final-value SD; a change-only study is excluded
    eff, excl = study_effects(frame(rows), EffectSettings(measure="SMD"))
    assert set(eff.smd_standardiser) == {"sd_final"}
    assert set(eff[eff.design == "parallel"].data_basis) == {"final"}
    assert "SYN-C" in set(excl.study_id)
    # change: the crossover SMD (between-person SD) is refused, the rest use SD of change
    pc = pair("SYN-P2", 20, 10.0, 3.0, 20, 8.0, 3.2, mean_change=2.0, sd_change=2.5)
    eff, excl = study_effects(frame(pc + change_only + xo),
                              EffectSettings(measure="SMD", data_basis="change"))
    assert set(eff.smd_standardiser) == {"sd_change"}
    assert "Handbook 10.5.2" in excl.set_index("study_id").reason["SYN-X"]
    # a mixture that reached pooling anyway stops the run instead of being pooled
    from metaproof import pipeline as pp
    real = pp.study_effects

    def mixed(*a, **k):
        e, x = real(*a, **k)
        e.loc[e.index[0], "smd_standardiser"] = "sd_change"
        return e, x
    monkeypatch.setattr(pp, "study_effects", mixed)
    with pytest.raises(ValueError, match="standardised by different SDs"):
        run_config(tmp_path, pair("SYN-P1", 20, 10.0, 3.0, 20, 8.0, 3.2)
                   + pair("SYN-P3", 22, 11.0, 3.0, 21, 8.5, 3.1),
                   {"SYNTHETIC_smd": {"outcome": "SYNTHETIC_x", "measure": "SMD"}})


# --- F09: crossover rows with within-period change only -------------------------------------

def test_f09_crossover_change_only_row_gets_a_clear_reason():
    rows = [arm("SYN-X", "NR", 0, 12, "crossover", mean_change=4.0, sd_change=6.0),
            arm("SYN-X", "PBO", 1, 12, "crossover", mean_change=0.5, sd_change=6.2)]
    eff, excl = study_effects(frame(rows), EffectSettings())
    assert len(eff) == 0 and "only within-period change" in excl.reason.iloc[0]


def test_f09_marginal_means_estimate_is_flagged():
    rows = [arm("SYN-X", "NR", 0, 12, "crossover", mean=45.0, sd=9.5),
            arm("SYN-X", "PBO", 1, 12, "crossover", mean=41.5, sd=9.0)]
    eff, _ = study_effects(frame(rows), EffectSettings())
    assert "assumes no period effect" in eff.method_notes.iloc[0]
    assert bool(eff.corr_within_assumed.iloc[0])


# --- F10: medians ---------------------------------------------------------------------------

def test_f10_range_only_excluded_by_default_and_flags_written():
    rows = [arm("SYN-M", "NR", 0, 8, median=20.0, min=5.0, max=300.0),
            arm("SYN-M", "PBO", 1, 8, median=15.0, min=4.0, max=40.0)]
    eff, excl = study_effects(frame(rows), EffectSettings(allow_range_only=False))
    assert len(eff) == 0 and "only a range" in excl.reason.iloc[0]
    eff, _ = study_effects(frame(rows), EffectSettings(allow_range_only=True))
    r = eff.iloc[0]
    assert bool(r.converted_from_median) and bool(r.sd_from_range_only)
    assert abs(r.median_asymmetry) > 0.8     # (5 + 300 - 40) / 295: heavy right skew


def test_f10_k5_and_k7_sensitivities_run_automatically(tmp_path):
    rows = [arm("SYN-M1", "NR", 0, 18, median=66.0, q1=52.0, q3=90.0),
            arm("SYN-M1", "PBO", 1, 17, median=31.0, q1=25.0, q3=38.0)]
    rows += pair("SYN-M2", 22, 64.0, 18.0, 22, 33.0, 9.0)
    rows += pair("SYN-M3", 20, 60.0, 15.0, 21, 35.0, 10.0)
    # K7 input: a change mean with only final and baseline SDs, so the SD of change is imputed
    rows += [arm("SYN-I", "NR", 0, 19, mean_change=31.0, sd=14.0, baseline_sd=12.0,
                 baseline_mean=30.0),
             arm("SYN-I", "PBO", 1, 19, mean_change=2.0, sd=11.0, baseline_sd=10.0,
                 baseline_mean=31.0)]
    rows += pair("SYN-M4", 25, 62.0, 16.0, 24, 34.0, 11.0)
    _, out = run_config(tmp_path, rows, {"SYNTHETIC_med": {"outcome": "SYNTHETIC_x"}})
    sens = pd.read_csv(out / "SYNTHETIC_med" / "sensitivity.csv").set_index("parameter")
    eff = pd.read_csv(out / "SYNTHETIC_med" / "study_effects.csv")
    assert eff.converted_from_median.tolist() == [True, False, False, False, False]
    assert eff.set_index("study_id").imputed_any.to_dict() == {
        "SYN-M1": False, "SYN-M2": False, "SYN-M3": False, "SYN-I": True, "SYN-M4": False}
    # K5b drops the converted study, K7 the imputed one; each row is the pool of what is left
    for param, drop in (("exclude_converted_from_median", "SYN-M1"),
                        ("exclude_imputed_or_digitised", "SYN-I")):
        row = sens.loc[param]
        keep = eff[eff.study_id != drop]
        ref, _, _ = pl.pooled_row(keep.yi.to_numpy(), keep.vi.to_numpy(), "REML",
                                  hksj_rule="revman")
        assert row.k_studies == 4 and drop in row.note
        assert row.estimate == pytest.approx(ref["estimate"], rel=1e-10)
        assert row.se == pytest.approx(ref["se"], rel=1e-10)
    assert "median_method" in sens.index


# --- F11: ratio of means ---------------------------------------------------------------------

def test_f11_geary_flag_and_exclude():
    rows = pair("SYN-G", 10, 2.0, 3.0, 10, 1.2, 2.0)
    eff, _ = study_effects(frame(rows), EffectSettings(measure="ROM"))
    assert eff.geary_ok.iloc[0] is False or not eff.geary_ok.iloc[0]
    assert "Geary" in eff.method_notes.iloc[0]
    eff, excl = study_effects(frame(rows), EffectSettings(measure="ROM", rom_geary="exclude"))
    assert len(eff) == 0 and "Geary" in excl.reason.iloc[0]


def test_f11_rom_bias_correction_runs_as_a_sensitivity(tmp_path):
    rows = []
    for i, (me, mc) in enumerate([(40.0, 30.0), (38.0, 29.0), (42.0, 31.0)]):
        rows += pair(f"SYN-R{i}", 12, me, 9.0, 12, mc, 8.0)
    _, out = run_config(tmp_path, rows, {"SYNTHETIC_rom": {"outcome": "SYNTHETIC_x",
                                                           "measure": "ROM"}})
    sens = pd.read_csv(out / "SYNTHETIC_rom" / "sensitivity.csv")
    assert "rom_bias_correction" in set(sens.parameter)


# --- F13: Benjamini-Hochberg --------------------------------------------------------------

def test_f13_bh_is_opt_in_and_one_primary_analysis_per_outcome(tmp_path):
    assert DEFAULTS["bh"] is False
    rows = []
    for i, (me, n) in enumerate([(14.0, 25), (11.0, 40), (16.5, 18)]):
        rows += pair(f"SYN-{i}", n, me, 4.0, n + 1, 10.0, 4.2)
    outcomes = {"SYNTHETIC_a": {"outcome": "SYNTHETIC_x"},
                "SYNTHETIC_b": {"outcome": "SYNTHETIC_x", "timepoint": "first"}}
    with pytest.raises(ValueError, match="both primary"):
        run_config(tmp_path, rows, outcomes, defaults={"bh": True})
    outcomes["SYNTHETIC_b"]["primary"] = False
    res, out = run_config(tmp_path, rows, outcomes, defaults={"bh": True})
    bh = pd.read_csv(out / "bh_across_outcomes.csv")
    assert bh.analysis.tolist() == ["SYNTHETIC_a"] and "significant_bh" not in bh.columns


# --- F14: small-study tests --------------------------------------------------------------

def test_f14_egger_gates(tmp_path):
    with pytest.raises(ValueError, match="egger_allow_small_k"):
        load_config_from(tmp_path, {"defaults": {"egger_min_k": 5}})
    y, v = np.linspace(0.1, 1.0, 12), np.linspace(0.01, 0.2, 12)
    assert dg.egger(y, v, 10, k_studies=9).p.isna().all()
    assert dg.egger(y, v, 10, k_studies=10).p.notna().all()
    prec = np.sqrt(1 / np.arange(10, 22) + 1 / np.arange(10, 22))
    assert "Pustejovsky" in dg.egger(y, v, 10, precision_se=prec).variant.iloc[0]
    flat = dg.egger(y, np.full(12, 0.05) * np.linspace(1, 1.2, 12), 10)
    assert "similar size" in flat.note.iloc[0]


def load_config_from(tmp_path, raw):
    p = tmp_path / "SYNTHETIC_cfg.yaml"
    p.write_text(yaml.safe_dump(raw))
    return load_config(p)


# --- F15: subgroups and meta-regression -----------------------------------------------------

def k_rows(k, pops=("a", "b")):
    rows = []
    for i in range(k):
        rows += pair(f"SYN-{i:02d}", 20 + i, 12.0 + 0.4 * i, 4.0, 20 + i, 10.0, 4.0,
                     population=f"SYNTHETIC_{pops[i % len(pops)]}", dose_mg_day=100 * (i + 1))
    return rows


def test_f15_subgroups_and_metareg_respect_protocol_gates(tmp_path):
    _, out = run_config(tmp_path, k_rows(6), {"SYNTHETIC_sg": {
        "outcome": "SYNTHETIC_x", "subgroups": ["population"], "moderators": ["dose_mg_day"]}})
    sg = pd.read_csv(out / "SYNTHETIC_sg" / "subgroups.csv")
    assert "not run" in sg.note.iloc[0]
    mr = pd.read_csv(out / "SYNTHETIC_sg" / "metareg.csv")
    assert "not run" in mr.note.iloc[0]
    assert DEFAULTS["subgroup_method"] == "common"


def test_f15_common_tau2_subgroups_and_permutation_p():
    eff, _ = study_effects(frame(k_rows(12)), EffectSettings(extra_columns=("population",)))
    sg = dg.subgroups(eff, "population", "REML", True, 0.95, tau2_mode="common", min_k=10,
                      min_per_level=2, hksj_rule="revman", study_col="source_study")
    assert sg.tau2.nunique() == 1 and (sg.tau2_mode == "common").all()
    # common tau2 is the tau2 of the meta-regression on the subgroup factor
    mr_pop = dg.metareg(eff, ["population"], "REML", True, 0.95, min_k=10,
                        hksj_rule="revman", study_col="source_study")
    assert sg.tau2.iloc[0] == pytest.approx(mr_pop.tau2.iloc[0], rel=1e-8, abs=1e-12)


def test_f15_permutation_p_matches_metafor_exact_permutest():
    """SYNTHETIC k = 6, so metafor can enumerate all 720 permutations: permutest(rma(yi, vi,
    mods = ~dose, method = 'REML', test = 'knha' | 'z'), exact = TRUE)$QMp, metafor 5.2.1,
    control threshold 1e-10, gives 68/720 (knha) and 72/720 (z). The Monte Carlo p must be a
    multiple of 1/(B + 1) within 3 Monte Carlo SEs of the exact value (Higgins & Thompson
    2004); the old test accepted any value in (0, 1]."""
    eff = pd.DataFrame({"study_id": [f"SYN-{i}" for i in range(6)],
                        "yi": [0.10, 0.55, 0.20, 0.90, 0.45, 1.20],
                        "vi": [0.04, 0.06, 0.05, 0.08, 0.03, 0.07],
                        "dose_mg_day": [100.0, 200, 300, 400, 500, 600]})
    B = 1999
    for hksj, exact in ((True, 68 / 720), (False, 72 / 720)):
        mr = dg.metareg(eff, ["dose_mg_day"], "REML", hksj, 0.95, min_k=0, permutations=B)
        p = mr.perm_p.iloc[0]
        assert p * (B + 1) == pytest.approx(round(p * (B + 1)), abs=1e-9)
        assert abs(p - exact) <= 3 * math.sqrt(exact * (1 - exact) / B)
        assert mr.QM.iloc[0] == pytest.approx(4.4696354622846242 if hksj else 4.2992808371723505,
                                              rel=1e-8)


# --- F17: SE of change ----------------------------------------------------------------------

def test_f17_se_change_is_converted_by_the_engine():
    rows = [arm("SYN-S", "NR", 0, 16, mean_change=-5.0, se_change=2.0),
            arm("SYN-S", "PBO", 1, 16, mean_change=-1.5, se_change=1.9)]
    eff, _ = study_effects(frame(rows), EffectSettings(data_basis="change"))
    assert eff.vi.iloc[0] == pytest.approx(2.0**2 + 1.9**2)
    assert not eff.sd_change_imputed.iloc[0]


# --- F19: engine-level properties -----------------------------------------------------------

@pytest.mark.parametrize("mode", ["split", "combine"])
def test_f19_control_participants_counted_once(mode):
    rng = np.random.default_rng(7)
    for _ in range(25):
        n_arms = int(rng.integers(2, 5))
        rows = [arm("SYN-R", f"NR{j}", 0, int(rng.integers(8, 40)),
                    mean=float(rng.uniform(20, 60)), sd=float(rng.uniform(4, 15)),
                    dose_mg_day=float(100 * (j + 1))) for j in range(n_arms)]
        n_c = int(rng.integers(8, 40))
        rows.append(arm("SYN-R", "PBO", 1, n_c, mean=30.0, sd=8.0, dose_mg_day=0.0))
        for measure in ("MD", "SMD", "ROM"):
            eff, excl = study_effects(frame(rows), EffectSettings(measure=measure,
                                                                  multi_arm=mode))
            assert len(excl) == 0
            assert eff.n_c.sum() == pytest.approx(n_c)
            assert eff.n_e.sum() == pytest.approx(sum(r["n"] for r in rows[:-1]))
            assert (eff.dose_min.fillna(1) > 0).all()


def test_f19_paired_only_path_never_reuses_the_control():
    rows = [arm("SYN-P", "NR", 0, 30, paired_md=-3.0, paired_se=1.2)]
    eff, _ = study_effects(frame(rows), EffectSettings())
    assert len(eff) == 1 and pd.isna(eff.n_c.iloc[0])
    rows.append(arm("SYN-P", "NR high", 0, 30, paired_md=-5.0, paired_se=1.3))
    eff, excl = study_effects(frame(rows), EffectSettings())
    assert len(eff) == 0 and len(excl) == 1


def test_f19_every_rule_keeps_one_timepoint_per_study():
    rng = np.random.default_rng(11)
    labels = ["wk2", "wk4", "day 21", "8 weeks", "3 months", "wk12"]
    rows = []
    for s in range(15):
        for tp in rng.choice(labels, size=int(rng.integers(1, 4)), replace=False):
            rows += pair(f"SYN-{s}", 20, 12.0, 4.0, 20, 10.0, 4.0, timepoint=str(tp))
    df = frame(rows)
    for rule in ("last", "first", ["wk4", "wk12", "8 weeks"], "wk12"):
        kept, _, excl = select_timepoint(df, rule)
        assert (kept.groupby("study_id").timepoint.nunique() <= 1).all()
        assert set(kept.study_id) | {e["study_id"] for e in excl} == set(df.study_id)


# === Second review, 2026-11-01: pipeline, config, CLI and plots (CR3-*, Astra) ==============

def three_pairs(**kw):
    rows = []
    for i, (me, n) in enumerate([(14.0, 25), (11.0, 40), (16.5, 18)]):
        rows += pair(f"SYN-{i}", n, me, 4.0, n + 1, 10.0, 4.2, **kw)
    return rows


def test_cr3_01_hyphenated_timepoint_labels():
    # '12-month' was read as 12 weeks and 'Wk-12' as -12 weeks, in both engines
    assert tp_weeks("12-month") == pytest.approx(365.25 / 7)
    assert tp_weeks("3-mo") == pytest.approx(365.25 / 28)
    assert tp_weeks("6-months") == pytest.approx(365.25 / 14)
    assert tp_weeks("Wk-12") == 12 and tp_weeks("week-12") == 12
    assert tp_weeks("day-21") == pytest.approx(3.0)
    for label, wk in [("wk6", 6), ("12 weeks", 12), ("day 21", 3), ("8", 8), ("12wk", 12),
                      ("W12", 12), ("Week 12 (end)", 12), ("month 6", 365.25 / 14),
                      ("3 months", 365.25 / 28), ("1.5 months", 1.5 * 365.25 / 84)]:
        assert tp_weeks(label) == pytest.approx(wk), label
    rows = (pair("SYN-H", 20, 30.0, 6.0, 20, 27.0, 6.0, timepoint="6 months")
            + pair("SYN-H", 20, 33.0, 6.0, 20, 27.0, 6.0, timepoint="12-month"))
    kept, _, _ = select_timepoint(frame(rows), "last")
    assert set(kept.timepoint) == {"12-month"}
    kept, _, excl = select_timepoint(frame(rows), 12)
    assert kept.empty and [e["study_id"] for e in excl] == ["SYN-H"]
    rows = (pair("SYN-W", 20, 30.0, 6.0, 20, 27.0, 6.0, timepoint="Wk-4")
            + pair("SYN-W", 20, 33.0, 6.0, 20, 27.0, 6.0, timepoint="Wk-12"))
    kept, _, _ = select_timepoint(frame(rows), "last")
    assert set(kept.timepoint) == {"Wk-12"}


def test_cr3_08_timepoint_rule_is_case_insensitive(tmp_path):
    for raw, want in [("Last", "last"), ("FIRST", "first"), ("Latest", "last"),
                      ("Earliest", "first"), ("wk12", "wk12"), ("Visit A", "Visit A")]:
        assert load_config_from(tmp_path, {"defaults": {"timepoint": raw}}
                                )["defaults"]["timepoint"] == want


def test_astra16_timepoint_weeks_conflicts_are_not_silent():
    # an arbitrary visit label resolved through timepoint_weeks, for numeric and last rules
    rows = (pair("SYN-V", 20, 30.0, 6.0, 20, 27.0, 6.0, timepoint="visit A", timepoint_weeks=12)
            + pair("SYN-V", 20, 31.0, 6.0, 20, 27.0, 6.0, timepoint="visit B", timepoint_weeks=4))
    kept, _, _ = select_timepoint(frame(rows), 12)
    assert set(kept.timepoint) == {"visit A"}
    kept, _, _ = select_timepoint(frame(rows), "last")
    assert set(kept.timepoint) == {"visit A"}
    # the label says 12 weeks, the column says 8: the column wins and the log says so
    rows = (pair("SYN-C", 20, 30.0, 6.0, 20, 27.0, 6.0, timepoint="12 weeks", timepoint_weeks=8)
            + pair("SYN-C", 20, 31.0, 6.0, 20, 27.0, 6.0, timepoint="wk10"))
    kept, notes, _ = select_timepoint(frame(rows), "last")
    assert set(kept.timepoint) == {"wk10"}
    assert any("12 weeks" in n and "timepoint_weeks 8" in n for n in notes)
    # rows of one label that disagree on timepoint_weeks cannot be placed
    rows = [arm("SYN-D", "NR", 0, 20, mean=30.0, sd=6.0, timepoint="end", timepoint_weeks=8),
            arm("SYN-D", "PBO", 1, 20, mean=27.0, sd=6.0, timepoint="end", timepoint_weeks=12)]
    rows += pair("SYN-D", 20, 31.0, 6.0, 20, 27.0, 6.0, timepoint="wk4")
    kept, _, excl = select_timepoint(frame(rows), "last")
    assert kept.empty and "timepoint_weeks" in excl[0]["reason"]


def test_astra13_quoted_booleans_are_parsed_not_truthy(tmp_path):
    d = load_config_from(tmp_path, {"defaults": {
        "hksj": "false", "rom_bias_correction": "FALSE", "prediction_interval": "False",
        "hksj_adhoc": "false", "metareg_joint": "true", "allow_range_only": 0}})["defaults"]
    assert (d["hksj"], d["rom_bias_correction"], d["prediction_interval"], d["hksj_adhoc"],
            d["metareg_joint"], d["allow_range_only"]) == (False, False, False, False, True,
                                                            False)
    assert d["hksj_rule"] == "revman"
    with pytest.raises(ValueError, match="hksj"):
        load_config_from(tmp_path, {"defaults": {"hksj": "maybe"}})


def test_cr3_13_config_gaps(tmp_path):
    with pytest.raises(ValueError, match="defaults.outcome"):
        load_config_from(tmp_path, {"defaults": {"outcome": "SYNTHETIC_x"}})
    with pytest.raises(ValueError, match="bh"):
        load_config_from(tmp_path, {"outcomes": {"a": {"bh": True}}})
    with pytest.raises(ValueError, match="hksj_adhoc"):
        load_config_from(tmp_path, {"defaults": {"hksj_adhoc": True},
                                    "outcomes": {"a": {"hksj_rule": "revman"}}})
    cfg = load_config_from(tmp_path, {"defaults": {"hksj_adhoc": True}, "outcomes": {"a": {}}})
    assert cfg["outcomes"]["a"]["hksj_rule"] == "truncated"
    with pytest.raises(ValueError, match="combine_arms"):
        load_config_from(tmp_path, {"defaults": {"combine_arms": ["x"]}})
    with pytest.raises(ValueError, match="corr_prepost"):
        load_config_from(tmp_path, {"defaults": {"corr_prepost": "half"}})
    with pytest.raises(ValueError, match="outcomes"):
        load_config_from(tmp_path, {"outcomes": [{"outcome": "a"}]})
    with pytest.raises(ValueError, match="ROM needs final values"):   # as run_meta.R refuses
        load_config_from(tmp_path, {"defaults": {"measure": "ROM", "data_basis": "change"}})


def test_subgroup_or_moderator_column_missing_from_the_csv_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="populaton"):
        run_config(tmp_path, three_pairs(), {"a": {"outcome": "SYNTHETIC_x",
                                                   "subgroups": ["populaton"]}})


def test_cr3_07_sensitivity_correlations_are_range_checked(tmp_path):
    with pytest.raises(ValueError, match="corr_prepost_sensitivity"):
        load_config_from(tmp_path, {"defaults": {"corr_prepost_sensitivity": [0.4, 80]}})
    with pytest.raises(ValueError, match="corr_within_sensitivity"):
        load_config_from(tmp_path, {"defaults": {"corr_within_sensitivity": ["high"]}})


def crossover_pairs(sds=(9.0, 9.0)):
    rows = []
    for i, (me, n) in enumerate([(45.0, 12), (47.0, 16), (44.0, 10)]):
        rows += [arm(f"SYN-X{i}", "NR", 0, n, "crossover", mean=me, sd=sds[0]),
                 arm(f"SYN-X{i}", "PBO", 1, n, "crossover", mean=41.0, sd=sds[1])]
    return rows


def test_astra14_assumed_crossover_r_gets_sensitivity_rows_by_default(tmp_path):
    _, out = run_config(tmp_path, crossover_pairs(), {"SYNTHETIC_xo": {"outcome": "SYNTHETIC_x"}})
    sens = pd.read_csv(out / "SYNTHETIC_xo" / "sensitivity.csv")
    assert sorted(sens[sens.parameter == "corr_within"].value.astype(float)) == [0.2, 0.5, 0.8]
    # parallel trials assume no crossover r, so the automatic rows are not added
    _, out = run_config(tmp_path, three_pairs(), {"SYNTHETIC_par": {"outcome": "SYNTHETIC_x"}})
    sens = pd.read_csv(out / "SYNTHETIC_par" / "sensitivity.csv")
    assert "corr_within" not in set(sens.parameter)
    # an explicit empty list is honoured, and the log says the assumption is untested
    _, out = run_config(tmp_path, crossover_pairs(), {"SYNTHETIC_xo": {
        "outcome": "SYNTHETIC_x", "corr_within_sensitivity": []}})
    assert "WARNING: corr_within" in (out / "SYNTHETIC_xo" / "log.txt").read_text()


def test_astra14_sensitivity_scenarios_never_disappear(tmp_path):
    # r = 1 with equal period SDs gives zero variance, so every study drops out
    _, out = run_config(tmp_path, crossover_pairs(), {"SYNTHETIC_xo": {
        "outcome": "SYNTHETIC_x", "corr_within_sensitivity": [0.2, 1.0]}})
    sens = pd.read_csv(out / "SYNTHETIC_xo" / "sensitivity.csv")
    sens = sens[sens.parameter == "corr_within"]
    r1 = sens[sens.value.astype(float) == 1.0]
    assert len(r1) == 1 and r1.k_studies.iloc[0] == 0
    assert "SYN-X0" in r1.note.iloc[0] and "drops" in r1.note.iloc[0]


def test_cr3_03_config_outcome_not_in_the_csv_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="SYNTHETIC_X.*SYNTHETIC_x"):
        run_config(tmp_path, three_pairs(), {"a": {"outcome": "SYNTHETIC_X"}})
    res, _ = run_config(tmp_path, three_pairs(), {"a": {"outcome": " SYNTHETIC_x "}})
    assert res.k.tolist() == [3]


def test_cr3_04_rerun_replaces_the_analysis_folder(tmp_path):
    outcomes = {"SYNTHETIC_a": {"outcome": "SYNTHETIC_x"},
                "SYNTHETIC_b": {"outcome": "SYNTHETIC_x", "timepoint": "first"}}
    _, out = run_config(tmp_path, three_pairs(), outcomes)
    assert (out / "SYNTHETIC_a" / "forest.png").exists()
    # same folder, now with a timepoint no study reports: k = 0, and SYNTHETIC_b is gone
    res, out = run_config(tmp_path, three_pairs(),
                          {"SYNTHETIC_a": {"outcome": "SYNTHETIC_x", "timepoint": "wk99"}})
    folder = out / "SYNTHETIC_a"
    assert not (folder / "forest.png").exists() and not (folder / "leave1out.csv").exists()
    pooled = pd.read_csv(folder / "pooled.csv")
    assert "no analysable studies" in pooled.note.iloc[0]
    assert not (out / "SYNTHETIC_b").exists()
    assert res.k.tolist() == [0] and not res.pooled.iloc[0]
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert len(manifest["extraction_sha256"]) == 64 and len(manifest["config_sha256"]) == 64
    assert [a["analysis"] for a in manifest["analyses"]] == ["SYNTHETIC_a"]
    assert manifest["removed_folders"] == ["SYNTHETIC_b"]
    assert not [p for p in out.iterdir() if p.name.startswith(".metaproof")]


def test_cr3_04_folders_metaproof_did_not_write_are_left_alone(tmp_path):
    (tmp_path / "out" / "my_notes").mkdir(parents=True)
    (tmp_path / "out" / "my_notes" / "n.txt").write_text("keep")
    _, out = run_config(tmp_path, three_pairs(), {"a": {"outcome": "SYNTHETIC_x"}})
    assert (out / "my_notes" / "n.txt").read_text() == "keep"


@pytest.mark.parametrize("names", [("NAD+", "NAD"), ("a/b", "a b"), ("SBP", "sbp")])
def test_cr3_05_analysis_names_that_share_a_folder_are_refused(tmp_path, names):
    outcomes = {n: {"outcome": "SYNTHETIC_x"} for n in names}
    with pytest.raises(ValueError, match="same folder"):
        run_config(tmp_path, three_pairs(), outcomes)


@pytest.mark.parametrize("name", ["..", ".", "summary_across_analyses.csv", "run_manifest.json"])
def test_astra23_analysis_names_cannot_escape_or_overwrite_run_files(tmp_path, name):
    with pytest.raises(ValueError, match="analysis name"):
        run_config(tmp_path, three_pairs(), {name: {"outcome": "SYNTHETIC_x"}})
    assert not (tmp_path / "study_effects.csv").exists()


def test_cr3_05_names_with_slashes_stay_inside_out(tmp_path):
    _, out = run_config(tmp_path, three_pairs(), {"a/b": {"outcome": "SYNTHETIC_x"}})
    assert (out / "a_b" / "pooled.csv").exists()
    with pytest.raises(ValueError, match="analysis name"):   # would be a hidden folder
        run_config(tmp_path, three_pairs(), {"../escape": {"outcome": "SYNTHETIC_x"}})
    assert not (tmp_path / "escape").exists()


@pytest.mark.parametrize("bad, words", [
    ({"n": None}, "n is required"),
    ({"n": "20 (ITT)"}, "n must be a number"),
    ({"mean": "12,5"}, "mean must be a number"),
    ({"is_control": "yes"}, "is_control must be 0 or 1"),
    ({"sd": -3.0}, "sd must be"),
    ({"corr_prepost": 1.1}, "corr_prepost must be <= 1"),
])
def test_cr3_02_06_bad_cells_stop_the_run_with_the_checker_message(tmp_path, bad, words):
    rows = three_pairs()
    rows[0] = {**rows[0], **bad}
    with pytest.raises(ValueError, match=r"row 2: " + words.replace("(", r"\(")):
        run_config(tmp_path, rows, {"a": {"outcome": "SYNTHETIC_x"}})
    assert "row 2" in (tmp_path / "out" / "extraction_check.txt").read_text()


def test_pipeline_reports_checker_findings_it_can_handle(tmp_path):
    # a study with two control rows is a checker error the engine handles per study
    rows = three_pairs() + [arm("SYN-0", "PBO2", 1, 20, mean=10.0, sd=4.0)]
    res, out = run_config(tmp_path, rows, {"a": {"outcome": "SYNTHETIC_x"}})
    text = (out / "extraction_check.txt").read_text()
    assert "needs exactly one control row" in text
    assert len(res) == 1


def test_cr3_09_bh_run_returns_every_analysis(tmp_path):
    rows = three_pairs()
    outcomes = {"SYNTHETIC_a": {"outcome": "SYNTHETIC_x"},
                "SYNTHETIC_b": {"outcome": "SYNTHETIC_x", "timepoint": "first",
                                "primary": False}}
    res, out = run_config(tmp_path, rows, outcomes, defaults={"bh": True})
    assert res.analysis.tolist() == ["SYNTHETIC_a", "SYNTHETIC_b"]
    assert res.p_bh.notna().tolist() == [True, False]


def test_cr3_10_modified_egger_precision_for_crossover_smd():
    from metaproof.pipeline import egger_precision
    rows = [arm("SYN-X", "NR", 0, 12, "crossover", mean=45.0, sd=9.0),
            arm("SYN-X", "PBO", 1, 12, "crossover", mean=41.0, sd=9.0)]
    rows += pair("SYN-P", 20, 14.0, 4.0, 22, 10.0, 4.2)
    eff, _ = study_effects(frame(rows), EffectSettings(measure="SMD", corr_within=0.5))
    prec = egger_precision(eff)
    xo = (eff.design == "crossover").to_numpy()
    assert prec[xo][0] == pytest.approx(math.sqrt(2 * (1 - 0.5) / 12))
    assert prec[~xo][0] == pytest.approx(math.sqrt(1 / 20 + 1 / 22))


def test_cr3_04_a_failure_part_way_leaves_the_earlier_run_whole(tmp_path, monkeypatch):
    from metaproof import pipeline
    outcomes = {"SYNTHETIC_a": {"outcome": "SYNTHETIC_x"},
                "SYNTHETIC_b": {"outcome": "SYNTHETIC_x", "timepoint": "first"}}
    _, out = run_config(tmp_path, three_pairs(), outcomes)
    before = (out / "SYNTHETIC_a" / "pooled.csv").read_text()
    real = pipeline.analyse

    def flaky(name, *a):
        if name == "SYNTHETIC_b":
            raise RuntimeError("SYNTHETIC failure")
        return real(name, *a)
    monkeypatch.setattr(pipeline, "analyse", flaky)
    with pytest.raises(RuntimeError):
        run_config(tmp_path, three_pairs(), {**outcomes, "SYNTHETIC_a": {
            "outcome": "SYNTHETIC_x", "estimator": "DL"}})
    assert (out / "SYNTHETIC_a" / "pooled.csv").read_text() == before
    assert not [p for p in out.iterdir() if p.name.startswith(".metaproof")]


def test_analysis_folders_get_the_normal_directory_mode(tmp_path):
    """Staging folders come from mkdtemp (mode 0700); the finished analysis folder must have
    the mode an ordinary mkdir gives, so results can be shared like any other folder."""
    import os
    import stat
    mask = os.umask(0)
    os.umask(mask)
    _, out = run_config(tmp_path, three_pairs(), {"SYNTHETIC_a": {"outcome": "SYNTHETIC_x"}})
    assert stat.S_IMODE((out / "SYNTHETIC_a").stat().st_mode) == 0o777 & ~mask


# --- Astra re-review, 2026-11-02 ---------------------------------------------------------

def test_astra2_p1_unowned_output_folder_is_not_deleted(tmp_path):
    """A folder under --out with an analysis's name that metaproof did not write is left
    alone, and the run stops with a clear message."""
    from pathlib import Path
    import metaproof.pipeline as pp
    root = Path(__file__).resolve().parents[1] / "examples" / "nicotinamide-riboside"
    data, cfg, out = (root / "SYNTHETIC_nr_like_extraction.csv",
                      root / "SYNTHETIC_nr_like_config.yaml", tmp_path / "out")
    run(data, cfg, out)
    victim = sorted(p for p in out.iterdir() if (p / pp.MARKER).is_file())[0]
    pp._remove(victim)
    victim.mkdir()
    (victim / "my_notes.txt").write_text("keep me")
    with pytest.raises(pp.UserError, match="did not write it"):
        run(data, cfg, out)
    assert (victim / "my_notes.txt").read_text() == "keep me"


def test_astra2_p2_zero_hksj_covariance_is_a_note_not_a_crash():
    rng = np.random.default_rng(1)
    n = 14
    d = pd.DataFrame({"study_id": [f"S{i}" for i in range(n)],
                      "yi": [10.0] * 12 + [3.0, 5.0],
                      "vi": rng.uniform(0.5, 1.5, n),
                      "dose": list(rng.uniform(100, 1000, 12)) + [np.nan, np.nan]})
    res = dg.metareg(d, ["dose"], "REML", True, 0.95, hksj_rule="always")
    assert "omnibus test undefined" in " ".join(res["note"].astype(str))


def test_astra2_p2_crossover_md_with_zero_period_sd_uses_paired_se():
    rows = [arm("X", "NR", 0, 25, "crossover", mean=12.0, sd=0.0, paired_se=0.4),
            arm("X", "PBO", 1, 25, "crossover", mean=10.0, sd=2.0)]
    eff, exc = study_effects(frame(rows), EffectSettings(measure="MD"))
    assert exc.empty
    assert (eff.iloc[0].yi, eff.iloc[0].vi) == pytest.approx((2.0, 0.16))
