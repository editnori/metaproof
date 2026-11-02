"""Regression tests for the differences found by validation/crosscheck.py (2026-11-01).

The expected values come from validation/reference_data/, written by
reference_r/make_reference_data.R from metafor 5.2-1 with tight convergence. Those files hold
real published methods datasets (dat.bcg, dat.normand1999, dat.curtis1998), not NR data.
Rows built inside these tests are SYNTHETIC.
"""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from metaproof import diagnostics as dg
from metaproof import effects as ef
from metaproof import pooling as pl
from metaproof.config import load_config
from metaproof.effects import EffectSettings, study_effects

REF = Path(__file__).resolve().parents[1] / "validation" / "reference_data"


def bcg():
    d = pd.read_csv(REF / "bcg_logrr_yi_vi.csv")
    return d.study_id.tolist(), d.yi.to_numpy(float), d.vi.to_numpy(float), d


def synth(rows):
    df = pd.DataFrame(rows)
    df.insert(1, "outcome", "SYNTHETIC_x")
    df.insert(2, "timepoint", "8 weeks")
    df["notes"] = "SYNTHETIC test row"
    return ef.coerce_extraction(df)


def arm(sid, label, ctl, n, design="parallel", **kw):
    return {"study_id": sid, "design": design, "arm": label, "is_control": ctl, "n": n, **kw}


def test_influence_matches_metafor_under_hksj():
    ids, y, v, _ = bcg()
    exp = pd.read_csv(REF / "bcg_expected_influence.csv")
    got = dg.influence(ids, y, v, "REML", hksj=True)
    for ours, theirs in (("rstudent", "rstudent"), ("dffits", "dffits"), ("cook_d", "cooks_d"),
                         ("cov_ratio", "cov_r"), ("dfbetas", "dfbetas"), ("hat", "hat"),
                         ("weight_pct", "weight"), ("tau2_del", "tau2_del")):
        assert got[ours].to_numpy(float) == pytest.approx(exp[theirs].to_numpy(float),
                                                          rel=1e-6, abs=1e-8), ours
    assert list(got.influential) == list(exp.influential)


def test_ml_takes_global_maximum_at_boundary():
    # SYNTHETIC: one dominant study (v = 0.016) and four small ones. Fisher scoring from the DL
    # value stops at a local maximum near tau2 = 0.765; the ML likelihood is higher at 0.
    y = np.array([0.8, 3.0, -1.5, 2.5, 4.0])
    v = np.array([0.016, 32 / 15, 8 / 3, 16 / 9, 16 / 7])
    X = np.ones((5, 1))
    assert pl.loglik(y, v, X, 0.0, False) > pl.loglik(y, v, X, 0.7647, False)
    f = pl.fit(y, v, None, "ML", hksj=False)
    assert f.tau2 == 0.0
    assert any("global maximum" in n for n in f.notes)


def test_pm_and_qprofile_upper_limit_are_exact():
    exp = pd.read_csv(REF / "normand1999_expected_pooled.csv")
    se = pd.read_csv(REF / "normand1999_expected_study_effects.csv")
    md = se[se.outcome == "los_days_md"]
    y, v = md.yi.to_numpy(float), md.vi.to_numpy(float)
    lo, hi = pl.tau2_qprofile_ci(y, v)
    # metafor's default tau2.max cap reported 10 x tau2_DL = 2054.09 for the DL row
    row = exp[(exp.estimator == "DL") & (exp.test == "z")].iloc[0]
    assert hi == pytest.approx(row.tau2_ci_ub, rel=1e-8)
    assert pl.q_gen(y, v, np.ones((len(y), 1)), hi) == pytest.approx(stats.chi2.ppf(0.025, 8))
    pm = exp[(exp.estimator == "PM") & (exp.test == "z")].iloc[0]
    assert pl.fit(y, v, None, "PM", hksj=False).tau2 == pytest.approx(pm.tau2, rel=1e-9)


def test_pi_type_metafor_matches_predict_default():
    _, y, v, _ = bcg()
    exp = pd.read_csv(REF / "bcg_expected_pooled.csv")
    for test in ("z", "knha"):
        row = exp[(exp.estimator == "REML") & (exp.test == test)].iloc[0]
        f = pl.fit(y, v, None, "REML", hksj=(test == "knha"))
        assert pl.prediction_interval(f, 0.95, "metafor") == pytest.approx(
            (row.pi_lb_metafor, row.pi_ub_metafor), rel=1e-7)
        assert pl.prediction_interval(f, 0.95, "riley") == pytest.approx(
            (row.pi_lb, row.pi_ub), rel=1e-7)


def test_q_between_uses_model_based_ses():
    ids, y, v, d = bcg()
    exp = pd.read_csv(REF / "bcg_expected_subgroups.csv")
    eff = pd.DataFrame({"yi": y, "vi": v, "population": d.alloc})
    out = dg.subgroups(eff, "population", "REML", True, 0.95)
    q = exp[exp.test == "Q_between_separate_tau2"].iloc[0]
    assert out.Q_between.iloc[0] == pytest.approx(q.test_stat, rel=1e-8)


def test_single_study_subgroup_has_no_tau2():
    eff = pd.DataFrame({"yi": [1.0, 1.2, 0.7, 3.0], "vi": [0.2, 0.3, 0.25, 0.4],
                        "population": ["SYNTHETIC_a"] * 3 + ["SYNTHETIC_b"]})
    out = dg.subgroups(eff, "population", "REML", True, 0.95).set_index("level")
    b = out.loc["SYNTHETIC_b"]
    assert b.estimate == 3.0 and math.isnan(b.tau2) and math.isnan(b.Q)
    assert b.ci_lb == pytest.approx(3.0 - stats.norm.ppf(0.975) * math.sqrt(0.4))


def test_r_engine_config_vocabulary(tmp_path):
    cfg = load_config(REF / "curtis1998_config.json")
    spec = cfg["outcomes"]["biomass_rom"]
    assert spec["pi_type"] == "riley" and spec["rom_bias_correction"] is False
    p = tmp_path / "c.json"
    p.write_text(json.dumps({"defaults": {"basis": "change", "median_method": "luo2018",
                                          "multi_arm": "split_control", "timepoint": "latest",
                                          "level": 95, "estimator": "EE"}}))
    d = load_config(p)["defaults"]
    assert (d["data_basis"], d["median_method"], d["multi_arm"], d["timepoint"], d["level"],
            d["estimator"]) == ("change", "luo", "split", "last", 0.95, "FE")
    p.write_text(json.dumps({"defaults": {"hksj_adhoc": True}}))
    assert load_config(p)["defaults"]["hksj_rule"] == "truncated"


def test_rom_bias_correction_matches_metafor():
    # metafor escalc("ROM", correct = TRUE) on dat.curtis1998 id 21: yi 0.5501921422790895,
    # vi unchanged 0.0384715432423457 (run in R 4.6.2, metafor 5.2-1)
    yi, vi = ef.rom_parallel(6.8169, 1.769982, 3, 3.945, 1.115797, 5, bias_correction=True)
    assert yi == pytest.approx(0.5501921422790895, rel=1e-12)
    assert vi == pytest.approx(0.0384715432423457, rel=1e-12)


def test_reported_between_group_md_without_control_row():
    _, _, _, d = bcg()
    ext = pd.read_csv(REF / "bcg_logrr_extraction.csv")
    eff, exc = study_effects(ef.coerce_extraction(ext),
                             EffectSettings(measure="MD", extra_columns=("ablat",)))
    assert exc.empty and len(eff) == 13
    assert eff.yi.to_numpy() == pytest.approx(d.yi.to_numpy())
    assert eff.vi.to_numpy() == pytest.approx(d.vi.to_numpy())
    assert eff.ablat.tolist() == d.ablat.tolist()   # extra moderator column carried through
    # but SMD cannot be passed through
    _, exc2 = study_effects(ef.coerce_extraction(ext), EffectSettings(measure="SMD"))
    assert len(exc2) == 13


def test_auto_basis_for_smd_uses_final_values():
    rows = [arm("S1", "pbo", 1, 20, mean=10, sd=3, mean_change=1, sd_change=2),
            arm("S1", "nr", 0, 20, mean=12, sd=3, mean_change=3, sd_change=2)]
    eff, _ = study_effects(synth(rows), EffectSettings(measure="SMD"))
    assert eff.iloc[0].data_basis == "final"
    eff_md, _ = study_effects(synth(rows), EffectSettings(measure="MD"))
    assert eff_md.iloc[0].data_basis == "change"      # reported change SDs: no imputation
    assert eff_md.iloc[0].vi == pytest.approx(4 / 20 + 4 / 20)


def test_sd_from_iqr_without_median():
    # Wan 2014 S2 SD needs only n and the quartiles
    rows = [arm("S1", "pbo", 1, 50, mean=15, q1=10, q3=22),
            arm("S1", "nr", 0, 50, mean=18, sd=9.0)]
    eff, exc = study_effects(synth(rows), EffectSettings(measure="MD"))
    assert exc.empty
    sd = 12 / (2 * stats.norm.ppf((0.75 * 50 - 0.125) / 50.25))
    assert eff.iloc[0].vi == pytest.approx(81 / 50 + sd**2 / 50)
