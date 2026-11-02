"""Leave-one-out, influence, subgroups, meta-regression, Egger and BH."""

import math

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
from scipy import stats

from metaproof import diagnostics as dg
from metaproof import pooling as pl
from test_pooling import ABLAT, bcg

Y = np.array([0.0, 1.0, 3.0, 0.5, -0.4, 0.8])
V = np.array([0.25, 0.5, 1.0, 0.3, 0.2, 0.4])
IDS = [f"SYNTH-{i}" for i in range(len(Y))]


def test_leave_one_out_equals_refit():
    out = dg.leave_one_out(IDS, Y, V, "REML", True, 0.95)
    for i, r in out.iterrows():
        keep = np.arange(len(Y)) != i
        f = pl.fit(Y[keep], V[keep], method="REML", hksj=True)
        assert r.omitted == IDS[i]
        assert r.estimate == pytest.approx(f.estimate)
        assert r.tau2 == pytest.approx(f.tau2)
        assert r.ci_ub == pytest.approx(f.ci_ub[0])


def test_influence_intercept_only_by_hand():
    out = dg.influence(IDS, Y, V, "REML")
    full = pl.fit(Y, V, method="REML", hksj=False)
    w = 1 / (V + full.tau2)
    # hat_i = w_i / sum w and they sum to p = 1
    assert out.hat.to_numpy() == pytest.approx(w / w.sum())
    assert out.hat.sum() == pytest.approx(1.0)
    for i, r in out.iterrows():
        keep = np.arange(len(Y)) != i
        d = pl.fit(Y[keep], V[keep], method="REML", hksj=False)
        assert r.rstudent == pytest.approx((Y[i] - d.estimate)
                                           / math.sqrt(V[i] + d.tau2 + d.se[0] ** 2))
        assert r.dffits == pytest.approx((full.estimate - d.estimate)
                                         / math.sqrt(r.hat * (d.tau2 + V[i])))
        # Cook: (mu - mu_(-i))^2 * sum w  (X'WX = sum w for intercept-only)
        assert r.cook_d == pytest.approx((full.estimate - d.estimate) ** 2 * w.sum())
        assert r.cov_ratio == pytest.approx(d.se[0] ** 2 / full.se[0] ** 2)
        assert r.tau2_del == pytest.approx(d.tau2)


def test_influence_flags_planted_outlier():
    y = np.array([0.1, 0.2, 0.15, 0.12, 0.18, 3.0, 0.11, 0.14])
    v = np.full(8, 0.01)
    out = dg.influence([str(i) for i in range(8)], y, v, "REML")
    assert out.rstudent.abs().idxmax() == 5
    assert bool(out.influential.iloc[5])


def test_subgroup_q_between_hand():
    # A: y = (1, 1), v = 1 -> mu 1, Var 1/2.  B: y = (3, 3) -> mu 3, Var 1/2.
    # mu_bar = 2; Q_between = 2(1-2)^2 + 2(3-2)^2 = 4 on 1 df, p = 0.0455
    eff = pd.DataFrame({"yi": [1, 1, 3, 3], "vi": [1.0] * 4,
                        "population": ["A", "A", "B", "B"]})
    out = dg.subgroups(eff, "population", "FE", False, 0.95)
    assert list(out.estimate) == pytest.approx([1, 3])
    assert out.Q_between.iloc[0] == pytest.approx(4.0)
    assert out.Q_between_df.iloc[0] == 1
    assert out.Q_between_p.iloc[0] == pytest.approx(0.045500263896358445)


def test_subgroup_single_level():
    eff = pd.DataFrame({"yi": [1, 2], "vi": [1.0, 1.0], "population": ["A", "A"]})
    out = dg.subgroups(eff, "population", "REML", True, 0.95)
    assert math.isnan(out.Q_between.iloc[0])


def test_metareg_bcg_wald_and_r2():
    y, v = bcg()
    eff = pd.DataFrame({"yi": y, "vi": v, "ablat": ABLAT})
    out = dg.metareg(eff, ["ablat"], "REML", False, 0.95)
    slope = out[out.term == "ablat"].iloc[0]
    assert slope.estimate == pytest.approx(-0.0291, abs=5e-5)
    assert slope.QM == pytest.approx(16.3571, abs=2e-3)  # see test_bcg_metaregression_...
    assert slope.QM_test == "chi2"
    # exact REML optimum: validation/reference_data/bcg_expected_metareg.csv (metafor 5.2.1 at
    # threshold 1e-10) gives R^2 75.626622; metafor prints 75.62% (75.6245) at its default
    # threshold 1e-5, where it stops at tau2 0.0763547 short of the optimum 0.076348
    assert slope.R2 == pytest.approx(75.626622, rel=1e-6)


def test_metareg_knapp_hartung():
    y, v = bcg()
    X = np.column_stack([np.ones(13), ABLAT])
    eff = pd.DataFrame({"yi": y, "vi": v, "ablat": ABLAT})
    out = dg.metareg(eff, ["ablat"], "REML", True, 0.95)
    slope = out[out.term == "ablat"].iloc[0]
    f = pl.fit(y, v, X, "REML", hksj=False)
    w = 1 / (v + f.tau2)
    r = y - X @ f.beta
    s2 = np.sum(w * r**2) / (13 - 2)
    assert slope.se == pytest.approx(math.sqrt(s2) * f.se[1])
    assert slope.df == 11
    # single moderator: F = t^2 on (1, k - p)
    assert slope.QM == pytest.approx(slope.stat**2)
    assert slope.QM_p == pytest.approx(slope.pval)


def test_metareg_categorical_dummies():
    eff = pd.DataFrame({"yi": Y, "vi": V, "rob": ["low", "high", "low", "high", "low", "some"]})
    out = dg.metareg(eff, ["rob"], "REML", True, 0.95)
    assert list(out.term) == ["intrcpt", "rob[low]", "rob[some]"]
    assert out.QM_df1.iloc[0] == 2


def test_metareg_not_estimable():
    eff = pd.DataFrame({"yi": [1.0, 2.0], "vi": [1.0, 1.0], "dose": [1, 2]})
    out = dg.metareg(eff, ["dose"], "REML", True, 0.95)
    assert "not estimable" in out.note.iloc[0]


def test_egger_matches_ols_and_k_rule():
    rng = np.random.default_rng(3)
    se = rng.uniform(0.1, 0.6, 12)
    y = 0.2 + 1.5 * se + rng.normal(0, se)
    out = dg.egger(y, se**2).iloc[0]
    ols = sm.OLS(y / se, sm.add_constant(1 / se)).fit()
    assert out.intercept == pytest.approx(ols.params[0])
    assert out.intercept_se == pytest.approx(ols.bse[0])
    assert out.p == pytest.approx(ols.pvalues[0])
    assert out.df == 10
    skipped = dg.egger(y[:9], se[:9] ** 2).iloc[0]
    assert math.isnan(skipped.p) and "9 studies < 10" in skipped.note


def test_egger_exact_line():
    # y/se = 2 + 0.5/se exactly except a small known perturbation; intercept ~ 2
    se = np.linspace(0.1, 1.0, 10)
    z = 2 + 0.5 / se + np.array([0.01, -0.01] * 5)
    out = dg.egger(z * se, se**2).iloc[0]
    assert out.intercept == pytest.approx(2, abs=0.05)
    assert out.slope == pytest.approx(0.5, abs=0.01)


def test_bh_hand():
    # sorted (0.005, 0.01, 0.03, 0.04) * 4/i = (0.02, 0.02, 0.04, 0.04)
    assert dg.bh_adjust([0.01, 0.04, 0.03, 0.005]) == pytest.approx([0.02, 0.04, 0.04, 0.02])
    # step-up monotonicity: (0.03*2/1, 0.031*2/2) = (0.06, 0.031) -> (0.031, 0.031)
    assert dg.bh_adjust([0.03, 0.031]) == pytest.approx([0.031, 0.031])
    out = dg.bh_adjust([0.2, np.nan, 0.9])
    assert out[0] == pytest.approx(0.4) and math.isnan(out[1]) and out[2] == pytest.approx(0.9)
    assert dg.bh_adjust([0.9, 0.95]) == pytest.approx([0.95, 0.95])


# --- fixes from the FX2 review round -------------------------------------------------

def _twelve(seed=1):
    rng = np.random.default_rng(seed)
    k = 12
    return pd.DataFrame({"yi": rng.normal(0, 0.6, k), "vi": rng.uniform(0.05, 0.3, k),
                         "source_study": [f"S{i}" for i in range(k)]})


@pytest.mark.parametrize("hksj", [True, False])
def test_metareg_one_level_moderator_is_a_note_not_a_crash(hksj):
    # CR2-1: every study with a non-missing value is 'adult'. Under the default config
    # (hksj, revman, tau2 > 0) this raised ZeroDivisionError; without HKSJ it wrote an
    # intercept-only row that looked like a meta-regression. run_meta.R writes this note.
    eff = _twelve().assign(pop="adult")
    eff.loc[0, "pop"] = None
    assert pl.fit(eff.yi, eff.vi, method="REML").tau2 > 0
    out = dg.metareg(eff, ["pop"], "REML", hksj, 0.95, min_k=10, hksj_rule="revman",
                     study_col="source_study")
    assert len(out) == 1 and "term" not in out
    assert out.note.iloc[0].startswith("moderator has fewer than 2 distinct values")
    # In a joint model the one-level moderator would otherwise drop out silently.
    eff["dose"] = np.arange(12.0)
    joint = dg.metareg(eff, ["dose", "pop"], "REML", hksj, 0.95, study_col="source_study")
    assert len(joint) == 1 and "(pop)" in joint.note.iloc[0]


def test_metareg_r2_is_zero_when_intercept_only_tau2_is_zero():
    # CR2-2: metafor rma.uni sets R2 = 0 when the intercept-only tau2 is 0 (and run_meta.R
    # copies f$R2); metaproof wrote NaN, a false cross-engine mismatch.
    y = np.array([0.1, 0.12, 0.09, 0.11, 0.1, 0.13, 0.08, 0.1, 0.11, 0.12])
    eff = pd.DataFrame({"yi": y, "vi": np.full(10, 0.1), "x": np.arange(10.0)})
    assert pl.fit(y, eff.vi, method="REML").tau2 == 0
    out = dg.metareg(eff, ["x"], "REML", True, 0.95)
    assert (out.R2 == 0.0).all()


def test_metareg_levels_sort_in_code_point_order():
    # CR2-3 contract: the reference level is the first in code-point (C-locale) order, which
    # run_meta.R has to match (R's default collation puts 'adults' first).
    eff = _twelve().assign(pop=["adults", "Older adults", "elderly"] * 4)
    out = dg.metareg(eff, ["pop"], "REML", True, 0.95, study_col="source_study")
    assert list(out.term) == ["intrcpt", "pop[adults]", "pop[elderly]"]


def test_metareg_gate_scales_with_number_of_moderators():
    # Astra 18: Handbook 10.11.5.1 asks for about ten studies per characteristic modelled.
    # Twelve studies passed the old gate for a joint model of two moderators.
    eff = _twelve().assign(dose=np.arange(12.0), weeks=np.tile([4.0, 8.0, 12.0], 4))
    joint = dg.metareg(eff, ["dose", "weeks"], "REML", True, 0.95, model_name="joint",
                       min_k=10, study_col="source_study")
    assert len(joint) == 1 and "term" not in joint
    assert "12 studies < metareg_min_k x 2 moderators = 10 x 2 = 20" in joint.note.iloc[0]
    single = dg.metareg(eff, ["dose"], "REML", True, 0.95, min_k=10, study_col="source_study")
    assert list(single.term) == ["intrcpt", "dose"]
    assert single.n_moderators.iloc[0] == 1 and single.k_per_param.iloc[0] == 12
    # With the gate lowered the joint model runs, and the note says it is under-supported.
    low = dg.metareg(eff, ["dose", "weeks"], "REML", True, 0.95, min_k=0,
                     study_col="source_study")
    assert low.k_per_param.iloc[0] == 6
    assert "fewer than 10 per characteristic" in low.note.iloc[0]


def test_metareg_reports_sparse_categories():
    eff = _twelve().assign(rob=["low"] * 6 + ["high"] * 5 + ["some"])
    out = dg.metareg(eff, ["rob"], "REML", True, 0.95, min_k=10, study_col="source_study")
    assert out.min_level_k.iloc[0] == 1
    assert out.k_per_param.iloc[0] == 6  # 12 studies, 2 dummies
    assert "a category of rob has 1 study" in out.note.iloc[0]


@pytest.mark.parametrize("vi", [0.25, 0.1])
def test_egger_constant_precision_not_estimable(vi):
    # Astra 12: equal SEs make the intercept and 1/se collinear. v = 0.25 raised
    # LinAlgError (stopping the pipeline); v = 0.1 returned intercept SE ~ 4e7 and p = 1.
    y = np.linspace(-1, 1, 10)
    out = dg.egger(y, np.full(10, vi)).iloc[0]
    assert math.isnan(out.p) and math.isnan(out.intercept)
    assert "constant precision" in out.note
    # Same for the modified SMD test when every trial has the same arm sizes.
    mod = dg.egger(y, np.linspace(0.1, 0.3, 10), precision_se=np.full(10, math.sqrt(2 / 20)))
    assert math.isnan(mod.iloc[0].p) and "constant precision" in mod.iloc[0].note


def test_egger_zero_residual_and_tiny_k_not_estimable():
    # Identical effects lie exactly on y/se = c/se: residual variance 0, intercept SE 0.
    se = np.linspace(0.1, 1.0, 10)
    out = dg.egger(np.full(10, 0.3), se**2).iloc[0]
    assert math.isnan(out.t) and "residual variance 0" in out.note
    two = dg.egger([0.1, 0.4], [0.04, 0.09], min_k=0).iloc[0]
    assert math.isnan(two.p) and "no residual df" in two.note
