"""Tests for the comparison rules in validation/crosscheck.py.

Every frame built here is SYNTHETIC except where a test reads validation/reference_data/
(metafor expected values for the published methods datasets dat.bcg and dat.curtis1998).
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
REF = HERE / "reference_data"
_spec = importlib.util.spec_from_file_location("crosscheck", HERE / "crosscheck.py")
cc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cc)

from metaproof import pooling as pl  # noqa: E402  (crosscheck puts the kit root on sys.path)
from metaproof.config import load_config  # noqa: E402


def run(left, right, cols, key="study_id", **kw):
    led = cc.Ledger()
    cc.compare(led, "t", left, right, key, key, cols, **kw)
    return led


def status(led):
    return led.lines[-1].split()[0]


def df(**cols):
    return pd.DataFrame(cols)


# ---------------------------------------------------------------------------------------
# column resolution, duplicate keys, malformed numbers, empty frames
# ---------------------------------------------------------------------------------------

def test_column_missing_on_one_side_is_not_compared_with_the_other_side():
    # Python lacks vi; the old merge fallback compared R's vi with itself and passed
    led = run(df(study_id=["s1"], yi=[1.0]), df(study_id=["s1"], yi=[1.0], vi=[123.0]),
              [("yi", "yi"), ("vi", "vi")])
    assert led.n_fail == 1 and "missing in Python" in "\n".join(led.lines)
    led = run(df(study_id=["s1"], yi=[1.0], vi=[123.0]), df(study_id=["s1"], yi=[1.0]),
              [("yi", "yi"), ("vi", "vi")])
    assert led.n_fail == 1 and "missing in R" in "\n".join(led.lines)


@pytest.mark.parametrize("dup_side", ["python", "R"])
def test_duplicate_study_rows_fail(dup_side):
    two, one = df(study_id=["s1", "s1"], yi=[1.0, 1.0]), df(study_id=["s1"], yi=[1.0])
    left, right = (two, one) if dup_side == "python" else (one, two)
    led = run(left, right, [("yi", "yi")])
    assert led.n_fail == 1 and "appears 2 times" in "\n".join(led.lines)


def test_malformed_numbers_do_not_agree_as_missing():
    led = run(df(study_id=["s1"], yi=["1.0x"]), df(study_id=["s1"], yi=["oops"]),
              [("yi", "yi")])
    assert led.n_fail == 1 and "not a number" in "\n".join(led.lines)
    # genuine missing values on both sides still agree
    led = run(df(study_id=["s1"], yi=[np.nan]), df(study_id=["s1"], yi=["NA"]), [("yi", "yi")])
    assert led.n_pass == 1


def test_empty_comparison_is_skipped_not_passed():
    led = run(df(study_id=[], yi=[]), df(study_id=[], yi=[]), [("yi", "yi")])
    assert (led.n_pass, led.n_fail, led.n_skip) == (0, 0, 1)
    assert status(led) == "SKIP"


def test_unclassified_column_is_a_failure():
    led = run(df(study_id=["s1"], mystery=[1.0]), df(study_id=["s1"], mystery=[1.0]),
              [("mystery", "mystery")])
    assert led.n_fail == 1 and "no tolerance class" in "\n".join(led.lines)


# ---------------------------------------------------------------------------------------
# tolerances: small values, p-values, unit invariance, real noise
# ---------------------------------------------------------------------------------------

def test_hundredfold_error_in_a_small_variance_fails():
    # Astra 7: Python vi 1e-6 against reference 1e-8 passed under max(1, |b|) scaling
    led = run(df(study_id=["s1"], vi=[1e-6]), df(study_id=["s1"], vi=[1e-8]), [("vi", "vi")])
    assert led.n_fail == 1


def test_doubled_small_pvalue_fails_against_metafor_values():
    ep = pd.read_csv(REF / "bcg_expected_pooled.csv")
    ep = ep[ep.test == "z"]
    got = ep.assign(pval=2 * ep.pval)
    for est in ep.estimator:
        led = run(got[got.estimator == est], ep[ep.estimator == est], [("pval", "pval")],
                  key="estimator")
        assert led.n_fail == 1, est           # REML/z p 7.05e-5 -> 1.41e-4


def test_scaled_leave1out_pvalue_and_small_hat_error_fail():
    el = pd.read_csv(REF / "curtis1998_expected_leave1out.csv")
    led = run(el.assign(pval=el.pval * 1000), el, [("pval", "pval")])
    assert led.n_fail == 1
    ei = pd.read_csv(REF / "curtis1998_expected_influence.csv")
    led = run(ei.assign(hat=ei.hat * 1.005), ei, [("hat", "hat")])
    assert led.n_fail == 1                    # +0.5% of a hat value below 0.015


def test_identical_values_pass_and_tiny_noise_passes():
    ep = pd.read_csv(REF / "bcg_expected_pooled.csv")
    ep = ep[ep.test == "knha"]
    cols = cc.POOLED + [("pval", "pval")]
    noisy = ep.copy()
    for _, c in cols:
        if c not in cc.EXACT_COLS:
            noisy[c] = noisy[c] * (1 + 1e-10)
    for est in ep.estimator:
        led = run(noisy[noisy.estimator == est], ep[ep.estimator == est], cols,
                  key="estimator", var_scale=0.001)
        assert led.n_pass == 1, "\n".join(led.lines)


@pytest.mark.parametrize("col, ref, got, extra", [
    ("I2", 1.87732452810868e-07, 0.0, {}),             # REML stopping noise near tau2 = 0
    ("Q", 0.0, 4.440892098500626e-16, {}),             # identical effects
    ("pval", 0.0, 6.16297582204e-33, {}),              # identical effects, SE 1.6e-16
    ("se", 0.0, 1.57009245868e-16, {"vs": 0.1}),        # HKSJ s2 = 0
    ("tau2_del", 1.53635065036e-11, 0.0, {"vs": 0.004}),
])
def test_convergence_noise_seen_in_clean_runs_passes(col, ref, got, extra):
    led = run(df(study_id=["s"], **{col: [got]}), df(study_id=["s"], **{col: [ref]}),
              [(col, col)], var_scale=extra.get("vs"))
    assert led.n_pass == 1, "\n".join(led.lines)


@pytest.mark.parametrize("scale", [1e-8, 1e-3, 1.0, 1e4])
def test_verdicts_do_not_change_with_measurement_units(scale):
    # one study-effects row and one pooled row in arbitrary units; location in units of
    # `scale`, variances in units of scale^2
    se_ref = df(study_id=["s"], yi=[0.3 * scale], vi=[0.04 * scale ** 2], sei=[0.2 * scale])
    po_ref = df(estimator=["REML"], estimate=[0.3 * scale], se=[0.05 * scale],
                ci_lb=[-0.002 * scale], ci_ub=[0.6 * scale], tau2=[0.01 * scale ** 2])
    cols_se = [("yi", "yi"), ("vi", "vi"), ("sei", "sei")]
    cols_po = [(c, c) for c in ("estimate", "se", "ci_lb", "ci_ub", "tau2")]
    vs = 0.04 * scale ** 2
    # agreement to 1e-9 relative passes at every scale
    ok = se_ref.copy()
    ok[["yi", "vi", "sei"]] *= 1 + 1e-9
    assert run(ok, se_ref, cols_se).n_pass == 1
    okp = po_ref.copy()
    okp[["estimate", "se", "ci_lb", "ci_ub", "tau2"]] *= 1 + 1e-9
    assert run(okp, po_ref, cols_po, key="estimator", var_scale=vs).n_pass == 1
    # the same relative errors fail at every scale
    for col, frame, cols, key in (("vi", se_ref, cols_se, "study_id"),
                                  ("yi", se_ref, cols_se, "study_id"),
                                  ("ci_lb", po_ref, cols_po, "estimator"),
                                  ("tau2", po_ref, cols_po, "estimator")):
        bad = frame.copy()
        bad[col] = bad[col] * 1.5 if col != "ci_lb" else bad[col] + 1e-3 * scale
        led = run(bad, frame, cols, key=key, var_scale=vs)
        assert led.n_fail == 1, (col, scale, led.lines)


# ---------------------------------------------------------------------------------------
# engine-to-engine folders: gates, sensitivity, pooled test, metareg terms
# ---------------------------------------------------------------------------------------

def _note(path, text):
    pd.DataFrame({"note": [text]}).to_csv(path, index=False)


def make_pair(tmp, k_studies=2, leave1out=False, egger=False, sens=True, r_pval=None,
              flags=None):
    """Write a minimal matching Python/R analysis folder pair (SYNTHETIC values)."""
    rng = np.random.default_rng(1)
    y = np.round(rng.normal(1, 0.5, k_studies), 3)
    v = np.round(rng.uniform(0.05, 0.2, k_studies), 4)
    ids = [f"SYNTHETIC-{i}" for i in range(k_studies)]
    spec = load_config(None)["defaults"]
    P, R = tmp / "py" / "a", tmp / "r" / "a"
    P.mkdir(parents=True)
    R.mkdir(parents=True)
    eff = pd.DataFrame({"study_id": ids, "yi": y, "vi": v, "sei": np.sqrt(v),
                        "source_study": ids, "design": "parallel"})
    for c in cc.K7_FLAGS + cc.K5B_FLAGS:
        eff[c] = False
    if flags is not None:
        eff.loc[flags, "sd_change_imputed"] = True
    eff.to_csv(P / "study_effects.csv", index=False)
    eff[["study_id", "yi", "vi", "sei", "source_study"]].to_csv(R / "study_effects.csv",
                                                                 index=False)
    for d in (P, R):
        pd.DataFrame(columns=["study_id", "reason"]).to_csv(d / "excluded.csv", index=False)
        (d / "forest.png").write_bytes(b"")
        (d / "funnel.png").write_bytes(b"")
    row, test, _ = pl.pooled_row(y, v, "REML", True, 0.95, True, "riley", spec["hksj_rule"])
    pooled = pd.DataFrame([{**row, "k_studies": k_studies, "note": ""}])
    pooled.to_csv(P / "pooled.csv", index=False)
    r_pooled = pooled.copy()
    if r_pval is not None:
        r_pooled["stat"], r_pooled["pval"] = test["stat"], r_pval
    r_pooled.to_csv(R / "pooled.csv", index=False)
    pd.DataFrame([test]).to_csv(P / "pooled_tests.csv", index=False)
    if leave1out:
        lo = pd.DataFrame([{"omitted": ids[0], **{c: row.get(c, 0.0) for c, _ in cc.LOO}}])
        lo.to_csv(P / "leave1out.csv", index=False)
        lo.rename(columns={"omitted": "study_id"}).to_csv(R / "leave1out.csv", index=False)
    else:
        _note(P / "leave1out.csv", "not run")
        pd.DataFrame(columns=["study_id"] + [c for c, _ in cc.LOO]).to_csv(
            R / "leave1out.csv", index=False)
    _note(P / "influence.csv", "not run")
    pd.DataFrame(columns=["study_id", "rstudent"]).to_csv(R / "influence.csv", index=False)
    _note(P / "subgroups.csv", "no subgroups configured")
    pd.DataFrame(columns=["moderator", "level", "estimate", "note"]).to_csv(
        R / "subgroups.csv", index=False)
    _note(P / "metareg.csv", "no moderators configured")
    pd.DataFrame(columns=["model", "term", "estimate", "note"]).to_csv(R / "metareg.csv",
                                                                        index=False)
    eg_p = {"test": "Egger regression", "k": k_studies, "k_studies": k_studies,
            "intercept": 0.5 if egger else np.nan, "intercept_se": 0.2 if egger else np.nan,
            "t": 2.5 if egger else np.nan, "df": k_studies - 2 if egger else np.nan,
            "p": 0.03 if egger else np.nan, "slope": 0.9 if egger else np.nan}
    pd.DataFrame([eg_p]).to_csv(P / "small_study.csv", index=False)
    pd.DataFrame([{"test": "egger", "k": k_studies, "k_studies": k_studies,
                   "coef_sei": eg_p["intercept"], "coef_se": eg_p["intercept_se"],
                   "statistic": eg_p["t"], "df": eg_p["df"], "pval": eg_p["p"],
                   "limit_estimate": eg_p["slope"]}]).to_csv(R / "small_study.csv", index=False)
    sens_cols = [c for c, _ in cc.SENS_COLS]
    if sens:
        sp = pd.DataFrame([{"parameter": "corr_prepost", "value": x, "estimator": "REML",
                            **{c: row[c] for c in sens_cols}} for x in (0.5, 0.2, 0.8)])
        sr = pd.DataFrame([{"analysis": a, "corr_prepost": x, "corr_within": 0.5,
                            **{c: row[c] for c in sens_cols}}
                           for a, x in (("primary", 0.5), ("corr_prepost=0.2", 0.2),
                                        ("corr_prepost=0.8", 0.8))])
        if flags is not None:
            keep = ~eff.sd_change_imputed.to_numpy()
            r2, _, _ = pl.pooled_row(y[keep], v[keep], "REML", True, 0.95, True, "riley",
                                     spec["hksj_rule"])
            sp = pd.concat([sp, pd.DataFrame([{"parameter": "exclude_imputed_or_digitised",
                                               "value": True, "estimator": "REML",
                                               **{c: r2[c] for c in sens_cols}}])])
        sp.to_csv(P / "sensitivity.csv", index=False)
        sr.to_csv(R / "sensitivity.csv", index=False)
    else:
        _note(P / "sensitivity.csv", "none")
        pd.DataFrame(columns=["analysis"]).to_csv(R / "sensitivity.csv", index=False)
    return P.parent, R.parent, spec


def engines(tmp, spec_override=None, **kw):
    py, r, spec = make_pair(tmp, **kw)
    spec = {**spec, **(spec_override or {})}
    led = cc.Ledger()
    cc.compare_engines(led, "T", py, r, "a", spec)
    return led, py, r, spec


def lines_for(led, word):
    return [ln for ln in led.lines if word in ln]


def test_matching_fixture_has_no_failures(tmp_path):
    led, *_ = engines(tmp_path)
    assert led.n_fail == 0, "\n".join(led.lines)


def test_missing_sensitivity_file_on_either_side_fails(tmp_path):
    py, r, spec = make_pair(tmp_path)
    (r / "a" / "sensitivity.csv").unlink()
    led = cc.Ledger()
    cc.compare_engines(led, "T", py, r, "a", spec)
    assert led.n_fail >= 1 and lines_for(led, "sensitivity.csv")


def test_correlation_rows_in_one_engine_only_fail(tmp_path):
    py, r, spec = make_pair(tmp_path)
    pd.DataFrame(columns=["analysis"]).to_csv(r / "a" / "sensitivity.csv", index=False)
    led = cc.Ledger()
    cc.compare_engines(led, "T", py, r, "a", spec)
    assert any(ln.startswith("FAIL") for ln in lines_for(led, "imputed correlations"))


def test_leave1out_skipped_by_both_is_checked_against_k(tmp_path):
    led, *_ = engines(tmp_path / "k12", k_studies=12, egger=True)
    assert any(ln.startswith("FAIL") for ln in lines_for(led, "leave1out"))
    led, *_ = engines(tmp_path / "k2", k_studies=2, spec_override={"min_k_pool": 3})
    loo = lines_for(led, "leave1out")
    assert loo and all(ln.startswith("SKIP") for ln in loo)


def test_egger_skipped_by_both_is_checked_against_egger_min_k(tmp_path):
    led, *_ = engines(tmp_path / "k12", k_studies=12, leave1out=True)
    assert any(ln.startswith("FAIL") for ln in lines_for(led, "Egger"))
    led, *_ = engines(tmp_path / "k4", k_studies=4, leave1out=True)
    assert all(ln.startswith("SKIP") for ln in lines_for(led, "Egger"))


def test_configured_rule_pvalue_is_skip_without_r_column_and_compared_with_it(tmp_path):
    led, *_ = engines(tmp_path / "none")
    ln = lines_for(led, "pooled test statistic")
    assert ln and ln[0].startswith("SKIP") and led.n_skip >= 1
    py, r, spec = make_pair(tmp_path / "ok", r_pval=0.0)
    p_ok = pd.read_csv(py / "a" / "pooled_tests.csv").pval.iloc[0]
    for pv, want in ((p_ok, "PASS"), (p_ok * 2, "FAIL")):
        rp = pd.read_csv(r / "a" / "pooled.csv").assign(pval=pv)
        rp.to_csv(r / "a" / "pooled.csv", index=False)
        led = cc.Ledger()
        cc.compare_engines(led, "T", py, r, "a", spec)
        assert lines_for(led, "pooled test statistic")[0].startswith(want)


def test_skips_are_not_counted_as_passes(tmp_path):
    led, *_ = engines(tmp_path)
    n_skip_lines = sum(ln.startswith("SKIP") for ln in led.lines)
    n_pass_lines = sum(ln.startswith("PASS") for ln in led.lines)
    assert n_skip_lines == led.n_skip > 0 and n_pass_lines == led.n_pass


def test_k7_sensitivity_row_is_re_pooled(tmp_path):
    py, r, spec = make_pair(tmp_path, k_studies=6, flags=[0, 1])
    led = cc.Ledger()
    cc.compare_engines(led, "T", py, r, "a", spec)
    assert lines_for(led, "exclude_imputed_or_digitised")[0].startswith("PASS")
    sp = pd.read_csv(py / "a" / "sensitivity.csv")
    sp.loc[sp.parameter == "exclude_imputed_or_digitised", "estimate"] *= 1.01
    sp.to_csv(py / "a" / "sensitivity.csv", index=False)
    led = cc.Ledger()
    cc.compare_engines(led, "T", py, r, "a", spec)
    assert lines_for(led, "exclude_imputed_or_digitised")[0].startswith("FAIL")


def test_metareg_rows_are_matched_by_term_not_position(tmp_path):
    py, r, spec = make_pair(tmp_path)
    mr = pd.DataFrame({"model": ["dose", "dose"], "term": ["intrcpt", "dose"], "k": [2, 2],
                       "estimate": [1.0, 0.002], "se": [0.3, 0.001], "stat": [3.3, 2.0],
                       "df": [np.nan, np.nan], "pval": [0.001, 0.05], "ci_lb": [0.4, 0.0],
                       "ci_ub": [1.6, 0.004], "tau2": [0.01, 0.01], "R2": [10.0, 10.0],
                       "QM": [4.0, 4.0], "QM_df1": [1, 1], "QM_df2": [np.nan, np.nan],
                       "QM_p": [0.05, 0.05], "QE": [3.0, 3.0], "QE_p": [0.2, 0.2]})
    mr.to_csv(py / "a" / "metareg.csv", index=False)
    # run_meta.R lists the same coefficients in the other order
    mr.iloc[::-1].rename(columns={"stat": "statistic"}).to_csv(r / "a" / "metareg.csv",
                                                               index=False)
    spec = {**spec, "moderators": ["dose"]}
    led = cc.Ledger()
    cc.compare_engines(led, "T", py, r, "a", spec)
    assert lines_for(led, "metareg")[0].startswith("PASS"), "\n".join(led.lines)
    mr.assign(term=["intrcpt", "dose_mg"]).to_csv(py / "a" / "metareg.csv", index=False)
    led = cc.Ledger()
    cc.compare_engines(led, "T", py, r, "a", spec)
    assert lines_for(led, "metareg")[0].startswith("FAIL")


def test_crosscheck_report_summary_counts_skips_separately():
    text = (HERE / "crosscheck_report.txt").read_text()
    summary = text.split("Summary")[-1]
    assert " skipped, " in summary
    assert "max(1," not in text
    n_skip = sum(ln.startswith("SKIP") for ln in text.splitlines())
    assert f"{n_skip} skipped" in summary and n_skip > 0


def test_lock_status_reports_versions_that_differ_from_the_locks(tmp_path, monkeypatch):
    """Astra 22: the report says whether the run used the verified dependency versions."""
    assert cc.lock_status().startswith("environment matches")
    root = tmp_path / "kit"
    (root / "reference_r" / "lib").mkdir(parents=True)
    lock = (cc.ROOT / "requirements-lock.txt").read_text()
    (root / "requirements-lock.txt").write_text(
        lock.replace(f"numpy=={np.__version__}", "numpy==0.0.1"))
    rl = (cc.ROOT / "reference_r" / "R_PACKAGES_LOCK.txt").read_text()
    (root / "reference_r" / "R_PACKAGES_LOCK.txt").write_text(rl)
    monkeypatch.setattr(cc, "ROOT", root)
    msg = cc.lock_status()
    assert msg.startswith("environment differs")
    assert f"numpy {np.__version__} (lock 0.0.1)" in msg
    assert "R metafor missing (lock 5.2-1)" in msg
