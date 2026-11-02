"""End-to-end run on the SYNTHETIC example extraction."""

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from metaproof.config import load_config
from metaproof.pipeline import run, select_timepoint
from metaproof.pooling import POOLED_COLUMNS

DATA = Path(__file__).parent / "data"
CSV = DATA / "SYNTHETIC_extraction_example.csv"
CFG = DATA / "SYNTHETIC_config_example.yaml"
FILES = ["study_effects.csv", "pooled.csv", "leave1out.csv", "influence.csv", "subgroups.csv",
         "metareg.csv", "small_study.csv", "forest.png", "funnel.png"]


def test_example_rows_are_labelled_synthetic():
    df = pd.read_csv(CSV)
    assert df["notes"].str.startswith("SYNTHETIC").all()
    assert df["study_id"].str.startswith("SYNTH").all()
    assert df["outcome"].str.startswith("SYNTHETIC").all()


def test_cli_end_to_end(tmp_path):
    out = tmp_path / "out"
    proc = subprocess.run([sys.executable, "-m", "metaproof", "run", str(CSV), "--config", str(CFG),
                           "--out", str(out)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    for name in ("SYNTHETIC_blood_marker", "SYNTHETIC_systolic_bp", "SYNTHETIC_grip_strength"):
        for f in FILES:
            assert (out / name / f).exists(), f"{name}/{f}"
        pooled = pd.read_csv(out / name / "pooled.csv")
        assert list(pooled.columns)[:len(POOLED_COLUMNS)] == POOLED_COLUMNS
        assert pooled.estimator.iloc[0] == "REML"
        eff = pd.read_csv(out / name / "study_effects.csv")
        assert list(eff.columns[:5]) == ["study_id", "yi", "vi", "sei", "method_notes"]
    summary = pd.read_csv(out / "summary_across_analyses.csv")
    assert len(summary) == 3 and summary.pooled.all()
    assert not (out / "bh_across_outcomes.csv").exists()  # BH is opt-in (protocol L)

    marker = out / "SYNTHETIC_blood_marker"
    egger = pd.read_csv(marker / "small_study.csv")
    # 12 comparisons from 11 studies (one split 3-arm); the k >= 10 rule counts studies
    assert egger.k.iloc[0] == 12 and egger.k_studies.iloc[0] == 11 and egger.p.notna().iloc[0]
    assert pd.read_csv(out / "SYNTHETIC_grip_strength" / "small_study.csv").p.isna().all()
    eff = pd.read_csv(marker / "study_effects.csv")
    assert {"SYNTH-A05|low_dose", "SYNTH-A05|high_dose"} <= set(eff.study_id)
    assert "SYNTH-A03: timepoints" in (marker / "log.txt").read_text()
    sens = pd.read_csv(out / "SYNTHETIC_systolic_bp" / "sensitivity.csv")
    assert "exclude_imputed_or_digitised" in set(sens.parameter)  # protocol K7, automatic
    sens = sens[sens.parameter == "corr_prepost"]
    assert sorted(sens.value.astype(float)) == [0.2, 0.5, 0.8]
    assert sens.estimate.nunique() == 3  # imputed change SDs move the estimate


def test_select_timepoint_rules():
    df = pd.DataFrame({"study_id": ["S", "S", "T"], "timepoint": ["4 weeks", "12 weeks", "8 wk"]})
    last, notes, _ = select_timepoint(df, "last")
    assert list(last.timepoint) == ["12 weeks", "8 wk"] and notes
    first, _, _ = select_timepoint(df, "first")
    assert list(first.timepoint) == ["4 weeks", "8 wk"]
    exact, _, excl = select_timepoint(df, "8 wk")
    assert list(exact.study_id) == ["T"]
    assert [e["study_id"] for e in excl] == ["S"]  # dropped studies are reported (F12c)
    # labels without a number no longer stop the run: the study is excluded with a reason
    kept, _, excl = select_timepoint(
        pd.DataFrame({"study_id": ["S", "S"], "timepoint": ["end", "start"]}), "last")
    assert kept.empty and "no number of weeks" in excl[0]["reason"]


def test_config_json_and_validation(tmp_path):
    p = tmp_path / "c.json"
    p.write_text(json.dumps({"defaults": {"estimator": "dl"},
                             "outcomes": {"a_12wk": {"outcome": "a", "measure": "smd"}}}))
    cfg = load_config(p)
    spec = cfg["outcomes"]["a_12wk"]
    assert (spec["estimator"], spec["measure"], spec["outcome"]) == ("DL", "SMD", "a")
    assert spec["corr_prepost_sensitivity"] == [0.2, 0.8]
    p.write_text(json.dumps({"defaults": {"estimater": "DL"}}))
    with pytest.raises(ValueError, match="unknown keys"):
        load_config(p)
    p.write_text(json.dumps({"defaults": {"measure": "OR"}}))
    with pytest.raises(ValueError, match="measure"):
        load_config(p)


def test_run_without_config_uses_all_outcomes(tmp_path):
    res = run(CSV, None, tmp_path)
    assert set(res.analysis) == {"SYNTHETIC_blood_marker", "SYNTHETIC_systolic_bp",
                                 "SYNTHETIC_grip_strength"}


# --- second review: user errors, exit codes and the forest plot -----------------------------

def _cli(*args):
    return subprocess.run([sys.executable, "-m", "metaproof", "run", *map(str, args)],
                          capture_output=True, text=True)


@pytest.mark.parametrize("case", ["missing_csv", "missing_config", "bad_value", "bad_yaml",
                                  "blank_n", "outcome_typo", "missing_column"])
def test_cli_user_errors_are_one_plain_line(tmp_path, case):
    cfg, csv = tmp_path / "c.yaml", tmp_path / "x.csv"
    df = pd.read_csv(CSV)
    df.to_csv(csv, index=False)
    cfg.write_text("defaults: {estimator: REML}\n")
    if case == "missing_csv":
        csv = tmp_path / "nope.csv"
    elif case == "missing_config":
        cfg = tmp_path / "nope.yaml"
    elif case == "bad_value":
        cfg.write_text("defaults: {estimator: BAYES}\n")
    elif case == "bad_yaml":
        cfg.write_text("defaults: {estimator: [REML\n")
    elif case == "blank_n":
        df.loc[3, "n"] = None
        df.to_csv(csv, index=False)
    elif case == "outcome_typo":
        cfg.write_text("outcomes:\n  a: {outcome: SYNTHETIC_blood_markr}\n")
    elif case == "missing_column":
        df.drop(columns=["n"]).to_csv(csv, index=False)
    proc = _cli(csv, "--config", cfg, "--out", tmp_path / "out")
    assert proc.returncode == 2, proc.stderr
    lines = [ln for ln in proc.stderr.splitlines() if ln.strip()]
    assert len(lines) == 1 and lines[0].startswith("metaproof: error: "), proc.stderr
    assert "Traceback" not in proc.stderr


def test_cli_prints_every_analysis_under_bh(tmp_path):
    cfg = tmp_path / "c.yaml"
    cfg.write_text("defaults: {bh: true, min_k_pool: 50}\n")
    proc = _cli(CSV, "--config", cfg, "--out", tmp_path / "out")
    assert proc.returncode == 0, proc.stderr
    assert "SYNTHETIC_blood_marker" in proc.stdout and "no analysable" not in proc.stderr


def _forest_rows(k, with_pi=True):
    eff = pd.DataFrame({"study_id": [f"SYNTH-{i}" for i in range(k)],
                        "yi": [0.1 * i for i in range(k)], "vi": [0.04] * k})
    pooled = {"estimate": 0.3, "ci_lb": 0.1, "ci_ub": 0.5, "tau2": 0.01, "I2": 20.0,
              "Q": 5.0, "Q_p": 0.3, "estimator": "REML", "ci_method": "HKSJ",
              "pi_lb": -0.1 if with_pi else float("nan"), "pi_ub": 0.7 if with_pi else float("nan")}
    return eff, pooled


@pytest.mark.parametrize("k", [1, 2, 4, 5, 6, 7, 10, 20])
def test_forest_footer_never_overlaps_the_axis_label(k):
    from metaproof import plots
    eff, pooled = _forest_rows(k)
    fig = plots.forest_figure(eff, pooled, "MD")
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    ax = fig.axes[0]
    foot = [t for t in fig.texts if t.get_text().startswith("k = ")][0].get_window_extent(r)
    xlab = ax.xaxis.label.get_window_extent(r)
    assert foot.y1 < xlab.y0, (k, foot, xlab)
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_forest_hides_the_prediction_interval_below_five_studies():
    from metaproof import plots
    eff, pooled = _forest_rows(4)
    fig = plots.forest_figure(eff, pooled, "MD", pi_note="prediction interval not shown: 4 "
                              "studies (< 5)")
    labels = [t.get_text() for t in fig.axes[0].get_yticklabels()]
    assert not any("prediction interval" in s for s in labels)
    assert any("prediction interval not shown" in t.get_text() for t in fig.texts)
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_pipeline_passes_the_pi_caveat_to_the_forest_plot(tmp_path, monkeypatch):
    from metaproof import plots
    seen = {}
    real = plots.forest

    def spy(path, eff, pooled, measure, level=0.95, title="", pi_note=None):
        seen[title.split(":")[0]] = pi_note
        return real(path, eff, pooled, measure, level, title, pi_note=pi_note)
    monkeypatch.setattr(plots, "forest", spy)
    run(CSV, CFG, tmp_path / "out")
    assert seen["SYNTHETIC_grip_strength"] is not None    # 4 studies
    assert seen["SYNTHETIC_blood_marker"] is None          # 11 studies
