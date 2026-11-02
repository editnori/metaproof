"""Tests for the review templates: schema, workbook, checker, and PRISMA drawing."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pandas as pd
import pytest
from openpyxl import load_workbook

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import build_extraction_sheet  # noqa: E402
import check_extraction  # noqa: E402
import prisma_flow  # noqa: E402
from extraction_schema import COLUMN_NAMES, EXTENSION_NAMES  # noqa: E402

INTERFACE = ("study_id outcome timepoint design arm is_control n mean sd se median q1 q3 "
             "min max mean_change sd_change baseline_mean baseline_sd corr_prepost corr_within "
             "paired_md paired_se dose_mg_day duration_weeks population risk_of_bias notes").split()


def test_schema_matches_interface_exactly():
    assert COLUMN_NAMES == INTERFACE


EXTENSION = ["se_change", "timepoint_weeks"]


def test_header_csv_matches_interface():
    assert EXTENSION_NAMES == EXTENSION
    with open(HERE / "extraction_template.csv", newline="") as f:
        assert next(csv.reader(f)) == INTERFACE + EXTENSION


def test_workbook_built_from_schema(tmp_path):
    paths = build_extraction_sheet.build(tmp_path)
    wb = load_workbook(paths["xlsx"])
    assert {"extraction", "data_dictionary", "validation_rules", "data_patterns"} <= set(wb.sheetnames)
    header = [c.value for c in wb["extraction"][1]]
    assert header == INTERFACE + EXTENSION
    dict_cols = [r[0].value for r in wb["data_dictionary"].iter_rows(min_row=2)]
    assert dict_cols == INTERFACE + EXTENSION
    sc = [c.value for c in wb["study_characteristics"][1]]
    assert {"washout_weeks", "n_per_sequence", "crossover_analysis", "result_selection"} <= set(sc)
    covered = set()
    for dv in wb["extraction"].data_validations.dataValidation:
        for rng in dv.sqref.ranges:
            covered.add(wb["extraction"].cell(row=1, column=rng.min_col).value)
    for col in ("design", "is_control", "risk_of_bias", "n", "sd", "corr_prepost", "corr_within"):
        assert col in covered


def _sheet_values(path):
    wb = load_workbook(path)
    return {ws.title: [[c.value for c in row] for row in ws.iter_rows()] for ws in wb}


def test_committed_workbook_is_current(tmp_path):
    fresh = build_extraction_sheet.build(tmp_path)
    for kind in ("dictionary", "header"):
        assert fresh[kind].read_text() == (HERE / fresh[kind].name).read_text()
    # The committed workbook itself, every sheet and cell (the validation_rules and
    # data_patterns sheets are what extractors read).
    assert _sheet_values(HERE / "extraction_sheet.xlsx") == _sheet_values(fresh["xlsx"])


def test_synthetic_example_passes_and_is_labelled():
    df = check_extraction.load(HERE / "extraction_example_SYNTHETIC.csv")
    errors, warnings, _ = check_extraction.check(df)
    assert errors == []
    assert df["study_id"].str.contains("SYNTHETIC").all()
    assert df["notes"].str.contains("SYNTHETIC").all()


def _example():
    return check_extraction.load(HERE / "extraction_example_SYNTHETIC.csv")


@pytest.mark.parametrize("mutate, expected", [
    (lambda d: d.drop(columns=["notes"]), "missing"),
    (lambda d: d[list(reversed(d.columns))], "different order"),
    (lambda d: d.assign(design=d["design"].where(d.index != 0, "cross-over")), "design must be"),
    (lambda d: d.assign(is_control=d["is_control"].where(d.index != 0, "2")), "is_control"),
    (lambda d: d.assign(n=d["n"].where(d.index != 0, "20.5")), "whole number"),
    (lambda d: d.assign(sd=d["sd"].where(d.index != 0, "-1")), "sd must be >= 0"),
    (lambda d: d.assign(mean=d["mean"].where(d.index != 0, "NR")), "must be a number"),
    (lambda d: d.assign(corr_within=d["corr_within"].where(d.index != 9, "1.5")), "<= 1"),
    (lambda d: d.assign(q1=d["q1"].where(d.index != 5, "60")), "min <= q1"),
    (lambda d: d.assign(corr_within=d["corr_within"].where(d.index != 0, "0.5")), "crossover trials only"),
    (lambda d: d.assign(is_control=d["is_control"].where(d.index != 0, "1")), "exactly one control"),
    (lambda d: pd.concat([d, d.iloc[[0]]], ignore_index=True), "duplicate"),
    (lambda d: d.assign(population=d["population"].where(d.index != 0, "other")), "population differs"),
    (lambda d: d.assign(mean=d["mean"].where(d.index != 0, None)), "no complete data pattern"),
    (lambda d: d.assign(paired_se=d["paired_se"].where(d.index != 7, None)), "given together"),
])
def test_checker_catches(mutate, expected):
    errors, _, _ = check_extraction.check(mutate(_example()))
    assert any(expected in e for e in errors), errors


def test_crossover_control_row_without_data_is_ok_when_paired_given():
    d = _example()
    d.loc[8, ["mean", "sd"]] = None
    errors, _, _ = check_extraction.check(d)
    assert errors == []


def test_checker_reads_xlsx_and_writes_csv(tmp_path):
    paths = build_extraction_sheet.build(tmp_path)
    wb = load_workbook(paths["xlsx"])
    ws = wb["extraction"]
    ex = _example()
    for _, row in ex.iterrows():
        ws.append([None if pd.isna(v) else v for v in row.tolist()])
    wb.save(paths["xlsx"])
    out = tmp_path / "out.csv"
    assert check_extraction.main([str(paths["xlsx"]), "--write-csv", str(out)]) == 0
    back = pd.read_csv(out)
    assert list(back.columns) == INTERFACE + EXTENSION
    assert len(back) == len(ex)
    assert (back.timepoint_weeks == back.duration_weeks).all()   # filled when blank


def test_prisma_template_and_example_render(tmp_path):
    for name in ("prisma_counts_template.csv", "prisma_counts_example_SYNTHETIC.csv"):
        counts = prisma_flow.read_counts(HERE / name)
        out = prisma_flow.draw(counts, tmp_path / (name + ".png"), title="SYNTHETIC")
        assert out.stat().st_size > 10_000


def test_prisma_example_adds_up():
    counts = prisma_flow.read_counts(HERE / "prisma_counts_example_SYNTHETIC.csv")
    assert prisma_flow.check_counts(counts) == []
    assert counts.has_other()


def test_prisma_flags_bad_arithmetic(tmp_path):
    p = tmp_path / "c.csv"
    p.write_text("key,label,n\ndatabases,,100\nregisters,,0\nduplicates_removed,,10\n"
                 "records_screened,,95\n")
    warnings = prisma_flow.check_counts(prisma_flow.read_counts(p))
    assert any("records_screened" in w for w in warnings)


def test_prisma_rejects_unknown_key(tmp_path):
    p = tmp_path / "c.csv"
    p.write_text("key,label,n\nrecords_found,,3\n")
    with pytest.raises(ValueError, match="unknown key"):
        prisma_flow.read_counts(p)


def test_prisma_without_other_methods(tmp_path):
    p = tmp_path / "c.csv"
    p.write_text("key,label,n\ndatabases,,10\nreports_excluded_reason,Wrong design,2\n")
    counts = prisma_flow.read_counts(p)
    assert not counts.has_other()
    assert prisma_flow.draw(counts, tmp_path / "x.pdf").exists()


def test_checklists_complete():
    with open(HERE / "prisma2020_checklist.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 42
    assert rows[-1]["item"] == "27"
    assert "location_reported_page" in rows[0]
    with open(HERE / "prisma2020_abstract_checklist.csv", newline="") as f:
        assert len(list(csv.DictReader(f))) == 12


def test_rob2_templates_have_domains():
    for name, domains in [("rob2_parallel_template.csv", {"D1", "D2", "D3", "D4", "D5", "Overall"}),
                          ("rob2_crossover_template.csv", {"D1", "DS", "D2", "D3", "D4", "D5", "Overall"})]:
        with open(HERE / name, newline="") as f:
            rows = list(csv.DictReader(f))
        assert domains <= {r["domain_id"] for r in rows}


# --- review findings F03, F09a, F11, F17, F18 (2026-11-01) -----------------------------------

def _rows(*specs):
    base = {c: None for c in INTERFACE}
    out = []
    for sid, arm, ctl, design, kw in specs:
        out.append({**base, "study_id": sid, "outcome": "SYNTHETIC_x", "timepoint": "wk8",
                    "design": design, "arm": arm, "is_control": ctl, "n": 20,
                    "dose_mg_day": 0 if ctl else 1000, "duration_weeks": 8,
                    "population": "SYNTHETIC_pop", "risk_of_bias": "low",
                    "notes": "SYNTHETIC test row", **kw})
    return pd.DataFrame(out, columns=INTERFACE, dtype=object)


def test_parallel_reported_md_without_control_row_is_accepted():
    d = _rows(("SYNTHETIC-P", "NR", 0, "parallel", {"paired_md": -3.0, "paired_se": 1.2}))
    errors, warnings, _ = check_extraction.check(d)
    assert errors == [] and any("MD analyses only" in w for w in warnings)


def test_several_reported_mds_against_one_control_are_rejected():
    d = _rows(("SYNTHETIC-P", "NR low", 0, "parallel", {"paired_md": -3.0, "paired_se": 1.2}),
              ("SYNTHETIC-P", "NR high", 0, "parallel", {"paired_md": -5.0, "paired_se": 1.3}),
              ("SYNTHETIC-P", "placebo", 1, "parallel", {"mean": 30.0, "sd": 8.0}))
    errors, _, _ = check_extraction.check(d)
    assert any("against one shared control" in e for e in errors), errors
    d.loc[0, ["mean", "sd"]] = [41.0, 9.0]
    d.loc[1, ["mean", "sd"]] = [55.0, 12.0]
    errors, warnings, _ = check_extraction.check(d)
    assert errors == [] and any("combine the arms" in w for w in warnings)


def test_crossover_change_only_row_is_rejected():
    d = _rows(("SYNTHETIC-X", "NR period", 0, "crossover", {"mean_change": 4.0, "sd_change": 6.0}),
              ("SYNTHETIC-X", "placebo period", 1, "crossover",
               {"mean_change": 0.5, "sd_change": 6.2}))
    errors, _, _ = check_extraction.check(d)
    assert any("only within-period change" in e for e in errors), errors


def test_fold_change_values_are_rejected():
    d = _rows(("SYNTHETIC-F", "NR", 0, "parallel", {"mean": 2.3, "sd": 0.6,
                                                    "notes": "SYNTHETIC: fold change vs baseline"}),
              ("SYNTHETIC-F", "placebo", 1, "parallel", {"mean": 1.02, "sd": 0.2}))
    errors, _, _ = check_extraction.check(d)
    assert any("fold or percent change" in e for e in errors), errors


def test_sem_entered_as_sd_is_flagged():
    specs = []
    for i, sd in enumerate([12.0, 13.0, 11.0, 12.5, 2.7]):   # last one: SEM = 12 / sqrt(20)
        sid = f"SYNTHETIC-S{i}"
        specs += [(sid, "NR", 0, "parallel", {"mean": 50.0, "sd": sd}),
                  (sid, "placebo", 1, "parallel", {"mean": 40.0, "sd": sd})]
    errors, warnings, _ = check_extraction.check(_rows(*specs))
    assert errors == []
    flagged = [w for w in warnings if "SEM" in w]
    assert len(flagged) == 1 and "SYNTHETIC-S4" in flagged[0]


def test_extension_columns_accepted_and_typed():
    d = _rows(("SYNTHETIC-C", "NR", 0, "parallel", {"mean_change": -5.0}),
              ("SYNTHETIC-C", "placebo", 1, "parallel", {"mean_change": -1.5}))
    d["se_change"] = [2.0, 1.9]
    d["timepoint_weeks"] = [8, 8]
    errors, _, clean = check_extraction.check(d)
    assert errors == [] and list(clean.columns) == INTERFACE + EXTENSION
    d["se_change"] = [2.0, -1.0]
    assert any("se_change must be >= 0" in e for e in check_extraction.check(d)[0])
    d["unknown"] = 1
    assert any("unexpected" in e for e in check_extraction.check(d)[0])


# --- review CR7 and second-review items on the checker and templates (2026-11-01) --------

ROOT = HERE.parent


def test_fold_note_on_absolute_values_near_100_is_only_a_warning():
    # CR7-1: fasting glucose ~100 mg/dL whose notes say a % change was reported but not used.
    note = "SYNTHETIC: Table 2; paper also reports % change from baseline which was not used"
    d = _rows(("SYNTHETIC-G", "NR", 0, "parallel", {"mean": 98.0, "sd": 10.0, "notes": note}),
              ("SYNTHETIC-G", "placebo", 1, "parallel", {"mean": 101.0, "sd": 11.0, "notes": note}))
    errors, warnings, _ = check_extraction.check(d)
    assert errors == [], errors
    assert any("near 100" in w for w in warnings), warnings


def test_fold_note_with_baseline_in_outcome_unit_is_only_a_warning():
    # Creatinine ~0.9 mg/dL is in the ratio band, but a baseline in the outcome unit shows the
    # values are not ratios to baseline.
    note = "SYNTHETIC: percent change also reported"
    d = _rows(("SYNTHETIC-K", "NR", 0, "parallel", {"mean": 0.95, "sd": 0.2, "notes": note}),
              ("SYNTHETIC-K", "placebo", 1, "parallel",
               {"mean": 0.9, "sd": 0.2, "baseline_mean": 0.92}))
    errors, warnings, _ = check_extraction.check(d)
    assert errors == [] and any("fold or percent" in w for w in warnings)


@pytest.mark.parametrize("nr, ctl", [
    # change mean + final SD + baseline SD: the engines impute the SD of change (6.5.2.8)
    ({"mean_change": -6.0, "sd": 10.0, "baseline_sd": 11.0},
     {"mean_change": -1.0, "sd": 9.0, "baseline_sd": 10.0}),
    ({"mean_change": -6.0, "se": 2.2, "baseline_sd": 11.0},
     {"mean_change": -1.0, "se": 2.0, "baseline_sd": 10.0}),
    # final mean + IQR with no median: Wan 2014 SD from the IQR
    ({"mean": 51.0, "q1": 40.0, "q3": 60.0}, {"mean": 45.0, "q1": 38.0, "q3": 57.0}),
    # final mean + baseline mean + SD of change: change mean = final - baseline
    ({"mean": 50.0, "baseline_mean": 55.0, "sd_change": 6.0},
     {"mean": 45.0, "baseline_mean": 46.0, "sd_change": 5.0}),
])
def test_engine_supported_patterns_pass_the_checker(nr, ctl):
    d = _rows(("SYNTHETIC-E", "NR", 0, "parallel", nr),
              ("SYNTHETIC-E", "placebo", 1, "parallel", ctl))
    errors, _, _ = check_extraction.check(d)
    assert errors == [], errors


def test_every_data_pattern_gives_the_engine_usable_arm_statistics():
    """Each checker pattern is one metaproof turns into a mean and SD (final or change) or a
    paired MD and SE, so the checker cannot pass a row the engine cannot use."""
    sys.path.insert(0, str(ROOT))
    from extraction_schema import DATA_PATTERNS
    import inspect
    from metaproof.effects import prepare_arm
    # Range-only rows are a sensitivity analysis (allow_range_only); older engine versions
    # refuse them inside prepare_arm unless told otherwise.
    kw = ({"allow_range_only": True}
          if "allow_range_only" in inspect.signature(prepare_arm).parameters else {})
    vals = {"mean": 50.0, "sd": 10.0, "se": 2.0, "median": 48.0, "q1": 40.0, "q3": 58.0,
            "min": 20.0, "max": 80.0, "mean_change": -4.0, "sd_change": 6.0, "se_change": 1.4,
            "baseline_mean": 54.0, "baseline_sd": 11.0, "paired_md": -3.0, "paired_se": 1.1}
    for label, cols in DATA_PATTERNS:
        row = {"arm": "NR", "is_control": 0, "n": 20, **{c: vals[c] for c in cols}}
        arm = prepare_arm(row, "luo", 0.5, **kw)
        usable = ((arm.mean is not None and arm.sd is not None)
                  or (arm.mean_change is not None and arm.sd_change is not None)
                  or (arm.paired_md is not None and arm.paired_se is not None))
        assert usable, label


def test_declared_moderator_columns_are_accepted(tmp_path):
    d = _example()
    d["baseline_bmi"] = "27.5"
    errors, _, _ = check_extraction.check(d)
    assert any("unexpected ['baseline_bmi']" in e and "--extra-columns" in e for e in errors)
    errors, _, clean = check_extraction.check(d, extra_columns=("baseline_bmi",))
    assert errors == [] and list(clean.columns) == INTERFACE + EXTENSION + ["baseline_bmi"]
    src, out = tmp_path / "in.csv", tmp_path / "out.csv"
    d.to_csv(src, index=False)
    assert check_extraction.main([str(src), "--write-csv", str(out),
                                  "--extra-columns", "baseline_bmi"]) == 0
    assert list(pd.read_csv(out).columns) == INTERFACE + EXTENSION + ["baseline_bmi"]
    errors, _, _ = check_extraction.check(d, extra_columns=("baseline_bmi", "sd"))
    assert any("interface columns already" in e for e in errors)


def test_zero_sd_is_judged_per_comparison():
    # A placebo arm at the assay floor (all values 0, SD 0) against an NR arm with spread:
    # the MD variance is SD_NR^2 / n_NR > 0, and both engines use it (as escalc does).
    d = _rows(("SYNTHETIC-Z", "NR", 0, "parallel", {"mean": 3.0, "sd": 1.5}),
              ("SYNTHETIC-Z", "placebo", 1, "parallel", {"mean": 0.0, "sd": 0.0}))
    errors, warnings, _ = check_extraction.check(d)
    assert errors == [] and any("sd is 0" in w for w in warnings)
    d.loc[0, "sd"] = 0.0
    errors, _, _ = check_extraction.check(d)
    assert any("zero variance" in e for e in errors), errors


def test_grade_sof_offers_three_level_imprecision_downgrade():
    # CR7-6: GRADE 34 (Zeng 2022) allows rating down three levels for imprecision.
    with open(HERE / "grade_sof_template.csv", newline="") as f:
        row = next(csv.DictReader(f))
    assert "extremely serious (-3)" in row["imprecision"]
    for dom in ("risk_of_bias", "inconsistency", "indirectness"):
        assert "(-3)" not in row[dom]


def test_prisma_blank_other_assessed_does_not_give_false_warning(tmp_path):
    # CR7-7: other-methods exclusions filled but other_reports_assessed blank: the
    # reports_included sum cannot be checked, so no warning.
    p = tmp_path / "c.csv"
    p.write_text("key,label,n\nreports_assessed,,10\nreports_excluded_reason,Wrong design,5\n"
                 "other_reports_excluded_reason,No NR arm,1\nreports_included,,6\n")
    warnings = prisma_flow.check_counts(prisma_flow.read_counts(p))
    assert not any("reports_included" in w for w in warnings), warnings
    p.write_text("key,label,n\nreports_assessed,,10\nreports_excluded_reason,Wrong design,5\n"
                 "other_reports_assessed,,3\nother_reports_excluded_reason,No NR arm,1\n"
                 "reports_included,,6\n")
    assert any("reports_included" in w for w in prisma_flow.check_counts(prisma_flow.read_counts(p)))


RULE = "\n----------------------------------------------------------------------\n"


def _protocol_section(letter):
    text = (HERE / "protocol_template.txt").read_text()
    start = text.index(f"{RULE}{letter}. ") + len(RULE)
    ends = [i for i in (text.find(RULE, start + len(RULE) + 4), text.find("\n=====", start))
            if i > 0]
    return text[start:min(ends)]


def test_protocol_adverse_event_method_exists():
    # CR7-2: Mantel-Haenszel is fixed-effect only (Handbook 10.4.1); the random-effects
    # method for dichotomous data is inverse variance (Handbook 10.4).
    o = _protocol_section("O")
    assert "Mantel-Haenszel random" not in o
    assert "inverse-variance" in o and "10.4.1" in o and "10.4.4.3" in o


def test_protocol_says_which_decisions_the_config_carries():
    # Second review: not every Part 2 decision is a config field (units, first-period
    # crossover data, eligibility, result selection and MIDs are applied by people).
    text = (HERE / "protocol_template.txt").read_text()
    assert "The analysis config (one block per outcome) implements Part 2." not in text
    assert "corr_within_sensitivity" in _protocol_section("K")


def test_multi_period_crossover_is_an_error():
    """The engines exclude a crossover with more than one NR period (Handbook 23.2.7, 23.3.4),
    so the checker says so before the run instead of after it."""
    d = _rows(("SYNTHETIC-X", "NR low period", 0, "crossover", {"mean": 44.0, "sd": 10.0}),
              ("SYNTHETIC-X", "NR high period", 0, "crossover", {"mean": 49.0, "sd": 11.0}),
              ("SYNTHETIC-X", "placebo period", 1, "crossover", {"mean": 40.0, "sd": 10.5}))
    errors, _, _ = check_extraction.check(d)
    assert any("2 NR periods in the same participants" in e for e in errors), errors
    d = _rows(("SYNTHETIC-X", "NR period", 0, "crossover", {"mean": 44.0, "sd": 10.0}),
              ("SYNTHETIC-X", "placebo A", 1, "crossover", {"mean": 40.0, "sd": 10.5}),
              ("SYNTHETIC-X", "placebo B", 1, "crossover", {"mean": 41.0, "sd": 10.0}))
    errors, _, _ = check_extraction.check(d)
    assert any("needs exactly one control row, found 2" in e for e in errors), errors


def test_crossover_unequal_n_is_an_error_unless_a_paired_md_is_reported():
    """Complete-pair N is not guessed (Handbook 23.2.7.1): unequal period n is an error, and
    only a warning when the NR row carries the trial's paired MD and SE (usable for MD)."""
    rows = [("SYNTHETIC-X", "NR period", 0, "crossover", {"mean": 46.0, "sd": 10.0, "n": 15}),
            ("SYNTHETIC-X", "placebo period", 1, "crossover", {"mean": 42.0, "sd": 9.5, "n": 13})]
    errors, _, _ = check_extraction.check(_rows(*rows))
    assert any("number of participants with both periods is not known" in e for e in errors)
    rows[0][4].update(paired_md=4.0, paired_se=1.5)
    errors, warnings, _ = check_extraction.check(_rows(*rows))
    assert errors == [] and any("complete-pair N" in w for w in warnings), (errors, warnings)
