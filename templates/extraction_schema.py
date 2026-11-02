"""Column definitions for the NR meta-analysis extraction CSV.

This module is the single source for the shared interface columns. The Excel
extraction sheet, the data dictionary CSV, the header-only CSV template and the
checker (check_extraction.py) are all generated from COLUMNS, so they cannot
drift apart.

One row = one study arm x one outcome x one timepoint.

The column order below is the interface order. Both analysis engines read a
CSV with these columns in this order. EXTENSION_COLUMNS may follow notes; both
engines read them when present and ignore their absence.
"""

from __future__ import annotations

DESIGNS = ("parallel", "crossover")
IS_CONTROL = (0, 1)
ROB_LEVELS = ("low", "some concerns", "high")

# kind: text | int | number | choice
# required: "always" | "conditional" | "optional"
# lo / hi: inclusive numeric bounds (None = unbounded)
# gt0: value must be strictly greater than zero
COLUMNS: list[dict] = [
    dict(name="study_id", kind="text", required="always", unit="",
         definition="Short stable study label, first author and year plus a letter "
                    "if needed (e.g. Smith2021a). Same value on every row of the study. "
                    "Use the same label in the RoB 2 and GRADE files."),
    dict(name="outcome", kind="text", required="always", unit="",
         definition="Outcome key from the protocol outcome list (lower_snake_case, "
                    "e.g. whole_blood_nad, sbp, fasting_glucose). One key per construct "
                    "and matrix; different assay matrices are different outcomes."),
    dict(name="timepoint", kind="text", required="always", unit="",
         definition="Timepoint label for this result (e.g. wk6, wk12, day 21). Every arm "
                    "of one comparison must use the same label. Put a number and unit in the "
                    "label; a label with no number (e.g. 'end') cannot be ordered, so fill "
                    "timepoint_weeks for it. The protocol says which timepoint enters the "
                    "primary analysis when a study reports several."),
    dict(name="design", kind="choice", choices=DESIGNS, required="always", unit="",
         definition="parallel or crossover. Same value on every row of the study."),
    dict(name="arm", kind="text", required="always", unit="",
         definition="Arm label as reported (e.g. NR 1000 mg/day, placebo). For a "
                    "crossover trial the arm is the treatment condition (NR period or "
                    "control period), not the sequence."),
    dict(name="is_control", kind="choice", choices=IS_CONTROL, required="always", unit="",
         definition="1 for the comparator arm (placebo or no treatment), 0 for an NR arm. "
                    "Exactly one control row per study x outcome x timepoint."),
    dict(name="n", kind="int", required="always", unit="participants", lo=1,
         definition="Participants analysed for this outcome in this arm at this timepoint "
                    "(not the number randomised, unless they are the same). Crossover: "
                    "participants contributing data to the paired analysis; same n on both rows."),
    dict(name="mean", kind="number", required="conditional", unit="outcome unit",
         definition="Arm mean at the timepoint (final value, post-intervention)."),
    dict(name="sd", kind="number", required="conditional", unit="outcome unit", lo=0,
         definition="SD of the final value. Leave blank if only an SE or CI is reported; "
                    "do not hand-convert SE to SD here. Many NR papers report mean +/- SEM: "
                    "check the table footnote and enter an SEM in se, never in sd. An SD of 0 is "
                    "accepted only if another arm of the comparison has an SD above 0 (the "
                    "comparison variance must be > 0); check that it was not rounded to 0."),
    dict(name="se", kind="number", required="conditional", unit="outcome unit", lo=0,
         definition="SE of the arm mean when no SD is reported. If only a 95% CI for the "
                    "arm mean is reported, enter se = (upper - lower) / (2 * t(0.975, n-1)) "
                    "and record the original CI in notes (Cochrane Handbook 6.5.2.2)."),
    dict(name="median", kind="number", required="conditional", unit="outcome unit",
         definition="Arm median at the timepoint, when the mean is not reported."),
    dict(name="q1", kind="number", required="conditional", unit="outcome unit",
         definition="First quartile (25th percentile) at the timepoint."),
    dict(name="q3", kind="number", required="conditional", unit="outcome unit",
         definition="Third quartile (75th percentile) at the timepoint."),
    dict(name="min", kind="number", required="conditional", unit="outcome unit",
         definition="Minimum at the timepoint."),
    dict(name="max", kind="number", required="conditional", unit="outcome unit",
         definition="Maximum at the timepoint."),
    dict(name="mean_change", kind="number", required="conditional", unit="outcome unit",
         definition="Mean change from baseline (timepoint minus baseline), as reported. "
                    "Sign convention: later minus earlier."),
    dict(name="sd_change", kind="number", required="conditional", unit="outcome unit",
         lo=0,
         definition="SD of the change from baseline, as reported. If only the SE of the "
                    "change is reported, leave this blank and fill se_change."),
    dict(name="baseline_mean", kind="number", required="optional", unit="outcome unit",
         definition="Arm mean at baseline. Needed to impute an SD of change or to "
                    "estimate corr_prepost from a study that reports all three SDs."),
    dict(name="baseline_sd", kind="number", required="optional", unit="outcome unit",
         lo=0,
         definition="Arm SD at baseline."),
    dict(name="corr_prepost", kind="number", required="optional", unit="", lo=-1, hi=1,
         definition="Baseline to follow-up correlation, only if the study reports it or "
                    "it can be computed from baseline, final and change SDs. Blank means "
                    "the config default (0.5) is used."),
    dict(name="corr_within", kind="number", required="optional", unit="", lo=-1, hi=1,
         definition="Crossover only: within-person correlation between the NR-period and "
                    "control-period measurements, if reported or derivable. Blank means "
                    "the config default (0.5) is used. Leave blank for parallel trials."),
    dict(name="paired_md", kind="number", required="optional", unit="outcome unit",
         definition="Reported NR-minus-control mean difference, entered on the NR row "
                    "(is_control = 0), blank on the control row. Crossover: the within-person "
                    "difference from a paired or period-adjusted analysis (preferred to "
                    "period means). Parallel: a reported between-group MD (e.g. ANCOVA-"
                    "adjusted), used as it stands for MD analyses only (generic inverse "
                    "variance); then the control row may be omitted. At most one such "
                    "parallel comparison per study: several MDs against one shared control "
                    "cannot be entered (Handbook 23.3.4)."),
    dict(name="paired_se", kind="number", required="optional", unit="outcome unit",
         lo=0, gt0=True,
         definition="SE of paired_md, on the same row. From a 95% CI in a small trial use "
                    "(upper - lower) / (2 t(0.975, df)) with the trial's df, not 3.92 "
                    "(Handbook 6.5.2.3)."),
    dict(name="dose_mg_day", kind="number", required="always", unit="mg/day", lo=0,
         definition="NR dose in mg per day for this arm (0 for the control arm), as "
                    "reported by the trial. Most trials report NR chloride; record the "
                    "salt form in notes."),
    dict(name="duration_weeks", kind="number", required="always", unit="weeks", lo=0, gt0=True,
         definition="Time from start of supplementation to this timepoint, in weeks "
                    "(days / 7). Crossover: length of each treatment period."),
    dict(name="population", kind="text", required="always", unit="",
         definition="Population category from the protocol list (e.g. healthy_older, "
                    "obesity_insulin_resistance, disease_<name>). Same value on every "
                    "row of the study."),
    dict(name="risk_of_bias", kind="choice", choices=ROB_LEVELS, required="always", unit="",
         definition="Overall RoB 2 judgement for this study result (this outcome and "
                    "timepoint), copied from rob2_summary_template.csv after consensus."),
    dict(name="notes", kind="text", required="optional", unit="",
         definition="Source of each number (table, figure, page, digitised, author "
                    "contact), unit conversions applied, CI or P values converted, and "
                    "anything a second extractor needs to reproduce the row."),
]

COLUMN_NAMES = [c["name"] for c in COLUMNS]

# A row needs at least one of these complete patterns to yield an arm-level
# estimate. Checked by check_extraction.py; listed in the data dictionary. Each pattern is
# one the engines turn into a mean and SD (final or change) or a paired MD and SE
# (metaproof.effects.prepare_arm and reference_r/run_meta.R final_stats / change_stats);
# templates/test_templates.py checks every pattern against metaproof.
DATA_PATTERNS = [
    ("final mean + SD", ("mean", "sd")),
    ("final mean + SE", ("mean", "se")),
    ("change mean + SD of change", ("mean_change", "sd_change")),
    ("median + IQR", ("median", "q1", "q3")),
    ("median + range", ("median", "min", "max")),
    ("median + IQR + range", ("median", "q1", "q3", "min", "max")),
    ("final mean + IQR, no median (SD from the IQR, Wan 2014)", ("mean", "q1", "q3")),
    ("final mean + range, no median (SD from the range, Wan 2014; excluded unless "
     "allow_range_only)", ("mean", "min", "max")),
    ("change mean + SE of change (extension column se_change)", ("mean_change", "se_change")),
    ("change mean + final SD + baseline SD (SD of change imputed with corr_prepost, "
     "Handbook 6.5.2.8)", ("mean_change", "sd", "baseline_sd")),
    ("change mean + final SE + baseline SD (SD of change imputed with corr_prepost, "
     "Handbook 6.5.2.8)", ("mean_change", "se", "baseline_sd")),
    ("final mean + baseline mean + SD of change (change mean = final - baseline)",
     ("mean", "baseline_mean", "sd_change")),
    ("final mean + baseline mean + SE of change (extension column se_change)",
     ("mean", "baseline_mean", "se_change")),
    ("crossover paired MD + SE (NR row only)", ("paired_md", "paired_se")),
    ("parallel reported between-group MD + SE (NR row only; MD analyses only; one per study)",
     ("paired_md", "paired_se")),
]

# Optional columns after notes. Both engines read them when present.
EXTENSION_COLUMNS: list[dict] = [
    dict(name="se_change", kind="number", required="optional", unit="outcome unit", lo=0,
         definition="SE of the change from baseline when no SD of change is reported. The "
                    "engines convert it (SD = SE x sqrt(n), Handbook 6.5.2.2); do not "
                    "hand-convert."),
    dict(name="timepoint_weeks", kind="number", required="optional", unit="weeks", lo=0,
         definition="The timepoint in weeks, used by the engines to order timepoints instead "
                    "of parsing the label. check_extraction.py --write-csv fills it from "
                    "duration_weeks when blank."),
]
EXTENSION_NAMES = [c["name"] for c in EXTENSION_COLUMNS]
ALL_COLUMNS = COLUMNS + EXTENSION_COLUMNS
ALL_NAMES = COLUMN_NAMES + EXTENSION_NAMES

VALIDATION_RULES = [
    "Column names and order match the interface exactly; only the extension columns "
    "(se_change, timepoint_weeks) and moderator or subgroup columns declared with "
    "check_extraction.py --extra-columns may follow notes.",
    "study_id, outcome, timepoint, design, arm, is_control, n, dose_mg_day, "
    "duration_weeks, population and risk_of_bias are filled on every row.",
    "design is parallel or crossover; is_control is 0 or 1; risk_of_bias is low, "
    "some concerns or high.",
    "n is a whole number >= 1. sd, se, sd_change, se_change and baseline_sd are >= 0; "
    "paired_se is > 0. corr_prepost and corr_within are between -1 and 1.",
    "An SD or SE of 0 gives a warning. It is an error only when every arm of a study x "
    "outcome x timepoint has SD (or SE) 0 for the same kind of value (final or change), "
    "because that comparison has zero variance and cannot be weighted.",
    "Every row has at least one complete data pattern (see the patterns list), except "
    "a crossover control row when the NR row carries paired_md and paired_se.",
    "Crossover rows need period means and SDs (or SEs) or paired_md/paired_se; a row with "
    "only within-period change is rejected because the engines cannot use it.",
    "Parallel paired_md/paired_se: on NR rows only, at most one such row per study x outcome "
    "x timepoint unless every NR row also has arm-level mean and SD (then the arms are "
    "combined); with no control row the study is usable for MD only.",
    "Fold change or percent change: when the notes mention one and mean is filled, the "
    "checker warns. It is an error only when the control mean is between 0.8 and 1.25 and "
    "the control baseline_mean is blank or exactly 1, which looks like a ratio to baseline; "
    "enter final values in the outcome unit, or fill baseline_mean if the values already "
    "are in the outcome unit.",
    "Order checks: min <= q1 <= median <= q3 <= max for whichever of these are filled.",
    "Each study x outcome x timepoint has exactly one control row and at least one NR row "
    "(or, parallel only, one NR row with paired_md/paired_se and no control row).",
    "design, population and study-level fields are constant within a study.",
    "Crossover: both rows of a comparison should have the same n (warning). corr_within is "
    "blank on parallel rows; paired_md/paired_se on a parallel row are allowed only as a "
    "reported between-group MD (see the rule above).",
    "No duplicate rows for the same study_id, outcome, timepoint and arm.",
    "Warnings: an SD that is small for its mean compared with the other studies of the "
    "outcome (coefficient of variation below half the median), which suggests an SEM "
    "entered as an SD; an SD more than twice the mean (SE/SD mix-up or skewed data); "
    "timepoint_weeks that differs from duration_weeks.",
]
