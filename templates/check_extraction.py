"""Check an extraction CSV or workbook against the shared interface rules.

Usage:
    .venv/bin/python templates/check_extraction.py extraction.csv
    .venv/bin/python templates/check_extraction.py extraction_sheet.xlsx --write-csv extraction.csv
    .venv/bin/python templates/check_extraction.py extraction.csv --extra-columns baseline_bmi,sponsor

Errors (exit code 1) break the interface or make a row unusable. Warnings are
printed for things a person should look at but that the engines can handle.
With --write-csv the extraction sheet (or CSV) is written back out in the
interface column order, followed by the extension columns (se_change,
timepoint_weeks), only when there are no errors. A blank timepoint_weeks is
filled from duration_weeks.

Columns that the analysis config names as subgroups or moderators but that are not in the
interface (e.g. baseline_bmi) are declared with --extra-columns. They may follow notes, are
kept as text, and are written after the extension columns. Any other column is an error.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import pandas as pd

from extraction_schema import (ALL_COLUMNS, ALL_NAMES, COLUMN_NAMES, DATA_PATTERNS,
                               EXTENSION_NAMES)

STUDY_CONSTANT = ("design", "population")
CROSSOVER_ONLY = ("corr_within",)
# Words in notes that mean the numbers are a ratio to baseline, not values in the outcome unit.
FOLD_WORDS = re.compile(r"fold|% ?change|percent(age)? change|relative to baseline|"
                        r"ratio to baseline|normali[sz]ed to baseline", re.I)
CV_RATIO = 0.5   # SEM-as-SD warning when CV < CV_RATIO x median CV of the other studies
# Spread columns that may be 0 in one arm, by kind of value: (SD column, SE column).
SPREADS = {"final": ("sd", "se"), "change": ("sd_change", "se_change")}


def load(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in (".xlsx", ".xlsm"):
        df = pd.read_excel(path, sheet_name="extraction", dtype=object)
    else:
        df = pd.read_csv(path, dtype=object, keep_default_na=False, na_values=[""])
    df = df.dropna(how="all")
    return df


def _blank(v) -> bool:
    return v is None or (isinstance(v, float) and math.isnan(v)) or str(v).strip() == ""


def _check_columns(cols: list[str], errors: list[str], declared: tuple = ()) -> bool:
    """Interface columns in order, then optionally the extension columns and the declared
    extra columns. Returns False when interface columns are missing (nothing else can be
    checked)."""
    head, tail = cols[:len(COLUMN_NAMES)], cols[len(COLUMN_NAMES):]
    missing = [c for c in COLUMN_NAMES if c not in cols]
    extra = [c for c in cols if c not in ALL_NAMES and c not in declared]
    if head == COLUMN_NAMES and not extra and len(set(tail)) == len(tail):
        return True
    msg = "columns do not match the interface"
    if missing:
        msg += f"; missing {missing}"
    if extra:
        msg += (f"; unexpected {extra} (only {EXTENSION_NAMES} may follow notes; declare a "
                "subgroup or moderator column named in the analysis config with "
                "--extra-columns)")
    if not missing and not extra:
        msg += "; same names but different order"
    errors.append(msg)
    return not missing


def check(df: pd.DataFrame, extra_columns=()) -> tuple[list[str], list[str], pd.DataFrame]:
    errors: list[str] = []
    warnings: list[str] = []
    declared = tuple(c for c in dict.fromkeys(extra_columns) if c)
    clash = [c for c in declared if c in ALL_NAMES]
    if clash:
        errors.append(f"--extra-columns {clash} are interface columns already")
        declared = tuple(c for c in declared if c not in ALL_NAMES)

    if not _check_columns(list(df.columns), errors, declared):
        return errors, warnings, df

    clean = df.copy()
    for name in EXTENSION_NAMES + [c for c in declared if c not in clean]:
        clean[name] = clean[name] if name in clean else None
    clean = clean[ALL_NAMES + list(declared)].copy()
    for name in declared:
        clean[name] = [None if _blank(v) else str(v).strip() for v in clean[name]]
    for c in ALL_COLUMNS:
        name = c["name"]
        col = clean[name]
        for i, v in col.items():
            row = f"row {i + 2}"  # spreadsheet row number (header is row 1)
            if _blank(v):
                if c["required"] == "always":
                    errors.append(f"{row}: {name} is required")
                clean.at[i, name] = None
                continue
            if c["kind"] == "text":
                clean.at[i, name] = str(v).strip()
                continue
            if c["kind"] == "choice":
                s = str(v).strip()
                if name == "is_control":
                    try:
                        x = float(s)
                    except ValueError:
                        x = None
                    if x not in (0.0, 1.0):
                        errors.append(f"{row}: is_control must be 0 or 1, got {v!r}")
                        clean.at[i, name] = None
                        continue
                    clean.at[i, name] = int(x)
                else:
                    s = s.lower()
                    if s not in c["choices"]:
                        errors.append(f"{row}: {name} must be one of {list(c['choices'])}, got {v!r}")
                        clean.at[i, name] = None
                        continue
                    clean.at[i, name] = s
                continue
            try:
                x = float(str(v).strip())
            except ValueError:
                errors.append(f"{row}: {name} must be a number, got {v!r}")
                clean.at[i, name] = None
                continue
            if not math.isfinite(x):
                errors.append(f"{row}: {name} must be finite, got {v!r}")
                clean.at[i, name] = None
                continue
            if c["kind"] == "int" and x != int(x):
                errors.append(f"{row}: {name} must be a whole number, got {v!r}")
            if c.get("gt0") and x <= 0:
                errors.append(f"{row}: {name} must be > 0, got {v!r}")
            elif c.get("lo") is not None and x < c["lo"]:
                errors.append(f"{row}: {name} must be >= {c['lo']}, got {v!r}")
            if c.get("hi") is not None and x > c["hi"]:
                errors.append(f"{row}: {name} must be <= {c['hi']}, got {v!r}")
            clean.at[i, name] = x

    def filled(r, names):
        return all(not _blank(r[n]) for n in names)

    paired_keys = set()
    for i, r in clean.iterrows():
        if r["design"] == "crossover" and r["is_control"] == 0 and filled(r, ("paired_md", "paired_se")):
            paired_keys.add((r["study_id"], r["outcome"], r["timepoint"]))

    for i, r in clean.iterrows():
        row = f"row {i + 2}"
        key = (r["study_id"], r["outcome"], r["timepoint"])
        if not any(filled(r, cols) for _, cols in DATA_PATTERNS):
            if not (r["design"] == "crossover" and r["is_control"] == 1 and key in paired_keys):
                errors.append(f"{row}: no complete data pattern (see data_patterns)")
        if r["design"] == "crossover":
            period = filled(r, ("mean", "sd")) or filled(r, ("mean", "se"))
            paired = filled(r, ("paired_md", "paired_se")) or key in paired_keys
            if not period and not paired and not _blank(r["mean_change"]):
                errors.append(f"{row}: crossover row has only within-period change; the "
                              "engines need the period means and SDs or paired_md/paired_se "
                              "(Handbook 23.2.7). Enter those, or the first period only as a "
                              "parallel trial if the protocol planned it (Handbook 23.2.4)")
        order = [n for n in ("min", "q1", "median", "q3", "max") if not _blank(r[n])]
        vals = [r[n] for n in order]
        if any(a > b for a, b in zip(vals, vals[1:])):
            errors.append(f"{row}: min <= q1 <= median <= q3 <= max is violated")
        if r["design"] == "parallel":
            for n in CROSSOVER_ONLY:
                if not _blank(r[n]):
                    errors.append(f"{row}: {n} is for crossover trials only")
        if r["is_control"] == 1 and not _blank(r["paired_md"]):
            errors.append(f"{row}: paired_md belongs on the NR row, not the control row")
        if filled(r, ("paired_md",)) != filled(r, ("paired_se",)):
            errors.append(f"{row}: paired_md and paired_se must be given together")
        if r["is_control"] == 1 and not _blank(r["dose_mg_day"]) and r["dose_mg_day"] != 0:
            warnings.append(f"{row}: control arm has dose_mg_day {r['dose_mg_day']}")
        if r["is_control"] == 0 and not _blank(r["dose_mg_day"]) and r["dose_mg_day"] == 0:
            warnings.append(f"{row}: NR arm has dose_mg_day 0")
        for n in ("sd", "se", "sd_change", "se_change", "baseline_sd"):
            if not _blank(r[n]) and r[n] == 0:
                warnings.append(f"{row}: {n} is 0; usable only if another arm of the "
                                "comparison has a value above 0; check it was not rounded to 0")
        if filled(r, ("mean", "sd")) and r["mean"] > 0 and r["sd"] > 2 * r["mean"]:
            warnings.append(f"{row}: sd is more than twice the mean; check for SE/SD mix-up "
                            "or skewed data")
        if filled(r, ("timepoint_weeks", "duration_weeks")) \
                and abs(r["timepoint_weeks"] - r["duration_weeks"]) > 1e-9:
            warnings.append(f"{row}: timepoint_weeks {r['timepoint_weeks']:g} differs from "
                            f"duration_weeks {r['duration_weeks']:g}")

    dup = clean.duplicated(subset=["study_id", "outcome", "timepoint", "arm"], keep=False)
    for i in clean.index[dup]:
        errors.append(f"row {i + 2}: duplicate study_id/outcome/timepoint/arm")

    for key, g in clean.groupby(["study_id", "outcome", "timepoint"], dropna=False):
        label = "/".join(str(k) for k in key)
        ctrl, nr = g[g["is_control"] == 1], g[g["is_control"] == 0]
        parallel = (g["design"] == "parallel").all()
        nr_paired = nr[nr["paired_md"].notna() & nr["paired_se"].notna()]
        paired_only = parallel and len(ctrl) == 0 and len(nr) == 1 and len(nr_paired) == 1
        if len(ctrl) != 1 and not paired_only:
            errors.append(f"{label}: needs exactly one control row, found {len(ctrl)}"
                          + (" (with no control row, only one NR row with paired_md/"
                             "paired_se is allowed)" if len(ctrl) == 0 and parallel else ""))
        if len(nr) < 1:
            errors.append(f"{label}: needs at least one NR row")
        if paired_only:
            warnings.append(f"{label}: reported between-group MD with no control row; usable "
                            "for MD analyses only (not SMD or ROM)")
        if parallel and len(nr_paired) > 1:
            arm_level = nr["mean"].notna() & (nr["sd"].notna() | nr["se"].notna())
            if len(ctrl) == 1 and arm_level.all():
                warnings.append(f"{label}: {len(nr_paired)} reported between-group MDs share "
                                "one control. Under multi_arm: combine the engines ignore them "
                                "and combine the arms from arm-level data; under split the "
                                "study is refused (Handbook 23.3.4)")
            else:
                errors.append(f"{label}: {len(nr_paired)} reported between-group MDs "
                              "(paired_md/paired_se) against one shared control. Entering "
                              "each would count the control more than once (Handbook "
                              "23.3.4): keep one prespecified comparison, or add arm-level "
                              "mean and SD to every arm so the arms can be combined")
        crossover = (g["design"] == "crossover").all()
        if crossover and len(nr) > 1:
            # The engines refuse it: the periods share participants and the covariance of
            # their contrasts is not in summary data (Handbook 23.2.7, 23.3.4).
            errors.append(f"{label}: crossover with {len(nr)} NR periods in the same "
                          "participants; the engines exclude it (the contrasts share people and "
                          "their covariance is not reported). Keep one prespecified NR-vs-placebo "
                          "pair of rows, or one paired estimate the trial reports")
        if crossover and g["n"].nunique(dropna=True) > 1:
            ns = sorted(g["n"].dropna().unique())
            if len(nr_paired):
                warnings.append(f"{label}: crossover rows have different n {ns}; the reported "
                                "paired_md/paired_se is used for MD, but SMD and ROM need the "
                                "complete-pair N and exclude the study")
            else:
                errors.append(f"{label}: crossover rows have different n {ns}; the number of "
                              "participants with both periods is not known, so the engines "
                              "exclude the study. Enter the complete-pair N on both rows with the "
                              "period summaries of those participants, or the trial's paired MD "
                              "and SE on the NR row (Handbook 23.2.7.1)")
        for kind, (sd_col, se_col) in SPREADS.items():
            spread = [r[sd_col] if not _blank(r[sd_col]) else r[se_col] for _, r in g.iterrows()]
            if len(spread) > 1 and len(nr_paired) == 0 \
                    and all(not _blank(v) and v == 0 for v in spread):
                errors.append(f"{label}: every arm has {sd_col} (or {se_col}) 0, so the "
                              f"comparison of {kind} values has zero variance and cannot be "
                              "weighted; check the source, or leave these blank")
        fold = _fold_check(label, g, ctrl)
        if fold:
            (errors if fold[0] else warnings).append(fold[1])

    for study, g in clean.groupby("study_id", dropna=False):
        for n in STUDY_CONSTANT:
            if g[n].nunique(dropna=True) > 1:
                errors.append(f"{study}: {n} differs between rows {sorted(g[n].dropna().unique())}")

    warnings += _sem_as_sd(clean)
    return errors, warnings, clean


def _fold_check(label: str, g: pd.DataFrame, ctrl: pd.DataFrame):
    """Notes that mention fold or percent change, with mean filled. Returns (is_error, message)
    or None. Many outcomes have values near 100 in their own unit (glucose in mg/dL, systolic
    BP, heart rate), and notes often say a percent change was reported but not used, so the
    wording alone is only a warning. It is an error when the control mean is in the ratio band
    0.8-1.25 and the control baseline is blank or exactly 1, the shape of a series expressed
    as a ratio to baseline (protocol Part 2 section A: fold change is not entered as a mean)."""
    text = " ".join(str(x) for x in g["notes"].dropna())
    if not (FOLD_WORDS.search(text) and g["mean"].notna().any()):
        return None
    msg = (f"{label}: notes describe fold or percent change and mean is filled; ROM of mean "
           "fold changes is a different estimand from a ratio of final means (enter final "
           "values in the outcome unit)")
    c = ctrl[ctrl["mean"].notna()]
    if not len(c):
        return False, msg
    cm, cb = float(c["mean"].iloc[0]), c["baseline_mean"].iloc[0]
    ratio_base = _blank(cb) or abs(float(cb) - 1) < 1e-9
    pct_base = _blank(cb) or abs(float(cb) - 100) < 1e-9
    if 0.8 <= cm <= 1.25 and ratio_base:
        return True, msg + (" -- control mean is near 1 with no baseline in the outcome unit, "
                            "so these look like ratios to baseline. If they are values in the "
                            "outcome unit, fill baseline_mean on the control row")
    if 80 <= cm <= 125 and pct_base:
        return False, msg + (" -- control mean is near 100: check these are not percentages "
                             "of baseline")
    return False, msg


def _sem_as_sd(clean: pd.DataFrame) -> list[str]:
    """Flag arms whose coefficient of variation (SD / mean) is far below that of the other
    studies of the same outcome: an SEM entered as an SD is about sqrt(n) times too small and
    makes the trial look far too precise (Handbook 6.5.2: SDs and SEs are occasionally
    confused). Needs final means > 0 and at least 3 other studies."""
    out = []
    d = clean[clean["mean"].notna() & clean["sd"].notna()].copy()
    d = d[d["mean"].astype(float) > 0]
    if not len(d):
        return out
    d["cv"] = d["sd"].astype(float) / d["mean"].astype(float)
    for outcome, g in d.groupby("outcome"):
        per_study = g.groupby("study_id")["cv"].median()
        if len(per_study) < 4:
            continue
        for study, cv in per_study.items():
            others = per_study.drop(study)
            ref = float(others.median())
            if cv < CV_RATIO * ref:
                rows = g[g["study_id"] == study]
                n = rows["n"].astype(float).median()
                hint = (f"; SD x sqrt(n) = {cv * math.sqrt(n):.3g} x mean, close to the others"
                        if n and abs(cv * math.sqrt(n) - ref) < 0.5 * ref else "")
                out.append(f"{study}/{outcome}: SD/mean {cv:.3g} is below half the median of "
                           f"the other studies ({ref:.3g}); check that an SEM was not entered "
                           f"as an SD{hint}")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("path", type=Path)
    ap.add_argument("--write-csv", type=Path, default=None)
    ap.add_argument("--extra-columns", default="",
                    help="comma-separated subgroup or moderator columns named in the analysis "
                         "config that are not interface columns")
    args = ap.parse_args(argv)
    df = load(args.path)
    extra = tuple(c.strip() for c in args.extra_columns.split(",") if c.strip())
    errors, warnings, clean = check(df, extra)
    for w in warnings:
        print(f"WARNING {w}")
    for e in errors:
        print(f"ERROR {e}")
    print(f"{len(df)} rows, {len(errors)} errors, {len(warnings)} warnings")
    if errors:
        return 1
    if args.write_csv:
        out = clean.copy()
        out["is_control"] = out["is_control"].astype(int)
        out["n"] = out["n"].astype(int)
        blank = out["timepoint_weeks"].isna()
        out.loc[blank, "timepoint_weeks"] = out.loc[blank, "duration_weeks"]
        out.to_csv(args.write_csv, index=False)
        print(f"wrote {args.write_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
