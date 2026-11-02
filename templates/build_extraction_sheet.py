"""Build the extraction workbook, the data dictionary CSV and the header-only CSV.

Usage:
    .venv/bin/python templates/build_extraction_sheet.py [--outdir templates]

Writes:
    extraction_sheet.xlsx        sheets: README, extraction, data_dictionary,
                                 data_patterns, validation_rules, lists,
                                 study_characteristics
    extraction_data_dictionary.csv
    extraction_template.csv      header row only: interface columns, then the
                                 extension columns (se_change, timepoint_weeks)

Everything comes from extraction_schema.py, so rerun this script after any
change to the schema.
"""

from __future__ import annotations

import argparse
import csv
from itertools import zip_longest
from pathlib import Path

from openpyxl import Workbook
from openpyxl.comments import Comment
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.datavalidation import DataValidation

from extraction_schema import (ALL_COLUMNS, ALL_NAMES, DATA_PATTERNS, DESIGNS,
                               IS_CONTROL, ROB_LEVELS, VALIDATION_RULES)

# The workbook carries the interface columns followed by the optional extension columns.
COLUMNS, COLUMN_NAMES = ALL_COLUMNS, ALL_NAMES

MAX_ROWS = 2000  # rows of the extraction sheet that carry validation rules

STUDY_CHAR_COLUMNS = [
    ("study_id", "Same label as the extraction sheet."),
    ("citation", "Full reference of the main report."),
    ("other_reports", "Other reports of the same study (registry record, abstracts, secondary papers)."),
    ("registration_id", "Trial registry ID (e.g. NCT number) or 'none found'."),
    ("country", "Country or countries of recruitment."),
    ("setting", "Single or multi-centre; outpatient, inpatient, community."),
    ("design", "parallel or crossover."),
    ("washout_weeks", "Crossover only: washout length in weeks between periods; 'none' if the arms switched with no washout, 'not reported' if the report does not say (RoB 2 crossover domain S)."),
    ("period_length_weeks", "Crossover only: length of each treatment period in weeks."),
    ("n_per_sequence", "Crossover only: participants randomised and analysed in each sequence (e.g. AB 10/9, BA 10/10); needed to judge period-effect balance (RoB 2 domain S)."),
    ("crossover_analysis", "Crossover only: how the paired estimate was obtained: period-adjusted model | paired test without period term | built by the kit from period means (assumes no period effect) | first period only (allowed only if planned in the protocol, Handbook 23.2.4)."),
    ("result_selection", "When the study reports several eligible results for one outcome (matrix, assay, metabolite, timepoint, adjusted or unadjusted): which one was extracted and the protocol hierarchy rule that chose it; any deviation feeds RoB 2 question 5.2."),
    ("n_randomised", "Total randomised."),
    ("population_description", "Inclusion criteria in brief, age, sex, BMI, condition."),
    ("intervention", "Product, salt form (NR chloride), dose, frequency, formulation."),
    ("co_intervention", "Anything given with NR (e.g. pterostilbene). Combination arms are eligible only if the comparator receives the same co-intervention."),
    ("comparator", "Placebo, no treatment, or other."),
    ("adherence_method", "How adherence was measured and the reported adherence."),
    ("outcomes_reported", "Outcome keys reported, and outcomes registered but not reported."),
    ("funding_source", "Funder(s); industry supply of product."),
    ("author_coi", "Declared conflicts of interest."),
    ("contacted_authors", "Date contacted, what was asked, reply."),
    ("extractor_1", "Initials."),
    ("extractor_2", "Initials."),
    ("notes", "Anything else."),
]

HEADER_FILL = PatternFill("solid", fgColor="DDEBF7")
BOLD = Font(bold=True)
WRAP = Alignment(wrap_text=True, vertical="top")


def _write_table(ws, header, rows, widths=None):
    ws.append(header)
    for cell in ws[1]:
        cell.font = BOLD
        cell.fill = HEADER_FILL
        cell.alignment = WRAP
    for r in rows:
        ws.append(list(r))
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = WRAP
    ws.freeze_panes = "A2"
    if widths:
        for i, w in enumerate(widths, start=1):
            ws.column_dimensions[get_column_letter(i)].width = w


def dictionary_rows():
    rows = []
    for c in COLUMNS:
        if c["kind"] == "choice":
            allowed = " | ".join(str(x) for x in c["choices"])
        elif c["kind"] in ("int", "number"):
            lo, hi = c.get("lo"), c.get("hi")
            parts = ["whole number" if c["kind"] == "int" else "number"]
            if c.get("gt0"):
                parts.append("> 0")
            elif lo is not None:
                parts.append(f">= {lo}")
            if hi is not None:
                parts.append(f"<= {hi}")
            allowed = ", ".join(parts)
        else:
            allowed = "text"
        rows.append((c["name"], c["required"], allowed, c["unit"], c["definition"]))
    return rows


def _add_validations(ws):
    col = {name: get_column_letter(i) for i, name in enumerate(COLUMN_NAMES, start=1)}
    last = MAX_ROWS + 1

    def rng(name):
        return f"{col[name]}2:{col[name]}{last}"

    def add(dv, *names):
        dv.showErrorMessage = True
        dv.errorStyle = "stop"
        for n in names:
            dv.add(rng(n))
        ws.add_data_validation(dv)

    add(DataValidation(type="list", formula1='"' + ",".join(DESIGNS) + '"', allow_blank=True,
                       errorTitle="design", error="Use parallel or crossover."), "design")
    add(DataValidation(type="list", formula1='"0,1"', allow_blank=True,
                       errorTitle="is_control", error="Use 1 for control, 0 for NR arm."),
        "is_control")
    add(DataValidation(type="list", formula1='"' + ",".join(ROB_LEVELS) + '"', allow_blank=True,
                       errorTitle="risk_of_bias", error="Use low, some concerns or high."),
        "risk_of_bias")
    add(DataValidation(type="whole", operator="greaterThanOrEqual", formula1="1",
                       allow_blank=True, errorTitle="n", error="Whole number, at least 1."), "n")

    positive = [c["name"] for c in COLUMNS if c.get("gt0")]
    add(DataValidation(type="decimal", operator="greaterThan", formula1="0", allow_blank=True,
                       errorTitle="Must be > 0", error="paired_se must be greater than 0."),
        *[p for p in positive if p != "duration_weeks"])
    add(DataValidation(type="decimal", operator="greaterThan", formula1="0", allow_blank=True,
                       errorTitle="duration_weeks", error="Weeks, greater than 0."),
        "duration_weeks")
    nonneg = [c["name"] for c in COLUMNS if c.get("lo") == 0 and not c.get("gt0")
              and c["name"] not in ("timepoint_weeks", "dose_mg_day")]
    add(DataValidation(type="decimal", operator="greaterThanOrEqual", formula1="0",
                       allow_blank=True, errorTitle="Must be >= 0",
                       error="SDs and SEs cannot be negative. An SD of 0 is allowed only if "
                             "another arm of the comparison has an SD above 0."),
        *nonneg)
    add(DataValidation(type="decimal", operator="greaterThanOrEqual", formula1="0",
                       allow_blank=True, errorTitle="timepoint_weeks",
                       error="Weeks, 0 or more."), "timepoint_weeks")
    add(DataValidation(type="decimal", operator="greaterThanOrEqual", formula1="0",
                       allow_blank=True, errorTitle="dose_mg_day",
                       error="mg/day, 0 for the control arm."), "dose_mg_day")
    add(DataValidation(type="decimal", operator="between", formula1="-1", formula2="1",
                       allow_blank=True, errorTitle="Correlation",
                       error="Correlation must be between -1 and 1."),
        "corr_prepost", "corr_within")
    others = [c["name"] for c in COLUMNS
              if c["kind"] == "number" and not c.get("gt0") and c.get("lo") is None
              and c.get("hi") is None]
    add(DataValidation(type="decimal", operator="between", formula1="-1E+307",
                       formula2="1E+307", allow_blank=True, errorTitle="Number",
                       error="Enter a number only. Put units and comments in notes."),
        *others)


def build(outdir: Path) -> dict[str, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    wb = Workbook()

    ws = wb.active
    ws.title = "README"
    readme = [
        "NR meta-analysis extraction workbook (template, contains no study data).",
        "One row per study arm x outcome x timepoint on the 'extraction' sheet.",
        "Column names and order are the shared interface; do not rename or reorder them.",
        "se_change and timepoint_weeks, after notes, are optional extension columns that both engines read.",
        "A subgroup or moderator column named in the analysis config may follow them; declare it to the",
        "checker with --extra-columns <name>. Any other added column is an error.",
        "Enter an SEM in se (or se_change), never in sd; many NR papers report mean +/- SEM.",
        "Leave a cell blank when the value is not reported. Never type NR, NA, '-' or units into number cells.",
        "Enter numbers in the protocol unit for the outcome; record any conversion in notes (see unit_conversions.csv).",
        "Two extractors fill separate copies; disagreements are resolved before analysis.",
        "Export: save the 'extraction' sheet as CSV UTF-8, or run",
        "  .venv/bin/python templates/check_extraction.py extraction_sheet.xlsx --write-csv extraction.csv",
        "which checks the rules on the 'validation_rules' sheet and writes the CSV the engines read.",
        "Study-level details go on 'study_characteristics' (not read by the engines).",
    ]
    for line in readme:
        ws.append([line])
    ws.column_dimensions["A"].width = 120

    ws = wb.create_sheet("extraction")
    ws.append(COLUMN_NAMES)
    for i, cell in enumerate(ws[1], start=1):
        cell.font = BOLD
        cell.fill = HEADER_FILL
        ws.column_dimensions[get_column_letter(i)].width = max(10, len(cell.value) + 2)
    ws.column_dimensions[get_column_letter(COLUMN_NAMES.index("notes") + 1)].width = 60
    ws.freeze_panes = "C2"
    _add_validations(ws)

    _write_table(wb.create_sheet("data_dictionary"),
                 ["column", "required", "allowed values", "unit", "definition"],
                 dictionary_rows(), widths=[18, 12, 26, 16, 100])
    _write_table(wb.create_sheet("data_patterns"),
                 ["pattern", "columns that must all be filled"],
                 [(label, ", ".join(cols)) for label, cols in DATA_PATTERNS],
                 widths=[70, 60])
    _write_table(wb.create_sheet("validation_rules"), ["rule"],
                 [(r,) for r in VALIDATION_RULES], widths=[130])
    _write_table(wb.create_sheet("lists"), ["design", "is_control", "risk_of_bias"],
                 zip_longest(DESIGNS, IS_CONTROL, ROB_LEVELS), widths=[14, 12, 16])
    sc = wb.create_sheet("study_characteristics")
    _write_table(sc, [c for c, _ in STUDY_CHAR_COLUMNS], [], widths=[22] * len(STUDY_CHAR_COLUMNS))
    for i, (_, desc) in enumerate(STUDY_CHAR_COLUMNS, start=1):
        sc.cell(row=1, column=i).comment = Comment(desc, "template")

    paths = {
        "xlsx": outdir / "extraction_sheet.xlsx",
        "dictionary": outdir / "extraction_data_dictionary.csv",
        "header": outdir / "extraction_template.csv",
    }
    wb.save(paths["xlsx"])
    with open(paths["dictionary"], "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["column", "required", "allowed_values", "unit", "definition"])
        w.writerows(dictionary_rows())
    with open(paths["header"], "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerow(COLUMN_NAMES)
    return paths


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent)
    args = ap.parse_args()
    for kind, p in build(args.outdir).items():
        print(f"{kind}: {p}")


if __name__ == "__main__":
    main()
