# metaproof

metaproof does the math for a meta-analysis. You enter each study's results in one spreadsheet, and it combines them into one overall result per outcome, with the standard plots for a paper.

It follows Cochrane's official guide for systematic reviews. The same analysis also runs a second time in a separate statistics program as a check, and the two give the same numbers.

Made by Dr. Layth Qassem, PharmD, MSACI. Free to use under the MIT license. See [REFERENCES.md](REFERENCES.md) for how to cite it.

## What you get for each outcome

- The overall effect across studies with its 95% confidence interval. It's a mean difference in the outcome's own unit (for example mmHg), a standardized difference when studies used different scales, or a ratio when the outcome is a lab level like NAD+
- A prediction interval, the range where a future study's result would likely fall
- How much the studies disagree with each other
- Each study's result and how much it counts toward the total
- The result with each study left out one at a time, and which studies change it the most
- A forest plot and a funnel plot

Any study it can't use goes in `excluded.csv` with the reason.

## Setup on a Mac

Open Terminal and run:

```
xcode-select --install
curl -LsSf https://astral.sh/uv/install.sh | sh
git clone https://github.com/editnori/metaproof
cd metaproof
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python --require-hashes -r requirements.txt
```

## Setup on Windows

Install Git from https://git-scm.com/download/win, then open PowerShell and run:

```
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
git clone https://github.com/editnori/metaproof
cd metaproof
uv venv --python 3.12 .venv
uv pip install --python .venv\Scripts\python.exe --require-hashes -r requirements.txt
```

If PowerShell says `uv` isn't found, close it and open it again. On Windows, use `.venv\Scripts\python` wherever the steps below say `.venv/bin/python`.

## Running it

1. Fill in `templates/extraction_sheet.xlsx`, one row per study group, outcome and time point. The `data_dictionary` tab explains each column.
2. Check the sheet for mistakes and save it as a CSV:
   `.venv/bin/python templates/check_extraction.py extraction_sheet.xlsx --write-csv extraction.csv`
3. Write a short settings file that lists your outcomes. `examples/nicotinamide-riboside/SYNTHETIC_nr_like_config.yaml` is a full example.
4. Run it:
   `.venv/bin/python -m metaproof run extraction.csv --config config.yaml --out results`

## How it combines studies

Every choice and the reason for it is written out in `templates/protocol_template.txt`, which you can use as your review protocol.

- It assumes the true effect can differ from study to study (a random-effects model), since studies differ in dose, length and population
- It widens the confidence interval when the studies disagree, so the result doesn't look more certain than it is
- It won't combine fewer than 3 studies
- When a study only reports the change from baseline without its spread, it fills that in with a standard assumption and reruns with two other assumptions so you can see if it matters
- Crossover trials, where the same people get NR and placebo, are analyzed as paired data
- When a trial has several NR doses against one placebo group, the placebo group is counted once
- Medians and IQRs are converted to means and SDs with the standard published formulas
- Subgroup analyses and the funnel plot test only run with 10 or more studies, since with fewer they are not reliable

The `templates` folder also has the PRISMA 2020 flow diagram and checklist, the risk of bias tool (RoB 2) and a GRADE summary of findings table.

## How it was checked

- Every result was run twice, once in Python and once in R with metafor, the statistics package most published meta-analyses use. On three published datasets, the NR example and a set of tricky made-up cases, 696 checks matched and 0 didn't.
- 346 automated tests check each calculation against worked examples.
- Two separate AI reviewers (Claude and OpenAI Codex) read every file, and every issue they found was fixed and tested.
- The core formulas were proved mathematically, for example that the overall result always falls between the smallest and largest study. `validation/proofs/RECEIPTS.txt` lists each proof.

## NR example

`examples/nicotinamide-riboside/` has a search strategy, a list of candidate NR trials and a sample spreadsheet. Every number in the sample spreadsheet is made up and labeled SYNTHETIC.

## Limits

- It handles continuous outcomes like blood pressure or lab levels. Yes/no outcomes like events or side effects aren't in yet.
- In the NR example, only the PubMed and ClinicalTrials.gov searches have been run. Trial list rows marked UNVERIFIED still need checking against the paper.

Your spreadsheet stays on your computer. metaproof doesn't send it anywhere.
