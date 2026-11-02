# metaproof

metaproof does the math for a meta-analysis. You put each study's results in one spreadsheet and it gives you the pooled result for each outcome, the forest plot and the funnel plot. It follows Cochrane's guide for systematic reviews, and every result was also run in R as a second check.

For each outcome you get the pooled effect with its 95% confidence interval, in the outcome's own unit (like mmHg) or as a ratio for lab levels like NAD+. You also get the range a new study would likely fall in, how much the studies disagree, how much each study counts, and what happens to the result when each study is left out. Any study it can't use goes in `excluded.csv` with the reason.

Made by Dr. Layth Qassem, PharmD, MSACI. Free to use under the MIT license. How to cite it is in [REFERENCES.md](REFERENCES.md).

## Setup

On a Mac, open Terminal and run:

```
xcode-select --install
curl -LsSf https://astral.sh/uv/install.sh | sh
git clone https://github.com/editnori/metaproof
cd metaproof
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python --require-hashes -r requirements.txt
```

On Windows, install Git from https://git-scm.com/download/win, then open PowerShell and run:

```
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
git clone https://github.com/editnori/metaproof
cd metaproof
uv venv --python 3.12 .venv
uv pip install --python .venv\Scripts\python.exe --require-hashes -r requirements.txt
```

If PowerShell can't find `uv`, close it and open it again. On Windows use `.venv\Scripts\python` wherever the commands below say `.venv/bin/python`.

## Use

Fill in `templates/extraction_sheet.xlsx`, one row per study group, outcome and time point. The `data_dictionary` tab says what goes in each column. Then check it for mistakes and save it as a CSV:

```
.venv/bin/python templates/check_extraction.py extraction_sheet.xlsx --write-csv extraction.csv
```

Write a short settings file listing your outcomes (`examples/nicotinamide-riboside/SYNTHETIC_nr_like_config.yaml` is a full example), then run it:

```
.venv/bin/python -m metaproof run extraction.csv --config config.yaml --out results
```

`examples/nicotinamide-riboside/` has a search strategy, a list of candidate NR trials and a sample spreadsheet. The numbers in the sample spreadsheet are made up and labeled SYNTHETIC. Only the PubMed and ClinicalTrials.gov searches have been run, and trial list rows marked UNVERIFIED still need checking against the paper.

## Methods

It uses a random-effects model, since studies differ in dose, length and population. The confidence interval gets wider when the studies disagree, and it won't pool fewer than 3 studies. If a study reports change from baseline without its SD, it fills one in with a standard assumption and reruns at two others so you can see if that choice matters. Crossover trials are analyzed as paired data, and when a trial has several doses against one placebo group, that placebo group is only counted once. Medians and IQRs get converted to means and SDs. Subgroups and the funnel plot test only run with 10 or more studies because they aren't reliable with fewer.

Every choice and its reason is in `templates/protocol_template.txt`, which you can use as your protocol. The `templates` folder also has the PRISMA 2020 flow diagram, the RoB 2 risk of bias tool and a GRADE summary of findings table.

## Checks

Every result was run in Python and again in R with metafor, the R package most meta-analyses are done in. Across three published datasets, the NR example and a set of edge cases, 696 checks matched and 0 didn't. 346 automated tests check each calculation against worked examples. The core formulas were also proved mathematically, like the pooled result always falling between the smallest and largest study (`validation/proofs/RECEIPTS.txt`).

## Limits

It only handles continuous outcomes like blood pressure or lab levels for now. Yes/no outcomes like adverse events aren't in yet. Your spreadsheet stays on your computer and metaproof doesn't send it anywhere.
