# metaproof

metaproof runs the stats for a systematic review and meta-analysis. You put each study's numbers into one spreadsheet. It pools them into one effect estimate per outcome, with forest and funnel plots, following the Cochrane Handbook. A second copy of the same analysis runs in R on metafor, the package most meta-analyses use, so you can check the two match.

Made by Dr. Layth Qassem, PharmD, MSACI. MIT license. See [REFERENCES.md](REFERENCES.md) for how to cite it and the methods it follows.

## What it measures

For each outcome you list, metaproof gives:

- the pooled effect and its 95% confidence interval, as a mean difference in the outcome's own unit, Hedges g when studies used different scales, or a ratio of means for lab values like NAD+
- the prediction interval, which is the range a new study's effect would likely fall in
- how much the studies disagree with each other (I² and tau²)
- each study's effect and weight, what happens when each study is left out, and which studies pull the result the most
- a forest plot and a funnel plot

Any study it leaves out goes in `excluded.csv` with the reason.

## Setup on a Mac

Open Terminal and run:

```
xcode-select --install
curl -LsSf https://astral.sh/uv/install.sh | sh
git clone https://github.com/editnori/metaproof
cd metaproof
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python --require-hashes -r requirements-lock.txt
```

The R check is optional. To add it, run `brew install r` and then `Rscript reference_r/install_packages.R`.

## Setup on Windows

Install Git from https://git-scm.com/download/win. Then open PowerShell and run:

```
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
git clone https://github.com/editnori/metaproof
cd metaproof
uv venv --python 3.12 .venv
uv pip install --python .venv\Scripts\python.exe --require-hashes -r requirements-lock.txt
```

Close and reopen PowerShell if `uv` isn't found after installing it.

The R check is optional. To add it, install R from https://cran.r-project.org/bin/windows/base/ and run `Rscript reference_r\install_packages.R`.

On Windows, write `.venv\Scripts\python` wherever the steps below say `.venv/bin/python`.

## Running it

1. Fill in `templates/extraction_sheet.xlsx`, one row per study arm, outcome and timepoint. The `data_dictionary` tab explains each column, and `templates/HOW_TO_USE_TEMPLATES.txt` covers every template.
2. Check the sheet and save it as a CSV:
   `.venv/bin/python templates/check_extraction.py extraction_sheet.xlsx --write-csv extraction.csv`
3. Write a config file that lists your outcomes. `examples/nicotinamide-riboside/SYNTHETIC_nr_like_config.yaml` is a full example.
4. Run it:
   `.venv/bin/python -m metaproof run extraction.csv --config config.yaml --out results`
5. To run the R copy on the same files:
   `Rscript reference_r/run_meta.R --data extraction.csv --config config.yaml --out results_r`

## The methods

Each default and the reason for it is in `templates/protocol_template.txt`.

- Random-effects pooling with REML, which estimates the between-study spread better than DerSimonian-Laird when there are few studies (Cochrane Handbook 10.10.4)
- Hartung-Knapp adjustment to the confidence interval when the studies disagree and there are at least 3 (Handbook 10.10.4.4)
- No pooled number from fewer than 3 studies
- Prediction interval, flagged when there are fewer than 5 studies (Handbook 10.10.4.3)
- A missing change-from-baseline SD is filled from a correlation of 0.5, and the analysis reruns at 0.2 and 0.8 to show if that matters (Handbook 6.5.2.8)
- Crossover trials analyzed as paired data (Handbook 23.2)
- Multiple dose arms combined, so the placebo group is counted once (Handbook 23.3.4)
- Medians and IQRs converted to means and SDs (Wan 2014, Luo 2018)
- Subgroups, meta-regression and the funnel plot test only with 10 or more studies (Handbook 10.11.5.1 and 13.3.4.4)
- Templates for PRISMA 2020, RoB 2 and GRADE

## How it was checked

- The Python and R copies ran on three published datasets that ship with metafor, the NR example and a set of made-up edge cases. 696 checks passed and 0 failed, across 28,904 values. Results are in `validation/crosscheck_report.txt`.
- 346 automated tests check each formula against worked numbers.
- The R copy was checked against metafor's own results: 597 passed, 0 failed.
- Claude and OpenAI Codex each reviewed every file. Every issue they found was fixed with a test or written up with why it was left.
- The formulas were proved with SymPy, for example that the pooled estimate is always between the smallest and largest study and that I² is always between 0 and 1. The counting rules, like splitting a placebo group across dose arms without changing the total, were proved with Bend2. `validation/proofs/RECEIPTS.txt` lists what each proof covers.

## Example

`examples/nicotinamide-riboside/` is a worked example for a review of nicotinamide riboside (NR). It has a search strategy, a list of candidate NR trials and a sample extraction sheet. Every number in the sample sheet is made up and labeled SYNTHETIC.

## Limits

- Continuous outcomes only. Risk ratios and odds ratios aren't in yet.
- No skewness test for medians (Shi 2023), so check that by hand.
- In the NR example, the Embase, CENTRAL and ICTRP searches are written but haven't been run, and rows marked UNVERIFIED in the trial list still need checking.
- Rerunning the Bend2 proofs needs a Linux machine. See `validation/proofs/README.txt`.

Your spreadsheet stays on your computer. metaproof doesn't send it anywhere.
