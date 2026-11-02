# metaproof

Meta-analysis you can check. You fill in an extraction sheet, metaproof pools the studies the way the Cochrane Handbook says to, and a second engine built on metafor (the R package most reviewers already trust) runs the same sheet so you can see both give the same answer.

Built by Dr. Layth Qassem, PharmD, MSACI. MIT license. How to cite it is in [REFERENCES.md](REFERENCES.md).

## What's in it

- `metaproof/` the main engine, in Python
- `reference_r/` the second engine, on metafor
- `templates/` protocol with every PROSPERO field, extraction workbook and checker, RoB 2, GRADE and PRISMA 2020
- `validation/` the Python vs R cross-check and the proofs
- `examples/nicotinamide-riboside/` a worked example for a nicotinamide riboside (NR) review, with a search strategy, a candidate trial list and a made-up extraction sheet. Every example number is labeled SYNTHETIC and none of it is trial data.

## Setup on a Mac

```
xcode-select --install
curl -LsSf https://astral.sh/uv/install.sh | sh
git clone https://github.com/editnori/metaproof
cd metaproof
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python --require-hashes -r requirements-lock.txt
```

You only need that much to run it. The R engine is optional:

```
brew install r
Rscript reference_r/install_packages.R
```

It installs metafor 5.2-1 from a dated CRAN snapshot, so you get the exact versions this was checked with.

## Running it on your review

1. Fill in `templates/extraction_sheet.xlsx`, one row per arm, outcome and timepoint. The `data_dictionary` sheet explains each column. `templates/HOW_TO_USE_TEMPLATES.txt` walks through each template.
2. Check the sheet and export it:
   `.venv/bin/python templates/check_extraction.py extraction_sheet.xlsx --write-csv extraction.csv`
3. Write a config listing your outcomes. `examples/nicotinamide-riboside/SYNTHETIC_nr_like_config.yaml` is a full example.
4. Run it:
   `.venv/bin/python -m metaproof run extraction.csv --config config.yaml --out out`
5. Optional, run the R engine on the same files and compare:
   `Rscript reference_r/run_meta.R --data extraction.csv --config config.yaml --out out_r`

Each outcome gets its own output with the pooled result, every study's effect, leave-one-out, influence checks, a forest plot and a funnel plot. Anything left out is listed in `excluded.csv` with the reason.

## What it does and the standard behind each step

- Pools studies with random effects (REML), which handles a small number of studies better than the older DerSimonian-Laird method. Cochrane Handbook 10.10.4, Veroniki 2016
- Widens the confidence interval with Hartung-Knapp when studies disagree and there are at least 3. Handbook 10.10.4.4
- Doesn't pool fewer than 3 studies
- Gives a prediction interval, the range where a new study's effect would likely land, and flags it when there are fewer than 5 studies. Handbook 10.10.4.3, IntHout 2016
- Reports heterogeneity (I², tau²) with a confidence interval. Handbook 10.10.2, Viechtbauer 2007
- Mean difference keeps the clinical unit. Hedges g when studies use different scales. Ratio of means for outcomes like NAD+ levels. Handbook 6.5.1, 10.5.2
- Fills a missing change-score SD from a correlation of 0.5, and reruns at 0.2 and 0.8 to show whether that choice matters. Handbook 6.5.2.8
- Handles crossover trials as paired data. Handbook 23.2, Elbourne 2002
- Combines multiple NR dose arms so the placebo group is only counted once. Handbook 23.3.4
- Converts medians and IQRs to means and SDs. Wan 2014, Luo 2018
- Runs subgroups and meta-regression only with 10 or more studies (Handbook 10.11.5.1), and the funnel plot test the same way (Handbook 13.3.4.4, Sterne 2011)
- Templates follow PRISMA 2020 (Page 2021), RoB 2 (Sterne 2019) and GRADE (Guyatt 2008)

Every default, its reason and its sensitivity check are in `templates/protocol_template.txt` (section F onward). Full formulas are in `metaproof/METHODS.txt`.

## How it's checked

Both engines run three published metafor datasets, the NR example and a set of made-up edge cases, and every value has to match within a relative 1e-6 (1e-4 for intervals and tests). On the last run 696 checks passed and 0 failed, comparing 28,904 values, and checks that don't apply to an analysis are listed as skipped instead of passed (`validation/crosscheck_report.txt`). On top of that, 346 automated tests cover every formula with worked numbers, and the R engine checks itself against metafor (597 passed, 0 failed). Claude and OpenAI Codex each reviewed every file on their own, and every finding was fixed with a test or rejected with a written reason.

The formulas are also proved. SymPy proves they behave for any number of studies: weights add to 1, the pooled estimate always sits between the smallest and largest study, I² stays between 0 and 1, and the prediction interval always contains the confidence interval. Bend2 proves the counting rules on whole numbers, like splitting a placebo group across dose arms keeping the total, and five deliberately broken versions get rejected. `validation/proofs/RECEIPTS.txt` says exactly what each proof covers and what it doesn't.

## Limits

- Continuous outcomes only. No risk ratios or odds ratios yet.
- When a study has several dose arms split for a dose analysis, Hartung-Knapp still counts comparisons, not studies.
- No skewness test for medians (Shi 2023). Check that by hand.
- The Embase, CENTRAL and ICTRP searches in the NR example are written but not run. Rows in the NR trial list marked UNVERIFIED still need checking against the paper.
- Rerunning the Bend2 proofs needs a Linux machine. See `validation/proofs/README.txt`.

## Your data

Your extraction sheet stays on your computer. metaproof doesn't send anything anywhere.
