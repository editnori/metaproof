"""Write SYNTHETIC_extraction_example.csv: invented numbers that exercise every input path.

SYNTHETIC. These rows are generated from a fixed random seed. They are not results of any
real trial and must never be presented as nicotinamide riboside evidence.

Run: .venv/bin/python tests/data/make_SYNTHETIC_example.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

COLUMNS = ["study_id", "outcome", "timepoint", "design", "arm", "is_control", "n", "mean", "sd",
           "se", "median", "q1", "q3", "min", "max", "mean_change", "sd_change",
           "baseline_mean", "baseline_sd", "corr_prepost", "corr_within", "paired_md",
           "paired_se", "dose_mg_day", "duration_weeks", "population", "risk_of_bias", "notes"]
TAG = "SYNTHETIC example row; invented numbers, not real trial data"


def main(path: Path) -> None:
    rng = np.random.default_rng(20261101)
    rows = []

    def add(**kw):
        rows.append({c: kw.get(c) for c in COLUMNS} | {"notes": f"{TAG}. {kw.get('notes', '')}".strip()})

    # Outcome 1: a positive biomarker for ROM; 11 studies so Egger runs.
    out = "SYNTHETIC_blood_marker"
    pops = ["SYNTHETIC_healthy", "SYNTHETIC_older"]
    for i in range(1, 12):
        sid = f"SYNTH-A{i:02d}"
        n = int(rng.integers(10, 60))
        base = 25 + rng.normal(0, 3)
        ratio = np.exp(rng.normal(0.35, 0.15))
        dose = [300, 500, 1000][i % 3]
        common = dict(study_id=sid, outcome=out, timepoint="12 weeks", dose_mg_day=dose,
                      duration_weeks=12, population=pops[i % 2],
                      risk_of_bias=["low", "some concerns"][i % 2])
        sd_c, sd_e = base * 0.25, base * ratio * 0.28
        if i == 5:   # three-arm trial: two doses share one control
            add(**common, design="parallel", arm="placebo", is_control=1, n=n, mean=base, sd=sd_c)
            add(**{**common, "dose_mg_day": 300}, design="parallel", arm="low_dose", is_control=0,
                n=n, mean=base * ratio * 0.9, sd=sd_e)
            add(**{**common, "dose_mg_day": 1000}, design="parallel", arm="high_dose",
                is_control=0, n=n, mean=base * ratio * 1.1, sd=sd_e)
        elif i == 7:  # crossover reported as period means and SDs
            add(**common, design="crossover", arm="placebo", is_control=1, n=n, mean=base, sd=sd_c)
            add(**common, design="crossover", arm="active", is_control=0, n=n,
                mean=base * ratio, sd=sd_e)
        elif i == 9:  # medians and IQRs only
            add(**common, design="parallel", arm="placebo", is_control=1, n=n, median=base,
                q1=base - 0.67 * sd_c, q3=base + 0.7 * sd_c, notes="median/IQR only")
            add(**common, design="parallel", arm="active", is_control=0, n=n,
                median=base * ratio, q1=base * ratio - 0.6 * sd_e, q3=base * ratio + 0.75 * sd_e,
                notes="median/IQR only")
        elif i == 10:  # SEs instead of SDs
            add(**common, design="parallel", arm="placebo", is_control=1, n=n, mean=base,
                se=sd_c / np.sqrt(n), notes="SE reported")
            add(**common, design="parallel", arm="active", is_control=0, n=n, mean=base * ratio,
                se=sd_e / np.sqrt(n), notes="SE reported")
        else:
            if i == 3:  # an extra, earlier timepoint the 'last' rule must drop
                add(**{**common, "timepoint": "4 weeks"}, design="parallel", arm="placebo",
                    is_control=1, n=n, mean=base, sd=sd_c)
                add(**{**common, "timepoint": "4 weeks"}, design="parallel", arm="active",
                    is_control=0, n=n, mean=base * 1.1, sd=sd_e)
            add(**common, design="parallel", arm="placebo", is_control=1, n=n, mean=base, sd=sd_c)
            add(**common, design="parallel", arm="active", is_control=0, n=n,
                mean=base * ratio, sd=sd_e)

    # Outcome 2: change from baseline, MD.
    out = "SYNTHETIC_systolic_bp"
    for i in range(1, 7):
        sid = f"SYNTH-B{i:02d}"
        n = int(rng.integers(12, 40))
        b = 130 + rng.normal(0, 4)
        eff = rng.normal(-3, 2)
        common = dict(study_id=sid, outcome=out, timepoint="8 weeks", dose_mg_day=1000,
                      duration_weeks=8, population=pops[i % 2], risk_of_bias="low")
        if i == 6:  # crossover with a reported paired difference
            add(**common, design="crossover", arm="active", is_control=0, n=n, mean=b + eff,
                sd=12, paired_md=eff, paired_se=8 / np.sqrt(n))
            add(**common, design="crossover", arm="placebo", is_control=1, n=n, mean=b, sd=12)
            continue
        for arm, ctl, delta in (("placebo", 1, rng.normal(0, 1)), ("active", 0, eff)):
            r = dict(common, design="parallel", arm=arm, is_control=ctl, n=n,
                     baseline_mean=b, baseline_sd=12, mean=b + delta, sd=12.5)
            if i in (1, 2):  # change stats reported directly
                r.update(mean_change=delta, sd_change=9.0)
            if i == 3:
                r.update(corr_prepost=0.7, notes="correlation reported")
            add(**r)

    # Outcome 3: SMD with one study that has two placebo arms.
    out = "SYNTHETIC_grip_strength"
    for i in range(1, 5):
        sid = f"SYNTH-C{i:02d}"
        n = int(rng.integers(15, 45))
        common = dict(study_id=sid, outcome=out, timepoint="6 weeks", dose_mg_day=500,
                      duration_weeks=6, population="SYNTHETIC_older", risk_of_bias="high")
        add(**common, design="parallel", arm="placebo", is_control=1, n=n, mean=30, sd=6)
        if i == 2:
            add(**common, design="parallel", arm="placebo_b", is_control=1, n=n, mean=31, sd=7)
        add(**common, design="parallel", arm="active", is_control=0, n=n,
            mean=30 + rng.normal(1.5, 1), sd=6.5)

    df = pd.DataFrame(rows, columns=COLUMNS)
    df.to_csv(path, index=False, float_format="%.4f")


if __name__ == "__main__":
    main(Path(__file__).with_name("SYNTHETIC_extraction_example.csv"))
