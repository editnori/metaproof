"""Engine conformance: the formulas proven in sympy_proofs.py, checked against the metaproof code.

sympy_proofs.py proves properties of formulas written inside that file. That says nothing
about the engine unless the engine computes the same formulas. This file closes the gap: each
check builds the symbolic form with sympy, evaluates it on random SYNTHETIC inputs, and
asserts it equals what the metaproof function returns. It also runs engine-level properties that
no formula proof can cover (unit of analysis: control participants counted once, one
timepoint per study, dose taken from experimental rows only).

Run from the repo root:
    .venv/bin/python validation/proofs/engine_conformance.py
Writes receipt_engine.txt beside this file; exits 0 only when every check passes.

Every input generated here is SYNTHETIC (random numbers from a fixed seed), not trial data.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from metaproof import conversions as cv  # noqa: E402
from metaproof import effects as ef  # noqa: E402
from metaproof import pooling as pl  # noqa: E402
from metaproof.pipeline import select_timepoint  # noqa: E402

SEED = 20261101
TRIALS = 200
RTOL = 1e-9


@dataclass
class Check:
    cid: str
    prop: str
    engine: str
    n: int
    max_rel: float
    ok: bool
    detail: str = ""


CHECKS = []


def check(fn):
    CHECKS.append(fn)
    return fn


def rel(a, b):
    return abs(a - b) / max(1.0, abs(b))


# --- symbolic forms (sympy expressions, lambdified for numeric comparison) ---------------

def _sym_pooled(k):
    y = sp.symbols(f"y1:{k + 1}")
    v = sp.symbols(f"v1:{k + 1}", positive=True)
    t = sp.Symbol("t", nonnegative=True)
    w = [1 / (vi + t) for vi in v]
    mu = sum(wi * yi for wi, yi in zip(w, y)) / sum(w)
    var = 1 / sum(w)
    return sp.lambdify((y, v, t), (mu, var), "math")


def _sym_dl(k):
    y = sp.symbols(f"y1:{k + 1}")
    v = sp.symbols(f"v1:{k + 1}", positive=True)
    w = [1 / vi for vi in v]
    S1, S2 = sum(w), sum(wi**2 for wi in w)
    mu = sum(wi * yi for wi, yi in zip(w, y)) / S1
    Q = sum(wi * (yi - mu) ** 2 for wi, yi in zip(w, y))
    C = S1 - S2 / S1
    tau2 = sp.Max(0, (Q - (k - 1)) / C)
    i2_q = sp.Max(0, (Q - (k - 1)) / Q)
    return sp.lambdify((y, v), (tau2, Q, C, i2_q), "math")


_DF = sp.Symbol("m", positive=True)
# mpmath: the gamma function overflows double precision beyond m of about 340
SYM_J_EXACT = sp.lambdify(_DF, sp.gamma(_DF / 2) / (sp.sqrt(_DF / 2) * sp.gamma((_DF - 1) / 2)),
                          "mpmath")
SYM_J_APPROX = sp.lambdify(_DF, 1 - sp.Integer(3) / (4 * _DF - 1), "math")
_sb, _sf, _r, _n = sp.symbols("sb sf r n", positive=True)
SYM_CHANGE_SD = sp.lambdify((_sb, _sf, _r), sp.sqrt(_sb**2 + _sf**2 - 2 * _r * _sb * _sf), "math")
SYM_PAIRED_VAR = sp.lambdify((_sb, _sf, _r, _n), (_sb**2 + _sf**2 - 2 * _r * _sb * _sf) / _n,
                             "math")


def _draw(rng, k):
    y = rng.normal(0.4, 0.5, k)
    v = rng.uniform(0.01, 0.3, k)
    return y, v


# --- formula conformance --------------------------------------------------------------------

@check
def e01_pooled_estimate():
    """P01/P02: inverse-variance pooled estimate and its model variance, at the engine's tau2."""
    rng = np.random.default_rng(SEED)
    worst, n = 0.0, 0
    for k in range(2, 9):
        f_sym = _sym_pooled(k)
        for _ in range(TRIALS // 7):
            y, v = _draw(rng, k)
            for method in ("FE", "DL", "REML"):
                f = pl.fit(y, v, method=method, hksj=False)
                mu, var = f_sym(tuple(y), tuple(v), f.tau2)
                worst = max(worst, rel(f.estimate, mu), rel(float(f.vb_model[0, 0]), var))
                n += 2
    return Check("E01", "pooled estimate and model variance = sum(w y)/sum(w), 1/sum(w)",
                 "pooling.fit", n, worst, worst < RTOL)


@check
def e02_dl_tau2_and_i2():
    """P03/P04: DL tau2, Cochran Q and I^2. The engine's random-effects rows report
    I^2 = tau2 / (tau2 + s^2), s^2 = df / C; under DL this equals the Q-based
    max(0, (Q - df)/Q) that P04 proves (the identity is step 5 of P04). FE rows use the
    Q-based form directly."""
    rng = np.random.default_rng(SEED + 1)
    worst, n = 0.0, 0
    for k in range(2, 9):
        f_sym = _sym_dl(k)
        for _ in range(TRIALS // 7):
            y, v = _draw(rng, k)
            tau2, Q, C, i2_q = (float(x) for x in f_sym(tuple(y), tuple(v)))
            dl = pl.fit(y, v, method="DL", hksj=False)
            fe = pl.fit(y, v, method="FE")
            worst = max(worst, rel(dl.tau2, tau2), rel(dl.QE, Q), rel(dl.I2, 100 * i2_q),
                        rel(fe.I2, 100 * i2_q))
            reml = pl.fit(y, v, method="REML", hksj=False)
            s2 = (k - 1) / C
            worst = max(worst, rel(reml.I2, 100 * reml.tau2 / (reml.tau2 + s2)))
            n += 5
    return Check("E02", "DL tau2, Q, I2 (Q form for FE and DL rows; tau2/(tau2+s2) for REML)",
                 "pooling.fit", n, worst, worst < 1e-8)


@check
def e03_hedges_j_exact():
    """P07 proves the approximate J = 1 - 3/(4m - 1). The engine uses the exact gamma form,
    so check the engine against sympy's gamma form and bound its distance from the
    approximation."""
    worst, gap, n = 0.0, 0.0, 0
    for m in list(range(2, 200)) + [250, 500, 1000]:
        j = ef.hedges_j(m)
        worst = max(worst, rel(j, float(SYM_J_EXACT(m))))
        gap = max(gap, abs(j - float(SYM_J_APPROX(m))))
        n += 1
    ok = worst < 1e-10 and gap < 0.01 and 0 < ef.hedges_j(2) < 1
    return Check("E03", "Hedges J exact gamma form (engine) = sympy gamma form",
                 "effects.hedges_j", n, worst, ok,
                 f"max |exact - approx| = {gap:.5f}; P07's approximate form is within 0.01")


@check
def e04_change_sd():
    """P08: SD of change = sqrt(SDb^2 + SDf^2 - 2 r SDb SDf) (Handbook 6.5.2.8)."""
    rng = np.random.default_rng(SEED + 3)
    worst = 0.0
    for _ in range(TRIALS):
        sb, sf, r = rng.uniform(1, 30), rng.uniform(1, 30), rng.uniform(-0.95, 0.95)
        worst = max(worst, rel(cv.sd_change_imputed(sb, sf, r), float(SYM_CHANGE_SD(sb, sf, r))))
    return Check("E04", "imputed change SD", "conversions.sd_change_imputed", TRIALS, worst,
                 worst < RTOL)


@check
def e05_crossover_md():
    """P09: paired crossover variance (SDe^2 + SDc^2 - 2 r SDe SDc) / N (Handbook 23.2.7)."""
    rng = np.random.default_rng(SEED + 4)
    worst = 0.0
    for _ in range(TRIALS):
        se_, sc, r, n = rng.uniform(1, 20), rng.uniform(1, 20), rng.uniform(-0.9, 0.9), \
            int(rng.integers(6, 60))
        _, vi = ef.crossover_md(10.0, se_, 8.0, sc, n, r)
        worst = max(worst, rel(vi, float(SYM_PAIRED_VAR(se_, sc, r, n))))
    return Check("E05", "crossover paired MD variance", "effects.crossover_md", TRIALS, worst,
                 worst < RTOL)


@check
def e06_prediction_interval():
    """P06: PI = mu +/- c_PI sqrt(tau2 + SE^2) with the critical value of the configured
    pi_type, for DL and REML, with HKSJ off and under every hksj_rule. riley uses t(k-2),
    revman t(k-1), and metafor the CI's own critical value (t(k-1) under HKSJ, z otherwise:
    metafor predict()). Also checks that the engine's PI contains its CI in every fit (P06:
    same SE, c_PI >= c_CI), and that the PI is NaN for FE fits and for riley at k = 2."""
    from scipy import stats
    rng = np.random.default_rng(SEED + 5)
    worst, n, bad = 0.0, 0, []
    configs = [(False, "always")] + [(True, rule) for rule in pl.HKSJ_RULES]
    for trial in range(TRIALS // 2):
        k = int(rng.integers(3, 12))
        y, v = _draw(rng, k)
        if trial % 3 == 0:
            y = rng.normal(0.4, 0.01, k)       # low heterogeneity: HK s2 < 1
        for method in ("DL", "REML"):
            for hk, rule in configs:
                f = pl.fit(y, v, method=method, hksj=hk, hksj_rule=rule)
                c_ci = (float(f.ci_ub[0]) - f.estimate) / f.se[0]
                for pit, crit in (("riley", stats.t.ppf(0.975, k - 2)),
                                  ("revman", stats.t.ppf(0.975, k - 1)), ("metafor", c_ci)):
                    lb, ub = pl.prediction_interval(f, 0.95, pit)
                    half = crit * math.sqrt(f.tau2 + f.se[0] ** 2)
                    worst = max(worst, rel(lb, f.estimate - half), rel(ub, f.estimate + half))
                    if not (lb <= f.ci_lb[0] + 1e-12 and ub >= f.ci_ub[0] - 1e-12):
                        bad.append(f"{method}/{rule if hk else 'no HKSJ'}/{pit}: PI inside CI")
                    n += 2
    y, v = _draw(rng, 2)
    if not all(math.isnan(x) for x in pl.prediction_interval(pl.fit(y, v, method="REML"), 0.95, "riley")):
        bad.append("riley PI at k = 2 is not NaN")
    if not all(math.isnan(x) for x in pl.prediction_interval(pl.fit(*_draw(rng, 5), method="FE"), 0.95, "riley")):
        bad.append("FE fit has a PI")
    return Check("E06", "prediction interval (riley t(k-2), revman t(k-1), metafor = CI crit; "
                 "DL/REML x HKSJ off and every hksj_rule); PI contains CI; NaN for FE and riley k=2",
                 "pooling.prediction_interval", n, worst, worst < RTOL and not bad, "; ".join(bad[:3]))


# --- engine-level properties (unit of analysis) -------------------------------------------

def _rows(rng, n_arms, sid="SYNTHETIC-R", n_controls=1, design="parallel"):
    rows = [{"study_id": sid, "outcome": "SYNTHETIC_x", "timepoint": "8 weeks",
             "design": design, "arm": f"NR{j}", "is_control": 0,
             "n": int(rng.integers(8, 40)), "mean": float(rng.uniform(20, 60)),
             "sd": float(rng.uniform(4, 15)), "dose_mg_day": 100.0 * (j + 1),
             "notes": "SYNTHETIC"} for j in range(n_arms)]
    for j in range(n_controls):
        rows.append({"study_id": sid, "outcome": "SYNTHETIC_x", "timepoint": "8 weeks",
                     "design": design, "arm": f"PBO{j}", "is_control": 1,
                     "n": int(rng.integers(8, 40)), "mean": float(rng.uniform(25, 35)),
                     "sd": float(rng.uniform(6, 10)), "dose_mg_day": 0.0, "notes": "SYNTHETIC"})
    return rows


def _pool_arms(arms):
    """Handbook v6 Table 6.5.a, applied pairwise: (n, mean, sd) of the merged group."""
    n, m, sd = arms[0]
    for n2, m2, sd2 in arms[1:]:
        tot = n + n2
        var = ((n - 1) * sd**2 + (n2 - 1) * sd2**2 + n * n2 / tot * (m - m2) ** 2) / (tot - 1)
        n, m, sd = tot, (n * m + n2 * m2) / tot, math.sqrt(var)
    return n, m, sd


def _effect(measure, e, c):
    """Independent effect and variance from (n, mean, sd) of the two groups.
    MD: Handbook 6.5.1.2. SMD: Hedges g with exact J, metafor vtype LS (Hedges 1981).
    ROM: log ratio of means, delta-method variance (Hedges, Gurevitch & Curtis 1999)."""
    (n1, m1, s1), (n2, m2, s2) = e, c
    if measure == "MD":
        return m1 - m2, s1**2 / n1 + s2**2 / n2
    if measure == "SMD":
        df = n1 + n2 - 2
        g = float(SYM_J_EXACT(df)) * (m1 - m2) / math.sqrt(((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / df)
        return g, 1 / n1 + 1 / n2 + g**2 / (2 * (n1 + n2))
    return math.log(m1 / m2), s1**2 / (n1 * m1**2) + s2**2 / (n2 * m2**2)


def _tie(bad, tag, row, measure, e, c):
    yi, vi = _effect(measure, e, c)
    d = max(rel(row.yi, yi), abs(row.vi - vi) / vi)     # vi: purely relative (vi > 0)
    if d > 1e-9:
        bad.append(f"{tag}: yi/vi {row.yi:.6g}/{row.vi:.6g} != {yi:.6g}/{vi:.6g} "
                   f"from n_e={e[0]:g}, n_c={c[0]:g}")
    return d


@check
def e07_control_counted_once():
    """Handbook v6 23.3.4: a shared control must not enter twice. On every arm-level path the
    engine accepts (split, combine, several control rows, explicit combine_arms) this checks
    that the control participants enter once in total, that each comparison's effect AND
    variance are those of the independent formula at the split n (so a split n_c with a
    variance from the full control fails), and that each row's dose is its own
    experimental arm's (split: that arm's dose; combined: min and max of the merged arms).
    A multi-period crossover must either be excluded or count the placebo period once.
    Several reported between-group MDs against one control must be refused. Splitting
    still leaves the comparisons correlated, which no count can show (Handbook 23.3.4)."""
    rng = np.random.default_rng(SEED + 6)
    bad, n, worst = [], 0, 0.0
    for trial in range(60):
        n_ctl = 2 if trial % 4 == 3 else 1
        rows = _rows(rng, int(rng.integers(2, 5)), n_controls=n_ctl)
        exps = [r for r in rows if not r["is_control"]]
        ctl = _pool_arms([(r["n"], r["mean"], r["sd"]) for r in rows if r["is_control"]])
        df = ef.coerce_extraction(pd.DataFrame(rows))
        k = len(exps)
        for measure in ("MD", "SMD", "ROM"):
            for mode in ("split", "combine"):
                eff, excl = ef.study_effects(df, ef.EffectSettings(measure=measure, multi_arm=mode))
                n += 1
                tag = f"{mode}/{measure}/{n_ctl} control rows"
                if len(excl):
                    bad.append(f"{tag}: excluded ({excl.reason.iloc[0][:60]})")
                    continue
                if abs(eff.n_c.sum() - ctl[0]) > 1e-9:
                    bad.append(f"{tag}: control n {eff.n_c.sum()} != {ctl[0]}")
                if mode == "split":
                    if len(eff) != k:
                        bad.append(f"{tag}: {len(eff)} rows for {k} arms")
                    for _, row in eff.iterrows():
                        arm = next((r for r in exps if r["arm"] == row.arm_e), None)
                        if arm is None:
                            bad.append(f"{tag}: arm_e {row.arm_e} is not an experimental arm")
                            continue
                        if not (row.dose_min == row.dose_max == arm["dose_mg_day"]):
                            bad.append(f"{tag}: {row.arm_e} dose {row.dose_min}-{row.dose_max} "
                                       f"!= its arm's {arm['dose_mg_day']}")
                        worst = max(worst, _tie(bad, tag, row, measure,
                                                (arm["n"], arm["mean"], arm["sd"]),
                                                (ctl[0] / k, ctl[1], ctl[2])))
                else:
                    if len(eff) != 1:
                        bad.append(f"{tag}: {len(eff)} rows")
                        continue
                    row = eff.iloc[0]
                    doses = [r["dose_mg_day"] for r in exps]
                    if (row.dose_min, row.dose_max) != (min(doses), max(doses)):
                        bad.append(f"{tag}: dose {row.dose_min}-{row.dose_max} != "
                                   f"{min(doses)}-{max(doses)}")
                    worst = max(worst, _tie(bad, tag, row, measure,
                                            _pool_arms([(r["n"], r["mean"], r["sd"]) for r in exps]),
                                            ctl))
        # explicit combine_arms of the first two NR arms, then split against the rest
        merged = [r["arm"] for r in exps[:2]]
        s = ef.EffectSettings(measure="MD", multi_arm="split",
                              combine={"SYNTHETIC-R": {"NRmerged": merged}})
        eff, excl = ef.study_effects(df, s)
        n += 1
        rest = [r for r in exps if r["arm"] not in merged]
        if len(excl) or len(eff) != 1 + len(rest) or abs(eff.n_c.sum() - ctl[0]) > 1e-9:
            bad.append(f"combine_arms: {len(eff)} rows, control n {eff.n_c.sum()} != {ctl[0]}")
        else:
            m_row = eff[eff.arm_e == "NRmerged"].iloc[0]
            worst = max(worst, _tie(bad, "combine_arms", m_row, "MD",
                                    _pool_arms([(r["n"], r["mean"], r["sd"]) for r in exps[:2]]),
                                    (ctl[0] / (1 + len(rest)), ctl[1], ctl[2])))
        # multi-period crossover: two NR periods and one placebo period in the same people
        xo = _rows(rng, 2, sid="SYNTHETIC-X", design="crossover")
        for r in xo:
            r["n"] = xo[0]["n"]
        eff, excl = ef.study_effects(ef.coerce_extraction(pd.DataFrame(xo)), ef.EffectSettings())
        n += 1
        if len(eff) and eff.n_c.sum() - xo[0]["n"] > 1e-9:
            bad.append(f"multi-period crossover: placebo period counted {eff.n_c.sum()} > {xo[0]['n']}")
        if not len(eff) and not len(excl):
            bad.append("multi-period crossover vanished without an exclusion row")
        paired = [dict(r, paired_md=1.0, paired_se=0.5) for r in rows if not r["is_control"]]
        eff, excl = ef.study_effects(ef.coerce_extraction(pd.DataFrame(paired)),
                                     ef.EffectSettings())
        n += 1
        if len(eff):
            bad.append("several paired-only MDs from one study were accepted")
    return Check("E07", "control participants counted once on every arm-level path; yi and vi "
                 "equal the independent formula at the split n; dose from the comparison's own "
                 "experimental arm(s)", "effects.study_effects", n, worst, not bad,
                 "; ".join(bad[:3]))


@check
def e08_one_timepoint_per_study():
    """Handbook 6.2.4: results from more than one timepoint of a study cannot enter one
    meta-analysis. Every timepoint rule keeps at most one label per study, and a study it
    cannot place is listed as excluded, not dropped silently."""
    rng = np.random.default_rng(SEED + 7)
    labels = ["wk2", "wk4", "day 21", "8 weeks", "3 months", "wk12", "end"]
    rows = []
    for s in range(30):
        for tp in rng.choice(labels, size=int(rng.integers(1, 4)), replace=False):
            for arm, ctl in (("NR", 0), ("PBO", 1)):
                rows.append({"study_id": f"SYNTHETIC-{s}", "outcome": "SYNTHETIC_x",
                             "timepoint": str(tp), "design": "parallel", "arm": arm,
                             "is_control": ctl, "n": 20, "mean": 10.0, "sd": 3.0})
    df = ef.coerce_extraction(pd.DataFrame(rows))
    bad, n = [], 0
    for rule in ("last", "first", ["wk4", "wk12"], "wk12", "8"):
        kept, _, excl = select_timepoint(df, rule)
        n += 1
        if (kept.groupby("study_id").timepoint.nunique() > 1).any():
            bad.append(f"{rule}: a study kept two timepoints")
        if set(kept.study_id) | {e["study_id"] for e in excl} != set(df.study_id):
            bad.append(f"{rule}: a study vanished without an exclusion row")
    return Check("E08", "one timepoint per study; unplaced studies listed in excluded",
                 "pipeline.select_timepoint", n, 0.0, not bad, "; ".join(bad))


def run():
    return [c() for c in CHECKS]


def receipt(results) -> str:
    lines = ["engine conformance receipt for validation/proofs/engine_conformance.py",
             f"date_utc: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
             f"python: {sys.version.split()[0]}  sympy: {sp.__version__}  numpy: {np.__version__}",
             "Inputs: random SYNTHETIC values from a fixed seed; not trial data.", "",
             "id | property | engine function | comparisons | max scaled diff | result"]
    for c in results:
        lines.append(f"{c.cid} | {c.prop} | {c.engine} | {c.n} | {c.max_rel:.2e} | "
                     f"{'PASS' if c.ok else 'FAIL'}" + (f" | {c.detail}" if c.detail else ""))
    lines += ["", f"overall: {'PASS' if all(c.ok for c in results) else 'FAIL'}"]
    return "\n".join(lines) + "\n"


def main() -> int:
    res = run()
    text = receipt(res)
    (HERE / "receipt_engine.txt").write_text(text)
    print(text)
    return 0 if all(c.ok for c in res) else 1


if __name__ == "__main__":
    raise SystemExit(main())
