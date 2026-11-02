"""Symbolic proofs of the properties the meta-analysis formulas must have.

Run from the repo root:
    .venv/bin/python validation/proofs/sympy_proofs.py

Writes receipt_sympy.txt beside this file and exits 0 only when every proof
step holds and every mutant is refuted.

How a step counts as proven
---------------------------
A step passes only when sympy *decides* it: an identity simplifies to 0, or
sympy's assumption system returns True for is_positive / is_nonnegative /
is_negative given the declared symbol assumptions (variances positive,
tau^2 >= 0, and so on). sympy returning None ("cannot decide") is a failure.

Steps tagged NUMERIC are scans, not proofs. They are used where the claim
depends on a special function sympy cannot reason about (t and normal
quantiles, the exact gamma-function form of Hedges' J). The receipt labels
them so.

Steps tagged k=2..8 are proven separately for each number of studies from 2
to 8. Steps tagged "general k" use sympy Sum over indexed symbols with sign
assumptions (w_i > 0, d_i >= 0) and hold for every number of studies k >= 1.
The Bend2 laws in bend/ are integer (Nat) models, not proofs of these
real-valued statements; RECEIPTS.txt gives the exact domain of every claim.

Mutants
-------
Each mutant is a deliberately wrong version of a formula. It is refuted only
when (a) sympy cannot prove the property for it and (b) a concrete SYNTHETIC
counterexample is found and printed. A mutant that survives is a failure of
this file: it would mean the checks cannot tell the right formula from a
wrong one.

Sources for the formulas
------------------------
- Inverse-variance and random-effects pooling, tau^2, I^2, prediction
  intervals: Cochrane Handbook for Systematic Reviews of Interventions v6,
  Chapter 10 (sections 10.3 and 10.10).
- DerSimonian R, Laird N. Meta-analysis in clinical trials. Controlled
  Clinical Trials 1986;7:177-188.
- Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis.
  Statistics in Medicine 2002;21:1539-1558.
- Higgins JPT, Thompson SG, Spiegelhalter DJ. A re-evaluation of
  random-effects meta-analysis. J R Stat Soc A 2009;172:137-159.
- Hedges LV. Distribution theory for Glass's estimator of effect size and
  related estimators. J Educational Statistics 1981;6:107-128.
- Change-from-baseline SD: Cochrane Handbook v6, Chapter 6, section 6.5.2.8.
- Cross-over trials: Cochrane Handbook v6, Chapter 23, section 23.2; Elbourne
  DR et al. Meta-analyses involving cross-over trials: methodological issues.
  Int J Epidemiol 2002;31:140-149.
- Shared control groups: Cochrane Handbook v6, Chapter 23, section 23.3.4.
"""

from __future__ import annotations

import random
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import sympy as sp

HERE = Path(__file__).resolve().parent
KS = range(2, 9)


# ---------------------------------------------------------------------------
# deciding helpers
# ---------------------------------------------------------------------------

def _forms(expr):
    expr = sp.sympify(expr)
    yield expr
    yield sp.together(expr)
    yield sp.factor(sp.together(expr))
    yield sp.simplify(expr)


def is_zero(expr) -> bool:
    return any(f == 0 for f in _forms(expr))


def decide(expr, prop: str) -> bool:
    """True only if sympy proves expr has prop under the symbol assumptions."""
    return any(getattr(f, prop) is True for f in _forms(expr))


@dataclass
class Result:
    pid: str
    prop: str
    method: str
    steps: list[tuple[str, bool]] = field(default_factory=list)

    def step(self, text: str, ok: bool) -> None:
        self.steps.append((text, bool(ok)))

    @property
    def ok(self) -> bool:
        return bool(self.steps) and all(ok for _, ok in self.steps)


@dataclass
class Mutant:
    mid: str
    target: str
    wrong: str
    proven: bool          # sympy proved the property for the mutant (bad)
    counterexample: str   # SYNTHETIC values where the mutant breaks the property

    @property
    def refuted(self) -> bool:
        return (not self.proven) and bool(self.counterexample)


PROOFS: list[Callable[[], Result]] = []
MUTANTS: list[Callable[[], Mutant]] = []


def proof(fn):
    PROOFS.append(fn)
    return fn


def mutant(fn):
    MUTANTS.append(fn)
    return fn


def search(pred, sampler, tries=20000, seed=20261101):
    """Return the first SYNTHETIC sample where pred is False, else ''."""
    rng = random.Random(seed)
    for _ in range(tries):
        vals = sampler(rng)
        if not pred(vals):
            return "SYNTHETIC " + ", ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in vals.items())
    return ""


# ---------------------------------------------------------------------------
# the formulas under test
# ---------------------------------------------------------------------------

def fe_weights(v, tau2):
    return [1 / vi for vi in v]


def re_weights(v, tau2):
    return [1 / (vi + tau2) for vi in v]


def pooled(w, y):
    return sum(wi * yi for wi, yi in zip(w, y)) / sum(w)


def dl_C(w):
    S1 = sum(w)
    S2 = sum(wi**2 for wi in w)
    return S1 - S2 / S1


def tau2_dl(Q, df, C):
    return sp.Max(0, (Q - df) / C)


def i2(Q, df):
    return sp.Max(0, (Q - df) / Q)


def hedges_j(df):
    return 1 - sp.Integer(3) / (4 * df - 1)


def change_sd_sq(sd_base, sd_final, r):
    return sd_base**2 + sd_final**2 - 2 * r * sd_base * sd_final


def paired_se_sq(s1, s2, r, n):
    return (s1**2 + s2**2 - 2 * r * s1 * s2) / n


def parallel_se_sq(s1, s2, n):
    return s1**2 / n + s2**2 / n


# ---------------------------------------------------------------------------
# P1 inverse-variance weights
# ---------------------------------------------------------------------------

def weights_positive_and_normalize(weight_fn, k):
    v = sp.symbols(f"v1:{k + 1}", positive=True)
    tau2 = sp.Symbol("tau2", nonnegative=True)
    w = weight_fn(v, tau2)
    W = sum(w)
    pos = all(decide(wi, "is_positive") for wi in w)
    norm = is_zero(sum(wi / W for wi in w) - 1)
    share = all(decide(wi / W, "is_positive") and decide(1 - wi / W, "is_positive") for wi in w)
    return pos, norm, share


@proof
def p01_weights():
    r = Result("P01", "inverse-variance weights are positive and normalize to 1 (fixed and random effects)",
               "sympy assumptions, k=2..8; general-k sum identity")
    for name, fn in (("FE w=1/v", fe_weights), ("RE w=1/(v+tau2)", re_weights)):
        for k in KS:
            pos, norm, share = weights_positive_and_normalize(fn, k)
            r.step(f"{name} k={k}: every w_i > 0", pos)
            r.step(f"{name} k={k}: sum_i w_i/W == 1", norm)
            r.step(f"{name} k={k}: 0 < w_i/W < 1", share)
    i, k = sp.symbols("i k", integer=True, positive=True)
    w = sp.IndexedBase("w")
    Wsym = sp.Symbol("W", positive=True)
    s = sp.Sum(w[i] / Wsym, (i, 1, k)).doit().subs(Wsym, sp.Sum(w[i], (i, 1, k)))
    r.step("general k: Sum_i w_i / Sum_j w_j == 1", sp.simplify(s - 1) == 0)
    return r


@mutant
def m01_normalize_by_se():
    v = sp.symbols("v1:3", positive=True)
    w = [1 / vi for vi in v]
    norm_by = sum(1 / sp.sqrt(vi) for vi in v)
    proven = is_zero(sum(wi / norm_by for wi in w) - 1)

    def pred(x):
        ww = [1 / x["v1"], 1 / x["v2"]]
        return abs(sum(ww) / (x["v1"] ** -0.5 + x["v2"] ** -0.5) - 1) < 1e-9

    cx = search(pred, lambda g: {"v1": g.uniform(0.01, 2), "v2": g.uniform(0.01, 2)})
    return Mutant("M01", "P01", "weights divided by sum of 1/SE instead of sum of 1/variance", proven, cx)


@mutant
def m02_re_weight_minus_tau2():
    v = sp.symbols("v1:3", positive=True)
    tau2 = sp.Symbol("tau2", nonnegative=True)
    w = [1 / (vi - tau2) for vi in v]
    proven = all(decide(wi, "is_positive") for wi in w)
    cx = search(lambda x: 1 / (x["v1"] - x["tau2"]) > 0,
                lambda g: {"v1": g.uniform(0.01, 1), "tau2": g.uniform(0, 1)})
    return Mutant("M02", "P01", "random-effects weight 1/(v - tau2)", proven, cx)


# ---------------------------------------------------------------------------
# P2 pooled estimate lies between the smallest and largest study estimate
# ---------------------------------------------------------------------------

def between_bounds(weight_fn, k):
    v = sp.symbols(f"v1:{k + 1}", positive=True)
    tau2 = sp.Symbol("tau2", nonnegative=True)
    w = weight_fn(v, tau2)
    lo, hi = sp.symbols("lo hi", real=True)
    d = sp.symbols(f"d1:{k + 1}", nonnegative=True)
    above_lo = decide(pooled(w, [lo + di for di in d]) - lo, "is_nonnegative")
    below_hi = decide(hi - pooled(w, [hi - di for di in d]), "is_nonnegative")
    return above_lo, below_hi


@proof
def p02_between():
    r = Result("P02", "fixed-effect and random-effects pooled estimates lie in [min y_i, max y_i]",
               "sympy assumptions, k=2..8 and general k (Sum over indexed v_i > 0, d_i >= 0)")
    for name, fn in (("FE", fe_weights), ("RE", re_weights)):
        for k in KS:
            a, b = between_bounds(fn, k)
            r.step(f"{name} k={k}: pooled - min >= 0 (y_i = min + d_i, d_i >= 0)", a)
            r.step(f"{name} k={k}: max - pooled >= 0 (y_i = max - d_i, d_i >= 0)", b)
    i, k = sp.symbols("i k", integer=True, positive=True)
    w, y = sp.IndexedBase("w"), sp.IndexedBase("y")
    m = sp.Symbol("m", real=True)
    ident = sp.Sum(w[i] * (y[i] - m), (i, 1, k)) - (sp.Sum(w[i] * y[i], (i, 1, k)) - m * sp.Sum(w[i], (i, 1, k)))
    r.step("general k: sum w_i y_i - m sum w_i == sum w_i (y_i - m)", sp.simplify(ident) == 0)
    for name, gap_lo, gap_hi in general_k_between():
        r.step(f"general k, {name}: y_i = min + d_i, d_i >= 0: pooled - min == sum w d / sum w >= 0", gap_lo)
        r.step(f"general k, {name}: y_i = max - d_i, d_i >= 0: max - pooled == sum w d / sum w >= 0", gap_hi)
    return r


def _sum_expand(expr):
    return expr.replace(sp.Sum, lambda *a: sp.Sum(*a).expand())


def _sum_identity(expr) -> bool:
    """True only if sympy simplifies the Sum expression to 0; an error is 'cannot decide'."""
    try:
        return sp.simplify(_sum_expand(expr)) == 0
    except Exception:
        return False


def general_k_between(pool=None):
    """For every k: write each y_i as min + d_i (or max - d_i) with d_i >= 0. Then
    pooled - min (max - pooled) is identically sum w_i d_i / sum w_i, which sympy shows is
    >= 0 from w_i > 0 and d_i >= 0. Every finite set of reals has such a min and max."""
    pool = pool or (lambda w, y, i, k: sp.Sum(w * y, (i, 1, k)) / sp.Sum(w, (i, 1, k)))
    i, k = sp.symbols("i k", integer=True, positive=True)
    v = sp.IndexedBase("v", positive=True)
    d = sp.IndexedBase("d", nonnegative=True)
    t = sp.Symbol("tau2", nonnegative=True)
    lo, hi = sp.symbols("lo hi", real=True)
    out = []
    for name, w in (("FE w = 1/v", 1 / v[i]), ("RE w = 1/(v + tau2)", 1 / (v[i] + t))):
        gap = sp.Sum(w * d[i], (i, 1, k)) / sp.Sum(w, (i, 1, k))
        lo_ok = _sum_identity(pool(w, lo + d[i], i, k) - lo - gap) and gap.is_nonnegative is True
        hi_ok = _sum_identity(hi - pool(w, hi - d[i], i, k) - gap) and gap.is_nonnegative is True
        out.append((name, lo_ok, hi_ok))
    return out


@mutant
def m03_pooled_without_division():
    k = 2
    v = sp.symbols(f"v1:{k + 1}", positive=True)
    w = fe_weights(v, 0)
    hi = sp.Symbol("hi", real=True)
    d = sp.symbols(f"d1:{k + 1}", nonnegative=True)
    wrong = sum(wi * (hi - di) for wi, di in zip(w, d))
    proven = decide(hi - wrong, "is_nonnegative")

    def pred(x):
        ws = [1 / x["v1"], 1 / x["v2"]]
        ys = [x["y1"], x["y2"]]
        return sum(a * b for a, b in zip(ws, ys)) <= max(ys) + 1e-12

    cx = search(pred, lambda g: {"v1": g.uniform(0.01, 1), "v2": g.uniform(0.01, 1),
                                  "y1": g.uniform(0.1, 1), "y2": g.uniform(0.1, 1)})
    return Mutant("M03", "P02", "pooled = sum(w*y) without dividing by sum(w)", proven, cx)


# ---------------------------------------------------------------------------
# P3 DerSimonian-Laird tau^2
# ---------------------------------------------------------------------------

@proof
def p03_dl_tau2():
    r = Result("P03", "DL tau^2 >= 0 after truncation; tau^2 == 0 when Q <= k-1; C > 0; Q >= 0",
               "sympy identities and assumptions, k=2..8 for C and Q")
    for k in KS:
        w = sp.symbols(f"w1:{k + 1}", positive=True)
        y = sp.symbols(f"y1:{k + 1}", real=True)
        S1, S2 = sum(w), sum(wi**2 for wi in w)
        cross = 2 * sum(w[a] * w[b] for a in range(k) for b in range(a + 1, k))
        r.step(f"k={k}: S1^2 - S2 == 2 sum_(i<j) w_i w_j", is_zero(sp.expand(S1**2 - S2 - cross)))
        r.step(f"k={k}: C = S1 - S2/S1 > 0", decide(dl_C(w), "is_positive"))
        th = pooled(w, y)
        Qdef = sum(wi * (yi - th) ** 2 for wi, yi in zip(w, y))
        Qcomp = sum(wi * yi**2 for wi, yi in zip(w, y)) - sum(wi * yi for wi, yi in zip(w, y)) ** 2 / S1
        r.step(f"k={k}: Q = sum w(y - pooled)^2 is >= 0", decide(Qdef, "is_nonnegative"))
        r.step(f"k={k}: Q definition == computing form sum w y^2 - (sum w y)^2/S1", is_zero(Qdef - Qcomp))
    Q = sp.Symbol("Q", nonnegative=True)
    C = sp.Symbol("C", positive=True)
    df = sp.Symbol("df", positive=True)
    s = sp.Symbol("s", nonnegative=True)
    sp_ = sp.Symbol("s", positive=True)
    r.step("tau2 = max(0, (Q - df)/C) >= 0 for every Q", decide(tau2_dl(Q, df, C), "is_nonnegative"))
    r.step("Q = df - s, s >= 0  =>  tau2 == 0", is_zero(tau2_dl(df - s, df, C)))
    r.step("Q = df + s, s > 0  =>  tau2 == s/C (truncation inactive)", is_zero(tau2_dl(df + sp_, df, C) - sp_ / C))
    return r


@mutant
def m04_untruncated_tau2():
    Q = sp.Symbol("Q", nonnegative=True)
    C = sp.Symbol("C", positive=True)
    df = sp.Symbol("df", positive=True)
    proven = decide((Q - df) / C, "is_nonnegative")
    cx = search(lambda x: (x["Q"] - x["df"]) / x["C"] >= 0,
                lambda g: {"Q": g.uniform(0, 5), "df": g.randint(1, 6), "C": g.uniform(0.1, 10)})
    return Mutant("M04", "P03", "tau2 = (Q - df)/C with no max(0, .)", proven, cx)


@mutant
def m05_C_without_division():
    w = sp.symbols("w1:3", positive=True)
    wrong_C = sum(w) - sum(wi**2 for wi in w)
    proven = decide(wrong_C, "is_positive")
    cx = search(lambda x: (x["w1"] + x["w2"]) - (x["w1"] ** 2 + x["w2"] ** 2) > 0,
                lambda g: {"w1": g.uniform(0.1, 50), "w2": g.uniform(0.1, 50)})
    return Mutant("M05", "P03", "C = S1 - S2 (missing division by S1)", proven, cx)


# ---------------------------------------------------------------------------
# P4 I^2
# ---------------------------------------------------------------------------

@proof
def p04_i2():
    r = Result("P04", "I^2 = max(0,(Q-df)/Q) lies in [0,1); equals (H^2-1)/H^2 and, at the DL tau2, tau2/(tau2+s^2) (the form the engine reports for random-effects rows; engine match: E02)",
               "sympy identities and assumptions")
    Q = sp.Symbol("Q", positive=True)
    df = sp.Symbol("df", positive=True)
    s = sp.Symbol("s", nonnegative=True)
    sp_ = sp.Symbol("s", positive=True)
    C = sp.Symbol("C", positive=True)
    r.step("I2 >= 0", decide(i2(Q, df), "is_nonnegative"))
    r.step("argument (Q-df)/Q < 1 since (Q-df)/Q - 1 = -df/Q < 0", decide((Q - df) / Q - 1, "is_negative"))
    branches = [e for e, _ in i2(Q, df).rewrite(sp.Piecewise).args]
    r.step("I2 < 1: each branch of I2 written as a Piecewise (0, and (Q-df)/Q) is < 1",
           len(branches) == 2 and all(decide(e - 1, "is_negative") for e in branches))
    H2 = Q / df
    r.step("(H^2 - 1)/H^2 == (Q - df)/Q with H^2 = Q/df", is_zero((H2 - 1) / H2 - (Q - df) / Q))
    tau2 = sp_ / C          # DL tau2 when Q = df + s
    s2 = df / C             # Higgins-Thompson typical within-study variance
    r.step("Q = df + s: tau2/(tau2 + s^2) == (Q - df)/Q, s^2 = df/C",
           is_zero(tau2 / (tau2 + s2) - ((df + sp_) - df) / (df + sp_)))
    for k in KS:
        w = sp.symbols(f"w1:{k + 1}", positive=True)
        S1, S2 = sum(w), sum(wi**2 for wi in w)
        r.step(f"k={k}: H&T typical variance (k-1) S1/(S1^2 - S2) == df/C",
               is_zero((k - 1) * S1 / (S1**2 - S2) - (k - 1) / dl_C(w)))
    r.step("df = Q + s, s >= 0 (Q <= df)  =>  I2 == 0", is_zero(i2(Q, Q + s)))
    return r


@mutant
def m06_i2_over_df():
    Q = sp.Symbol("Q", positive=True)
    df = sp.Symbol("df", positive=True)
    proven = decide((Q - df) / df - 1, "is_negative")
    cx = search(lambda x: max(0.0, (x["Q"] - x["df"]) / x["df"]) < 1,
                lambda g: {"Q": g.uniform(0.1, 40), "df": g.randint(1, 10)})
    return Mutant("M06", "P04", "I2 = max(0, (Q - df)/df) (wrong denominator)", proven, cx)


@mutant
def m07_i2_untruncated():
    Q = sp.Symbol("Q", positive=True)
    df = sp.Symbol("df", positive=True)
    proven = decide((Q - df) / Q, "is_nonnegative")
    cx = search(lambda x: (x["Q"] - x["df"]) / x["Q"] >= 0,
                lambda g: {"Q": g.uniform(0.1, 10), "df": g.randint(1, 10)})
    return Mutant("M07", "P04", "I2 = (Q - df)/Q with no max(0, .)", proven, cx)


# ---------------------------------------------------------------------------
# P5 random-effects weights become more equal as tau^2 grows
# ---------------------------------------------------------------------------

def ratio_shrinks(weight_fn):
    vi = sp.Symbol("v_i", positive=True)
    delta = sp.Symbol("delta", positive=True)
    t = sp.Symbol("tau2", nonnegative=True)
    wi, wj = weight_fn([vi, vi + delta], t)
    ratio = sp.simplify(wi / wj)
    return ratio, t, sp.diff(ratio, t)


@proof
def p05_re_weights_equalize():
    r = Result("P05", "random-effects weights become more equal as tau^2 grows",
               "sympy derivative signs and limits")
    ratio, t, dr = ratio_shrinks(re_weights)
    r.step("v_i < v_j: weight ratio w_i/w_j = (v_j+tau2)/(v_i+tau2) > 1", decide(ratio - 1, "is_positive"))
    r.step("d(w_i/w_j)/d tau2 < 0 (ratio falls toward 1)", decide(dr, "is_negative"))
    r.step("lim tau2->oo w_i/w_j == 1", sp.limit(ratio, t, sp.oo) == 1)
    v1 = sp.Symbol("v1", positive=True)
    delta = sp.Symbol("delta", positive=True)
    w = re_weights([v1, v1 + delta], t)
    share = w[0] / sum(w)
    r.step("k=2: normalized weight of the more precise study falls as tau2 grows",
           decide(sp.diff(share, t), "is_negative"))
    v = sp.symbols("v1:4", positive=True)
    w3 = re_weights(v, t)
    r.step("k=3: every normalized weight tends to 1/3 as tau2 -> oo",
           all(sp.limit(wi / sum(w3), t, sp.oo) == sp.Rational(1, 3) for wi in w3))
    return r


@mutant
def m08_multiplicative_tau2():
    def wrong(v, t):
        return [1 / (vi * (1 + t)) for vi in v]
    ratio, t, dr = ratio_shrinks(wrong)
    proven = decide(dr, "is_negative")
    vals = {sp.Symbol("v_i", positive=True): sp.Rational(1, 10), sp.Symbol("delta", positive=True): sp.Rational(3, 10)}
    r0, r10 = ratio.subs(vals).subs(t, 0), ratio.subs(vals).subs(t, 10)
    cx = ""
    if is_zero(dr) and r0 == r10:
        cx = f"SYNTHETIC v_i=0.1, v_j=0.4: ratio {float(r0):.4g} at tau2=0 and {float(r10):.4g} at tau2=10"
    return Mutant("M08", "P05", "weight 1/(v*(1+tau2)): ratio does not move with tau2", proven, cx)


# ---------------------------------------------------------------------------
# P6 prediction interval contains the confidence interval
# ---------------------------------------------------------------------------

def quantile_scan():
    from scipy import stats
    z = stats.norm.ppf(0.975)
    bad = [df for df in list(range(1, 5001)) + [10**5, 10**7] if not stats.t.ppf(0.975, df) > z]
    mono = all(stats.t.ppf(0.975, d) > stats.t.ppf(0.975, d + 1) for d in range(1, 2000))
    return bad, mono


@proof
def p06_pi_contains_ci():
    r = Result("P06", "prediction interval mu +/- c_PI sqrt(tau2 + SE^2) contains the CI mu +/- c_CI SE when both use the same SE (Wald, or Hartung-Knapp-scaled under HKSJ) and c_PI >= c_CI (engine pairs: riley t_(k-2), revman t_(k-1), metafor = c_CI; CI t_(k-1) or z; engine check E06)",
               "sympy for the half-widths; NUMERIC scan for the quantile order")
    se = sp.Symbol("SE", positive=True)
    t = sp.Symbol("tau2", nonnegative=True)
    root = sp.sqrt(t + se**2)
    r.step("sqrt(tau2 + SE^2) is real and >= 0", decide(root, "is_nonnegative"))
    r.step("sqrt(tau2 + SE^2)^2 - SE^2 == tau2 >= 0, so sqrt(tau2 + SE^2) >= SE",
           is_zero(root**2 - se**2 - t) and decide(t, "is_nonnegative"))
    z, g, h = sp.symbols("z g h", nonnegative=True)
    s = sp.Symbol("s", nonnegative=True)
    r.step("c_PI = c_CI + g, sqrt(tau2 + SE^2) = SE + h with g, h >= 0: c_PI*(SE + h) - c_CI*SE >= 0",
           decide(sp.expand((z + g) * (s + h) - z * s), "is_nonnegative"))
    mu = sp.Symbol("mu", real=True)
    a, c = sp.symbols("a c", nonnegative=True)     # CI half-width a, PI half-width a + c
    r.step("CI half-width a, PI half-width a + c (c >= 0): PI lower <= CI lower and CI upper <= PI upper",
           decide((mu - a) - (mu - (a + c)), "is_nonnegative") and decide((mu + (a + c)) - (mu + a), "is_nonnegative"))
    bad, mono = quantile_scan()
    r.step("NUMERIC: t_(0.975, df) > z_0.975 for df = 1..5000, 1e5, 1e7", not bad)
    r.step("NUMERIC: t_(0.975, df) decreasing in df for df = 1..2000 (so t_(k-2) > t_(k-1))", mono)
    r.step("NUMERIC: engine pairs c_PI >= c_CI for k = 3..2000: riley t_(k-2) and revman t_(k-1) vs CI "
           "t_(k-1) (HKSJ) or z; metafor uses c_CI itself (g = 0)", engine_crit_pairs_ok())
    return r


def engine_crit_pairs_ok():
    from scipy import stats
    z = stats.norm.ppf(0.975)
    for k in range(3, 2001):
        ci = (stats.t.ppf(0.975, k - 1), z)
        for c_pi in (stats.t.ppf(0.975, k - 2), stats.t.ppf(0.975, k - 1)):
            if not all(c_pi >= c for c in ci):
                return False
    return True


@mutant
def m09_pi_without_se():
    se = sp.Symbol("SE", positive=True)
    t = sp.Symbol("tau2", nonnegative=True)
    z = sp.Symbol("z", positive=True)
    proven = decide(z * sp.sqrt(t) - z * se, "is_nonnegative")
    cx = search(lambda x: 1.96 * x["tau2"] ** 0.5 >= 1.96 * x["SE"],
                lambda g: {"tau2": g.uniform(0, 0.05), "SE": g.uniform(0.05, 0.5)})
    return Mutant("M09", "P06", "PI half-width z*sqrt(tau2), SE^2 left out", proven, cx)


# ---------------------------------------------------------------------------
# P7 Hedges' correction J
# ---------------------------------------------------------------------------

def exact_j_scan():
    import numpy as np
    from scipy.special import gammaln
    df = np.arange(2, 1001, dtype=float)
    j = np.exp(gammaln(df / 2) - gammaln((df - 1) / 2)) / np.sqrt(df / 2)
    approx = 1 - 3 / (4 * df - 1)
    return bool(np.all((j > 0) & (j < 1))), bool(np.all(np.diff(j) > 0)), float(np.max(np.abs(j - approx)))


@proof
def p07_hedges_j():
    r = Result("P07", "approximate Hedges J = 1 - 3/(4 df - 1) lies in (0,1) for df > 1 and increases with df (the engine uses the exact gamma form: NUMERIC steps here, engine match E03)",
               "sympy assumptions; NUMERIC scan for the exact gamma form")
    e = sp.Symbol("e", positive=True)          # df = 1 + e
    df = sp.Symbol("df", positive=True)
    J = hedges_j(1 + e)
    r.step("df = 1 + e, e > 0: J > 0", decide(J, "is_positive"))
    r.step("df = 1 + e, e > 0: 1 - J > 0", decide(1 - J, "is_positive"))
    r.step("dJ/d df = 12/(4 df - 1)^2 > 0 for df > 1", decide(sp.diff(hedges_j(df), df).subs(df, 1 + e), "is_positive"))
    r.step("lim df->oo J == 1", sp.limit(hedges_j(df), df, sp.oo) == 1)
    r.step("boundary: J(df = 1) == 0, so (0,1) needs df > 1 (two-arm trial: n1 + n2 >= 4)", hedges_j(sp.Integer(1)) == 0)
    inside, increasing, gap = exact_j_scan()
    r.step("NUMERIC: exact J = G(df/2)/(sqrt(df/2) G((df-1)/2)) in (0,1) for df = 2..1000", inside)
    r.step("NUMERIC: exact J increasing for df = 2..1000", increasing)
    r.step(f"NUMERIC: max |exact - approx| over df = 2..1000 is {gap:.5f} (< 0.01)", gap < 0.01)
    return r


@mutant
def m10_j_sign_flip():
    e = sp.Symbol("e", positive=True)
    wrong = 1 + sp.Integer(3) / (4 * (1 + e) - 1)
    proven = decide(1 - wrong, "is_positive")
    cx = search(lambda x: 1 + 3 / (4 * x["df"] - 1) < 1, lambda g: {"df": g.randint(2, 200)})
    return Mutant("M10", "P07", "J = 1 + 3/(4 df - 1) (sign flipped)", proven, cx)


@mutant
def m11_j_domain_df1():
    proven = decide(hedges_j(sp.Integer(1)), "is_positive")
    cx = "SYNTHETIC df=1: J = 1 - 3/3 = 0, not in (0,1)" if hedges_j(sp.Integer(1)) == 0 else ""
    return Mutant("M11", "P07", "claim J in (0,1) for df >= 1 (boundary df = 1 included)", proven, cx)


# ---------------------------------------------------------------------------
# P8 change-from-baseline SD
# ---------------------------------------------------------------------------

@proof
def p08_change_sd():
    r = Result("P08", "SD_change = sqrt(SDb^2 + SDf^2 - 2 r SDb SDf) is a valid SD for r in [-1,1] and falls as r rises",
               "sympy identities and assumptions")
    a, b = sp.symbols("SDb SDf", positive=True)
    s, u = sp.symbols("s u", nonnegative=True)      # r = 1 - s and r = -1 + u
    rr = sp.Symbol("r", real=True)
    R = change_sd_sq(a, b, rr)
    r.step("radicand - (SDb - SDf)^2 == 2 SDb SDf (1 - r)", is_zero(R - (a - b) ** 2 - 2 * a * b * (1 - rr)))
    r.step("r <= 1 (r = 1 - s): radicand - (SDb - SDf)^2 >= 0", decide(change_sd_sq(a, b, 1 - s) - (a - b) ** 2, "is_nonnegative"))
    r.step("so radicand >= (SDb - SDf)^2 >= 0 and the SD is real", decide((a - b) ** 2, "is_nonnegative"))
    r.step("r >= -1 (r = -1 + u): (SDb + SDf)^2 - radicand >= 0, so SD <= SDb + SDf",
           decide((a + b) ** 2 - change_sd_sq(a, b, -1 + u), "is_nonnegative"))
    r.step("d radicand / dr = -2 SDb SDf < 0", decide(sp.diff(R, rr), "is_negative"))
    Rpos = sp.Symbol("R", positive=True)
    dsd = sp.diff(sp.sqrt(R), rr).subs(R, Rpos)
    r.step("d SD / dr = -SDb SDf / SD < 0 wherever SD > 0", decide(dsd, "is_negative"))
    r.step("r = 0 gives sqrt(SDb^2 + SDf^2)", is_zero(change_sd_sq(a, b, 0) - (a**2 + b**2)))
    r.step("edge: r = 1 and SDb = SDf gives SD 0 (degenerate, flag in the kit)", is_zero(change_sd_sq(a, a, 1)))
    return r


@mutant
def m12_change_sd_plus():
    a, b = sp.symbols("SDb SDf", positive=True)
    rr = sp.Symbol("r", real=True)
    wrong = a**2 + b**2 + 2 * rr * a * b
    proven = decide(sp.diff(wrong, rr), "is_negative")
    cx = search(lambda x: (x["a"] ** 2 + x["b"] ** 2 + 2 * 0.9 * x["a"] * x["b"]) < (x["a"] ** 2 + x["b"] ** 2 + 2 * 0.1 * x["a"] * x["b"]),
                lambda g: {"a": g.uniform(0.5, 3), "b": g.uniform(0.5, 3)})
    return Mutant("M12", "P08", "radicand SDb^2 + SDf^2 + 2 r SDb SDf (sign flipped)", proven, cx)


@mutant
def m13_change_sd_domain():
    a, b = sp.symbols("SDb SDf", positive=True)
    s = sp.Symbol("s", nonnegative=True)
    proven = decide(change_sd_sq(a, b, 2 - s), "is_nonnegative")   # r up to 2
    cx = search(lambda x: change_sd_sq(x["a"], x["b"], x["r"]) >= 0,
                lambda g: {"a": 1.0, "b": 1.0, "r": g.uniform(1.0, 2.0)})
    return Mutant("M13", "P08", "allow correlation up to 2 (domain check removed)", proven, cx)


# ---------------------------------------------------------------------------
# P9 cross-over paired SE vs parallel-design SE
# ---------------------------------------------------------------------------

@proof
def p09_crossover():
    r = Result("P09", "cross-over paired SE < parallel-design SE when within-person r > 0",
               "sympy assumptions")
    s1, s2, n = sp.symbols("s1 s2 n", positive=True)
    rp = sp.Symbol("r", positive=True)
    rn = sp.Symbol("r", negative=True)
    diff = parallel_se_sq(s1, s2, n) - paired_se_sq(s1, s2, rp, n)
    r.step("SE_par^2 - SE_paired^2 == 2 r s1 s2 / n", is_zero(diff - 2 * rp * s1 * s2 / n))
    r.step("r > 0: SE_par^2 - SE_paired^2 > 0", decide(diff, "is_positive"))
    q = sp.Symbol("q", positive=True)               # r = 1 - q < 1
    r.step("r < 1 (r = 1 - q): SE_paired^2 = ((s1 - s2)^2 + 2 q s1 s2)/n > 0, so both roots are real and SE_paired < SE_par",
           is_zero(paired_se_sq(s1, s2, 1 - q, n) - ((s1 - s2) ** 2 + 2 * q * s1 * s2) / n)
           and decide(((s1 - s2) ** 2 + 2 * q * s1 * s2) / n, "is_positive"))
    r.step("r = 0: the two SEs are equal", is_zero(parallel_se_sq(s1, s2, n) - paired_se_sq(s1, s2, 0, n)))
    r.step("r < 0: paired SE is larger", decide(paired_se_sq(s1, s2, rn, n) - parallel_se_sq(s1, s2, n), "is_positive"))
    return r


@mutant
def m14_paired_plus():
    s1, s2, n = sp.symbols("s1 s2 n", positive=True)
    rp = sp.Symbol("r", positive=True)
    wrong = (s1**2 + s2**2 + 2 * rp * s1 * s2) / n
    proven = decide(parallel_se_sq(s1, s2, n) - wrong, "is_positive")
    cx = search(lambda x: (x["s1"] ** 2 + x["s2"] ** 2) > (x["s1"] ** 2 + x["s2"] ** 2 + 2 * x["r"] * x["s1"] * x["s2"]),
                lambda g: {"s1": g.uniform(0.5, 3), "s2": g.uniform(0.5, 3), "r": g.uniform(0.05, 0.95)})
    return Mutant("M14", "P09", "paired variance with + 2 r s1 s2", proven, cx)


# ---------------------------------------------------------------------------
# P10 splitting a shared control keeps total n
# ---------------------------------------------------------------------------

@proof
def p10_split_shared_control():
    r = Result("P10", "splitting a shared control across k comparisons keeps total n (and the risk); split comparisons stay correlated, which a count cannot show (engine paths: E07)",
               "sympy general k (the engine's real n/k split); the Bend2 split laws are an integer model the engine does not run")
    i, k = sp.symbols("i k", integer=True, positive=True)
    n, e = sp.symbols("n e", positive=True)
    r.step("Sum_(i=1..k) n/k == n", sp.simplify(sp.Sum(n / k, (i, 1, k)).doit() - n) == 0)
    r.step("Sum_(i=1..k) events/k == events", sp.simplify(sp.Sum(e / k, (i, 1, k)).doit() - e) == 0)
    r.step("each sub-arm keeps the risk: (e/k)/(n/k) == e/n", is_zero((e / k) / (n / k) - e / n))
    return r


@mutant
def m15_full_control_each():
    i = sp.Symbol("i", integer=True, positive=True)
    k = sp.Symbol("k", integer=True, positive=True)
    n = sp.Symbol("n", positive=True)
    total = sp.Sum(n, (i, 1, k)).doit()
    proven = is_zero(total - n)
    cx = "" if proven else f"SYNTHETIC n=40, k=2: total = {total.subs({n: 40, k: 2})} != 40"
    return Mutant("M15", "P10", "shared control entered in full in each comparison", proven, cx)


# ---------------------------------------------------------------------------
# runner
# ---------------------------------------------------------------------------

def run() -> tuple[list[Result], list[Mutant]]:
    return [p() for p in PROOFS], [m() for m in MUTANTS]


def receipt(results, mutants) -> str:
    lines = [
        "sympy receipt for validation/proofs/sympy_proofs.py",
        f"date_utc: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"python: {sys.version.split()[0]}  sympy: {sp.__version__}",
        "",
        "PROOFS: id | property | method | steps passed | result",
    ]
    for r in results:
        passed = sum(ok for _, ok in r.steps)
        lines.append(f"{r.pid} | {r.prop} | {r.method} | {passed}/{len(r.steps)} | {'PASS' if r.ok else 'FAIL'}")
    lines += ["", "MUTANTS: id | targets | wrong version | sympy proved it? | counterexample | result"]
    for m in mutants:
        lines.append(f"{m.mid} | {m.target} | {m.wrong} | {'yes' if m.proven else 'no'} | "
                     f"{m.counterexample or 'none found'} | {'REFUTED' if m.refuted else 'SURVIVED'}")
    lines += ["", "STEPS"]
    for r in results:
        lines.append(f"{r.pid} {r.prop}")
        for text, ok in r.steps:
            lines.append(f"  [{'ok' if ok else 'FAIL'}] {text}")
    ok = all(r.ok for r in results) and all(m.refuted for m in mutants)
    lines += ["", f"overall: {'PASS' if ok else 'FAIL'}"]
    return "\n".join(lines) + "\n"


def main() -> int:
    results, mutants = run()
    text = receipt(results, mutants)
    (HERE / "receipt_sympy.txt").write_text(text)
    print(text)
    return 0 if all(r.ok for r in results) and all(m.refuted for m in mutants) else 1


if __name__ == "__main__":
    raise SystemExit(main())
