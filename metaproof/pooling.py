"""Inverse-variance fixed-effect and random-effects (mixed-effects) models.

One fitting routine serves both pooling (X = column of ones) and meta-regression.

Sources
-------
* DerSimonian R, Laird N. Control Clin Trials 1986;7(3):177-188 (DL).
* Paule RC, Mandel J. J Res Natl Bur Stand 1982;87(5):377-385 (PM).
* Viechtbauer W. J Educ Behav Stat 2005;30(3):261-293 (ML, REML, Fisher scoring).
* Hartung J, Knapp G. Stat Med 2001;20(24):3875-3889; Sidik K, Jonkman JN. Stat Med
  2002;21(21):3153-3159; Knapp G, Hartung J. Stat Med 2003;22(17):2693-2710 (HKSJ / KH, and
  the ad hoc truncation of the scale at 1).
* When to use HKSJ: Cochrane Handbook v6.5 section 10.10.4.4; Roever C, Knapp G, Friede T.
  BMC Med Res Methodol 2015;15:99; Jackson D, Law M, Ruecker G, Schwarzer G. Stat Med
  2017;36(25):3923-3934 (HKSJ can be narrower than the Wald interval when tau2 is estimated
  as 0, and very wide when k is very small).
* Viechtbauer W. Stat Med 2007;26(1):37-52 (Q-profile CI for tau^2).
* Higgins JPT, Thompson SG. Stat Med 2002;21(11):1539-1558 (I^2, H^2).
* Higgins JPT, Thompson SG, Spiegelhalter DJ. J R Stat Soc A 2009;172(1):137-159;
  IntHout J, Ioannidis JPA, Rovers MM, Goeman JJ. BMJ Open 2016;6:e010247 (prediction interval).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy import optimize, stats

ESTIMATORS = ("FE", "DL", "PM", "ML", "REML")
HKSJ_RULES = ("always", "revman", "truncated")
# tau2 counts as 0 for the revman rule below this fraction of the mean within-study variance,
# so the two engines take the same branch when an estimate is 0 up to rounding.
TAU2_ZERO = 1e-8
TOL = 1e-10
MAXITER = 1000


@dataclass
class Fit:
    method: str
    k: int
    p: int
    beta: np.ndarray
    vb: np.ndarray          # variance-covariance used for inference (HK-scaled when hksj)
    vb_model: np.ndarray    # model-based (X'WX)^-1, never HK-scaled
    se: np.ndarray
    stat: np.ndarray
    pval: np.ndarray
    ci_lb: np.ndarray
    ci_ub: np.ndarray
    df: float               # inf for z tests
    tau2: float
    QE: float               # Cochran residual Q with FE weights, df k - p
    QE_df: int
    QE_p: float
    I2: float
    H2: float
    hksj: bool
    s2_hk: float
    hksj_rule: str = "always"
    converged: bool = True
    iterations: int = 0
    notes: list = field(default_factory=list)

    @property
    def estimate(self) -> float:
        return float(self.beta[0])


def _as_X(k, X):
    if X is None:
        return np.ones((k, 1))
    X = np.asarray(X, dtype=float)
    return X.reshape(k, -1)


def _wls(y, v, X, tau2):
    """Weighted least squares at a given tau2. Returns beta, (X'WX)^-1, W, P.

    beta comes from a least-squares solve of the sqrt(W)-scaled system rather than from
    (X'WX)^-1 X'Wy, which keeps the residuals accurate when the effects share a large
    common offset."""
    w = 1.0 / (v + tau2)
    sw = np.sqrt(w)
    XtW = X.T * w
    vb = np.linalg.inv(XtW @ X)
    beta = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    W = np.diag(w)
    P = W - (W @ X) @ vb @ (X.T @ W)
    return beta, vb, w, P


def _center(y, X):
    """Shift y by its median when the column space of X contains the constant vector (every
    model in this kit: an intercept, or one dummy per subgroup). Returns y - c, c and a with
    X a = 1, so beta(y) = beta(y - c) + c a. Q, tau2 and the standard errors do not depend
    on c; fitting the centred values keeps their accuracy when the effects share a large
    offset (a mean difference in raw units far from 0)."""
    if len(y) == 0:
        return y, 0.0, np.zeros(X.shape[1])
    a = np.linalg.lstsq(X, np.ones(len(y)), rcond=None)[0]
    if not np.allclose(X @ a, 1.0, rtol=0, atol=1e-10):
        return y, 0.0, np.zeros(X.shape[1])
    c = float(np.median(y))
    return y - c, c, a


def _wrss(y, X, beta, w) -> float:
    """sum w_i (y_i - x_i b)^2. Equal to y'Py, but y'Py subtracts two numbers of size
    sum w y^2 and loses every significant digit when the y share a large offset (for
    y = (0, 1, 2) + 1e8 it gives a negative Q). The residual form does not cancel; metafor
    centres y before forming y'Py for the same reason."""
    r = y - X @ beta
    return float(np.sum(w * r * r))


def q_gen(y, v, X, tau2) -> float:
    """Generalised Q(tau2) = y' P(tau2) y = sum w_i (y_i - x_i b(tau2))^2."""
    y, _, _ = _center(np.asarray(y, dtype=float), X)
    beta, _, w, _ = _wls(y, v, X, tau2)
    return _wrss(y, X, beta, w)


def _tau2_dl(y, v, X):
    k, p = X.shape
    beta, _, w, P0 = _wls(y, v, X, 0.0)
    QE = _wrss(y, X, beta, w)
    return max(0.0, (QE - (k - p)) / np.trace(P0))


def _upper_bracket(f, start):
    hi = max(start, 1.0)
    for _ in range(200):
        if f(hi) < 0:
            return hi
        hi *= 2
    raise RuntimeError("could not bracket root")


def _tau2_pm(y, v, X):
    k, p = X.shape
    f = lambda t2: q_gen(y, v, X, t2) - (k - p)
    if f(0.0) <= 0:
        return 0.0, True, 0
    hi = _upper_bracket(f, 10 * float(np.var(y) + np.max(v)))
    root, res = optimize.brentq(f, 0.0, hi, xtol=1e-14, rtol=1e-14, full_output=True)
    return root, res.converged, res.iterations


def loglik(y, v, X, tau2, reml: bool) -> float:
    w = 1.0 / (v + tau2)
    vb = np.linalg.inv((X.T * w) @ X)
    beta = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
    ll = -0.5 * (len(y) * math.log(2 * math.pi) + np.sum(np.log(v + tau2))
                 + _wrss(y, X, beta, w))
    if reml:
        p = X.shape[1]
        ll += 0.5 * p * math.log(2 * math.pi) + 0.5 * np.linalg.slogdet(vb)[1] \
            + 0.5 * np.linalg.slogdet(X.T @ X)[1]
    return float(ll)


def _score_info(y, v, X, tau2, reml):
    beta, vb, w, P = _wls(y, v, X, tau2)
    # P y = W (y - X b), so both scores use the residuals and avoid forming P y directly.
    wr = w * (y - X @ beta)
    if reml:
        return float(wr @ wr - np.trace(P)), float(np.trace(P @ P))
    return float(wr @ wr - np.sum(w)), float(np.sum(w**2))


def _tau2_fisher(y, v, X, reml):
    """Fisher scoring with step halving; returns tau2, converged, iterations, notes.

    Update: tau2 += U / I with U = y'PPy - tr(P), I = tr(PP) (REML), or U = r'W^2r - tr(W),
    I = tr(W^2) (ML). Truncated at 0. Converged when |delta| < TOL and the KKT condition
    holds: score ~ 0 for tau2 > 0, score <= 0 at tau2 = 0."""
    notes = []
    tau2 = _tau2_dl(y, v, X)
    ll = loglik(y, v, X, tau2, reml)
    it = 0
    for it in range(1, MAXITER + 1):
        U, I = _score_info(y, v, X, tau2, reml)
        step = U / I
        new = max(0.0, tau2 + step)
        new_ll = loglik(y, v, X, new, reml)
        halvings = 0
        while new_ll < ll - 1e-12 and halvings < 60:
            step /= 2
            new = max(0.0, tau2 + step)
            new_ll = loglik(y, v, X, new, reml)
            halvings += 1
        delta = abs(new - tau2)
        tau2, ll = new, new_ll
        if delta < TOL * max(1.0, tau2):
            break
    U, I = _score_info(y, v, X, tau2, reml)
    scale = max(1.0, abs(I) * max(tau2, 1e-8))
    kkt = (tau2 > 0 and abs(U) / I < 1e-6 * max(1.0, tau2)) or (tau2 == 0 and U <= 1e-8 * scale)
    converged = it < MAXITER and kkt
    if not converged:
        # Fallback: direct bounded maximisation of the (restricted) likelihood.
        hi = 100 * float(np.var(y) + np.max(v)) + 1.0
        res = optimize.minimize_scalar(lambda t: -loglik(y, v, X, t, reml), bounds=(0, hi),
                                       method="bounded", options={"xatol": 1e-12})
        cand = float(res.x)
        if loglik(y, v, X, cand, reml) >= loglik(y, v, X, tau2, reml):
            tau2 = cand
        if loglik(y, v, X, 0.0, reml) >= loglik(y, v, X, tau2, reml):
            tau2 = 0.0
        notes.append("Fisher scoring did not meet convergence checks; bounded optimiser used")
        converged = bool(res.success)
    tau2, moved = _global_max(y, v, X, reml, tau2)
    if moved:
        notes.append(moved)
    return tau2, converged, it, notes


def _global_max(y, v, X, reml, tau2):
    """The (restricted) likelihood in tau2 can have more than one local maximum, for example
    when one study is far more precise than the rest; Fisher scoring then stops at whichever
    it reaches first and still passes the score check. ML and REML are defined as the global
    maximiser on [0, inf), so scan a grid that includes the boundary tau2 = 0 and move to a
    better point when one exists (metafor 5.2 checks the boundary in the same situation)."""
    ll_hat = loglik(y, v, X, tau2, reml)
    scale = float(np.var(y) + np.max(v))
    grid = np.concatenate([[0.0], np.geomspace(1e-8 * scale, 1e3 * scale, 300)])
    lls = np.array([loglik(y, v, X, t, reml) for t in grid])
    j = int(np.argmax(lls))
    if not lls[j] > ll_hat + 1e-9 * max(1.0, abs(ll_hat)):
        return tau2, ""
    if j == 0:
        return 0.0, (f"likelihood higher at tau2 = 0 than at the Fisher-scoring solution "
                     f"{tau2:.6g}; tau2 set to 0 (global maximum)")
    lo, hi = grid[j - 1], grid[min(j + 1, len(grid) - 1)]
    res = optimize.minimize_scalar(lambda t: -loglik(y, v, X, t, reml), bounds=(lo, hi),
                                   method="bounded", options={"xatol": 1e-14 * max(1.0, hi)})
    new = float(res.x) if -res.fun >= lls[j] else float(grid[j])
    return new, (f"Fisher scoring reached a local maximum at tau2 = {tau2:.6g}; global maximum "
                 f"at {new:.6g} used")


def estimate_tau2(y, v, X, method):
    if method == "FE":
        return 0.0, True, 0, []
    if method == "DL":
        return _tau2_dl(y, v, X), True, 0, []
    if method == "PM":
        t2, ok, it = _tau2_pm(y, v, X)
        return t2, ok, it, []
    if method in ("ML", "REML"):
        return _tau2_fisher(y, v, X, reml=(method == "REML"))
    raise ValueError(f"unknown estimator {method!r}; choose from {ESTIMATORS}")


def hksj_applies(rule: str, method: str, k: int, p: int, tau2: float, v) -> bool:
    """Whether the HKSJ adjustment is used for this fit under the configured rule.

    always / truncated : every random-effects fit with k - p >= 1.
    revman             : only when tau2 > 0 and k - p >= 2 (for pooling, k >= 3); otherwise
                         the Wald z interval. This is the rule the review cites from Cochrane
                         Handbook v6.5 section 10.10.4.4 for RevMan. It avoids the two known
                         failures of HKSJ with few studies: an interval narrower than the
                         fixed-effect one when tau2 = 0, and t(1) = 12.71 at k = 2."""
    if rule not in HKSJ_RULES:
        raise ValueError(f"unknown hksj_rule {rule!r}; choose from {HKSJ_RULES}")
    if method == "FE" or k - p < 1:
        return False
    if rule == "revman":
        return k - p >= 2 and tau2 > TAU2_ZERO * float(np.mean(v))
    return True


def fit(y, v, X=None, method="REML", hksj=True, level=0.95, tau2=None,
        hksj_rule="always") -> Fit:
    """Fit a fixed-effect (method='FE') or random/mixed-effects model.

    hksj applies the Hartung-Knapp(-Sidik-Jonkman) adjustment to random-effects fits, as
    decided by hksj_rule (see hksj_applies): Var_HK(b) = s2 (X'WX)^-1, s2 = r'Wr / (k - p),
    t on k - p df; under 'truncated' s2 is replaced by max(1, s2). Never applied to FE.
    """
    y = np.asarray(y, dtype=float)
    v = np.asarray(v, dtype=float)
    k = len(y)
    X = _as_X(k, X)
    p = X.shape[1]
    if k < p:
        raise ValueError(f"k={k} < p={p}")
    y, c, a = _center(y, X)
    notes = []
    if tau2 is None:
        if k - p < 1 and method != "FE":
            tau2_hat, conv, it, n2 = 0.0, True, 0, []
            notes.append("k - p = 0: tau2 not estimable, set to 0")
        else:
            tau2_hat, conv, it, n2 = estimate_tau2(y, v, X, method)
        notes += n2
    else:
        tau2_hat, conv, it = float(tau2), True, 0

    beta_c, vb_model, w, P = _wls(y, v, X, tau2_hat)
    beta = beta_c + c * a
    use_hk = bool(hksj) and hksj_applies(hksj_rule, method, k, p, tau2_hat, v)
    if hksj and method != "FE" and not use_hk:
        if k - p < 1:
            notes.append("HKSJ not applicable with k - p = 0; Wald z used")
        else:
            notes.append("HKSJ not used (hksj_rule revman: tau2 = 0 or k - p < 2); Wald z used")
    alpha = 1 - level
    if use_hk:
        s2 = _wrss(y, X, beta_c, w) / (k - p)
        if hksj_rule == "truncated":
            s2 = max(1.0, s2)
        vb = s2 * vb_model
        df = float(k - p)
        crit = stats.t.ppf(1 - alpha / 2, df)
    else:
        s2 = float("nan")
        vb = vb_model
        df = math.inf
        crit = stats.norm.ppf(1 - alpha / 2)
    se = np.sqrt(np.diag(vb))
    stat = beta / se
    if use_hk:
        pval = 2 * stats.t.sf(np.abs(stat), df)
    else:
        pval = 2 * stats.norm.sf(np.abs(stat))

    # Cochran's residual Q always uses FE weights.
    beta0, _, w0, P0 = _wls(y, v, X, 0.0)
    QE = _wrss(y, X, beta0, w0)
    QE_df = k - p
    QE_p = float(stats.chi2.sf(QE, QE_df)) if QE_df > 0 else float("nan")
    if QE_df <= 0:
        I2 = H2 = float("nan")
    elif method == "FE":
        I2 = max(0.0, 100 * (QE - QE_df) / QE) if QE > 0 else 0.0
        H2 = QE / QE_df
    else:
        s2_typ = QE_df / float(np.trace(P0))  # "typical" within-study variance
        I2 = 100 * tau2_hat / (tau2_hat + s2_typ)
        H2 = (tau2_hat + s2_typ) / s2_typ
    return Fit(method=method, k=k, p=p, beta=beta, vb=vb, vb_model=vb_model, se=se, stat=stat,
               pval=pval, ci_lb=beta - crit * se, ci_ub=beta + crit * se, df=df, tau2=tau2_hat,
               QE=QE, QE_df=QE_df, QE_p=QE_p, I2=I2, H2=H2, hksj=use_hk, s2_hk=s2,
               converged=conv, iterations=it, notes=notes, hksj_rule=hksj_rule)


def prediction_interval(f: Fit, level=0.95, pi_type="riley"):
    """mu +/- crit sqrt(tau2 + SE(mu)^2). SE(mu) is the one used for inference (HK-scaled
    when HKSJ is on).

    pi_type 'riley'  : crit = t_{k-2} (Higgins 2009; Riley RD, Higgins JPT, Deeks JJ.
                       BMJ 2011;342:d549; IntHout 2016). Needs k >= 3.
    pi_type 'metafor': metafor predict() default, the CI's own multiplier: t_{k-1} when the
                       fit used HKSJ, z otherwise. Under hksj_rule revman this is what
                       RevMan prints: Handbook v6.5 10.10.4.4 says RevMan uses the same
                       multiplier for the PI as for the CI (t on k - 1 df with HKSJ, the
                       normal quantile with the Wald interval).
    pi_type 'revman' : crit = t_{k-1} whatever the CI method: the formula of Handbook v6.5
                       10.10.4.3 ("t_{k-1} is the 97.5th percentile of a t-distribution with
                       k-1 degrees of freedom"). It equals RevMan's output only when the
                       HKSJ interval is used; for RevMan's Wald-interval output use 'metafor'.
                       Needs k >= 2.

    With few studies the t quantile dominates the width (t_1 = 12.71), so the pipeline flags
    any PI from fewer than 5 studies as not interpretable (Handbook v6.5 10.10.4.3)."""
    if f.method == "FE" or f.p != 1:
        return float("nan"), float("nan")
    q = 1 - (1 - level) / 2
    if pi_type == "metafor":
        crit = stats.t.ppf(q, f.df) if f.hksj else stats.norm.ppf(q)
    elif pi_type == "revman":
        if f.k < 2:
            return float("nan"), float("nan")
        crit = stats.t.ppf(q, f.k - 1)
    elif f.k < 3:
        return float("nan"), float("nan")
    else:
        crit = stats.t.ppf(q, f.k - 2)
    half = crit * math.sqrt(f.tau2 + f.se[0] ** 2)
    return f.estimate - half, f.estimate + half


def tau2_qprofile_ci(y, v, X=None, level=0.95):
    """Q-profile CI (Viechtbauer 2007): solve Q_gen(tau2) = chi2_{df, 1-a/2} for the lower
    bound and = chi2_{df, a/2} for the upper; a bound is 0 when Q_gen(0) is below its target."""
    y = np.asarray(y, dtype=float)
    v = np.asarray(v, dtype=float)
    X = _as_X(len(y), X)
    y, _, _ = _center(y, X)
    df = len(y) - X.shape[1]
    if df < 1:
        return float("nan"), float("nan")
    a = 1 - level
    targets = (stats.chi2.ppf(1 - a / 2, df), stats.chi2.ppf(a / 2, df))
    bounds = []
    for target in targets:
        f = lambda t2: q_gen(y, v, X, t2) - target
        if f(0.0) <= 0:
            bounds.append(0.0)
            continue
        hi = _upper_bracket(f, 10 * float(np.var(y) + np.max(v)))
        bounds.append(optimize.brentq(f, 0.0, hi, xtol=1e-14, rtol=1e-12))
    return bounds[0], bounds[1]


POOLED_COLUMNS = ["estimator", "k", "estimate", "se", "ci_lb", "ci_ub", "pi_lb", "pi_ub",
                  "tau2", "tau2_ci_lb", "tau2_ci_ub", "I2", "H2", "Q", "Q_p"]
# Written after the interface columns. ci_lb/ci_ub are the interval of the configured rule
# (ci_method); the three intervals that rule chooses between are always written as well, so
# protocol sensitivity K6 (truncated HKSJ) and "report both when HKSJ is narrower" need no
# second run.
POOLED_EXTRA_COLUMNS = ["ci_method", "df", "wald_ci_lb", "wald_ci_ub", "hksj_ci_lb",
                        "hksj_ci_ub", "hksj_trunc_ci_lb", "hksj_trunc_ci_ub"]


def _interval(y, v, method, level, rule, hksj=True):
    f = fit(y, v, None, method, hksj, level, hksj_rule=rule)
    return float(f.ci_lb[0]), float(f.ci_ub[0]), f


def ci_method(f: Fit) -> str:
    if not f.hksj:
        return "Wald z"
    lab = "HKSJ truncated" if f.hksj_rule == "truncated" else "HKSJ"
    return f"{lab} t({f.df:g})"


def pooled_row(y, v, method, hksj=True, level=0.95, pi=True, pi_type="riley",
               hksj_rule="always"):
    """One row of pooled.csv plus the test statistics, for a given estimator."""
    f = fit(y, v, None, method, hksj, level, hksj_rule=hksj_rule)
    pi_lb, pi_ub = prediction_interval(f, level, pi_type) if pi else (float("nan"), float("nan"))
    if method == "FE":
        t_lb = t_ub = float("nan")
    else:
        t_lb, t_ub = tau2_qprofile_ci(y, v, None, level)
    row = {"estimator": method, "k": f.k, "estimate": f.estimate, "se": float(f.se[0]),
           "ci_lb": float(f.ci_lb[0]), "ci_ub": float(f.ci_ub[0]), "pi_lb": pi_lb, "pi_ub": pi_ub,
           "tau2": f.tau2, "tau2_ci_lb": t_lb, "tau2_ci_ub": t_ub, "I2": f.I2, "H2": f.H2,
           "Q": f.QE, "Q_p": f.QE_p}
    nan = float("nan")
    w_lb, w_ub, _ = _interval(y, v, method, level, "always", hksj=False)
    if method != "FE" and f.k >= 2:
        h_lb, h_ub, _ = _interval(y, v, method, level, "always")
        t2_lb, t2_ub, _ = _interval(y, v, method, level, "truncated")
    else:
        h_lb = h_ub = t2_lb = t2_ub = nan
    row.update({"ci_method": ci_method(f), "df": f.df if f.hksj else nan,
                "wald_ci_lb": w_lb, "wald_ci_ub": w_ub, "hksj_ci_lb": h_lb, "hksj_ci_ub": h_ub,
                "hksj_trunc_ci_lb": t2_lb, "hksj_trunc_ci_ub": t2_ub})
    test = {"estimator": method, "k": f.k, "test": "t" if f.hksj else "z",
            "stat": float(f.stat[0]), "df": f.df if f.hksj else float("nan"),
            "pval": float(f.pval[0]), "hksj": f.hksj, "converged": f.converged,
            "iterations": f.iterations, "notes": "; ".join(f.notes)}
    return row, test, f
