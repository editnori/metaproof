"""Leave-one-out, influence diagnostics, subgroup analysis, meta-regression,
Egger's test and Benjamini-Hochberg adjustment.

Sources
-------
* Viechtbauer W, Cheung MW-L. Outlier and influence diagnostics for meta-analysis.
  Res Synth Methods 2010;1(2):112-125.
* Egger M, Davey Smith G, Schneider M, Minder C. BMJ 1997;315(7109):629-634.
* Sterne JAC, Sutton AJ, Ioannidis JPA, et al. BMJ 2011;343:d4002 (asymmetry tests only
  with at least 10 studies, not when the studies are of similar size); Cochrane Handbook
  v6.5, chapter 13.
* Pustejovsky JE, Rodgers MA. Res Synth Methods 2019;10(1):57-71 (Egger test for SMD with
  sqrt(1/n1 + 1/n2) in place of the SE).
* Higgins JPT, Thompson SG. Stat Med 2004;23(11):1663-1682 (permutation test for
  meta-regression); Cochrane Handbook v6.5, sections 10.11.3.1 and 10.11.5.1 (subgroups
  and meta-regression need about ten studies per characteristic).
* Benjamini Y, Hochberg Y. J R Stat Soc B 1995;57(1):289-300.
* Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis.
  Wiley; 2009 (subgroup analyses, Q-between).
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats

from . import pooling as pl

NAN = float("nan")


def leave_one_out(ids, y, v, method, hksj, level, pi=True, pi_type="riley", groups=None,
                  hksj_rule="always") -> pd.DataFrame:
    """Refit without each study. With groups (the source study of each comparison), all
    comparisons of one study are dropped together, so a split multi-arm trial is left out
    as one trial."""
    rows = []
    y, v = np.asarray(y, float), np.asarray(v, float)
    groups = list(ids) if groups is None else list(groups)
    for g in dict.fromkeys(groups):
        keep = np.array([gi != g for gi in groups])
        if keep.sum() == 0:
            continue
        row, test, _ = pl.pooled_row(y[keep], v[keep], method, hksj, level, pi, pi_type,
                                     hksj_rule)
        rows.append({"omitted": g, **{c: row[c] for c in pl.POOLED_COLUMNS if c != "estimator"},
                     "ci_method": row["ci_method"], "pval": test["pval"]})
    return pd.DataFrame(rows)


def influence(ids, y, v, method, X=None, hksj=False, hksj_rule="always") -> pd.DataFrame:
    """Case-deletion diagnostics (Viechtbauer & Cheung 2010), computed the way metafor
    influence.rma.uni (version 5.2) does, so the two engines agree.

    For each i the model is refitted without study i (same estimator, same test). With
    s2_(-i) the Knapp-Hartung scale of that refit (1 when hksj is off), Var(.) the variance
    used for inference (HK-scaled when hksj is on) and h_i the full-model hat value:

    rstudent_i = (y_i - x_i b_(-i)) / sqrt(s2_(-i) (v_i + tau2_(-i)) + x_i Var(b_(-i)) x_i')
    dffits_i   = (x_i b - x_i b_(-i)) / sqrt(s2_(-i) h_i (tau2_(-i) + v_i))
    cook_d_i   = (b - b_(-i))' Var(b)^-1 (b - b_(-i))
    cov_ratio  = det Var(b_(-i)) / det Var(b)
    dfbetas    = (b - b_(-i)) / sqrt(s2_(-i) diag (X' W_(-i) X)^-1), W_(-i) = diag 1/(v_j +
                 tau2_(-i)) over all k studies
    hat        = diag(X (X'WX)^-1 X' W)
    A study is flagged when |dffits| > 3 sqrt(p/(k-p)), P(chi2_p <= cook_d) > 0.5,
    hat > 3p/k or any |dfbetas| > 1 (the metafor rules of thumb, not formal tests).

    hksj_rule decides whether the full model uses HKSJ; every deleted fit then uses the same
    test as the full model (as metafor does), so the measures are comparable across studies.
    """
    y, v = np.asarray(y, float), np.asarray(v, float)
    k = len(y)
    X = pl._as_X(k, X)
    p = X.shape[1]
    full = pl.fit(y, v, X, method, hksj=hksj, hksj_rule=hksj_rule)
    hksj = full.hksj
    del_rule = "truncated" if hksj_rule == "truncated" else "always"
    if full.hksj and not full.s2_hk > 0:
        # Identical fitted values give a Knapp-Hartung scale of 0, so Var(b) = 0 and the
        # measures are undefined (metafor's influence() stops with an error here too).
        return pd.DataFrame([{"note": "influence not defined: HKSJ scale s2 = 0 "
                                      "(all residuals are 0)"}])
    w = 1 / (v + full.tau2)
    H = X @ full.vb_model @ (X.T * w)
    hat = np.diag(H)
    svb = np.linalg.inv(full.vb)
    rows = []
    for i in range(k):
        keep = np.arange(k) != i
        row = {"study_id": ids[i], "weight_pct": 100 * w[i] / w.sum(), "hat": hat[i]}
        if keep.sum() < p + 1 and method != "FE":
            rows.append({**row, "note": "too few studies after deletion"})
            continue
        d = pl.fit(y[keep], v[keep], X[keep], method, hksj=hksj, hksj_rule=del_rule)
        s2 = d.s2_hk if d.hksj else 1.0
        xi = X[i]
        delpred = float(xi @ d.beta)
        resid = y[i] - delpred
        if abs(resid) < 100 * np.finfo(float).eps:
            resid = 0.0
        row["rstudent"] = resid / math.sqrt(s2 * (v[i] + d.tau2) + xi @ d.vb @ xi)
        row["dffits"] = (xi @ full.beta - delpred) / math.sqrt(s2 * hat[i] * (d.tau2 + v[i]))
        diff = full.beta - d.beta
        row["cook_d"] = float(diff @ svb @ diff)
        row["cov_ratio"] = float(np.linalg.det(d.vb) / np.linalg.det(full.vb))
        row["tau2_del"] = d.tau2
        row["QE_del"] = d.QE
        vb_del = np.linalg.inv((X.T / (v + d.tau2)) @ X)
        dfb = diff / np.sqrt(s2 * np.diag(vb_del))
        row["dfbetas"] = float(dfb[0]) if p == 1 else ";".join(f"{x:.6g}" for x in dfb)
        row["influential"] = bool(abs(row["dffits"]) > 3 * math.sqrt(p / (k - p))
                                  or stats.chi2.cdf(row["cook_d"], p) > 0.5
                                  or hat[i] > 3 * p / k
                                  or bool(np.any(np.abs(dfb) > 1)))
        rows.append(row)
    cols = ["study_id", "rstudent", "dffits", "cook_d", "cov_ratio", "tau2_del", "QE_del",
            "hat", "weight_pct", "dfbetas", "influential"]
    out = pd.DataFrame(rows)
    return out.reindex(columns=cols + [c for c in out.columns if c not in cols])


def _not_run(base: dict, note: str) -> pd.DataFrame:
    return pd.DataFrame([{**base, "note": note}])


def subgroups(eff: pd.DataFrame, moderator: str, method, hksj, level, pi=True,
              pi_type="riley", tau2_mode="separate", min_k=0, min_per_level=0,
              hksj_rule="always", study_col=None) -> pd.DataFrame:
    """Subgroup analysis with a test for subgroup differences.

    tau2_mode 'common' (protocol I): one mixed-effects model with a dummy per subgroup and no
      intercept, so every subgroup shares tau2; the test for differences is the omnibus
      moderator test of the same model with an intercept (F with Knapp-Hartung, chi2
      otherwise), as metafor rma(mods = ~ factor) reports it.
    tau2_mode 'separate': a model per subgroup (own tau2). Q_between = sum_j w_j (mu_j -
      mu_bar)^2, w_j = 1 / SE_model(mu_j)^2, df = J - 1, where SE_model is the model-based
      (non-HKSJ) SE of each subgroup estimate.
    Not run (one 'note' row) when the number of studies is below min_k or any subgroup has
    fewer than min_per_level studies (protocol I: 10 and 2). Studies with a missing moderator
    value are dropped from this analysis."""
    d = eff.dropna(subset=[moderator])
    sc = study_col if study_col in d else None
    k_st = d[sc].nunique() if sc else len(d)
    counts = d.groupby(moderator)[sc].nunique() if sc else d.groupby(moderator).size()
    base = {"moderator": moderator, "k": len(d), "k_studies": k_st}
    if k_st < min_k:
        return _not_run(base, f"not run: {k_st} studies < subgroup_min_k = {min_k} (protocol I)")
    if len(counts) and counts.min() < min_per_level:
        return _not_run(base, f"not run: a subgroup has {int(counts.min())} studies < "
                              f"subgroup_min_per_level = {min_per_level} (protocol I)")
    if tau2_mode == "common":
        return _subgroups_common(d, moderator, method, hksj, level, hksj_rule, counts)
    rows, mus, ses = [], [], []
    for level_val, g in d.groupby(moderator, sort=True):
        y, v = g["yi"].to_numpy(float), g["vi"].to_numpy(float)
        row, test, f = pl.pooled_row(y, v, method, hksj, level, pi, pi_type, hksj_rule)
        if len(y) == 1:
            # One study: tau2 and heterogeneity are not estimable; report them as missing
            # rather than 0 (the estimate and its z CI are the study's own).
            for c in ("tau2", "tau2_ci_lb", "tau2_ci_ub", "I2", "H2", "Q", "Q_p"):
                row[c] = NAN
        wald = pl.fit(y, v, None, method, hksj=False, level=level)
        mus.append(wald.estimate)
        ses.append(float(wald.se[0]))
        rows.append({"moderator": moderator, "level": level_val,
                     **{c: row[c] for c in pl.POOLED_COLUMNS if c != "estimator"},
                     "pval": test["pval"], "se_model": float(wald.se[0]),
                     "ci_method": row["ci_method"]})
    out = pd.DataFrame(rows)
    if len(mus) >= 2:
        wj = 1 / np.square(ses)
        mu_bar = np.sum(wj * mus) / np.sum(wj)
        qb = float(np.sum(wj * (np.asarray(mus) - mu_bar) ** 2))
        out["Q_between"] = qb
        out["Q_between_df"] = len(mus) - 1
        out["Q_between_p"] = stats.chi2.sf(qb, len(mus) - 1)
    else:
        out["Q_between"] = out["Q_between_df"] = out["Q_between_p"] = NAN
    out["Q_between_test"] = "chi2 (separate tau2, model-based SEs)"
    out["tau2_mode"] = "separate"
    out.insert(2, "estimator", method)
    return out


def _subgroups_common(d, moderator, method, hksj, level, hksj_rule, counts):
    levels = sorted(d[moderator].astype(str).unique())
    y, v = d["yi"].to_numpy(float), d["vi"].to_numpy(float)
    lab = d[moderator].astype(str).to_numpy()
    X0 = np.column_stack([(lab == lev).astype(float) for lev in levels])
    if method == "FE":
        f = pl.fit(y, v, X0, "FE", hksj=False, level=level)
        f_test = pl.fit(y, v, np.column_stack([np.ones(len(y)), X0[:, 1:]]), "FE",
                        hksj=False, level=level)
    else:
        f = pl.fit(y, v, X0, method, hksj, level, hksj_rule=hksj_rule)
        f_test = pl.fit(y, v, np.column_stack([np.ones(len(y)), X0[:, 1:]]), method, hksj,
                        level, hksj_rule=hksj_rule)
    rows = []
    for j, lev in enumerate(levels):
        g = d[lab == lev]
        rows.append({"moderator": moderator, "level": lev, "k": len(g),
                     "estimate": f.beta[j], "se": f.se[j], "ci_lb": f.ci_lb[j],
                     "ci_ub": f.ci_ub[j], "pi_lb": NAN, "pi_ub": NAN, "tau2": f.tau2,
                     "tau2_ci_lb": NAN, "tau2_ci_ub": NAN, "I2": NAN, "H2": NAN, "Q": NAN,
                     "Q_p": NAN, "pval": f.pval[j], "se_model": math.sqrt(f.vb_model[j, j]),
                     "ci_method": pl.ci_method(f)})
    out = pd.DataFrame(rows)
    m = len(levels) - 1
    if m >= 1 and len(y) > m + 1:
        b = f_test.beta[1:]
        QM = float(b @ np.linalg.inv(f_test.vb[1:, 1:]) @ b)
        if f_test.hksj:
            out["Q_between"], out["Q_between_df"] = QM / m, m
            out["Q_between_df2"] = f_test.df
            out["Q_between_p"] = float(stats.f.sf(QM / m, m, f_test.df))
            out["Q_between_test"] = "F (common tau2, Knapp-Hartung)"
        else:
            out["Q_between"], out["Q_between_df"], out["Q_between_df2"] = QM, m, NAN
            out["Q_between_p"] = float(stats.chi2.sf(QM, m))
            out["Q_between_test"] = "chi2 (common tau2)"
        out["QE_residual"], out["QE_residual_p"] = f.QE, f.QE_p
    else:
        out["Q_between"] = out["Q_between_df"] = out["Q_between_df2"] = out["Q_between_p"] = NAN
        out["Q_between_test"] = "not estimable"
    out["tau2_mode"] = "common"
    out.insert(2, "estimator", method)
    return out


def design_matrix(eff: pd.DataFrame, moderators: list[str]):
    """Intercept plus numeric columns as-is and categorical columns as treatment dummies
    (reference = first level in sorted order). Rows with missing moderators are dropped.

    Levels are sorted by code point (C-locale order: 'Older adults' sorts before 'adults').
    R's factor() sorts by the session collation instead, so run_meta.R must sort levels in
    C-locale order too or the two engines choose different reference levels."""
    d = eff.dropna(subset=moderators)
    parts, names = [np.ones((len(d), 1))], ["intrcpt"]
    for m in moderators:
        col = d[m]
        num = pd.to_numeric(col, errors="coerce")
        if num.notna().all():
            parts.append(num.to_numpy(float)[:, None])
            names.append(m)
        else:
            levels = sorted(col.astype(str).unique())
            for lev in levels[1:]:
                parts.append((col.astype(str) == lev).to_numpy(float)[:, None])
                names.append(f"{m}[{lev}]")
    return d, np.hstack(parts), names


def _level_counts(d: pd.DataFrame, moderators, study_col):
    """Studies per category of each categorical moderator ({moderator: Series})."""
    out = {}
    for m in moderators:
        col = d[m]
        if pd.to_numeric(col, errors="coerce").notna().all():
            continue
        lab = col.astype(str)
        out[m] = (d.groupby(lab)[study_col].nunique() if study_col in d
                  else lab.value_counts())
    return out


def metareg(eff: pd.DataFrame, moderators: list[str], method, hksj, level,
            model_name=None, min_k=0, hksj_rule="always", permutations=0, seed=20261101,
            study_col=None) -> pd.DataFrame:
    """Mixed-effects meta-regression with Knapp-Hartung inference when hksj is on.

    Omnibus test of all non-intercept terms: QM = b' Var(b)^-1 b, reported as
    F = QM / m on (m, k - p) df with KH, chi2_m otherwise.
    Pseudo-R2 = max(0, (tau2_0 - tau2_1) / tau2_0) with tau2_0 from the intercept-only model
    on the same studies and estimator; R2 = 0 when tau2_0 = 0 (there is no heterogeneity
    to explain), as metafor rma.uni reports it, with 0 meaning below pooling.TAU2_ZERO times
    the mean within-study variance.

    Gate: min_k studies per moderator (Handbook v6.5 10.11.5.1: about ten studies for each
    characteristic modelled), so a joint model of 3 moderators with metareg_min_k = 10 needs
    30 studies. A categorical moderator is one characteristic however many levels it has;
    the 'k_per_param' column (studies per non-intercept coefficient) and 'min_level_k'
    (fewest studies in any category) show how thinly the model is supported, and the note
    says when a category has a single study (its coefficient then rests on that study)."""
    name = model_name or "+".join(moderators)
    d, X, names = design_matrix(eff, moderators)
    k, p = X.shape
    k_st = d[study_col].nunique() if study_col in d else k
    n_mod = len(moderators)
    counts = _level_counts(d, moderators, study_col)
    min_level = min((int(c.min()) for c in counts.values() if len(c)), default=NAN)
    base = {"model": name, "k": k, "p": p}
    need = min_k * n_mod
    if k_st < need:
        gate = (f"metareg_min_k = {min_k}" if n_mod == 1
                else f"metareg_min_k x {n_mod} moderators = {min_k} x {n_mod} = {need}")
        return pd.DataFrame([{**base, "note": f"not run: {k_st} studies < {gate} (Handbook "
                                              "10.11.5.1: about ten studies per characteristic; "
                                              "protocol I)"}])
    if method == "FE":
        return pd.DataFrame([{**base, "note": "meta-regression needs a random-effects estimator"}])
    one_level = [mod for mod, c in counts.items() if len(c) < 2]
    if p == 1 or one_level:
        # Same wording as run_meta.R. Fitting would give an intercept-only model labelled as
        # a meta-regression, and QM / m with m = 0 under HKSJ.
        return pd.DataFrame([{**base, "note": "moderator has fewer than 2 distinct values"
                              + (f" ({', '.join(one_level)})" if one_level else "")}])
    if k <= p or np.linalg.matrix_rank(X) < p:
        return pd.DataFrame([{**base, "note": "not estimable: k <= p or collinear moderators"}])
    y, v = d["yi"].to_numpy(float), d["vi"].to_numpy(float)
    f = pl.fit(y, v, X, method, hksj, level, hksj_rule=hksj_rule)
    f0 = pl.fit(y, v, None, method, hksj=False)
    # tau2_0 counts as 0 below the shared TAU2_ZERO fraction of the mean variance: metafor
    # stops a boundary REML fit at tau2 ~ 1e-11 to 1e-6 rather than 0, and dividing by that
    # gives an R2 of 80-100% for a model with nothing to explain.
    R2 = (max(0.0, 100 * (f0.tau2 - f.tau2) / f0.tau2)
          if f0.tau2 > pl.TAU2_ZERO * float(np.mean(v)) else 0.0)
    m = p - 1
    bm = f.beta[1:]
    vb = f.vb[1:, 1:]
    if not np.all(np.isfinite(vb)) or np.linalg.matrix_rank(vb) < m:
        # Under HKSJ the scale is 0 when the model fits every study exactly, so the moderator
        # covariance is 0 and the omnibus test is undefined (as in influence()).
        return pd.DataFrame([{**base, "note": "omnibus test undefined: moderator covariance "
                              "is zero (the model fits every study exactly under HKSJ)"}])
    QM = float(bm @ np.linalg.inv(vb) @ bm)
    if f.hksj:
        QM_stat, QM_df2 = QM / m, f.df
        QM_p = float(stats.f.sf(QM_stat, m, f.df))
    else:
        QM_stat, QM_df2 = QM, NAN
        QM_p = float(stats.chi2.sf(QM, m))
    notes = []
    if k_st < 10 * n_mod:
        notes.append(f"{k_st} studies for {n_mod} moderator(s), fewer than 10 per "
                     "characteristic: meta-regression is unreliable (Handbook 10.11.5.1)")
    sparse = [mod for mod, c in counts.items() if int(c.min()) < 2]
    if sparse:
        notes.append(f"a category of {', '.join(sparse)} has 1 study: its coefficient rests on "
                     "that study alone")
    note = "; ".join(notes)
    perm_p = NAN
    if permutations:
        # Higgins & Thompson 2004: permute the moderator rows against (y, v) and refit;
        # p = (1 + #{QM_perm >= QM}) / (1 + B).
        rng = np.random.default_rng(seed)
        hits = 0
        for _ in range(int(permutations)):
            Xp = X[rng.permutation(k)]
            if np.linalg.matrix_rank(Xp) < p:
                continue
            fp = pl.fit(y, v, Xp, method, hksj, level, hksj_rule=hksj_rule)
            bp = fp.beta[1:]
            hits += float(bp @ np.linalg.inv(fp.vb[1:, 1:]) @ bp) >= QM - 1e-12
        perm_p = (1 + hits) / (1 + int(permutations))
    rows = []
    for j, term in enumerate(names):
        rows.append({**base, "term": term, "estimate": f.beta[j], "se": f.se[j],
                     "stat": f.stat[j], "test": "t" if f.hksj else "z",
                     "df": f.df if f.hksj else NAN, "pval": f.pval[j],
                     "ci_lb": f.ci_lb[j], "ci_ub": f.ci_ub[j], "tau2": f.tau2, "R2": R2,
                     "QM": QM_stat, "QM_test": "F" if f.hksj else "chi2", "QM_df1": m,
                     "QM_df2": QM_df2, "QM_p": QM_p, "QE": f.QE, "QE_df": f.QE_df,
                     "QE_p": f.QE_p, "estimator": method, "perm_p": perm_p,
                     "ci_method": pl.ci_method(f), "n_moderators": n_mod,
                     "k_per_param": k_st / m, "min_level_k": min_level, "note": note})
    return pd.DataFrame(rows)


def egger(y, v, min_k=10, precision_se=None, k_studies=None, measure=None) -> pd.DataFrame:
    """Egger 1997: OLS of y_i / se_i on 1 / se_i; the intercept measures asymmetry,
    tested with t on k - 2 df. Not run when the number of studies is below min_k
    (Sterne 2011; Handbook chapter 13).

    precision_se replaces se_i in the regression. For SMD the pipeline passes
    sqrt(1/n1 + 1/n2), which does not depend on g (Pustejovsky & Rodgers 2019); the classic
    test on SMD is biased because g and its SE are correlated.

    Not estimable (NaN row with a note) when every study has the same precision: the
    intercept and 1/se columns are then collinear, so asymmetry cannot be separated from the
    pooled effect (Sterne 2011: not when studies are of similar size). Also not estimable
    with fewer than 3 studies (no residual df) or when the line fits exactly (residual
    variance 0, so the intercept SE is 0 and t is undefined). metafor regtest() stops with
    an error in the constant-precision case."""
    y, v = np.asarray(y, float), np.asarray(v, float)
    k = len(y)
    k_st = k if k_studies is None else int(k_studies)
    se = np.sqrt(v) if precision_se is None else np.asarray(precision_se, float)
    variant = ("modified, precision sqrt(1/n1 + 1/n2) (Pustejovsky & Rodgers 2019)"
               if precision_se is not None else "classic (Egger 1997)")
    base = {"test": "Egger regression", "k": k, "k_studies": k_st, "variant": variant}
    if k_st < min_k:
        return pd.DataFrame([{**base, "intercept": NAN, "intercept_se": NAN, "t": NAN,
                              "df": NAN, "p": NAN, "slope": NAN,
                              "note": f"not run: {k_st} studies < {min_k}"}])
    def not_estimable(why):
        return pd.DataFrame([{**base, "intercept": NAN, "intercept_se": NAN, "t": NAN,
                              "df": NAN, "p": NAN, "slope": NAN, "se_ratio": ratio,
                              "note": f"not estimable: {why}"}])

    ratio = float(se.max() / se.min()) if k else NAN
    if k < 3:
        return not_estimable(f"{k} studies leave no residual df for the regression")
    z, prec = y / se, 1 / se
    if np.ptp(prec) <= 1e-8 * float(np.max(np.abs(prec))):
        return not_estimable("every study has the same precision, so the intercept and the "
                             "precision slope cannot be separated (constant precision)")
    A = np.column_stack([np.ones(k), prec])
    coef, *_ = np.linalg.lstsq(A, z, rcond=None)
    resid = z - A @ coef
    rss = float(resid @ resid)
    if rss <= (64 * np.finfo(float).eps) ** 2 * float(z @ z):
        return not_estimable("y/se lies exactly on a line in 1/se (residual variance 0), so "
                             "the intercept SE is 0 and the t test is undefined")
    sigma2 = rss / (k - 2)
    cov = sigma2 * np.linalg.inv(A.T @ A)
    se_int = math.sqrt(cov[0, 0])
    t = coef[0] / se_int
    notes = []
    if ratio < 2:
        notes.append(f"SEs span only {ratio:.2f}-fold: studies are of similar size and the "
                     "test has little power (Sterne 2011); interpret with caution")
    if measure in ("SMD", "ROM") and precision_se is None:
        notes.append(f"classic test on {measure}: the estimate and its SE are correlated, "
                     "which can produce asymmetry by itself")
    if k != k_st:
        notes.append(f"{k} comparisons from {k_st} studies: comparisons sharing a control "
                     "are not independent")
    return pd.DataFrame([{**base, "intercept": coef[0], "intercept_se": se_int, "t": t,
                          "df": k - 2, "p": 2 * stats.t.sf(abs(t), k - 2), "slope": coef[1],
                          "se_ratio": ratio, "note": "; ".join(notes)}])


def bh_adjust(p) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values: min over j >= i of p_(j) m / j, capped at 1.
    NaN entries are left as NaN and not counted in m."""
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = len(pv)
    if m == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(m)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out
