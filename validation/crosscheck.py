#!/usr/bin/env python
"""Cross-check the Python engine (metaproof) against the R/metafor reference (run_meta.R).

One command reruns everything and rewrites validation/crosscheck_report.txt:

    .venv/bin/python validation/crosscheck.py            # exit 0 when no check fails
    .venv/bin/python validation/crosscheck.py --keep DIR # also keep both engines' outputs
    .venv/bin/python validation/crosscheck.py --reuse DIR  # compare outputs kept earlier

Sections
  A  reference datasets (dat.bcg, dat.normand1999, dat.curtis1998: real published methods
     datasets from metafor/metadat, not NR data) against the metafor expected values in
     validation/reference_data/:
       A1 Python, called directly, for every estimator x test combination in the expected files
       A2 run_meta.R through its command line
       A3 Python CLI against R CLI, the numeric CSV outputs both engines write
  B  examples/nicotinamide-riboside/SYNTHETIC_nr_like_extraction.csv with its config, Python CLI against R CLI
  C  constructed SYNTHETIC edge cases (written to validation/edge_cases/), Python CLI against
     R CLI, plus hand calculations from the cited formulas for the data-handling paths

What engine-to-engine comparison covers: study_effects, excluded study IDs, pooled.csv,
leave1out, influence, subgroups, metareg (matched by model and term), Egger, the imputed-
correlation sensitivity rows, and the pooled test statistic and p-value under the configured
hksj_rule when run_meta.R writes them. Not covered: plot content (forest.png and funnel.png
are checked for existence only), log.txt, and the Python-only sensitivity scenarios other
than the K5b/K7 exclusions, which are re-pooled here from the flagged studies. Every
comparison that does not run is reported as SKIP with its reason, never as PASS. When both
engines skip an analysis step, the reason is checked against the config and the study count;
a skip whose reason does not hold is a FAIL.

Tolerances: a value passes when |a - b| <= tol * |b| + floor (b = reference value).
  tol (relative)
    1e-6   estimates, SEs and tau2 (yi, vi, sei, pooled and leave-one-out estimates, meta-
           regression and Egger coefficients and SEs)
    1e-4   CIs, prediction intervals, Q-profile tau2 limits, influence measures, and every
           other statistic (I2, H2, Q, test statistics, p-values, R2)
  floor (absolute, by quantity)
    vi, sei        none. Per-study variances are positive and enter as inverse-variance
                   weights, so only their ratio matters; a fixed floor would let a variance in
                   small units (1e-8) be off by a factor of 100 and still pass.
    se             tol * sqrt(vmin), vmin = the smallest reference vi of the analysis. Only
                   matters when a pooled SE is 0 (identical effects under HKSJ).
    tau2           tol * vmin (tau2, tau2_del, Q-profile limits). A tau2 error d changes no
                   weight 1/(vi + tau2) by more than d / vmin, so this bounds the weight error
                   at tol. Without vmin, tol * SE^2 of the row.
    locations      tol * the largest |location| or SE in the same row (yi, estimates, CI and
                   PI limits, coefficients), so a CI limit near zero is judged on the scale of
                   its own estimate and SE. Both floors above and this one scale with the
                   measurement units, so rescaling the data does not change a verdict.
    p-values       1e-30, so p is compared relative down to 1e-30: a tail p from a statistic
                   with relative error e moves by about z^2 e relative (1e-7 at p = 1e-20),
                   far inside tol, while a floor such as 1e-15 would pass a p of 1e-20 that
                   is off by a factor of 1000. Below 1e-30 p is set by rounding in a near-
                   zero SE (identical effects under HKSJ: 6e-33 in Python, 0 in R).
    I2, R2, weight 1e-5 percentage points. At tau2 near 0 the REML stopping rule leaves I2
                   noise of about 2e-7 points (B whole_blood_nad leave-one-out); I2 is
                   reported to 0.1 points.
    other unit-free statistics  1e-8 (Q, QM, test statistics, H2, rstudent, dffits, dfbetas,
                   cov_r, Cook's d, hat, Egger intercept and SE). tau2 convergence noise gives
                   about 1e-9 in dffits; identical effects give Q = 4e-16 against 0.
  exact    k, degrees of freedom, influential flags, excluded-study sets
Both engines missing a value (NA/NaN) counts as agreement; one side missing is a mismatch, and
so is a value that is not a number. A key that appears twice on either side is a mismatch.
"""

from __future__ import annotations

import argparse
import math
import platform
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from metaproof import diagnostics as dg  # noqa: E402
from metaproof import pooling as pl  # noqa: E402
from metaproof.config import load_config  # noqa: E402

REF = ROOT / "validation" / "reference_data"
EDGE = ROOT / "validation" / "edge_cases"
REPORT = ROOT / "validation" / "crosscheck_report.txt"
RUN_META = ROOT / "reference_r" / "run_meta.R"
PY = sys.executable
REUSE = False            # --reuse: compare outputs already on disk instead of running engines

TOL_EST, TOL_CI = 1e-6, 1e-4
EST_COLS = {"estimate", "se", "tau2", "yi", "vi", "sei", "intercept", "intercept_se",
            "coef_sei", "coef_se", "slope", "limit_estimate"}
EXACT_COLS = {"k", "df", "QM_df1", "QM_df2", "test_df1", "test_df2", "influential", "k_studies"}
STRING_COLS = {"ci_method"}

# Quantity classes for the absolute floor (module docstring). Classified by the reference
# column name; a column in no class is reported as a failure, not compared loosely.
SE_COLS = {"se"}                   # SE of a pooled estimate or coefficient
STUDY_COLS = {"vi", "sei"}         # per-study sampling variance and SE
TAU2_COLS = {"tau2", "tau2_ci_lb", "tau2_ci_ub", "tau2_del"}
LOC_COLS = {"yi", "estimate", "ci_lb", "ci_ub", "pi_lb", "pi_ub", "pi_lb_metafor",
            "pi_ub_metafor", "wald_ci_lb", "wald_ci_ub", "hksj_ci_lb", "hksj_ci_ub",
            "hksj_trunc_ci_lb", "hksj_trunc_ci_ub", "intercept", "slope", "limit_estimate",
            "limit_ci_lb", "limit_ci_ub"}
P_COLS = {"pval", "p", "Q_p", "QM_p", "QE_p", "test_p", "Q_between_p", "perm_p"}
PCT_COLS = {"I2", "R2", "weight", "weight_pct"}
STAT_COLS = {"Q", "QE", "QM", "Q_between", "test_stat", "statistic", "stat", "t", "zval", "H2",
             "QE_del", "rstudent", "dffits", "dfbetas", "cov_r", "cooks_d", "hat", "coef_sei",
             "coef_se", "intercept_se", "geary_min"}
FLOOR_P, FLOOR_PCT, FLOOR_STAT = 1e-30, 1e-5, 1e-8
NA_STRINGS = {"", "NA", "NAN", "NONE", "NULL", "<NA>"}


# ---------------------------------------------------------------------------------------
# Comparison bookkeeping
# ---------------------------------------------------------------------------------------

class Ledger:
    def __init__(self):
        self.lines: list[str] = []
        self.n_pass = self.n_fail = self.n_skip = self.n_values = 0
        self.failures: list[str] = []
        self.skipped: list[str] = []

    def section(self, title):
        self.lines += ["", title, "-" * len(title)]

    def note(self, text):
        self.lines.append(f"      {text}")

    def record(self, name, n, max_diff, problems):
        self.n_values += n
        status = "PASS" if not problems else "FAIL"
        if problems:
            self.n_fail += 1
            self.failures.append(name)
        else:
            self.n_pass += 1
        md = "-" if max_diff is None else f"{max_diff:.1e}"
        self.lines.append(f"{status}  {name:<66} {n:>5} values  max diff/allowed {md}")
        for p in problems[:12]:
            self.lines.append(f"        {p}")
        if len(problems) > 12:
            self.lines.append(f"        ... {len(problems) - 12} more")

    def skip(self, name, reason):
        """A comparison that did not run. Never counted as a pass."""
        self.n_skip += 1
        self.skipped.append(name)
        self.lines.append(f"SKIP  {name:<66}     0 values  {reason}")


def _num(x):
    """(value, ok): NA-like inputs give (nan, True); text that is not a number gives
    (nan, False), so two malformed cells never agree as two missing values."""
    if x is None:
        return math.nan, True
    if isinstance(x, (bool, np.bool_)):
        return float(x), True
    if isinstance(x, str):
        s = x.strip().upper()
        if s in NA_STRINGS:
            return math.nan, True
        if s in ("TRUE", "FALSE"):
            return (1.0 if s == "TRUE" else 0.0), True
        try:
            return float(s), True
        except ValueError:
            return math.nan, False
    try:
        return float(x), True
    except (TypeError, ValueError):
        return math.nan, False


def tol_for(col):
    if col in EXACT_COLS:
        return 0.0
    return TOL_EST if col in EST_COLS else TOL_CI


def kind_of(col):
    for kind, cols in (("exact", EXACT_COLS), ("se", SE_COLS), ("study", STUDY_COLS),
                       ("tau2", TAU2_COLS), ("loc", LOC_COLS), ("p", P_COLS),
                       ("pct", PCT_COLS), ("stat", STAT_COLS)):
        if col in cols:
            return kind
    return None


def allowed(col, b, se_scale, loc_scale, var_scale=None):
    """Largest |a - b| accepted for reference value b (module docstring). var_scale is the
    smallest sampling variance vi of the reference study effects, when the caller has it."""
    tol, kind = tol_for(col), kind_of(col)
    if kind == "exact":
        return 0.0
    if kind == "study":
        return tol * abs(b)
    if kind == "se":
        return tol * (abs(b) + (math.sqrt(var_scale) if var_scale else 0.0))
    if kind == "tau2":
        return tol * (abs(b) + (var_scale if var_scale else se_scale ** 2))
    if kind == "loc":
        return tol * max(abs(b), loc_scale)
    floor = {"p": FLOOR_P, "pct": FLOOR_PCT, "stat": FLOOR_STAT}[kind]
    return tol * abs(b) + floor


def _row_scales(vals):
    """SE scale and location scale of one reference row, from its finite values."""
    def fin(cols):
        return [abs(v) for c, v in vals.items() if c in cols and math.isfinite(v)]
    se_scale = max(fin({"se", "sei"}) + [math.sqrt(v) for v in fin({"vi"})], default=0.0)
    return se_scale, max(fin(LOC_COLS) + [se_scale])


def _keys(df, key):
    cols = [key] if isinstance(key, str) else list(key)
    return ["/".join(_kstr(v) for v in vals) for vals in zip(*(df[c] for c in cols))]


def _kstr(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "NA"
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


def compare(led, name, left, right, lkey, rkey, cols, *, by_position=False, var_scale=None):
    """cols: list of (left_col, right_col); tolerance and floor are chosen by right_col.
    var_scale: smallest reference vi of the analysis, the scale for tau2 and pooled-SE floors.

    Columns are resolved on each side before the merge, so a column missing on one side is
    reported instead of being compared with the other side's copy, and keys must be unique
    on both sides (validate = one_to_one)."""
    if left is None or right is None:
        led.record(name, 0, None, ["missing output: " + ("python" if left is None else "R")])
        return
    problems, n, worst = [], 0, 0.0
    for side, df, key in (("Python", left, lkey), ("R/reference", right, rkey)):
        if not by_position:
            for c in [key] if isinstance(key, str) else list(key):
                if c not in df.columns:
                    problems.append(f"key column {c!r} missing in {side} output")
    for lc, rc in cols:
        if lc not in left.columns:
            problems.append(f"column {lc!r} missing in Python output")
        if rc not in right.columns:
            problems.append(f"column {rc!r} missing in R/reference output")
        if rc not in STRING_COLS and kind_of(rc) is None:
            problems.append(f"column {rc!r} has no tolerance class")
    if problems:
        led.record(name, 0, None, problems)
        return
    L = pd.DataFrame({"_key": range(len(left)) if by_position else _keys(left, lkey)})
    R = pd.DataFrame({"_key": range(len(right)) if by_position else _keys(right, rkey)})
    for i, (lc, rc) in enumerate(cols):
        L[f"_l{i}"] = left[lc].to_numpy()
        R[f"_r{i}"] = right[rc].to_numpy()
    for side, df in (("Python", L), ("R/reference", R)):
        dup = df["_key"][df["_key"].duplicated()]
        for k in sorted(set(dup)):
            problems.append(f"row {k!r} appears {int((df['_key'] == k).sum())} times in "
                            f"{side} output")
    L, R = L.drop_duplicates("_key"), R.drop_duplicates("_key")
    if len(L) == 0 and len(R) == 0 and not problems:
        led.skip(name, "no rows in either output")
        return
    lk, rk = set(L["_key"]), set(R["_key"])
    for k in sorted(lk - rk):
        problems.append(f"row {k!r} only in Python output")
    for k in sorted(rk - lk):
        problems.append(f"row {k!r} only in R/reference output")
    M = L.merge(R, on="_key", validate="one_to_one")
    for _, m in M.iterrows():
        ref = {rc: _num(m[f"_r{i}"])[0] for i, (_, rc) in enumerate(cols)
               if rc not in STRING_COLS}
        se_scale, loc_scale = _row_scales(ref)
        key = m["_key"]
        for i, (lc, rc) in enumerate(cols):
            a, b = m[f"_l{i}"], m[f"_r{i}"]
            n += 1
            if rc in STRING_COLS:
                if str(a) != str(b):
                    problems.append(f"{key}: {rc} python={a!r} R={b!r}")
                continue
            (a, oka), (b, okb) = _num(a), _num(b)
            if not (oka and okb):
                problems.append(f"{key}: {rc} is not a number: python={m[f'_l{i}']!r} "
                                f"R={m[f'_r{i}']!r}")
                continue
            if math.isnan(a) and math.isnan(b):
                continue
            if math.isnan(a) != math.isnan(b) or ((math.isinf(a) or math.isinf(b)) and a != b):
                problems.append(f"{key}: {rc} python={a!r} R={b!r}")
                continue
            if math.isinf(a) and a == b:
                continue
            d, allow = abs(a - b), allowed(rc, b, se_scale, loc_scale, var_scale)
            ratio = 0.0 if d == 0 else (d / allow if allow > 0 else math.inf)
            worst = max(worst, ratio)
            if ratio > 1:
                problems.append(f"{key}: {rc} python={a:.12g} R={b:.12g} |diff| {d:.2e} > "
                                f"allowed {allow:.2e}")
    led.record(name, n, worst, problems)


def read(path):
    p = Path(path)
    if not p.exists():
        return None
    try:
        return pd.read_csv(p)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def has_rows(df, col):
    return df is not None and col in df.columns and df[col].notna().any()


# ---------------------------------------------------------------------------------------
# Running the engines
# ---------------------------------------------------------------------------------------

def rscript():
    for c in (shutil.which("Rscript"), "/opt/homebrew/bin/Rscript", "/usr/local/bin/Rscript"):
        if c and Path(c).exists():
            return c
    raise SystemExit("Rscript not found")


def run_python(csv, config, out):
    if REUSE:
        return None
    p = subprocess.run([PY, "-m", "metaproof", "run", str(csv), "--config", str(config),
                        "--out", str(out)], cwd=ROOT, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"metaproof failed on {csv}:\n{p.stderr[-3000:]}")
    return p


def run_r(csv, config, out):
    if REUSE:
        return None
    p = subprocess.run([rscript(), str(RUN_META), "--data", str(csv), "--config", str(config),
                        "--out", str(out)], cwd=ROOT, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"run_meta.R failed on {csv}:\n{p.stderr[-3000:]}")
    return p


def r_versions():
    code = 'cat(R.version.string, "| metafor", as.character(packageVersion("metafor")))'
    lib = ROOT / "reference_r" / "lib"
    p = subprocess.run([rscript(), "-e", f'.libPaths(c("{lib}", .libPaths())); {code}'],
                       capture_output=True, text=True)
    return p.stdout.strip() or p.stderr.strip()


def lock_status():
    """Whether the running Python environment and reference_r/lib match the versions the kit
    was verified with (requirements-lock.txt, reference_r/R_PACKAGES_LOCK.txt): Astra 22."""
    import importlib.metadata as md
    import re as _re
    norm = lambda n: _re.sub(r"[-_.]+", "-", n).lower()
    want = {}
    for n, v, marker in _re.findall(r"^([A-Za-z0-9_.-]+)==([^ \\;]+)\s*(?:;([^\\]*))?",
                                    (ROOT / "requirements-lock.txt").read_text(), _re.M):
        if "win32" in (marker or "") and sys.platform != "win32":
            continue                      # Windows-only dependency (colorama for pytest)
        want[n] = v
    have = {norm(d.metadata["Name"]): d.version for d in md.distributions()}
    py = [f"{n} {v} (lock {want[n]})" if n in want else f"{n} {v} (not in lock)"
          for n, v in sorted(have.items()) if want.get(n) != v]
    py += [f"{n} missing (lock {v})" for n, v in sorted(want.items()) if norm(n) not in have]
    rl = (ROOT / "reference_r" / "R_PACKAGES_LOCK.txt").read_text()
    rwant = dict(_re.findall(r"^(\w+)\s+(\d[\w.-]*)\s+\S", rl.split("Package")[1], _re.M))
    code = ('ip <- installed.packages(lib.loc = "%s")[, "Version"]; '
            'cat(paste(names(ip), ip), sep = "\\n")' % (ROOT / "reference_r" / "lib"))
    pr = subprocess.run([rscript(), "-e", code], capture_output=True, text=True)
    rhave = dict(line.split(" ", 1) for line in pr.stdout.splitlines() if " " in line)
    rr = [f"R {n} {rhave.get(n, 'missing')} (lock {v})" for n, v in sorted(rwant.items())
          if rhave.get(n) != v]
    diff = py + rr
    return ("environment matches requirements-lock.txt and reference_r/R_PACKAGES_LOCK.txt"
            if not diff else "environment differs from the locks: " + "; ".join(diff))


# ---------------------------------------------------------------------------------------
# Engine-to-engine comparison of one analysis folder
# ---------------------------------------------------------------------------------------

POOLED = [(c, c) for c in pl.POOLED_COLUMNS if c != "estimator"]
POOLED_EXTRA = [(c, c) for c in ["k_studies"] + pl.POOLED_EXTRA_COLUMNS]
LOO = [(c, c) for c in ("k", "estimate", "se", "ci_lb", "ci_ub", "pval", "tau2", "I2", "H2",
                        "Q", "Q_p")]
INFL = [("rstudent", "rstudent"), ("dffits", "dffits"), ("cook_d", "cooks_d"),
        ("cov_ratio", "cov_r"), ("tau2_del", "tau2_del"), ("QE_del", "QE_del"), ("hat", "hat"),
        ("weight_pct", "weight"), ("dfbetas", "dfbetas"), ("influential", "influential")]
SUBG = [(c, c) for c in ("k", "estimate", "se", "ci_lb", "ci_ub", "pi_lb", "pi_ub", "tau2",
                         "I2", "Q", "Q_p")]
MREG = [("k", "k"), ("estimate", "estimate"), ("se", "se"), ("stat", "statistic"),
        ("df", "df"), ("pval", "pval"), ("ci_lb", "ci_lb"), ("ci_ub", "ci_ub"),
        ("tau2", "tau2"), ("R2", "R2"), ("QM", "QM"), ("QM_df1", "QM_df1"),
        ("QM_df2", "QM_df2"), ("QM_p", "QM_p"), ("QE", "QE"), ("QE_p", "QE_p")]
EGGER = [("k", "k"), ("intercept", "coef_sei"), ("intercept_se", "coef_se"),
         ("t", "statistic"), ("df", "df"), ("p", "pval"), ("slope", "limit_estimate")]
SENS_PARAMS = ("corr_prepost", "corr_within")   # the sensitivity rows both engines write
SENS_COLS = [(c, c) for c in ("k", "estimate", "se", "ci_lb", "ci_ub", "tau2", "pi_lb",
                              "pi_ub")]
# protocol K7: drops imputed change SDs, assumed crossover r, digitised values and
# range-only SDs; K5b: drops effects converted from medians
K7_FLAGS = ("sd_change_imputed", "corr_within_assumed", "digitised", "sd_from_range_only")
K5B_FLAGS = ("converted_from_median",)
OUTPUT_FILES = ("study_effects.csv", "excluded.csv", "pooled.csv", "leave1out.csv",
                "influence.csv", "subgroups.csv", "metareg.csv", "small_study.csv",
                "sensitivity.csv", "forest.png", "funnel.png")


def min_vi(df):
    """Smallest positive finite vi of a study-effects table (the reference side), or None."""
    if df is None or "vi" not in df.columns:
        return None
    v = pd.to_numeric(df["vi"], errors="coerce")
    v = v[np.isfinite(v) & (v > 0)]
    return float(v.min()) if len(v) else None


def _flag(df, cols):
    """Rows where any of cols is true; None if a column is missing."""
    if df is None or any(c not in df.columns for c in cols):
        return None
    out = np.zeros(len(df), dtype=bool)
    for c in cols:
        out |= np.array([_num(x)[0] == 1.0 for x in df[c]])
    return out


def _n_studies(eff, col=None):
    """Distinct source studies, optionally among rows where col is present."""
    if eff is None or not len(eff):
        return 0
    d = eff if col is None else eff[eff[col].notna()]
    sid = "source_study" if "source_study" in d.columns else "study_id"
    return int(d[sid].astype(str).str.split("|").str[0].nunique())


def _notes(df, by="moderator"):
    """{model or moderator: note} for the note rows of a subgroups/metareg file."""
    if df is None or "note" not in df.columns:
        return {}
    out = {}
    for _, r in df.iterrows():
        if isinstance(r.note, str) and r.note.strip():
            out[str(r.get(by, "")) if by in df.columns else ""] = r.note.strip()
    return out


def _metareg_estimable(eff, mods):
    """True when an intercept + mods design has full column rank and k > p."""
    d = eff.dropna(subset=list(mods))
    if not len(d):
        return False
    parts = [np.ones((len(d), 1))]
    for m in mods:
        col = pd.to_numeric(d[m], errors="coerce")
        if col.notna().all():
            parts.append(col.to_numpy(float)[:, None])
        else:
            parts.append(pd.get_dummies(d[m].astype(str), drop_first=True).to_numpy(float))
    X = np.hstack(parts)
    return len(d) > X.shape[1] and np.linalg.matrix_rank(X) == X.shape[1]


def compare_engines(led, tag, py_dir: Path, r_dir: Path, analysis: str, spec: dict):
    """Compare one analysis folder of the two engines. spec is the resolved config of the
    analysis (metaproof.config.load_config); it is used to check the reason for every step
    that both engines skip."""
    P, R = py_dir / analysis, r_dir / analysis

    def nm(s):
        return f"{tag} {analysis}: {s}"

    if not P.exists() or not R.exists():
        led.record(nm("output folders"), 0, None,
                   [f"missing: python={P.exists()} R={R.exists()}"])
        return
    for f in OUTPUT_FILES:
        if not (P / f).exists() or not (R / f).exists():
            led.record(nm(f"{f} written"), 0, None,
                       [f"python={(P / f).exists()} R={(R / f).exists()}"])
    pe, re_ = read(P / "study_effects.csv"), read(R / "study_effects.csv")

    # excluded studies (compare the source study IDs)
    ex_p, ex_r = read(P / "excluded.csv"), read(R / "excluded.csv")
    sp = set() if ex_p is None or not len(ex_p) else {str(s).split("|")[0] for s in ex_p.study_id}
    sr = set() if ex_r is None or not len(ex_r) else {str(s).split("|")[0] for s in ex_r.study_id}
    led.record(nm("excluded studies"), max(1, len(sp | sr)), None,
               [] if sp == sr else [f"python={sorted(sp)} R={sorted(sr)}"])

    if pe is not None and re_ is not None and len(pe) == 0 and len(re_) == 0:
        led.record(nm("no analysable studies in either engine"), 1, None, [])
        led.skip(nm("pooled and later outputs"), "no analysable studies")
        return
    compare(led, nm("study_effects"), pe, re_, "study_id", "study_id",
            [("yi", "yi"), ("vi", "vi"), ("sei", "sei")])
    if pe is None or re_ is None:
        return                      # recorded as a failure just above
    if spec["measure"] == "ROM":    # Geary's check: R writes geary_min to study_inputs.csv
        compare(led, nm("Geary min mean sqrt(n)/SD"), pe, read(R / "study_inputs.csv"),
                "study_id", "study_id", [("geary_min", "geary_min")])
    vs = min_vi(re_)

    po_p, po_r = read(P / "pooled.csv"), read(R / "pooled.csv")
    compare(led, nm("pooled (all estimators)"), po_p, po_r, "estimator", "estimator", POOLED,
            var_scale=vs)
    compare(led, nm("pooled Wald/HKSJ/truncated-HKSJ intervals and rule"), po_p, po_r,
            "estimator", "estimator", POOLED_EXTRA, var_scale=vs)
    if po_p is None or not len(po_p):
        led.record(nm("pooled.csv has rows (python)"), 0, None, ["no rows"])
        return
    prim = po_p.iloc[0]
    k, k_st = int(prim["k"]), int(prim["k_studies"])
    pooled_ok = bool(pd.notna(prim["estimate"]))
    method = str(prim.get("ci_method", ""))

    _pooled_test(led, nm, P, R, po_r, pooled_ok, spec)
    _leave1out(led, nm, P, R, pooled_ok, k_st, vs)
    _influence(led, nm, P, R, pooled_ok, k, method, _num(prim.get("Q"))[0], vs)
    _subgroups(led, nm, P, R, spec, pe, pooled_ok, vs)
    _metareg(led, nm, P, R, spec, pe, vs)
    _egger(led, nm, P, R, spec, k_st)
    _sensitivity(led, nm, P, R, spec, pe, pooled_ok, vs)


def _pooled_test(led, nm, P, R, po_r, pooled_ok, spec):
    """Test statistic and p-value of each estimator under the configured hksj_rule (the p
    reported and fed to Benjamini-Hochberg). Python writes them to pooled_tests.csv."""
    name = nm(f"pooled test statistic and p (hksj_rule {spec['hksj_rule']})")
    if not pooled_ok:
        led.skip(name, "not pooled (k_studies < min_k_pool, checked in pooled.csv)")
        return
    stat_r = next((c for c in ("stat", "statistic", "zval") if po_r is not None
                   and c in po_r.columns), None)
    if po_r is None or "pval" not in po_r.columns or stat_r is None:
        led.skip(name, "run_meta.R pooled.csv has no test statistic/pval columns; A1 checks "
                       "the p-value against metafor for hksj_rule always only")
        return
    compare(led, name, read(P / "pooled_tests.csv"), po_r, "estimator", "estimator",
            [("stat", stat_r), ("pval", "pval")])


def _leave1out(led, nm, P, R, pooled_ok, k_st, vs):
    lp, lr = read(P / "leave1out.csv"), read(R / "leave1out.csv")
    if has_rows(lp, "estimate") or has_rows(lr, "estimate"):
        compare(led, nm("leave1out"), lp, lr, "omitted", "study_id", LOO, var_scale=vs)
    elif not pooled_ok or k_st < 3:
        led.skip(nm("leave1out"), f"not run by either engine; checked: k_studies = {k_st}, "
                                  f"pooled = {pooled_ok} (runs only when pooled and k >= 3)")
    else:
        led.record(nm("leave1out"), 0, None,
                   [f"neither engine ran it, but k_studies = {k_st} >= 3 and pooled"])


def _influence(led, nm, P, R, pooled_ok, k, method, Q, vs):
    ip, ir = read(P / "influence.csv"), read(R / "influence.csv")
    if has_rows(ip, "rstudent") or has_rows(ir, "rstudent"):
        compare(led, nm("influence"), ip, ir, "study_id", "study_id", INFL, var_scale=vs)
        return
    s2_zero = method.startswith("HKSJ") and "trunc" not in method.lower() \
        and math.isfinite(Q) and abs(Q) <= FLOOR_STAT
    if not pooled_ok or k < 3 or s2_zero:
        why = ("HKSJ scale s2 = 0 (Q = 0)" if s2_zero and pooled_ok and k >= 3
               else f"k = {k}, pooled = {pooled_ok} (runs only when pooled and k >= 3)")
        led.skip(nm("influence"), f"not computed by either engine; checked: {why}")
    else:
        led.record(nm("influence"), 0, None,
                   [f"neither engine computed it, but k = {k} >= 3, pooled, Q = {Q:.3g}"])


def _subgroups(led, nm, P, R, spec, pe, pooled_ok, vs):
    sp_, sr_ = read(P / "subgroups.csv"), read(R / "subgroups.csv")
    if has_rows(sp_, "estimate") or has_rows(sr_, "estimate"):
        lv_r = sr_[sr_.level.astype(str) != "ALL"] if sr_ is not None else None
        sp_lv = sp_[sp_.estimate.notna()] if sp_ is not None and "estimate" in sp_ else sp_
        lv_r = lv_r[lv_r.estimate.notna()] if lv_r is not None else None
        compare(led, nm("subgroups (per level)"), sp_lv, lv_r,
                ("moderator", "level"), ("moderator", "level"), SUBG, var_scale=vs)
        common = "tau2_mode" in sp_ and (sp_.tau2_mode == "common").any()
        qcols = ["moderator", "Q_between", "Q_between_df", "Q_between_p"] + (
            ["Q_between_df2"] if common else [])
        qb_p = sp_lv.groupby("moderator", as_index=False).first()[qcols]
        want = "QM_common_tau2" if common else "Q_between_separate_tau2"
        qb_r = sr_[sr_.test.astype(str) == want] if sr_ is not None else None
        pairs = [("Q_between", "test_stat"), ("Q_between_df", "test_df1"),
                 ("Q_between_p", "test_p")] + ([("Q_between_df2", "test_df2")] if common else [])
        compare(led, nm(f"subgroups test for differences ({want})"), qb_p, qb_r,
                "moderator", "moderator", pairs)
        return
    if not spec["subgroups"]:
        led.skip(nm("subgroups"), "no subgroups configured (checked in the config)")
        return
    np_, nr_ = _notes(sp_), _notes(sr_)
    problems, why = [], []
    for m in spec["subgroups"]:
        a, b = np_.get(str(m), ""), nr_.get(str(m), "")
        if not (a.startswith("not run") and b.startswith("not run")):
            problems.append(f"{m}: no estimates, notes python={a!r} R={b!r}")
            continue
        d = pe[pe[m].notna()] if pe is not None and m in pe.columns else None
        kst = _n_studies(d)
        counts = pd.Series(dtype=float)
        if d is not None and len(d):
            sid = d["source_study"] if "source_study" in d.columns else d["study_id"]
            counts = d.assign(_s=sid.astype(str)).groupby(m)["_s"].nunique()
        lvl_min = int(counts.min()) if len(counts) else 0
        if not pooled_ok or kst < spec["subgroup_min_k"] or lvl_min < spec[
                "subgroup_min_per_level"]:
            why.append(f"{m}: {kst} studies (subgroup_min_k {spec['subgroup_min_k']}), "
                       f"smallest level {lvl_min} (min {spec['subgroup_min_per_level']})")
        else:
            problems.append(f"{m}: both engines skipped, but {kst} studies >= "
                            f"{spec['subgroup_min_k']} and every level has >= "
                            f"{spec['subgroup_min_per_level']}")
    if problems:
        led.record(nm("subgroups"), 0, None, problems)
    else:
        led.skip(nm("subgroups"), "not run by either engine (protocol gate); checked: "
                 + "; ".join(why))


def _metareg(led, nm, P, R, spec, pe, vs):
    mp, mr = read(P / "metareg.csv"), read(R / "metareg.csv")
    if has_rows(mp, "estimate") or has_rows(mr, "estimate"):
        mp2 = mp[mp.estimate.notna()] if has_rows(mp, "estimate") else mp
        mr2 = mr[mr.estimate.notna()] if has_rows(mr, "estimate") else mr
        compare(led, nm("metareg"), mp2, mr2, ("model", "term"), ("model", "term"), MREG,
                var_scale=vs)
        return
    mods = list(spec["moderators"])
    if not mods:
        led.skip(nm("metareg"), "no moderators configured (checked in the config)")
        return
    models = {str(m): [m] for m in mods}
    if spec.get("metareg_joint") and len(mods) > 1:
        models["joint"] = mods
    np_, nr_ = _notes(mp, "model"), _notes(mr, "model")
    problems, why = [], []
    for model, ms in models.items():
        a = np_.get(model, np_.get("", ""))
        b = nr_.get(model, nr_.get("", ""))
        if not a or not b:
            problems.append(f"{model}: no estimates, notes python={a!r} R={b!r}")
            continue
        gate_a, gate_b = a.startswith("not run"), b.startswith("not run")
        d = pe.dropna(subset=[m for m in ms if m in pe.columns]) if pe is not None else None
        kst = _n_studies(d)
        need = spec["metareg_min_k"] * len(ms)
        if gate_a != gate_b:
            problems.append(f"{model}: different reasons, python={a!r} R={b!r}")
        elif gate_a and kst < need:
            why.append(f"{model}: {kst} studies < {need}")
        elif not gate_a and (any(m not in pe.columns for m in ms)
                             or not _metareg_estimable(pe, ms)
                             or spec["estimator"] == "FE"):
            why.append(f"{model}: not estimable (k <= p, collinear or constant moderator)")
        else:
            problems.append(f"{model}: both engines skipped ({a!r}), but {kst} studies and "
                            "the design has full rank")
    if problems:
        led.record(nm("metareg"), 0, None, problems)
    else:
        led.skip(nm("metareg"), "not fitted by either engine; checked: " + "; ".join(why))


def _egger(led, nm, P, R, spec, k_st):
    ep, er = read(P / "small_study.csv"), read(R / "small_study.csv")
    er_e = er[er.test == "egger"] if er is not None and "test" in er.columns else er
    if has_rows(ep, "p") or has_rows(er_e, "pval"):
        compare(led, nm("Egger (k >= egger_min_k)"), ep, er_e, None, None,
                EGGER + [("k_studies", "k_studies")], by_position=True)
    elif k_st < spec["egger_min_k"]:
        led.skip(nm("Egger"), f"not run by either engine; checked: k_studies = {k_st} < "
                              f"egger_min_k = {spec['egger_min_k']}")
    else:
        led.record(nm("Egger"), 0, None, [f"neither engine ran it, but k_studies = {k_st} "
                                          f">= egger_min_k = {spec['egger_min_k']}"])


def _sensitivity(led, nm, P, R, spec, pe, pooled_ok, vs):
    sens_p, sens_r = read(P / "sensitivity.csv"), read(R / "sensitivity.csv")
    if sens_p is None or sens_r is None:
        led.record(nm("sensitivity.csv"), 0, None,
                   [f"missing: python={sens_p is not None} R={sens_r is not None}"])
        return
    py_rows = sens_p if "parameter" in sens_p.columns else sens_p.iloc[0:0].assign(parameter=[])
    corr_p = py_rows[py_rows.parameter.isin(SENS_PARAMS)].copy()
    r_rows = sens_r if "analysis" in sens_r.columns else pd.DataFrame(columns=["analysis"])
    corr_r = r_rows[r_rows.analysis.astype(str) != "primary"]
    expected = pooled_ok and bool(spec["corr_prepost_sensitivity"]
                                  or spec["corr_within_sensitivity"])
    if len(corr_p) or len(corr_r):
        rows = []
        prim = r_rows[r_rows.analysis.astype(str) == "primary"]
        for _, r in corr_r.iterrows():
            par, val = str(r.analysis).split("=")
            rows.append({**r.to_dict(), "_k": f"{par}={float(val):g}"})
        # the primary row stands for the configured value unless R also wrote that value
        for par in SENS_PARAMS:
            if (corr_p.parameter == par).any() and len(prim):
                p0 = prim.iloc[0]
                key = f"{par}={float(p0[par]):g}"
                if key not in {r["_k"] for r in rows}:
                    rows.append({**p0.to_dict(), "_k": key})
        sr2 = pd.DataFrame(rows, columns=list(r_rows.columns) + ["_k"])
        corr_p["_k"] = [f"{a}={float(b):g}" for a, b in zip(corr_p.parameter, corr_p.value)]
        compare(led, nm("sensitivity (imputed correlations)"), corr_p, sr2, "_k", "_k",
                SENS_COLS, var_scale=vs)
    elif expected:
        led.record(nm("sensitivity (imputed correlations)"), 0, None,
                   ["configured and pooled, but neither engine wrote correlation rows"])
    else:
        led.skip(nm("sensitivity (imputed correlations)"),
                 f"none in either engine; checked: pooled = {pooled_ok}, correlation "
                 "sensitivity lists empty" if pooled_ok else
                 "none in either engine; checked: not pooled")
    for _, r in py_rows[~py_rows.parameter.isin(SENS_PARAMS)].iterrows():
        par = str(r.parameter)
        flags = {"exclude_imputed_or_digitised": K7_FLAGS,
                 "exclude_converted_from_median": K5B_FLAGS}.get(par)
        name = nm(f"sensitivity {par} (Python only)")
        if flags is None:
            led.skip(name, "no R counterpart and no independent re-pool")
            continue
        drop = _flag(pe, flags)
        if drop is None:
            led.record(name, 0, None, [f"study_effects.csv lacks one of {list(flags)}"])
            continue
        sub = pe[~drop]
        exp = {c: math.nan for c, _ in SENS_COLS}
        exp["k"] = len(sub)
        if _n_studies(sub) >= spec["min_k_pool"] and len(sub):
            row, _, _ = pl.pooled_row(sub.yi.to_numpy(float), sub.vi.to_numpy(float),
                                      str(r.estimator), spec["hksj"], spec["level"],
                                      spec["prediction_interval"], spec["pi_type"],
                                      spec["hksj_rule"])
            exp = {c: row[c] for c, _ in SENS_COLS}
        compare(led, nm(f"sensitivity {par}: re-pooled here from the unflagged studies"),
                pd.DataFrame([r]), pd.DataFrame([exp]), None, None, SENS_COLS,
                by_position=True, var_scale=min_vi(sub))


# ---------------------------------------------------------------------------------------
# Section A: reference datasets against the metafor expected values
# ---------------------------------------------------------------------------------------

DATASETS = {
    "bcg": ("bcg_logrr_extraction.csv", "bcg_config.json", ["tb_log_rr"]),
    "normand1999": ("normand1999_extraction.csv", "normand1999_config.json",
                    ["los_days_md", "los_days_smd"]),
    "curtis1998": ("curtis1998_extraction.csv", "curtis1998_config.json", ["biomass_rom"]),
}


def section_a(led, work: Path):
    for ds, (csv, cfg, outcomes) in DATASETS.items():
        py_out, r_out = work / "A" / ds / "py", work / "A" / ds / "r"
        run_python(REF / csv, REF / cfg, py_out)
        run_r(REF / csv, REF / cfg, r_out)
        exp_pooled = read(REF / f"{ds}_expected_pooled.csv")
        exp_se = read(REF / f"{ds}_expected_study_effects.csv")
        if ds == "bcg":
            yv = read(REF / "bcg_logrr_yi_vi.csv")
            exp_se = yv.assign(outcome="tb_log_rr")

        led.section(f"A1/A2 {ds}: Python (direct calls) and run_meta.R against metafor "
                    "expected values")
        for oc in outcomes:
            eff = read(py_out / oc / "study_effects.csv")
            reff = read(r_out / oc / "study_effects.csv")
            exp_oc = exp_se[exp_se.outcome == oc]
            vs = min_vi(exp_oc)
            compare(led, f"A1 {ds} {oc}: study effects (python)", eff, exp_oc,
                    "study_id", "study_id", [("yi", "yi"), ("vi", "vi"), ("sei", "sei")])
            compare(led, f"A2 {ds} {oc}: study effects (R engine)", reff, exp_oc,
                    "study_id", "study_id", [("yi", "yi"), ("vi", "vi"), ("sei", "sei")])
            y, v = eff.yi.to_numpy(float), eff.vi.to_numpy(float)
            ids = list(eff.study_id)

            # pooled, every estimator x test
            ep = exp_pooled[exp_pooled.outcome == oc]
            rows = []
            for _, e in ep.iterrows():
                hk = e.test == "knha"
                row, test, f = pl.pooled_row(y, v, e.estimator, hk, 0.95, True, "riley")
                rowm, _, _ = pl.pooled_row(y, v, e.estimator, hk, 0.95, True, "metafor")
                rows.append({**row, "test": e.test, "pval": test["pval"],
                             "pi_lb_metafor": rowm["pi_lb"], "pi_ub_metafor": rowm["pi_ub"]})
            got = pd.DataFrame(rows)
            for (est, tst), g in got.groupby(["estimator", "test"], sort=False):
                compare(led, f"A1 {ds} {oc}: pooled {est}/{tst}", g,
                        ep[(ep.estimator == est) & (ep.test == tst)], "estimator", "estimator",
                        POOLED + [("pval", "pval"), ("pi_lb_metafor", "pi_lb_metafor"),
                                  ("pi_ub_metafor", "pi_ub_metafor")], var_scale=vs)
            # R engine pooled.csv: config test (knha) for random effects, z for FE
            rp = read(r_out / oc / "pooled.csv")
            ep_cfg = ep[((ep.estimator != "FE") & (ep.test == "knha"))
                        | ((ep.estimator == "FE") & (ep.test == "z"))]
            compare(led, f"A2 {ds} {oc}: pooled.csv (R engine, HKSJ config)", rp, ep_cfg,
                    "estimator", "estimator", POOLED, var_scale=vs)

            # expected files a dataset does not have are reported, not silently passed over
            em = read(REF / f"{ds}_expected_metareg.csv")
            for what, exp_df in (("metareg", em),
                                 ("Egger", read(REF / f"{ds}_expected_small_study.csv")),
                                 ("leave1out", read(REF / f"{ds}_expected_leave1out.csv")),
                                 ("influence", read(REF / f"{ds}_expected_influence.csv")),
                                 ("subgroups", read(REF / f"{ds}_expected_subgroups.csv"))):
                if exp_df is None:
                    led.skip(f"A1/A2 {ds} {oc}: {what}",
                             "no metafor expected values for this dataset")

            # meta-regression
            if em is not None:
                mod = em.model.iloc[0]
                eff_m = eff.copy()
                for (est, tst), g in em.groupby(["estimator", "test"], sort=False):
                    m = dg.metareg(eff_m, [mod], est, tst == "knha", 0.95)
                    compare(led, f"A1 {ds} {oc}: metareg ~{mod} {est}/{tst}", m, g, "term",
                            "term", MREG, var_scale=vs)
                rm = read(r_out / oc / "metareg.csv")
                compare(led, f"A2 {ds} {oc}: metareg.csv (R engine)", rm,
                        em[(em.estimator == "REML") & (em.test == "knha")], "term", "term",
                        [(r, r) for _, r in MREG], var_scale=vs)
                led.note("metareg I2 is not compared: metaproof does not report a residual I2.")

            # Egger
            es = read(REF / f"{ds}_expected_small_study.csv")
            if es is not None:
                eg = dg.egger(y, v)
                compare(led, f"A1 {ds} {oc}: Egger regression test", eg,
                        es[es.test == "egger"], None, None, EGGER, by_position=True, var_scale=vs)
                rs = read(r_out / oc / "small_study.csv")
                compare(led, f"A2 {ds} {oc}: small_study.csv Egger (R engine)",
                        rs[rs.test == "egger"], es[es.test == "egger"], None, None,
                        [(r, r) for _, r in EGGER] + [("limit_ci_lb", "limit_ci_lb"),
                                                      ("limit_ci_ub", "limit_ci_ub")],
                        by_position=True, var_scale=vs)
                led.note("egger_rma rows (regtest with model = 'rma') are R-only; metaproof "
                         "implements the classic Egger 1997 test.")

            # leave-one-out and influence (REML + HKSJ)
            el = read(REF / f"{ds}_expected_leave1out.csv")
            if el is not None:
                lo = dg.leave_one_out(ids, y, v, "REML", True, 0.95)
                compare(led, f"A1 {ds} {oc}: leave1out REML/knha", lo, el, "omitted",
                        "study_id", LOO, var_scale=vs)
                compare(led, f"A2 {ds} {oc}: leave1out.csv (R engine)",
                        read(r_out / oc / "leave1out.csv"), el, "study_id", "study_id",
                        [(r, r) for _, r in LOO], var_scale=vs)
            ei = read(REF / f"{ds}_expected_influence.csv")
            if ei is not None:
                inf = dg.influence(ids, y, v, "REML", hksj=True)
                compare(led, f"A1 {ds} {oc}: influence REML/knha", inf, ei, "study_id",
                        "study_id", INFL, var_scale=vs)
                compare(led, f"A2 {ds} {oc}: influence.csv (R engine)",
                        read(r_out / oc / "influence.csv"), ei, "study_id", "study_id",
                        [(r, r) for _, r in INFL], var_scale=vs)

            # subgroups
            eg_ = read(REF / f"{ds}_expected_subgroups.csv")
            if eg_ is not None:
                sg = dg.subgroups(eff, "population", "REML", True, 0.95, True, "riley")
                lv = eg_[eg_.level != "ALL"]
                compare(led, f"A1 {ds} {oc}: subgroups per level", sg, lv, "level", "level",
                        SUBG, var_scale=vs)
                qb = sg.iloc[[0]][["Q_between", "Q_between_df", "Q_between_p"]]
                compare(led, f"A1 {ds} {oc}: subgroups Q_between (model-based SEs)", qb,
                        eg_[eg_.test == "Q_between_separate_tau2"], None, None,
                        [("Q_between", "test_stat"), ("Q_between_df", "test_df1"),
                         ("Q_between_p", "test_p")], by_position=True, var_scale=vs)
                mq = dg.metareg(eff, ["population"], "REML", True, 0.95).iloc[[0]]
                compare(led, f"A1 {ds} {oc}: subgroup test with common tau2 (QM)", mq,
                        eg_[eg_.test == "QM_common_tau2"], None, None,
                        [("QM", "test_stat"), ("QM_df1", "test_df1"), ("QM_df2", "test_df2"),
                         ("QM_p", "test_p"), ("tau2", "tau2"), ("QE", "Q"), ("QE_p", "Q_p")],
                        by_position=True, var_scale=vs)
                rsg = read(r_out / oc / "subgroups.csv")
                compare(led, f"A2 {ds} {oc}: subgroups.csv (R engine)", rsg, eg_,
                        ("level", "test"), ("level", "test"),
                        [(c, c) for c in ("k", "estimate", "se", "ci_lb", "ci_ub", "pi_lb",
                                          "pi_ub", "tau2", "I2", "Q", "Q_p", "test_stat",
                                          "test_df1", "test_df2", "test_p")], var_scale=vs)

        led.section(f"A3 {ds}: Python CLI against R CLI, the CSV outputs both engines write")
        specs = load_config(REF / cfg)["outcomes"]
        for oc in outcomes:
            compare_engines(led, "A3", py_out, r_out, oc, specs[oc])


# ---------------------------------------------------------------------------------------
# Section B: the SYNTHETIC NR-like example
# ---------------------------------------------------------------------------------------

def section_b(led, work: Path):
    csv = ROOT / "examples" / "nicotinamide-riboside" / "SYNTHETIC_nr_like_extraction.csv"
    cfg = ROOT / "examples" / "nicotinamide-riboside" / "SYNTHETIC_nr_like_config.yaml"
    py_out, r_out = work / "B" / "py", work / "B" / "r"
    run_python(csv, cfg, py_out)
    pr = run_r(csv, cfg, r_out)
    led.section("B  examples/nicotinamide-riboside/SYNTHETIC_nr_like_extraction.csv (SYNTHETIC data): Python CLI "
                "against R CLI")
    for name, spec in load_config(cfg)["outcomes"].items():
        compare_engines(led, "B", py_out, r_out, name, spec)
    for line in (pr.stderr.splitlines() if pr is not None else []):
        if line.startswith("note:"):
            led.note(f"R engine: {line}")


# ---------------------------------------------------------------------------------------
# Section C: constructed SYNTHETIC edge cases
# ---------------------------------------------------------------------------------------

NOTE = "SYNTHETIC edge case for validation/crosscheck.py; invented numbers, not trial data"
COLS = ["study_id", "outcome", "timepoint", "design", "arm", "is_control", "n", "mean", "sd",
        "se", "median", "q1", "q3", "min", "max", "mean_change", "sd_change", "baseline_mean",
        "baseline_sd", "corr_prepost", "corr_within", "paired_md", "paired_se", "dose_mg_day",
        "duration_weeks", "population", "risk_of_bias", "notes",
        "se_change", "timepoint_weeks"]   # optional extension columns


def _row(outcome, sid, arm, ctl, n, design="parallel", **kw):
    r = {c: None for c in COLS}
    r.update(study_id=f"SYNTHETIC-{sid}", outcome=f"SYNTHETIC_edge_{outcome}", timepoint="wk8",
             design=design, arm=f"SYNTHETIC {arm}", is_control=ctl, n=n, notes=NOTE,
             dose_mg_day=0 if ctl else 1000, duration_weeks=8, population="SYNTHETIC_pop_a",
             risk_of_bias="low")
    r.update(kw)
    return r


def _pair(outcome, sid, ne, me, sde, nc, mc, sdc, **kw):
    e_kw = {k[2:]: v for k, v in kw.items() if k.startswith("e_")}
    c_kw = {k[2:]: v for k, v in kw.items() if k.startswith("c_")}
    both = {k: v for k, v in kw.items() if not k.startswith(("e_", "c_"))}
    return [_row(outcome, sid, "NR", 0, ne, mean=me, sd=sde, **both, **e_kw),
            _row(outcome, sid, "placebo", 1, nc, mean=mc, sd=sdc, **both, **c_kw)]


def edge_rows():
    rows = []
    # k = 2 and k = 3
    rows += _pair("k2", "K2-1", 20, 10.0, 3.0, 20, 8.0, 3.2)
    rows += _pair("k2", "K2-2", 30, 9.0, 4.0, 31, 8.5, 3.6)
    for i, (me, mc, n) in enumerate([(14.0, 10.0, 25), (11.0, 10.5, 40), (16.5, 10.2, 18)]):
        rows += _pair("k3", f"K3-{i + 1}", n, me, 4.0 + i / 2, n + 1, mc, 4.2,
                      population=f"SYNTHETIC_pop_{'ab'[i % 2]}", dose_mg_day=300 * (i + 1))
    # tau2 = 0 exactly: Q far below k - 1 (also below the lower chi-square quantile)
    for i, (me, n) in enumerate([(12.00, 20), (12.05, 25), (11.97, 30), (12.02, 22)]):
        rows += _pair("tau2_zero", f"T0-{i + 1}", n, me, 5.0, n, 10.0, 5.0)
    # one dominant study (n = 2000 per arm) among small ones
    for i, (me, n, pop) in enumerate([(10.8, 2000, "a"), (13.0, 15, "b"), (8.5, 12, "a"),
                                      (12.5, 18, "b"), (14.0, 14, "a")]):
        rows += _pair("dominant", f"DOM-{i + 1}", n, me, 4.0, n, 10.0, 4.0,
                      population=f"SYNTHETIC_pop_{pop}", e_dose_mg_day=250 * (i + 1))
    # identical estimates: every MD is exactly 2.0 -> Q = 0, HKSJ scale 0
    for i, n in enumerate([15, 25, 40, 60]):
        rows += _pair("identical", f"ID-{i + 1}", n, 12.0, 3.0, n, 10.0, 3.0)
    # crossover-only outcome: four ways a crossover trial is reported
    for oc in ("crossover_md", "crossover_smd", "crossover_rom"):
        rows += [_row(oc, "X-1", "NR period", 0, 14, "crossover", mean=48.0, sd=11.0,
                      corr_within=0.6),
                 _row(oc, "X-1", "placebo period", 1, 14, "crossover", mean=40.0, sd=10.0,
                      corr_within=0.6)]
        rows += [_row(oc, "X-2", "NR period", 0, 12, "crossover", mean=45.0, sd=9.5),
                 _row(oc, "X-2", "placebo period", 1, 12, "crossover", mean=41.5, sd=9.0)]
        rows += [_row(oc, "X-3", "NR period", 0, 16, "crossover", mean=50.0, sd=12.0,
                      paired_md=7.5, paired_se=2.1),
                 _row(oc, "X-3", "placebo period", 1, 16, "crossover", mean=42.5, sd=11.5)]
        rows += [_row(oc, "X-4", "NR period", 0, 10, "crossover", mean=47.0, se=3.1,
                      corr_within=0.4),
                 _row(oc, "X-4", "placebo period", 1, 10, "crossover", mean=43.0, se=2.9,
                      corr_within=0.4)]
        # three-period crossover: two NR doses and placebo in the same participants
        rows += [_row(oc, "X-5", "NR low period", 0, 18, "crossover", mean=44.0, sd=10.0,
                      corr_within=0.5),
                 _row(oc, "X-5", "NR high period", 0, 18, "crossover", mean=49.0, sd=11.0,
                      corr_within=0.5),
                 _row(oc, "X-5", "placebo period", 1, 18, "crossover", mean=40.0, sd=10.5,
                      corr_within=0.5)]
        # period n differ (15 vs 13): the number of complete pairs is not known, so both
        # engines exclude it rather than take min(n) (Handbook 23.2.7.1)
        rows += [_row(oc, "X-6", "NR period", 0, 15, "crossover", mean=46.0, sd=10.0,
                      corr_within=0.5),
                 _row(oc, "X-6", "placebo period", 1, 13, "crossover", mean=42.0, sd=9.5,
                      corr_within=0.5)]
    # change from baseline: reported change SD; imputed SD; imputed with row r; change only
    oc = "change"
    rows += [_row(oc, "C-1", "NR", 0, 22, mean_change=-4.0, sd_change=8.0),
             _row(oc, "C-1", "placebo", 1, 21, mean_change=-1.0, sd_change=8.5)]
    rows += [_row(oc, "C-2", "NR", 0, 25, mean=128.0, sd=13.0, baseline_mean=134.0,
                  baseline_sd=14.0),
             _row(oc, "C-2", "placebo", 1, 24, mean=132.0, sd=12.5, baseline_mean=133.5,
                  baseline_sd=13.0)]
    rows += [_row(oc, "C-3", "NR", 0, 30, mean=126.0, sd=12.0, baseline_mean=131.0,
                  baseline_sd=12.5, corr_prepost=0.65),
             _row(oc, "C-3", "placebo", 1, 30, mean=130.5, sd=12.2, baseline_mean=131.5,
                  baseline_sd=12.8, corr_prepost=0.65)]
    rows += [_row(oc, "C-4", "NR", 0, 18, mean_change=-6.0, sd=14.0, baseline_sd=15.0),
             _row(oc, "C-4", "placebo", 1, 19, mean_change=-2.0, sd=13.5, baseline_sd=14.5)]
    # change reported with the SE of the change (se_change extension column)
    rows += [_row(oc, "C-5", "NR", 0, 16, mean_change=-5.0, se_change=2.0),
             _row(oc, "C-5", "placebo", 1, 16, mean_change=-1.5, se_change=1.9)]
    # three-arm parallel trial among two-arm trials
    oc = "three_arm"
    rows += [_row(oc, "TA-1", "NR 300 mg/day", 0, 20, mean=41.0, sd=9.0, dose_mg_day=300),
             _row(oc, "TA-1", "NR 1000 mg/day", 0, 21, mean=55.0, sd=12.0, dose_mg_day=1000),
             _row(oc, "TA-1", "placebo", 1, 24, mean=30.0, sd=8.0, dose_mg_day=0)]
    rows += _pair(oc, "TA-2", 18, 52.0, 13.0, 17, 31.0, 9.0, e_dose_mg_day=500)
    rows += _pair(oc, "TA-3", 30, 46.0, 11.0, 29, 33.0, 10.0, e_dose_mg_day=250)
    rows += _pair(oc, "TA-4", 15, 60.0, 15.0, 16, 29.0, 8.5, e_dose_mg_day=2000)
    # median conversions: S1 (range), S2 (IQR), S3 (both), plus a mean/SD study
    oc = "median"
    rows += [_row(oc, "M-1", "NR", 0, 20, median=70.0, min=40.0, max=120.0),
             _row(oc, "M-1", "placebo", 1, 20, median=32.0, min=20.0, max=55.0)]
    rows += [_row(oc, "M-2", "NR", 0, 18, median=66.0, q1=52.0, q3=90.0),
             _row(oc, "M-2", "placebo", 1, 17, median=31.0, q1=25.0, q3=38.0)]
    rows += [_row(oc, "M-3", "NR", 0, 26, median=72.0, q1=55.0, q3=93.0, min=35.0, max=140.0),
             _row(oc, "M-3", "placebo", 1, 25, median=30.0, q1=24.0, q3=37.0, min=15.0,
                  max=60.0)]
    rows += _pair(oc, "M-4", 22, 64.0, 18.0, 22, 33.0, 9.0)
    # crossover reported only as change within each period: not usable by either engine
    oc = "crossover_change_only"
    rows += [_row(oc, "XC-1", "NR period", 0, 12, "crossover", mean_change=4.0, sd_change=6.0),
             _row(oc, "XC-1", "placebo period", 1, 12, "crossover", mean_change=0.5,
                  sd_change=6.2)]
    for i, (me, mc) in enumerate([(45.0, 41.5), (48.0, 40.0), (50.0, 44.0)]):
        rows += [_row(oc, f"XC-{i + 2}", "NR period", 0, 14, "crossover", mean=me, sd=9.5,
                      corr_within=0.5),
                 _row(oc, f"XC-{i + 2}", "placebo period", 1, 14, "crossover", mean=mc,
                      sd=9.0, corr_within=0.5)]
    # several reported between-group MDs against one shared control (Handbook 23.3.4)
    oc = "paired_multi"
    rows += [_row(oc, "PM-1", "NR 300", 0, 20, mean=44.0, sd=9.0, paired_md=11.0, paired_se=2.6,
                  dose_mg_day=300),
             _row(oc, "PM-1", "NR 1000", 0, 21, mean=52.0, sd=10.0, paired_md=19.0, paired_se=2.8,
                  dose_mg_day=1000),
             _row(oc, "PM-1", "placebo", 1, 22, mean=33.0, sd=8.5, dose_mg_day=0)]
    rows += _pair(oc, "PM-2", 18, 50.0, 12.0, 18, 34.0, 9.0)
    rows += _pair(oc, "PM-3", 25, 47.0, 11.0, 24, 35.0, 10.0)
    rows += _pair(oc, "PM-4", 15, 55.0, 14.0, 16, 31.0, 9.5)
    # one timepoint per study under a list rule (priority order), and labels in mixed units
    oc = "two_timepoints"
    for tp, me in (("wk6", 15.0), ("wk12", 18.0)):
        rows += [_row(oc, "TT-1", "NR", 0, 20, mean=me, sd=4.0, timepoint=tp),
                 _row(oc, "TT-1", "placebo", 1, 20, mean=11.0, sd=4.2, timepoint=tp)]
    rows += [_row(oc, "TT-2", "NR", 0, 22, mean=16.0, sd=4.5, timepoint="wk12"),
             _row(oc, "TT-2", "placebo", 1, 22, mean=12.0, sd=4.0, timepoint="wk12")]
    rows += [_row(oc, "TT-3", "NR", 0, 18, mean=14.0, sd=3.8, timepoint="wk6"),
             _row(oc, "TT-3", "placebo", 1, 18, mean=11.5, sd=4.1, timepoint="wk6")]
    rows += [_row(oc, "TT-4", "NR", 0, 24, mean=15.5, sd=4.0, timepoint="wk6"),
             _row(oc, "TT-4", "placebo", 1, 24, mean=11.0, sd=4.4, timepoint="wk6")]
    oc = "tp_units"
    for tp, me in (("day 21", 30.0), ("wk6", 33.0)):          # day 21 = 3 weeks < 6 weeks
        rows += [_row(oc, "U-1", "NR", 0, 20, mean=me, sd=6.0, timepoint=tp),
                 _row(oc, "U-1", "placebo", 1, 20, mean=27.0, sd=6.0, timepoint=tp)]
    for tp, me in (("end", 31.0), ("wk6", 32.0)):              # 'end' has no number of weeks
        rows += [_row(oc, "U-2", "NR", 0, 18, mean=me, sd=6.5, timepoint=tp),
                 _row(oc, "U-2", "placebo", 1, 18, mean=28.0, sd=6.0, timepoint=tp)]
    rows += [_row(oc, "U-3", "NR", 0, 22, mean=34.0, sd=7.0, timepoint="end"),
             _row(oc, "U-3", "placebo", 1, 22, mean=28.5, sd=6.5, timepoint="end")]
    for tp, me in (("3 months", 35.0), ("wk8", 31.0)):         # 3 months = 13 weeks
        rows += [_row(oc, "U-4", "NR", 0, 25, mean=me, sd=6.0, timepoint=tp),
                 _row(oc, "U-4", "placebo", 1, 25, mean=28.0, sd=6.2, timepoint=tp)]
    for tp, wk, me in (("visit 1", 4, 30.5), ("visit 2", 12, 33.5)):  # timepoint_weeks given
        rows += [_row(oc, "U-5", "NR", 0, 19, mean=me, sd=6.0, timepoint=tp, timepoint_weeks=wk),
                 _row(oc, "U-5", "placebo", 1, 19, mean=28.0, sd=6.0, timepoint=tp,
                      timepoint_weeks=wk)]
    # ROM where one trial fails Geary's check (mean sqrt(n) / SD < 3)
    oc = "geary"
    rows += _pair(oc, "G-1", 20, 40.0, 9.0, 20, 30.0, 8.0)
    rows += _pair(oc, "G-2", 24, 38.0, 10.0, 23, 29.0, 8.5)
    rows += _pair(oc, "G-3", 10, 2.0, 3.0, 10, 1.2, 2.0)
    rows += _pair(oc, "G-4", 16, 42.0, 11.0, 17, 31.0, 9.0)
    # SMD with parallel trials reporting final and change values plus one crossover: under
    # data_basis: change the crossover (between-person SD) is refused (Handbook 10.5.2)
    oc = "smd_mixed"
    for i, (me, mc, ce, cc, n) in enumerate([(52.0, 44.0, 6.0, 1.0, 20), (50.0, 45.0, 4.5, 0.5, 24),
                                             (55.0, 46.0, 7.0, 1.5, 18)]):
        rows += _pair(oc, f"SM-{i + 1}", n, me, 10.0, n + 1, mc, 9.5, e_mean_change=ce,
                      e_sd_change=7.0, c_mean_change=cc, c_sd_change=6.5)
    rows += [_row(oc, "SM-4", "NR period", 0, 14, "crossover", mean=49.0, sd=10.0,
                  corr_within=0.6),
             _row(oc, "SM-4", "placebo period", 1, 14, "crossover", mean=44.0, sd=9.0,
                  corr_within=0.6)]
    # adjusted MD reported on the NR row next to a control row holding only mean/SD (Astra 4)
    oc = "adjusted_md"
    rows += [_row(oc, "AM-1", "NR", 0, 30, paired_md=2.0, paired_se=0.5),
             _row(oc, "AM-1", "placebo", 1, 30, mean=40.0, sd=8.0)]
    rows += _pair(oc, "AM-2", 25, 43.0, 8.5, 25, 40.5, 8.0)
    rows += _pair(oc, "AM-3", 28, 42.0, 9.0, 27, 40.0, 8.8)
    # reported change statistics next to an unused median/range of the final values: the
    # range is not used, so allow_range_only: false does not exclude the study (Astra 10)
    oc = "range_unused"
    rows += [_row(oc, "RU-1", "NR", 0, 20, mean_change=-5.0, sd_change=6.0, median=120.0,
                  min=95.0, max=150.0),
             _row(oc, "RU-1", "placebo", 1, 21, mean_change=-1.0, sd_change=6.5, median=126.0,
                  min=100.0, max=155.0)]
    for sid, n_, ce, se_, cc, sc in (("RU-2", 22, -4.0, 6.0, -1.5, 6.2),
                                     ("RU-3", 18, -6.0, 7.0, -2.0, 6.8)):
        rows += [_row(oc, sid, "NR", 0, n_, mean_change=ce, sd_change=se_),
                 _row(oc, sid, "placebo", 1, n_ + 1, mean_change=cc, sd_change=sc)]
    # Egger boundary: 9 and 10 studies
    rng = np.random.default_rng(20261101)
    for oc, k in (("k9", 9), ("k10", 10)):
        for i in range(k):
            n = int(rng.integers(12, 80))
            rows += _pair(oc, f"{oc.upper()}-{i + 1}", n, float(round(10 + rng.normal(2, 1.5), 3)),
                          float(round(rng.uniform(3, 6), 3)), n, 10.0,
                          float(round(rng.uniform(3, 6), 3)))
    return pd.DataFrame(rows, columns=COLS)


EDGE_CONFIG = """# SYNTHETIC edge-case config for validation/crosscheck.py (regenerated on every run).
# Every row of SYNTHETIC_edge_cases_extraction.csv is invented; none is trial data.
defaults:
  estimator: REML
  hksj: true
  prediction_interval: true
  level: 0.95
  median_method: luo
  multi_arm: split
  timepoint: last
  corr_prepost: 0.5
  corr_prepost_sensitivity: [0.2, 0.8]
  corr_within: 0.5
  corr_within_sensitivity: [0.2, 0.8]
outcomes:
  # protocol defaults: min_k_pool 3 (k = 2 is not pooled), hksj_rule revman
  SYNTHETIC_edge_k2: {measure: MD}
  SYNTHETIC_edge_k2_pool2: {outcome: SYNTHETIC_edge_k2, measure: MD, min_k_pool: 2}
  SYNTHETIC_edge_k2_pool2_always: {outcome: SYNTHETIC_edge_k2, measure: MD, min_k_pool: 2, hksj_rule: always}
  SYNTHETIC_edge_k2_nohksj: {outcome: SYNTHETIC_edge_k2, measure: MD, hksj: false, min_k_pool: 2}
  SYNTHETIC_edge_k3: {measure: MD, subgroups: [population], moderators: [dose_mg_day]}
  SYNTHETIC_edge_k3_gates_off: {outcome: SYNTHETIC_edge_k3, measure: MD, subgroups: [population], moderators: [dose_mg_day], subgroup_min_k: 0, subgroup_min_per_level: 0, metareg_min_k: 0, metareg_allow_small_k: true, hksj_rule: always}
  SYNTHETIC_edge_k3_separate: {outcome: SYNTHETIC_edge_k3, measure: MD, subgroups: [population], subgroup_method: separate, subgroup_min_k: 0, subgroup_min_per_level: 0, hksj_rule: always}
  SYNTHETIC_edge_k3_truncated: {outcome: SYNTHETIC_edge_k3, measure: MD, hksj_rule: truncated}
  SYNTHETIC_edge_k3_ml_nohksj: {outcome: SYNTHETIC_edge_k3, measure: MD, estimator: ML, hksj: false}
  SYNTHETIC_edge_k3_pi_metafor: {outcome: SYNTHETIC_edge_k3, measure: MD, pi_type: metafor}
  SYNTHETIC_edge_k3_pi_revman: {outcome: SYNTHETIC_edge_k3, measure: MD, pi_type: revman}
  SYNTHETIC_edge_dominant_pi_revman: {outcome: SYNTHETIC_edge_dominant, measure: MD, pi_type: revman}
  SYNTHETIC_edge_k3_smd_ub: {outcome: SYNTHETIC_edge_k3, measure: SMD, smd_variance: UB}
  SYNTHETIC_edge_k3_smd_ls2: {outcome: SYNTHETIC_edge_k3, measure: SMD, smd_variance: LS2}
  SYNTHETIC_edge_k3_rom_bc: {outcome: SYNTHETIC_edge_k3, measure: ROM, rom_bias_correction: true}
  SYNTHETIC_edge_tau2_zero: {measure: MD}
  SYNTHETIC_edge_tau2_zero_always: {outcome: SYNTHETIC_edge_tau2_zero, measure: MD, hksj_rule: always}
  SYNTHETIC_edge_tau2_zero_dl: {outcome: SYNTHETIC_edge_tau2_zero, measure: MD, estimator: DL}
  SYNTHETIC_edge_dominant: {measure: MD, subgroups: [population], moderators: [dose_mg_day]}
  SYNTHETIC_edge_dominant_gates_off: {outcome: SYNTHETIC_edge_dominant, measure: MD, subgroups: [population], moderators: [dose_mg_day], subgroup_min_k: 0, metareg_min_k: 0, metareg_allow_small_k: true}
  SYNTHETIC_edge_dominant_separate: {outcome: SYNTHETIC_edge_dominant, measure: MD, subgroups: [population], subgroup_method: separate, subgroup_min_k: 0}
  SYNTHETIC_edge_dominant_truncated: {outcome: SYNTHETIC_edge_dominant, measure: MD, hksj_rule: truncated}
  SYNTHETIC_edge_dominant_pm: {outcome: SYNTHETIC_edge_dominant, measure: MD, estimator: PM}
  SYNTHETIC_edge_dominant_fe: {outcome: SYNTHETIC_edge_dominant, measure: MD, estimator: FE}
  SYNTHETIC_edge_identical: {measure: MD}
  SYNTHETIC_edge_identical_always: {outcome: SYNTHETIC_edge_identical, measure: MD, hksj_rule: always}
  SYNTHETIC_edge_crossover_md: {measure: MD}
  SYNTHETIC_edge_crossover_smd: {measure: SMD}
  SYNTHETIC_edge_crossover_smd_noj: {outcome: SYNTHETIC_edge_crossover_smd, measure: SMD, crossover_smd_j: false}
  SYNTHETIC_edge_crossover_rom: {measure: ROM}
  SYNTHETIC_edge_crossover_rom_bc: {outcome: SYNTHETIC_edge_crossover_rom, measure: ROM, rom_bias_correction: true}
  SYNTHETIC_edge_smd_mixed_auto: {outcome: SYNTHETIC_edge_smd_mixed, measure: SMD}
  SYNTHETIC_edge_smd_mixed_change: {outcome: SYNTHETIC_edge_smd_mixed, measure: SMD, data_basis: change}
  SYNTHETIC_edge_adjusted_md: {measure: MD}
  SYNTHETIC_edge_range_unused: {measure: MD}
  SYNTHETIC_edge_change_r02: {outcome: SYNTHETIC_edge_change, measure: MD, data_basis: change, corr_prepost: 0.2}
  SYNTHETIC_edge_change_r05: {outcome: SYNTHETIC_edge_change, measure: MD, data_basis: change, corr_prepost: 0.5}
  SYNTHETIC_edge_change_r08: {outcome: SYNTHETIC_edge_change, measure: MD, data_basis: change, corr_prepost: 0.8}
  SYNTHETIC_edge_change_auto: {outcome: SYNTHETIC_edge_change, measure: MD, data_basis: auto}
  SYNTHETIC_edge_change_smd: {outcome: SYNTHETIC_edge_change, measure: SMD, data_basis: change, min_k_pool: 2}
  SYNTHETIC_edge_three_arm: {measure: MD, moderators: [dose_mg_day], metareg_min_k: 0, metareg_allow_small_k: true}
  SYNTHETIC_edge_three_arm_combine: {outcome: SYNTHETIC_edge_three_arm, measure: MD, multi_arm: combine}
  SYNTHETIC_edge_three_arm_combine_dose: {outcome: SYNTHETIC_edge_three_arm, measure: MD, multi_arm: combine, moderators: [dose_mg_day], metareg_min_k: 0, metareg_allow_small_k: true}
  SYNTHETIC_edge_crossover_change_only: {measure: MD}
  SYNTHETIC_edge_const_dose: {outcome: SYNTHETIC_edge_tau2_zero, measure: MD, moderators: [dose_mg_day], metareg_min_k: 0, metareg_allow_small_k: true}
  SYNTHETIC_edge_three_arm_smd: {outcome: SYNTHETIC_edge_three_arm, measure: SMD}
  SYNTHETIC_edge_three_arm_rom: {outcome: SYNTHETIC_edge_three_arm, measure: ROM}
  SYNTHETIC_edge_paired_multi_split: {outcome: SYNTHETIC_edge_paired_multi, measure: MD}
  SYNTHETIC_edge_paired_multi_combine: {outcome: SYNTHETIC_edge_paired_multi, measure: MD, multi_arm: combine}
  SYNTHETIC_edge_two_timepoints: {measure: MD, timepoint: [wk6, wk12]}
  SYNTHETIC_edge_tp_units_last: {outcome: SYNTHETIC_edge_tp_units, measure: MD, timepoint: last}
  SYNTHETIC_edge_tp_units_wk6: {outcome: SYNTHETIC_edge_tp_units, measure: MD, timepoint: wk6, min_k_pool: 2}
  SYNTHETIC_edge_geary_flag: {outcome: SYNTHETIC_edge_geary, measure: ROM}
  SYNTHETIC_edge_geary_exclude: {outcome: SYNTHETIC_edge_geary, measure: ROM, rom_geary: exclude}
  SYNTHETIC_edge_median: {measure: MD, median_method: luo, allow_range_only: true}
  SYNTHETIC_edge_median_strict: {outcome: SYNTHETIC_edge_median, measure: MD, median_method: luo}
  SYNTHETIC_edge_median_wan: {outcome: SYNTHETIC_edge_median, measure: MD, median_method: wan, allow_range_only: true}
  SYNTHETIC_edge_k9: {measure: MD}
  SYNTHETIC_edge_k10: {measure: MD}
"""


def hand_checks(led, py_out: Path, r_out: Path):
    """Independent hand calculations from the cited formulas for the data-handling paths."""
    from scipy.stats import norm

    def eff(analysis, engine_out):
        d = read(engine_out / analysis / "study_effects.csv")
        return d.set_index("study_id")

    rows_py, rows_r, exp = [], [], []

    def add(analysis, sid, yi, vi):
        for out, sink in ((py_out, rows_py), (r_out, rows_r)):
            e = eff(analysis, out)
            sink.append({"key": f"{analysis}:{sid}", "yi": e.loc[sid, "yi"] if sid in e.index
                         else math.nan, "vi": e.loc[sid, "vi"] if sid in e.index else math.nan})
        exp.append({"key": f"{analysis}:{sid}", "yi": yi, "vi": vi})

    # Handbook 6.5.2.8 change SD at r = 0.2, 0.5, 0.8 (C-2 assumed r; C-3 row r = 0.65 wins)
    for tag, r in (("r02", 0.2), ("r05", 0.5), ("r08", 0.8)):
        sd = lambda b, f, rr: math.sqrt(b * b + f * f - 2 * rr * b * f)
        add(f"SYNTHETIC_edge_change_{tag}", "SYNTHETIC-C-2", (128 - 134) - (132 - 133.5),
            sd(14, 13, r) ** 2 / 25 + sd(13, 12.5, r) ** 2 / 24)
        add(f"SYNTHETIC_edge_change_{tag}", "SYNTHETIC-C-3", (126 - 131) - (130.5 - 131.5),
            sd(12.5, 12, .65) ** 2 / 30 + sd(12.8, 12.2, .65) ** 2 / 30)
        add(f"SYNTHETIC_edge_change_{tag}", "SYNTHETIC-C-4", -6 - -2,
            sd(15, 14, r) ** 2 / 18 + sd(14.5, 13.5, r) ** 2 / 19)
    # auto: C-2 and C-3 have complete final values -> final MD, no imputation
    add("SYNTHETIC_edge_change_auto", "SYNTHETIC-C-2", 128 - 132, 13**2 / 25 + 12.5**2 / 24)
    # Handbook 23.3.4: shared placebo n 24 split over two comparisons
    add("SYNTHETIC_edge_three_arm", "SYNTHETIC-TA-1|SYNTHETIC NR 300 mg/day", 41 - 30,
        81 / 20 + 64 / 12)
    add("SYNTHETIC_edge_three_arm", "SYNTHETIC-TA-1|SYNTHETIC NR 1000 mg/day", 55 - 30,
        144 / 21 + 64 / 12)
    # Handbook Table 6.5.a: combine the two NR arms
    n, m = 41, (20 * 41 + 21 * 55) / 41
    s2 = (19 * 81 + 20 * 144 + 20 * 21 / 41 * (41**2 + 55**2 - 2 * 41 * 55)) / 40
    add("SYNTHETIC_edge_three_arm_combine", "SYNTHETIC-TA-1", m - 30, s2 / n + 64 / 24)
    # Crossover MD (Handbook 23.2.7; Elbourne 2002): period means with r, default r, paired
    add("SYNTHETIC_edge_crossover_md", "SYNTHETIC-X-1", 8.0,
        (11**2 + 10**2 - 2 * .6 * 11 * 10) / 14)
    add("SYNTHETIC_edge_crossover_md", "SYNTHETIC-X-2", 3.5,
        (9.5**2 + 9**2 - 2 * .5 * 9.5 * 9) / 12)
    add("SYNTHETIC_edge_crossover_md", "SYNTHETIC-X-3", 7.5, 2.1**2)
    sde, sdc = 3.1 * math.sqrt(10), 2.9 * math.sqrt(10)
    add("SYNTHETIC_edge_crossover_md", "SYNTHETIC-X-4", 4.0,
        (sde**2 + sdc**2 - 2 * .4 * sde * sdc) / 10)
    # three-period crossover X-5: its two NR-vs-placebo contrasts share participants and their
    # covariance is not available from summary data, so it is excluded, not split (Handbook
    # 23.2.7, 23.3.4); checked in section_c, "multi-period crossover excluded"
    # Crossover SMD with r derived from the paired SE (X-3); default J on N - 1 = 15 df
    from scipy.special import gammaln
    J = lambda m: math.exp(gammaln(m / 2) - 0.5 * math.log(m / 2) - gammaln((m - 1) / 2))
    r3 = (12**2 + 11.5**2 - 16 * 2.1**2) / (2 * 12 * 11.5)
    d3 = 7.5 / math.sqrt((12**2 + 11.5**2) / 2)
    add("SYNTHETIC_edge_crossover_smd_noj", "SYNTHETIC-X-3", d3,
        (1 / 16 + d3**2 / 32) * 2 * (1 - r3))
    g3 = J(15) * d3
    add("SYNTHETIC_edge_crossover_smd", "SYNTHETIC-X-3", g3, (1 / 16 + g3**2 / 32) * 2 * (1 - r3))
    # crossover log ROM with the Lajeunesse 2015 term for each period mean (X-1, r = 0.6)
    v1 = 11**2 / (14 * 48**2) + 10**2 / (14 * 40**2) - 2 * .6 * 11 * 10 / (14 * 48 * 40)
    add("SYNTHETIC_edge_crossover_rom_bc", "SYNTHETIC-X-1",
        math.log(48 / 40) + .5 * (11**2 / (14 * 48**2) - 10**2 / (14 * 40**2)), v1)
    add("SYNTHETIC_edge_crossover_rom", "SYNTHETIC-X-1", math.log(48 / 40), v1)
    # NR-row adjusted MD used directly although the control row has arm data (Astra 4)
    add("SYNTHETIC_edge_adjusted_md", "SYNTHETIC-AM-1", 2.0, 0.25)
    # reported change MD; the unused median/range does not exclude it (Astra 10)
    add("SYNTHETIC_edge_range_unused", "SYNTHETIC-RU-1", -5.0 - -1.0, 36 / 20 + 6.5**2 / 21)
    # SD of change from its SE (Handbook 6.5.2.2): SD = SE sqrt(n)
    add("SYNTHETIC_edge_change_r05", "SYNTHETIC-C-5", -5.0 - -1.5, 2.0**2 + 1.9**2)
    # list timepoint rule: TT-1 reports wk6 and wk12 and keeps wk6 (first listed)
    add("SYNTHETIC_edge_two_timepoints", "SYNTHETIC-TT-1", 15.0 - 11.0, 16 / 20 + 4.2**2 / 20)
    # 'day 21' is 3 weeks, so wk6 is the latest; U-4 keeps '3 months'; U-5 keeps visit 2
    add("SYNTHETIC_edge_tp_units_last", "SYNTHETIC-U-1", 33.0 - 27.0, 36 / 20 + 36 / 20)
    add("SYNTHETIC_edge_tp_units_last", "SYNTHETIC-U-4", 35.0 - 28.0, 36 / 25 + 6.2**2 / 25)
    add("SYNTHETIC_edge_tp_units_last", "SYNTHETIC-U-5", 33.5 - 28.0, 36 / 19 + 36 / 19)
    # combine the two NR arms of PM-1 (Table 6.5.a) instead of using the two reported MDs
    n, m = 41, (20 * 44 + 21 * 52) / 41
    s2 = (19 * 81 + 20 * 100 + 20 * 21 / 41 * (44**2 + 52**2 - 2 * 44 * 52)) / 40
    add("SYNTHETIC_edge_paired_multi_combine", "SYNTHETIC-PM-1", m - 33, s2 / n + 8.5**2 / 22)
    # Wan 2014 SD and Luo 2018 / Wan 2014 means
    xi = lambda n_: 2 * norm.ppf((n_ - .375) / (n_ + .25))
    eta = lambda n_: 2 * norm.ppf((.75 * n_ - .125) / (n_ + .25))
    sd_e, sd_c = (90 - 52) / eta(18), (38 - 25) / eta(17)
    w_e, w_c = .7 + .39 / 18, .7 + .39 / 17
    add("SYNTHETIC_edge_median", "SYNTHETIC-M-2",
        (w_e * (52 + 90) / 2 + (1 - w_e) * 66) - (w_c * (25 + 38) / 2 + (1 - w_c) * 31),
        sd_e**2 / 18 + sd_c**2 / 17)
    add("SYNTHETIC_edge_median_wan", "SYNTHETIC-M-2", (52 + 66 + 90) / 3 - (25 + 31 + 38) / 3,
        sd_e**2 / 18 + sd_c**2 / 17)
    sd1e, sd1c = (120 - 40) / xi(20), (55 - 20) / xi(20)
    add("SYNTHETIC_edge_median_wan", "SYNTHETIC-M-1", (40 + 140 + 120) / 4 - (20 + 64 + 55) / 4,
        sd1e**2 / 20 + sd1c**2 / 20)
    sd3e = (140 - 35) / (2 * xi(26)) + (93 - 55) / (2 * eta(26))
    sd3c = (60 - 15) / (2 * xi(25)) + (37 - 24) / (2 * eta(25))
    add("SYNTHETIC_edge_median_wan", "SYNTHETIC-M-3",
        (35 + 110 + 144 + 186 + 140) / 8 - (15 + 48 + 60 + 74 + 60) / 8,
        sd3e**2 / 26 + sd3c**2 / 25)

    E = pd.DataFrame(exp)
    compare(led, "C  hand formulas: Python study effects", pd.DataFrame(rows_py), E, "key", "key",
            [("yi", "yi"), ("vi", "vi")])
    compare(led, "C  hand formulas: R study effects", pd.DataFrame(rows_r), E, "key", "key",
            [("yi", "yi"), ("vi", "vi")])


def section_c(led, work: Path):
    EDGE.mkdir(parents=True, exist_ok=True)
    csv = EDGE / "SYNTHETIC_edge_cases_extraction.csv"
    cfg = EDGE / "SYNTHETIC_edge_cases_config.yaml"
    edge_rows().to_csv(csv, index=False)
    cfg.write_text(EDGE_CONFIG)
    py_out, r_out = work / "C" / "py", work / "C" / "r"
    run_python(csv, cfg, py_out)
    pr = run_r(csv, cfg, r_out)
    led.section("C  constructed SYNTHETIC edge cases (validation/edge_cases/): Python CLI "
                "against R CLI")
    for name, spec in load_config(cfg)["outcomes"].items():
        compare_engines(led, "C", py_out, r_out, name, spec)
    for line in (pr.stderr.splitlines() if pr is not None else []):
        if line.startswith("note:"):
            led.note(f"R engine: {line}")
    led.section("C  hand calculations from the cited formulas (Handbook 6.5.2.8, 23.2.7, "
                "23.3.4, Table 6.5.a; Wan 2014; Luo 2018; Elbourne 2002)")
    hand_checks(led, py_out, r_out)
    # the boundary conditions the edge cases exist for
    po = read(py_out / "SYNTHETIC_edge_tau2_zero" / "pooled.csv")
    led.record("C  tau2_zero: every random-effects tau2 is exactly 0 (python)", 4, None,
               [] if (po[po.estimator != "FE"].tau2 == 0).all() else [po.to_string()])
    led.record("C  tau2_zero: hksj_rule revman falls back to the Wald interval (python)", 4,
               None, [] if (po.ci_method == "Wald z").all() else [po[["estimator", "ci_method"]]
                                                                .to_string()])
    p2 = read(py_out / "SYNTHETIC_edge_k2" / "pooled.csv")
    led.record("C  k = 2: no pooled estimate at the default min_k_pool = 3 (python)", 5, None,
               [] if p2.estimate.isna().all() else ["k = 2 was pooled"])
    se = read(py_out / "SYNTHETIC_edge_three_arm_combine_dose" / "study_effects.csv")
    ta = se[se.study_id == "SYNTHETIC-TA-1"]
    led.record("C  combined three-arm trial: dose blank, range 300-1000, never the placebo 0 "
               "(python)", 1, None,
               [] if len(ta) == 1 and ta.dose_mg_day.isna().all() and ta.dose_min.iloc[0] == 300
               and ta.dose_max.iloc[0] == 1000 else [ta.to_string()])
    for a in ("SYNTHETIC_edge_crossover_change_only", "SYNTHETIC_edge_paired_multi_split",
              "SYNTHETIC_edge_tp_units_last",
              "SYNTHETIC_edge_geary_exclude", "SYNTHETIC_edge_median_strict",
              "SYNTHETIC_edge_change_smd"):
        ex = read(py_out / a / "excluded.csv")
        led.record(f"C  {a}: expected exclusion present (python)", 1, None,
                   [] if ex is not None and len(ex) else ["no exclusion"])
    for a in ("SYNTHETIC_edge_crossover_md", "SYNTHETIC_edge_crossover_smd",
              "SYNTHETIC_edge_crossover_rom"):
        for eng, out in (("python", py_out), ("R", r_out)):
            ex, se_ = read(out / a / "excluded.csv"), read(out / a / "study_effects.csv")
            gone = ex is not None and len(ex) and "SYNTHETIC-X-5" in {
                str(x).split("|")[0] for x in ex.study_id}
            kept = se_ is not None and len(se_) and se_.study_id.astype(str).str.startswith(
                "SYNTHETIC-X-5").any()
            led.record(f"C  {a}: multi-period crossover X-5 excluded ({eng})", 1, None,
                       [] if gone and not kept else [f"excluded={bool(gone)} "
                                                     f"in study_effects={bool(kept)}"])
    # exclusions both engines must make, with the reason they must give
    for a, sid, why in (("SYNTHETIC_edge_crossover_md", "SYNTHETIC-X-6", "complete pairs"),
                        ("SYNTHETIC_edge_crossover_smd", "SYNTHETIC-X-6", "complete pairs"),
                        ("SYNTHETIC_edge_crossover_rom", "SYNTHETIC-X-6", "complete pairs"),
                        ("SYNTHETIC_edge_smd_mixed_change", "SYNTHETIC-SM-4", "10.5.2")):
        for eng, out in (("python", py_out), ("R", r_out)):
            ex, se_ = read(out / a / "excluded.csv"), read(out / a / "study_effects.csv")
            hit = ex[ex.study_id.astype(str) == sid] if ex is not None and len(ex) else ex
            ok_ = hit is not None and len(hit) and hit.reason.astype(str).str.contains(
                why, regex=False).all()
            kept = se_ is not None and len(se_) and (se_.study_id.astype(str) == sid).any()
            led.record(f"C  {a}: {sid} excluded ('{why}') ({eng})", 1, None,
                       [] if ok_ and not kept else
                       [f"excluded={bool(ok_)} in study_effects={bool(kept)}"])
    for eng, out in (("python", py_out), ("R", r_out)):
        se_ = read(out / "SYNTHETIC_edge_smd_mixed_auto" / "study_effects.csv")
        led.record(f"C  SYNTHETIC_edge_smd_mixed_auto: crossover SM-4 kept under auto ({eng})",
                   1, None, [] if (se_.study_id.astype(str) == "SYNTHETIC-SM-4").any()
                   else ["SM-4 missing"])
    ss9 = read(py_out / "SYNTHETIC_edge_k9" / "small_study.csv")
    ss10 = read(py_out / "SYNTHETIC_edge_k10" / "small_study.csv")
    led.record("C  Egger runs at k = 10 and not at k = 9 (python)", 2, None,
               [] if ss9.p.isna().all() and ss10.p.notna().all() else ["Egger k rule broken"])


# ---------------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--keep", help="directory to keep both engines' outputs in")
    ap.add_argument("--reuse", help="compare the outputs kept earlier with --keep DIR instead "
                    "of running the engines (the edge-case inputs are still rewritten)")
    ap.add_argument("--report", default=str(REPORT))
    args = ap.parse_args(argv)
    import scipy
    led = Ledger()
    head = [
        "metaproof (Python) vs run_meta.R (R/metafor) cross-check",
        "=====================================================",
        f"run {datetime.now().strftime('%Y-%m-%d %H:%M')} by validation/crosscheck.py",
        f"Python {platform.python_version()} | numpy {np.__version__} | scipy {scipy.__version__}"
        f" | pandas {pd.__version__}",
        r_versions(),
        lock_status(),
        "Tolerance |a - b| <= tol * |reference| + floor. tol 1e-6 for estimates, SEs and tau2,",
        "1e-4 for CIs, prediction intervals, Q-profile limits, influence measures and all other",
        "statistics. Floor: none for study vi and sei; tol x the row's largest |location| or",
        "SE for locations; tol x smallest vi for tau2; 1e-30 for p-values; 1e-5 points for I2,",
        "R2 and weights; 1e-8 for other unit-free statistics (see crosscheck.py). k, df, flags",
        "exact. 'max diff/allowed' is the largest |a - b| / allowed in a check (<= 1 passes).",
        "A SKIP line is a comparison that did not run, with its reason; it never counts as a",
        "pass. Plots are checked for existence only. Reference data: real published metafor/",
        "metadat datasets (methods benchmarks, not NR data). Example and edge-case data are",
        "SYNTHETIC.",
    ]
    global REUSE
    if args.reuse:
        REUSE = True
        work = Path(args.reuse).resolve()
        ctx = None
        head.append(f"--reuse: outputs compared from {work}; the engines were not rerun")
    elif args.keep:
        work = Path(args.keep).resolve()
        work.mkdir(parents=True, exist_ok=True)
        ctx = None
    else:
        ctx = tempfile.TemporaryDirectory(prefix="metaproof_crosscheck_")
        work = Path(ctx.name)
    try:
        section_a(led, work)
        section_b(led, work)
        section_c(led, work)
    finally:
        if ctx is not None:
            ctx.cleanup()
    tail = ["", "Summary", "-------",
            f"{led.n_pass} checks passed, {led.n_fail} failed, {led.n_skip} skipped, "
            f"{led.n_values} values compared."]
    if led.failures:
        tail += ["Failed checks:"] + [f"  {f}" for f in led.failures]
    if led.skipped:
        tail += [f"Skipped ({led.n_skip}; reasons are on each SKIP line above)."]
    text = "\n".join(head + led.lines + tail) + "\n"
    Path(args.report).write_text(text)
    print(text)
    return 0 if led.n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
