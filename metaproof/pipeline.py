"""End-to-end run: extraction CSV + config -> out/<analysis>/ files."""

from __future__ import annotations

import difflib
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from . import diagnostics as dg
from . import plots
from . import pooling as pl
from .config import AUTO_CORR_WITHIN_SENSITIVITY, UserError, load_config
from .effects import REQUIRED_COLUMNS, EffectSettings, coerce_extraction, study_effects

NAN = float("nan")
# A timepoint label is a number with an optional unit; it is converted to weeks so labels in
# different units compare correctly ('day 21' is 3 weeks, earlier than 'wk6'). A bare number
# is read as weeks. A hyphen joins a number and its unit ('12-month', 'Wk-12'); it is never a
# minus sign, so a label's weeks are never negative. A timepoint_weeks column, when filled,
# is used instead of the label.
_TP = re.compile(r"(\d*\.?\d+)\s*-?\s*([a-z]*)|([a-z]+)\s*[-:]?\s*(\d*\.?\d+)")
TP_LABEL_TOLERANCE = 0.5   # weeks; a label and timepoint_weeks further apart are logged
_UNIT_WEEKS = {"": 1.0, "w": 1.0, "wk": 1.0, "wks": 1.0, "week": 1.0, "weeks": 1.0,
               "d": 1 / 7, "day": 1 / 7, "days": 1 / 7,
               "h": 1 / 168, "hr": 1 / 168, "hrs": 1 / 168, "hour": 1 / 168, "hours": 1 / 168,
               "m": 365.25 / 84, "mo": 365.25 / 84, "mon": 365.25 / 84, "month": 365.25 / 84,
               "months": 365.25 / 84,
               "y": 365.25 / 7, "yr": 365.25 / 7, "yrs": 365.25 / 7, "year": 365.25 / 7,
               "years": 365.25 / 7}
PI_MIN_K = 5   # Handbook v6.5 10.10.4.3: prediction intervals with fewer studies mislead


def tp_weeks(label) -> float:
    """Weeks implied by a timepoint label ('wk6', '12 weeks', 'day 21', '3 months', '8');
    NaN when the label has no number or an unknown unit."""
    s = str(label).strip().lower()
    m = _TP.search(s)
    if not m:
        return NAN
    if m.group(1) is not None:
        num, unit = float(m.group(1)), m.group(2)
    else:
        unit, num = m.group(3), float(m.group(4))
    unit = unit.rstrip(".")
    return num * _UNIT_WEEKS[unit] if unit in _UNIT_WEEKS else NAN


def _weeks(g: pd.DataFrame):
    """label -> weeks for one study; the timepoint_weeks column wins when filled.

    Returns (weeks, notes, conflicts). notes: labels whose own number of weeks differs from
    timepoint_weeks (the column is used, and the log says so). conflicts: labels whose rows
    give different timepoint_weeks, so the study cannot be placed."""
    out, notes, bad = {}, [], []
    for tp, rows in g.groupby("timepoint", sort=False):
        col = pd.to_numeric(rows["timepoint_weeks"], errors="coerce").dropna().unique() \
            if "timepoint_weeks" in rows else []
        label = tp_weeks(tp)
        if len(col) > 1:
            bad.append(f"rows labelled {tp!r} give timepoint_weeks "
                       f"{', '.join(f'{x:g}' for x in sorted(col))}")
            out[tp] = NAN
        elif len(col) == 1:
            out[tp] = float(col[0])
            if np.isfinite(label) and abs(label - out[tp]) > TP_LABEL_TOLERANCE:
                notes.append(f"label {tp!r} reads as {label:g} weeks; timepoint_weeks "
                             f"{out[tp]:g} is used")
        else:
            out[tp] = label
    return out, notes, bad


def select_timepoint(df: pd.DataFrame, rule):
    """Keep one timepoint per study. Returns (rows kept, log notes, excluded rows).

    rule 'last' / 'first': the latest / earliest timepoint of each study, in weeks.
    an exact label: that label (a bare number matches labels at that many weeks).
    a list of labels: the first label in list order that the study reports.
    A study the rule cannot place is written to excluded.csv with the reason (PRISMA 2020
    item 13a: record which studies enter each synthesis)."""
    notes, excluded, keep = [], [], []
    for sid, g in df.groupby("study_id", sort=False):
        tps = list(dict.fromkeys(g["timepoint"]))
        if isinstance(rule, list):
            pick = next((str(r) for r in rule if str(r) in tps), None)
            if pick is None:
                excluded.append({"study_id": sid, "reason": f"none of the timepoints {rule} "
                                 f"reported (has {tps})"})
                continue
            if len(tps) > 1:
                notes.append(f"{sid}: timepoints {tps}; kept {pick!r} (first in {rule})")
            keep.append(g[g["timepoint"] == pick])
            continue
        if rule not in ("last", "first"):
            want = str(rule)
            hit = g["timepoint"] == want
            if not hit.any() and np.isfinite(tp_weeks(want)) \
                    and re.fullmatch(r"[-+]?\d*\.?\d+", want.strip()):
                wk, wk_notes, conflict = _weeks(g)
                if conflict:
                    excluded.append({"study_id": sid, "reason": "; ".join(conflict)
                                     + "; fix timepoint_weeks"})
                    continue
                notes += [f"{sid}: {n}" for n in wk_notes]
                hit = g["timepoint"].map(lambda t: wk.get(t) == float(want))
            if not hit.any():
                excluded.append({"study_id": sid, "reason": f"no rows at timepoint {want} "
                                 f"(has {tps})"})
                continue
            keep.append(g[hit])
            continue
        if len(tps) == 1:
            keep.append(g)
            continue
        wk, wk_notes, conflict = _weeks(g)
        if conflict:
            excluded.append({"study_id": sid, "reason": "; ".join(conflict)
                             + "; fix timepoint_weeks"})
            continue
        notes += [f"{sid}: {n}" for n in wk_notes]
        bad = [t for t in tps if not np.isfinite(wk[t])]
        if bad:
            excluded.append({"study_id": sid, "reason": f"timepoint labels {bad} have no "
                             "number of weeks; fill timepoint_weeks or relabel"})
            continue
        best = (max if rule == "last" else min)(wk.values())
        picks = [t for t in tps if wk[t] == best]
        if len(picks) > 1:
            excluded.append({"study_id": sid, "reason": f"timepoints {picks} are both "
                             f"{best:g} weeks; keep one label"})
            continue
        notes.append(f"{sid}: timepoints {tps}; kept {picks[0]!r} ({rule}, {best:g} weeks)")
        keep.append(g[g["timepoint"] == picks[0]])
    kept = pd.concat(keep) if keep else df.iloc[0:0]
    return kept, notes, excluded


def _settings(spec, **over) -> EffectSettings:
    return EffectSettings(
        measure=spec["measure"], smd_variance=spec["smd_variance"],
        data_basis=spec["data_basis"],
        corr_prepost=float(over.get("corr_prepost", spec["corr_prepost"])),
        corr_within=float(over.get("corr_within", spec["corr_within"])),
        median_method=over.get("median_method", spec["median_method"]),
        multi_arm=spec["multi_arm"],
        rom_bias_correction=bool(over.get("rom_bias_correction", spec["rom_bias_correction"])),
        extra_columns=tuple(dict.fromkeys(list(spec["subgroups"]) + list(spec["moderators"]))),
        combine={str(k): v for k, v in (spec["combine_arms"] or {}).items()},
        allow_range_only=bool(over.get("allow_range_only", spec["allow_range_only"])),
        crossover_smd_j=bool(spec["crossover_smd_j"]), rom_geary=spec["rom_geary"])


def _safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_") or "outcome"


# Files the run writes next to the analysis folders; no analysis folder may take these names.
RUN_FILES = ("summary_across_analyses.csv", "bh_across_outcomes.csv", "extraction_check.txt",
             "run_manifest.json")
MARKER = ".metaproof_analysis"     # inside every folder metaproof writes; only these are ever removed
STAGING = ".metaproof-staging-"


def analysis_folders(names) -> dict:
    """analysis name -> folder name under --out. Refused: names that would leave --out or
    hide the folder, take a run file's name, or share a folder with another analysis (also
    when they differ only in case: macOS and Windows file systems ignore case)."""
    folders, seen = {}, {}
    reserved = {f.lower() for f in RUN_FILES}
    for name in names:
        f = _safe_name(name)
        if f.startswith(".") or f.lower() in reserved:
            raise UserError(f"analysis name {name!r} gives the folder {f!r}, which "
                            + ("would be hidden or outside --out" if f.startswith(".")
                               else "is a file the run writes") + "; rename the analysis")
        if f.lower() in seen:
            raise UserError(f"analyses {seen[f.lower()]!r} and {name!r} would write to the same "
                            f"folder {f!r}; rename one of them")
        seen[f.lower()] = name
        folders[name] = f
    return folders


def _note_frame(text: str) -> pd.DataFrame:
    return pd.DataFrame([{"note": text}])


def _k_studies(eff: pd.DataFrame) -> int:
    return int(eff["source_study"].nunique()) if len(eff) else 0


def _pool(eff, spec, method):
    """pooled.csv row and test for one estimator, or blanks below min_k_pool."""
    k, k_st = len(eff), _k_studies(eff)
    if k_st < spec["min_k_pool"]:
        row = {c: NAN for c in pl.POOLED_COLUMNS + pl.POOLED_EXTRA_COLUMNS}
        row.update(estimator=method, k=k, ci_method="not pooled")
        note = (f"not pooled: {k_st} studies < min_k_pool = {spec['min_k_pool']} "
                "(protocol F); study results are shown without a pooled estimate")
        return row, {"estimator": method, "k": k, "pval": NAN, "notes": note}, None, note
    y, v = eff["yi"].to_numpy(float), eff["vi"].to_numpy(float)
    row, test, f = pl.pooled_row(y, v, method, spec["hksj"], spec["level"],
                                 spec["prediction_interval"], spec["pi_type"], spec["hksj_rule"])
    notes = []
    if method != "FE" and k_st < PI_MIN_K:
        if spec["prediction_interval"] and np.isfinite(row["pi_lb"]):
            notes.append(f"prediction interval not interpretable with {k_st} studies (< 5; "
                         "Handbook 10.10.4.3); do not use it to judge inconsistency")
        notes.append(f"tau2 CI uninformative with {k_st} studies (< 5)")
    if k != k_st:
        notes.append(f"{k} comparisons from {k_st} studies (split multi-arm); comparisons that "
                     "share a control are correlated and k, df and tau2 treat them as "
                     "independent")
    return row, test, f, "; ".join(notes)


def _membership(base, alt, excl=None) -> str:
    """Which studies a sensitivity scenario drops or adds relative to the main analysis, with
    the exclusion reason when the scenario's own run excluded the study."""
    a = set(map(str, base["source_study"])) if len(base) else set()
    b = set(map(str, alt["source_study"])) if len(alt) else set()
    why = dict(zip(excl["study_id"].astype(str), excl["reason"])) if excl is not None \
        and len(excl) else {}
    parts = []
    if a - b:
        parts.append("drops " + ", ".join(f"{s} ({why[s]})" if s in why else s
                                          for s in sorted(a - b)))
    if b - a:
        parts.append("adds " + ", ".join(sorted(b - a)))
    return "; ".join(parts)


def _sens_row(param, value, eff, spec, method, note="", base=None, excl=None):
    row, test, _, n2 = _pool(eff, spec, method)
    member = _membership(base, eff, excl) if base is not None else ""
    return {"parameter": param, "value": value, "k_studies": _k_studies(eff),
            "n_change_sd_imputed": int(eff["sd_change_imputed"].sum()) if len(eff) else 0,
            "n_crossover": int((eff["design"] == "crossover").sum()) if len(eff) else 0,
            **{c: row[c] for c in pl.POOLED_COLUMNS}, "ci_method": row["ci_method"],
            "pval": test["pval"], "note": "; ".join(x for x in (note, member, n2) if x)}


def egger_precision(eff: pd.DataFrame):
    """Precision term of the modified Egger test for SMD (Pustejovsky & Rodgers 2019, Res
    Synth Methods 10:57-71): the part of SE(g) that does not depend on g. Parallel rows:
    sqrt(1/n_e + 1/n_c). Crossover rows: sqrt(2 (1 - r) / N), the g-free part of
    effects.crossover_smd's variance (1/N + g^2/(2N)) 2 (1 - r), recovered exactly as
    sqrt(vi / (1 + g^2 / 2)) so the r actually used (assumed, reported or derived) is kept.
    None when a parallel row lacks an arm n (the classic test is used)."""
    y, v = eff["yi"].to_numpy(float), eff["vi"].to_numpy(float)
    cross = (eff["design"] == "crossover").to_numpy()
    ne, nc = eff["n_e"].to_numpy(float), eff["n_c"].to_numpy(float)
    ok = np.isfinite(ne) & np.isfinite(nc) & (ne > 0) & (nc > 0)
    if not np.all(ok | cross):
        return None
    with np.errstate(divide="ignore", invalid="ignore"):
        par = np.sqrt(1 / ne + 1 / nc)
    return np.where(cross, np.sqrt(v / (1 + y**2 / 2)), par)


def _sens_values(spec, param, eff):
    vals = spec[f"{param}_sensitivity"]
    if vals is None:      # corr_within: automatic when a crossover comparison assumed r
        assumed = eff["corr_within_assumed"].astype(bool).any() if len(eff) else False
        vals = AUTO_CORR_WITHIN_SENSITIVITY if assumed else []
    return list(vals)


def analyse(name: str, spec: dict, df: pd.DataFrame, folder: Path) -> dict:
    """Write every output of one analysis into folder (run() passes a fresh staging folder,
    so nothing from an earlier run can survive). Returns its summary row."""
    folder.mkdir(parents=True, exist_ok=True)
    log = [f"analysis: {name}", f"outcome rows: {spec['outcome']}",
           *(f"config {k}: {spec[k]}" for k in sorted(spec) if k != "outcome")]
    sub = df[df["outcome"] == str(spec["outcome"]).strip()]
    sub, tp_notes, tp_excl = select_timepoint(sub, spec["timepoint"])
    log += tp_notes

    eff, excl = study_effects(sub, _settings(spec))
    excl = pd.concat([pd.DataFrame(tp_excl, columns=["study_id", "reason"]), excl],
                     ignore_index=True)
    eff.to_csv(folder / "study_effects.csv", index=False)
    excl.to_csv(folder / "excluded.csv", index=False)
    for _, r in excl.iterrows():
        log.append(f"excluded {r.study_id}: {r.reason}")
    # SMDs divided by different SDs (an SD of change, or the SD of final values or crossover
    # period measurements) are different scales and are not pooled (Handbook 10.5.2). The
    # effects layer already refuses the mixtures it can produce (crossover SMD under
    # data_basis: change); this stops any that would reach pooling.
    std = eff["smd_standardiser"].dropna().unique() if len(eff) else []
    if spec["measure"] == "SMD" and len(std) > 1:
        raise UserError(f"analysis {name!r}: SMD rows are standardised by different SDs "
                        f"({', '.join(sorted(map(str, std)))}); they are not pooled together "
                        "(Handbook 10.5.2). Set data_basis: final or change for this analysis")
    k, k_st = len(eff), _k_studies(eff)
    log.append(f"k = {k} comparisons from {k_st} studies")
    for m in dict.fromkeys(list(spec["subgroups"]) + list(spec["moderators"])):
        if m not in eff:
            continue
        miss = eff.loc[eff[m].isna(), "study_id"].tolist()
        if not miss:
            continue
        why = (" (combined arms with different doses, or no dose on the NR rows; see "
               "dose_min/dose_max)" if m == "dose_mg_day" else "")
        log.append(f"{m}: blank for {miss}{why}; these comparisons are left out of the "
                   f"{m} subgroup and meta-regression models")
    if k == 0:
        note = ("no analysable studies: none reported this outcome at the selected timepoint, "
                "or every study was excluded (see excluded.csv)")
        _note_frame(note).to_csv(folder / "pooled.csv", index=False)
        log.append(note)
        (folder / "log.txt").write_text("\n".join(log) + "\n")
        return {"analysis": name, "outcome": spec["outcome"], "measure": spec["measure"],
                "estimator": spec["estimator"], "k": 0, "k_studies": 0, "pooled": False,
                "primary": spec["primary"], "ci_method": "not pooled", "estimate": NAN,
                "ci_lb": NAN, "ci_ub": NAN, "pval": NAN}

    method, hksj, level, pi = spec["estimator"], spec["hksj"], spec["level"], spec["prediction_interval"]
    pit, rule = spec["pi_type"], spec["hksj_rule"]
    y, v = eff["yi"].to_numpy(float), eff["vi"].to_numpy(float)
    ids, groups = list(eff["study_id"]), list(eff["source_study"])
    pooled_ok = k_st >= spec["min_k_pool"]

    order = spec["report_estimators"]          # primary estimator first (config.py)
    rows, tests = [], []
    for m in order:
        row, test, f, note = _pool(eff, spec, m)
        rows.append({**row, "k_studies": k_st, "note": note})
        tests.append(test)
        if f is not None:
            if not f.converged:
                log.append(f"WARNING: {m} did not converge")
            log += [f"{m}: {n}" for n in f.notes]
    if not pooled_ok:
        log.append(rows[0]["note"])
    cols = pl.POOLED_COLUMNS + ["k_studies"] + pl.POOLED_EXTRA_COLUMNS + ["note"]
    pooled = pd.DataFrame(rows, columns=cols)
    pooled.to_csv(folder / "pooled.csv", index=False)
    pd.DataFrame(tests).to_csv(folder / "pooled_tests.csv", index=False)
    log.append(f"primary estimator {method} (first row of pooled.csv); HKSJ={hksj}, "
               f"hksj_rule={rule}; FE rows never use HKSJ; interval used: {rows[0]['ci_method']}")

    # Leave-one-out drops one study (all its comparisons); with 2 studies each "pool" is a
    # single study, so (as in run_meta.R) it is skipped.
    if pooled_ok and k_st >= 3:
        dg.leave_one_out(ids, y, v, method, hksj, level, pi, pit, groups=groups,
                         hksj_rule=rule).to_csv(folder / "leave1out.csv", index=False)
    else:
        _note_frame("not run: fewer than 3 studies or not pooled").to_csv(
            folder / "leave1out.csv", index=False)
    if pooled_ok and k >= 3:
        dg.influence(ids, y, v, method, hksj=hksj, hksj_rule=rule).to_csv(
            folder / "influence.csv", index=False)
    else:
        _note_frame("not run: fewer than 3 comparisons or not pooled").to_csv(
            folder / "influence.csv", index=False)

    if spec["subgroups"]:
        pd.concat([dg.subgroups(eff, m, method, hksj, level, pi, pit,
                                tau2_mode=spec["subgroup_method"],
                                min_k=spec["subgroup_min_k"],
                                min_per_level=spec["subgroup_min_per_level"], hksj_rule=rule,
                                study_col="source_study") for m in spec["subgroups"]],
                  ignore_index=True).to_csv(folder / "subgroups.csv", index=False)
    else:
        _note_frame("no subgroups configured").to_csv(folder / "subgroups.csv", index=False)

    mods = spec["moderators"]
    if mods:
        kw = dict(min_k=spec["metareg_min_k"], hksj_rule=rule,
                  permutations=spec["metareg_permutations"], study_col="source_study")
        parts = [dg.metareg(eff, [m], method, hksj, level, **kw) for m in mods]
        if spec["metareg_joint"] and len(mods) > 1:
            parts.append(dg.metareg(eff, mods, method, hksj, level, model_name="joint", **kw))
        pd.concat(parts, ignore_index=True).to_csv(folder / "metareg.csv", index=False)
    else:
        _note_frame("no moderators configured").to_csv(folder / "metareg.csv", index=False)

    prec = None
    if spec["measure"] == "SMD":
        prec = egger_precision(eff)
        if prec is None:
            log.append("Egger: some SMD rows lack arm n; classic test used")
    dg.egger(y, v, spec["egger_min_k"], precision_se=prec, k_studies=k_st,
             measure=spec["measure"]).to_csv(folder / "small_study.csv", index=False)

    # Imputed correlations: every requested scenario is written, also when it leaves fewer
    # studies (or none), with the studies it drops and why (Handbook 6.5.2.8, 23.2.7).
    sens = []
    used = {"corr_prepost": ("sd_change_imputed", "comparisons with an imputed change SD",
                             "6.5.2.8"),
            "corr_within": ("corr_within_assumed", "crossover comparisons with an assumed "
                            "within-person r", "23.2.7")}
    for param, (flag, what, sec) in used.items():
        values = _sens_values(spec, param, eff)
        n_used = int(eff[flag].astype(bool).sum())
        if n_used and not values:
            log.append(f"WARNING: {param} = {spec[param]:g} is assumed for {n_used} {what}, "
                       f"and {param}_sensitivity is empty: the assumption is not tested "
                       f"(Handbook {sec})")
        if not pooled_ok or not values:
            continue
        for val in [spec[param]] + [x for x in values if x != spec[param]]:
            e2, x2 = study_effects(sub, _settings(spec, **{param: val}))
            sens.append(_sens_row(param, val, e2, spec, method, base=eff, excl=x2))
    if pooled_ok:
        sens += _extra_sensitivities(sub, eff, excl, spec, method)
    if sens:
        pd.DataFrame(sens).to_csv(folder / "sensitivity.csv", index=False)
    else:
        _note_frame("no sensitivity analyses: not pooled, or none configured or applicable"
                    ).to_csv(folder / "sensitivity.csv", index=False)

    rule_lab = rows[0]["ci_method"]
    title = f"{name}: {spec['measure']}, {method}" + ("" if method == "FE" else f", {rule_lab}")
    primary = rows[0]
    # The figure is what gets published, so it carries the same caveat as pooled.csv: with
    # fewer than 5 studies the prediction interval is not drawn (Handbook 10.10.4.3).
    pi_note = None
    if pooled_ok and k_st < PI_MIN_K and np.isfinite(primary["pi_lb"]):
        pi_note = (f"prediction interval not shown: {k_st} studies (< {PI_MIN_K}; Handbook "
                   "10.10.4.3); pooled.csv keeps it with that caveat")
    plots.forest(folder / "forest.png", eff, primary if pooled_ok else None, spec["measure"],
                 level, title if pooled_ok else f"{name}: {spec['measure']}, not pooled "
                 f"({k_st} studies < {spec['min_k_pool']})", pi_note=pi_note)
    plots.funnel(folder / "funnel.png", eff, primary["estimate"] if pooled_ok else None,
                 spec["measure"], level, title)

    (folder / "log.txt").write_text("\n".join(log) + "\n")
    return {"analysis": name, "outcome": spec["outcome"], "measure": spec["measure"],
            "estimator": method, "k": k, "k_studies": k_st, "pooled": pooled_ok,
            "primary": spec["primary"], "ci_method": primary["ci_method"],
            "estimate": primary["estimate"], "ci_lb": primary["ci_lb"],
            "ci_ub": primary["ci_ub"], "pval": tests[0]["pval"]}


def _extra_sensitivities(sub, eff, excl, spec, method):
    """Protocol K5 (median conversion), K7 (imputed or digitised values) and the ROM bias
    correction, run automatically whenever they would change the input."""
    out = []
    conv = eff["converted_from_median"].astype(bool)
    if conv.any():
        other = "wan" if spec["median_method"] == "luo" else "luo"
        e2, x2 = study_effects(sub, _settings(spec, median_method=other))
        out.append(_sens_row("median_method", other, e2, spec, method, "K5a", eff, x2))
        out.append(_sens_row("exclude_converted_from_median", True, eff[~conv], spec, method,
                             "K5b", eff))
    flag = eff["imputed_any"].astype(bool) | eff["digitised"].astype(bool) \
        | eff["sd_from_range_only"].astype(bool)
    if flag.any():
        out.append(_sens_row("exclude_imputed_or_digitised", True, eff[~flag], spec, method,
                             "K7: drops imputed change SDs, assumed crossover r, digitised "
                             "values and range-only SDs", eff))
    if not spec["allow_range_only"] and excl["reason"].astype(str).str.contains(
            "only a range").any():
        e2, x2 = study_effects(sub, _settings(spec, allow_range_only=True))
        out.append(_sens_row("allow_range_only", True, e2, spec, method,
                             "adds studies whose SD comes from a range alone", eff, x2))
    if spec["measure"] == "ROM":
        val = not spec["rom_bias_correction"]
        e2, x2 = study_effects(sub, _settings(spec, rom_bias_correction=val))
        out.append(_sens_row("rom_bias_correction", val, e2, spec, method,
                             "Lajeunesse 2015 small-sample correction", eff, x2))
    return out


# --- the extraction checker ---------------------------------------------------------------

_TEMPLATES = Path(__file__).resolve().parents[1] / "templates"
# Checker errors that make a value unusable as entered: blank or unreadable identifiers,
# arm flags and n, non-numeric or out-of-range numbers, and quartiles out of order. The run
# stops on these. The other checker errors describe a study the engine excludes with a
# reason in excluded.csv (no usable data pattern, control rows, crossover change-only rows),
# or a reporting field the engines do not compute with; they go to extraction_check.txt.
_FATAL = re.compile(r"^row \d+: (?:(?:study_id|outcome|design|arm|is_control|n) is required"
                    r"|\w+ must be |min <= q1 <= median)")


def _checker():
    try:
        if str(_TEMPLATES) not in sys.path:
            sys.path.insert(0, str(_TEMPLATES))
        import check_extraction
    except ImportError as exc:
        raise UserError(f"the extraction checker {_TEMPLATES / 'check_extraction.py'} could "
                        f"not be loaded ({exc}); run metaproof from the kit folder") from None
    return check_extraction


def check_extraction(raw: pd.DataFrame, extra_columns=()):
    """Run templates/check_extraction.py check() on the rows as the engine reads them.
    extra_columns: subgroup or moderator columns the config names that are not interface
    columns; they are declared to the checker, so they are not reported as unexpected.
    Returns (errors that stop the run, other errors, warnings)."""
    ck = _checker()
    extra = tuple(c for c in dict.fromkeys(extra_columns) if c not in ck.ALL_NAMES)
    frame = raw.copy()
    for c in ck.ALL_NAMES:
        if c not in frame:
            frame[c] = None
    frame = frame[list(ck.ALL_NAMES) + [c for c in raw.columns if c not in ck.ALL_NAMES]]
    errors, warnings, _ = ck.check(frame.astype(object), extra)   # as check_extraction.load reads
    fatal = [e for e in errors if _FATAL.match(e)]
    return fatal, [e for e in errors if not _FATAL.match(e)], list(warnings)


def read_extraction(path) -> pd.DataFrame:
    path = Path(path)
    if not path.is_file():
        raise UserError(f"extraction file {path} not found")
    try:
        raw = pd.read_csv(path)
    except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError) as exc:
        first = str(exc).strip().splitlines()[0] if str(exc).strip() else type(exc).__name__
        raise UserError(f"extraction file {path} could not be read as CSV: {first}") from None
    raw = raw.dropna(how="all")
    missing = [c for c in REQUIRED_COLUMNS if c not in raw.columns]
    if missing:
        raise UserError(f"extraction file {path} is missing required columns {missing}")
    return raw


def _write_check(out: Path, path, fatal, other, warnings, unused):
    lines = [f"extraction check of {path}", "rules: templates/check_extraction.py, the checker "
             "extractors run; the run stops on the first group", ""]
    for head, items in ((f"{len(fatal)} errors that stop the run (values that cannot be used "
                         "as entered)", fatal),
                        (f"{len(other)} other checker errors (the engine excludes the study "
                         "with a reason in excluded.csv, or does not compute with the field)",
                         other),
                        (f"{len(warnings)} warnings", warnings),
                        (f"{len(unused)} outcomes in the extraction file that no analysis "
                         "uses", unused)):
        lines += [head + ":"] + [f"  {x}" for x in items] + [""]
    (out / "extraction_check.txt").write_text("\n".join(lines))


# --- output folders --------------------------------------------------------------------------

def _sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _remove(p: Path) -> None:
    if p.is_symlink() or p.is_file():
        p.unlink()
    elif p.exists():
        shutil.rmtree(p)


def _dir_mode() -> int:
    """The mode a plain mkdir would give under the current umask (mkdtemp always uses 0700,
    which would leave the analysis folders unreadable to anyone else sharing the results)."""
    mask = os.umask(0)
    os.umask(mask)
    return 0o777 & ~mask


def _check_owned(final: Path) -> None:
    """Refuse to replace anything metaproof did not write: a file, a link, or a folder that
    has contents but no marker. An empty folder is fine."""
    if final.is_symlink() or final.is_file():
        raise UserError(f"{final} is a file or link, not an analysis folder; move it and rerun")
    if final.is_dir() and any(final.iterdir()) and not (final / MARKER).is_file():
        raise UserError(f"{final} already exists and metaproof did not write it; "
                        "move it or choose another --out")


def _replace(stage: Path, final: Path) -> None:
    """Swap a finished staging folder in for the analysis folder of an earlier run."""
    if final.is_symlink() or final.exists():
        _check_owned(final)
        trash = Path(tempfile.mkdtemp(prefix=STAGING + "old-", dir=final.parent))
        final.rename(trash / "old")
        stage.rename(final)
        _remove(trash)
    else:
        stage.rename(final)


def run(extraction: str | Path, config: str | Path | None, out: str | Path) -> pd.DataFrame:
    """Run every configured analysis. Each analysis folder is written in full to a staging
    folder and then replaces the folder of any earlier run, so no file from an earlier run
    survives; folders of analyses no longer configured are removed (only folders metaproof
    wrote, marked with .metaproof_analysis). run_manifest.json records the input hashes.
    Problems in the inputs raise UserError, which the command line prints as one line."""
    cfg = load_config(config)
    raw = read_extraction(extraction)
    specs = [cfg["defaults"], *cfg["outcomes"].values()]
    declared = [c for sp in specs for c in list(sp["subgroups"]) + list(sp["moderators"])
                if c in raw.columns]
    fatal, other, warnings = check_extraction(raw, declared)
    df = coerce_extraction(raw) if not fatal else None
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for p in out.glob(STAGING + "*"):          # left by an interrupted run
        _remove(p)
    if fatal:
        _write_check(out, extraction, fatal, other, warnings, [])
        raise UserError(f"{Path(extraction).name}: {len(fatal)} value(s) cannot be used as "
                        f"entered; first: {fatal[0]}. Every finding is listed in "
                        f"{out / 'extraction_check.txt'}")
    analyses = cfg["outcomes"] or {
        o: {**cfg["defaults"], "outcome": o} for o in dict.fromkeys(df["outcome"])}
    known = list(dict.fromkeys(df["outcome"]))
    missing = []
    for n, s in analyses.items():
        o = str(s["outcome"]).strip()
        if o not in known:
            close = difflib.get_close_matches(o, known, n=1, cutoff=0.6) \
                or [k for k in known if k.lower() == o.lower()]
            missing.append(f"{n!r} selects outcome {o!r}"
                           + (f" (did you mean {close[0]!r}?)" if close else ""))
    if missing:
        raise UserError("outcome not in the extraction file: " + "; ".join(missing))
    for n, s in analyses.items():
        absent = [m for m in dict.fromkeys(list(s["subgroups"]) + list(s["moderators"]))
                  if m not in raw.columns]
        if absent:
            raise UserError(f"analysis {n!r}: subgroup or moderator columns {absent} are not "
                            "in the extraction file")
    unused = [o for o in known if o not in {str(s["outcome"]).strip()
                                           for s in analyses.values()}]
    folders = analysis_folders(analyses)
    for f in folders.values():
        _check_owned(out / f)
    bh = cfg["defaults"]["bh"]
    if bh:
        fam = [(n, s["outcome"]) for n, s in analyses.items() if s["primary"]]
        seen = {}
        for n, o in fam:
            if o in seen:
                raise UserError(f"bh: analyses {seen[o]!r} and {n!r} are both primary for "
                                f"outcome {o!r}; set primary: false on all but one")
            seen[o] = n
    _write_check(out, extraction, fatal, other, warnings, unused)

    # Every analysis is written to its own staging folder first; only when all of them have
    # finished do they replace the folders of the earlier run, so a failure part-way leaves
    # the earlier run whole.
    summary, stages = [], {}
    try:
        for name, spec in analyses.items():
            stages[name] = Path(tempfile.mkdtemp(prefix=STAGING, dir=out))
            stages[name].chmod(_dir_mode())   # mkdtemp makes it 0700; keep the umask's mode
            summary.append(analyse(name, spec, df, stages[name]))
            (stages[name] / MARKER).write_text(f"written by metaproof {__version__}; the whole "
                                               "folder is replaced on every run\n")
        for name, stage in stages.items():
            _replace(stage, out / folders[name])
    finally:
        for stage in stages.values():
            if stage.exists():
                _remove(stage)
    written = {f.lower() for f in folders.values()}
    removed = sorted(p.name for p in out.iterdir() if p.is_dir() and not p.is_symlink()
                     and (p / MARKER).is_file() and p.name.lower() not in written)
    for name in removed:
        _remove(out / name)

    res = pd.DataFrame(summary)
    stale = out / "bh_across_outcomes.csv"
    if bh and len(res):
        fam = res[res["primary"].astype(bool) & res["pooled"].astype(bool)].copy()
        fam["p_bh"] = dg.bh_adjust(fam["pval"].to_numpy(float))
        fam["bh_alpha"] = cfg["defaults"]["bh_alpha"]
        fam["family"] = "primary pooled analyses, one per outcome"
        fam.to_csv(stale, index=False)
        res = res.merge(fam[["analysis", "p_bh"]], on="analysis", how="left")
    elif stale.exists():
        stale.unlink()
    res.to_csv(out / "summary_across_analyses.csv", index=False)

    def _rel(p):
        import os
        return os.path.relpath(Path(p).resolve(), Path.cwd())

    manifest = {
        "metaproof_version": __version__,
        "run_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "extraction": _rel(extraction), "extraction_sha256": _sha256(extraction),
        "config": _rel(config) if config else None,
        "config_sha256": _sha256(config) if config else None,
        "analyses": [{"analysis": r["analysis"], "folder": folders[r["analysis"]],
                      "outcome": r["outcome"], "k": int(r["k"]),
                      "k_studies": int(r["k_studies"]), "pooled": bool(r["pooled"])}
                     for r in summary],
        "removed_folders": removed,
        "extraction_check": {"other_errors": len(other), "warnings": len(warnings),
                             "unused_outcomes": unused},
    }
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    notices = []
    if other or warnings:
        notices.append(f"extraction checker: {len(other)} findings the run can go on with "
                       "(fields not used in the calculation, or studies excluded with a "
                       f"reason) and {len(warnings)} warnings; see "
                       f"{out / 'extraction_check.txt'}")
    if unused:
        notices.append(f"outcomes no analysis uses: {unused}")
    if removed:
        notices.append(f"removed folders of analyses no longer configured: {removed}")
    res.attrs["notices"] = notices
    return res
