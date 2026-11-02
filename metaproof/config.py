"""Analysis configuration: YAML or JSON, a `defaults` block plus per-outcome overrides.

Example (YAML)::

    defaults:
      measure: MD
      estimator: REML
    outcomes:
      whole_blood_nad:
        measure: ROM
        subgroups: [population]
        moderators: [dose_mg_day]

Each key under `outcomes` names one analysis and its output folder. The `outcome` field
selects the extraction rows and defaults to the key, so one outcome can be analysed at two
timepoints under two keys. With no `outcomes` block, every outcome in the CSV is run with
the defaults.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from .pooling import ESTIMATORS


class UserError(ValueError):
    """A problem in the user's extraction file or config. The command line prints it as one
    line, without a traceback."""


DEFAULTS = {
    "measure": "MD",                 # MD | SMD | ROM (log ratio of means)
    "smd_variance": "LS",            # LS | LS2 | UB, see effects.smd_parallel
    "estimator": "REML",             # REML | DL | PM | ML | FE
    "hksj": True,
    # When HKSJ is used for a random-effects fit. revman: only when tau2 > 0 and k - p >= 2,
    # otherwise the Wald z interval (Handbook v6.5 10.10.4.4). always: every random-effects
    # fit with k - p >= 1 (metafor test = "knha"). truncated: always, with the scale floored
    # at 1 (Knapp & Hartung 2003 ad hoc; metafor test = "adhoc"). pooled.csv always carries
    # the Wald, HKSJ and truncated-HKSJ intervals side by side.
    "hksj_rule": "revman",
    "prediction_interval": True,
    "level": 0.95,
    "min_k_pool": 3,                 # protocol F: no pooled estimate below 3 studies
    "data_basis": "auto",            # auto | change | final, see effects.choose_basis
    "corr_prepost": 0.5,
    "corr_prepost_sensitivity": [0.2, 0.8],
    "corr_within": 0.5,
    # None (or "auto"): rerun at AUTO_CORR_WITHIN_SENSITIVITY whenever a crossover comparison
    # uses the assumed corr_within (Handbook 23.2.7); a list is always run; [] turns it off.
    "corr_within_sensitivity": None,
    "median_method": "luo",          # luo (Luo 2018 mean + Wan 2014 SD) | wan (Wan 2014 both)
    "allow_range_only": False,       # median/range without IQR (Wan S1) is excluded unless true
    "multi_arm": "combine",          # combine experimental arms (protocol H main analysis) |
                                     # split the shared control n (dose analyses only)
    "combine_arms": {},              # {study_id: {new_arm: [arm, arm]}}
    "timepoint": "last",             # last | first | <value> | [<value>, ...] (priority order)
    "crossover_smd_j": True,         # Hedges J on N - 1 df for crossover SMD
    "rom_geary": "flag",             # flag | exclude arms failing Geary's mean sqrt(n)/SD >= 3
    "subgroups": [],
    "subgroup_method": "common",     # common tau2 (protocol I) | separate tau2 per subgroup
    "subgroup_min_k": 10,            # protocol I: subgroup test only with >= 10 studies ...
    "subgroup_min_per_level": 2,     # ... and >= 2 studies in every subgroup
    "moderators": [],
    "metareg_joint": False,          # also fit one model with all moderators
    # Studies per moderator (Handbook 10.11.5.1: about ten per characteristic modelled), so a
    # joint model of m moderators needs metareg_min_k x m studies (protocol I).
    "metareg_min_k": 10,
    "metareg_allow_small_k": False,  # metareg_min_k < 10 is refused unless this is true
    "metareg_permutations": 0,       # > 0 adds a permutation p (Higgins & Thompson 2004)
    "egger_min_k": 10,
    "egger_allow_small_k": False,    # egger_min_k < 10 is refused unless this is true
    "bh": False,                     # Benjamini-Hochberg across primary analyses (opt-in)
    "bh_alpha": 0.05,
    "primary": True,                 # analysis belongs to the BH family when bh is on
    "pi_type": "riley",              # riley: t on k-2 df (Higgins 2009, Riley 2011) |
                                     # revman: t on k-1 df always (Handbook 10.10.4.3 formula) |
                                     # metafor: the CI's multiplier, t(k-1) under HKSJ, z
                                     # otherwise (what RevMan prints, Handbook 10.10.4.4)
    "rom_bias_correction": False,    # Lajeunesse 2015 bias-corrected log ratio of means
    "report_estimators": ["REML", "DL", "PM", "ML", "FE"],
    "hksj_adhoc": False,             # R-engine spelling of hksj_rule: truncated
}

AUTO_CORR_WITHIN_SENSITIVITY = [0.2, 0.8]
BOOL_KEYS = ("rom_bias_correction", "allow_range_only", "crossover_smd_j", "bh", "primary",
             "egger_allow_small_k", "metareg_allow_small_k", "metareg_joint", "hksj",
             "prediction_interval", "hksj_adhoc")
RUN_LEVEL_KEYS = ("bh", "bh_alpha")   # act across analyses, so only under defaults

# The R reference engine (reference_r/run_meta.R) uses a few different key names and values.
# Both vocabularies are accepted so one config file drives both engines.
KEY_ALIASES = {
    "effect_measure": "measure", "tau2_estimator": "estimator", "method": "estimator",
    "hartung_knapp": "hksj", "knha": "hksj", "hk": "hksj", "prediction": "prediction_interval",
    "median_conversion": "median_method", "median_iqr_method": "median_method",
    "multiarm": "multi_arm", "multi_arm_handling": "multi_arm", "subgroup": "subgroups",
    "metareg": "moderators", "meta_regression": "moderators", "estimators": "report_estimators",
    "basis": "data_basis",
}
VALUE_ALIASES = {
    # luo2018_wan2014 is the protocol's spelling (Luo 2018 mean, Wan 2014 SD), the luo method
    "median_method": {"wan2014": "wan", "luo2018": "luo", "luo2016": "luo", "luowan": "luo",
                      "luo2018_wan2014": "luo", "luo2018wan2014": "luo"},
    "multi_arm": {"split_control": "split", "combine_treatments": "combine"},
    "hksj_rule": {"knha": "always", "adhoc": "truncated", "trunc": "truncated"},
    "subgroup_method": {"common_tau2": "common", "separate_tau2": "separate"},
    "timepoint": {"latest": "last", "earliest": "first", "last": "last", "first": "first"},
    "estimator": {"EE": "FE", "CE": "FE", "FIXED": "FE"},
    "measure": {"G": "SMD", "HEDGESG": "SMD", "LOGROM": "ROM", "LNRR": "ROM"},
}


def _normalise(spec: dict) -> dict:
    out = {}
    for key, val in (spec or {}).items():
        key = KEY_ALIASES.get(key, key)
        if key in VALUE_ALIASES and isinstance(val, str):
            table = VALUE_ALIASES[key]
            norm = val.upper() if key in ("estimator", "measure") else val.lower()
            norm = re.sub(r"[^A-Za-z0-9_]", "", norm) if key in ("measure", "median_method") else norm
            val = table.get(norm, val)
        if key == "level" and val is not None:
            val = _number(val, "level")
            if val > 1:
                val = val / 100.0        # R side writes 95, Python side 0.95
        out[key] = val
    return out


def _number(val, where: str) -> float:
    if isinstance(val, bool):
        raise UserError(f"{where} must be a number, got {val!r}")
    try:
        return float(val)
    except (TypeError, ValueError):
        raise UserError(f"{where} must be a number, got {val!r}") from None


def _bool(val, where: str) -> bool:
    """true/false as YAML or JSON booleans, 0/1, or the strings true/false in any case (what
    R's as.logical reads). Anything else is refused: bool("false") would be True."""
    if isinstance(val, bool):
        return val
    if isinstance(val, (int, float)) and val in (0, 1):
        return bool(val)
    if isinstance(val, str) and val.strip().lower() in ("true", "false"):
        return val.strip().lower() == "true"
    raise UserError(f"{where} must be true or false, got {val!r}")


CHOICES = {
    "measure": {"MD", "SMD", "ROM"},
    "smd_variance": {"LS", "LS2", "UB"},
    "estimator": set(ESTIMATORS),
    "data_basis": {"auto", "change", "final"},
    "median_method": {"luo", "wan"},
    "multi_arm": {"split", "combine"},
    "pi_type": {"riley", "metafor", "revman"},
    "hksj_rule": {"revman", "always", "truncated"},
    "rom_geary": {"flag", "exclude"},
    "subgroup_method": {"common", "separate"},
}


def _read(path) -> dict:
    if path is None:
        return {}
    if not Path(path).is_file():
        raise UserError(f"config file {path} not found")
    text = Path(path).read_text()
    try:
        if str(path).lower().endswith(".json"):
            raw = json.loads(text)
        else:
            import yaml
            raw = yaml.safe_load(text) or {}
    except Exception as exc:          # yaml.YAMLError, json.JSONDecodeError
        first = str(exc).strip().splitlines()[0] if str(exc).strip() else type(exc).__name__
        raise UserError(f"config {path} could not be read: {first}") from None
    if not isinstance(raw, dict):
        raise UserError(f"config {path} must be a mapping with defaults and outcomes")
    return raw


def _section(raw, key: str) -> dict:
    val = raw.get(key)
    if val is None:
        return {}
    if not isinstance(val, dict):
        raise UserError(f"{key} must be a mapping of setting names to values"
                        + (" (one entry per analysis name)" if key == "outcomes" else ""))
    return val


def load_config(path: str | Path | None) -> dict:
    raw = _read(path)
    unknown_top = set(raw) - {"defaults", "outcomes"}
    if unknown_top:
        raise UserError(f"unknown top-level config keys: {sorted(unknown_top)}")
    raw_defaults = _normalise(_section(raw, "defaults"))
    if "outcome" in raw_defaults:
        raise UserError("defaults.outcome is not allowed: outcome selects the rows of one "
                        "analysis, so set it under outcomes.<name>")
    defaults = {**DEFAULTS, **raw_defaults}
    _validate(defaults, "defaults", raw_defaults)
    outcomes = {}
    for name, spec in _section(raw, "outcomes").items():
        if spec is not None and not isinstance(spec, dict):
            raise UserError(f"outcomes.{name} must be a mapping of settings")
        own = _normalise(spec)
        for key in RUN_LEVEL_KEYS:
            if key in own:
                raise UserError(f"outcomes.{name}.{key}: {key} applies across analyses; set "
                                "it under defaults")
        merged = {**defaults, **own}
        merged.setdefault("outcome", name)
        _validate(merged, f"outcomes.{name}", {**raw_defaults, **own})
        outcomes[str(name)] = merged
    return {"defaults": defaults, "outcomes": outcomes}


def _validate(cfg: dict, where: str, given: dict | None = None) -> None:
    """Check and canonicalise one resolved analysis config in place. given: the keys the
    user set (defaults block plus this analysis), to tell explicit settings from defaults."""
    given = given or {}
    allowed = set(DEFAULTS) | {"outcome"}
    unknown = set(cfg) - allowed
    if unknown:
        raise UserError(f"{where}: unknown keys {sorted(unknown)}")
    for key, opts in CHOICES.items():
        val = cfg[key]
        if isinstance(val, str):
            val = val.upper() if key in ("measure", "smd_variance", "estimator") else val.lower()
            cfg[key] = val
        if val not in opts:
            raise UserError(f"{where}.{key}={cfg[key]!r}; choose from {sorted(opts)}")
    if cfg["measure"] == "ROM" and cfg["data_basis"] == "change":
        raise UserError(f"{where}: measure ROM needs final values (a ratio of mean changes has "
                        "no defined log); use data_basis: final or auto")
    for key in ("corr_prepost", "corr_within"):
        cfg[key] = _number(cfg[key], f"{where}.{key}")
        if not -1 <= cfg[key] <= 1:
            raise UserError(f"{where}.{key} must be in [-1, 1], got {cfg[key]:g}")
    if isinstance(cfg["corr_within_sensitivity"], str) \
            and cfg["corr_within_sensitivity"].strip().lower() == "auto":
        cfg["corr_within_sensitivity"] = None
    for key in ("subgroups", "moderators", "corr_prepost_sensitivity", "corr_within_sensitivity",
                "report_estimators"):
        if key == "corr_within_sensitivity" and cfg[key] is None:
            continue
        if not isinstance(cfg[key], list):
            raise UserError(f"{where}.{key} must be a list")
    for key in ("corr_prepost_sensitivity", "corr_within_sensitivity"):
        if cfg[key] is None:
            continue
        vals = [_number(x, f"{where}.{key} value") for x in cfg[key]]
        bad = [x for x in vals if not -1 <= x <= 1]
        if bad:
            raise UserError(f"{where}.{key} values {bad} must be correlations in [-1, 1]")
        cfg[key] = vals
    est = [str(e).upper() for e in cfg["report_estimators"]]
    est = [VALUE_ALIASES["estimator"].get(e, e) for e in est]
    if set(est) - set(ESTIMATORS):
        raise UserError(f"{where}.report_estimators: unknown {sorted(set(est) - set(ESTIMATORS))}")
    cfg["report_estimators"] = list(dict.fromkeys([cfg["estimator"]] + est))
    for key in BOOL_KEYS:
        cfg[key] = _bool(cfg[key], f"{where}.{key}")
    if cfg["hksj_adhoc"]:
        # R spelling of hksj_rule: truncated. An explicit different rule is a contradiction
        # the two engines would resolve differently, so it is refused.
        if "hksj_rule" in given and cfg["hksj_rule"] != "truncated":
            raise UserError(f"{where}: hksj_adhoc: true conflicts with hksj_rule: "
                            f"{cfg['hksj_rule']}; remove one of them")
        cfg["hksj_rule"] = "truncated"
    cfg["level"] = _number(cfg["level"], f"{where}.level")
    if not 0 < cfg["level"] < 1:
        raise UserError(f"{where}.level must be in (0, 1)")
    cfg["bh_alpha"] = _number(cfg["bh_alpha"], f"{where}.bh_alpha")
    if not 0 < cfg["bh_alpha"] < 1:
        raise UserError(f"{where}.bh_alpha must be in (0, 1)")
    for key in ("min_k_pool", "subgroup_min_k", "subgroup_min_per_level", "metareg_min_k",
                "metareg_permutations", "egger_min_k"):
        x = _number(cfg[key], f"{where}.{key}")
        if x != int(x) or x < 0:
            raise UserError(f"{where}.{key} must be a whole number >= 0, got {cfg[key]!r}")
        cfg[key] = int(x)
    if cfg["min_k_pool"] < 1:
        raise UserError(f"{where}.min_k_pool must be >= 1")
    if cfg["egger_min_k"] < 10 and not cfg["egger_allow_small_k"]:
        raise UserError(f"{where}.egger_min_k = {cfg['egger_min_k']} < 10: tests for funnel "
                        "asymmetry need at least 10 studies (Handbook v6.5 chapter 13; Sterne "
                        "2011). Set egger_allow_small_k: true to override, and say so.")
    if cfg["metareg_min_k"] < 10 and not cfg["metareg_allow_small_k"]:
        raise UserError(f"{where}.metareg_min_k = {cfg['metareg_min_k']} < 10: meta-regression "
                        "needs about ten studies per characteristic modelled (Handbook v6.5 "
                        "10.11.5.1). Set metareg_allow_small_k: true to override, and say why.")
    ca = cfg["combine_arms"] or {}
    if not isinstance(ca, dict) or not all(
            isinstance(g, dict) and all(isinstance(m, list) and m for m in g.values())
            for g in ca.values()):
        raise UserError(f"{where}.combine_arms must map study_id to {{new_arm: [arm, arm]}}")
    cfg["combine_arms"] = ca
    tp = cfg["timepoint"]
    if isinstance(tp, bool) or not (isinstance(tp, (str, int, float))
                                    or (isinstance(tp, list) and tp)):
        raise UserError(f"{where}.timepoint must be last, first, a label or a non-empty list")
    if "outcome" in cfg:
        cfg["outcome"] = str(cfg["outcome"]).strip()
