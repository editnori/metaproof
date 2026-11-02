#!/usr/bin/env Rscript
# run_meta.R -- independent reference implementation of the NR meta-analysis
# kit, written in R on top of metafor (Viechtbauer 2010, J Stat Softw 36(3)).
#
# Usage:
#   Rscript reference_r/run_meta.R --data extraction.csv --config config.yaml --out OUTDIR
#   Rscript reference_r/run_meta.R extraction.csv config.yaml OUTDIR
#
# Writes OUTDIR/<outcome>/ with study_effects.csv, pooled.csv, leave1out.csv,
# influence.csv, subgroups.csv, metareg.csv, small_study.csv, forest.png and
# funnel.png (the shared interface), plus sensitivity.csv, study_inputs.csv,
# excluded.csv and config_resolved.json (extra, for cross-engine diffing).
#
# What comes from metafor and what is plain R:
#   metafor : escalc (MD, SMD = Hedges' g, ROM) for parallel arm-level data,
#             rma (all estimators, test = "knha" for HKSJ), predict (CI, PI),
#             confint (Q-profile CI for tau^2), leave1out, influence, regtest,
#             forest, funnel.
#   plain R : crossover effects (Cochrane Handbook ch. 23, section 23.2.7),
#             change-score SD imputation (Handbook 6.5.2.8), SD from SE
#             (Handbook 6.5.2.2), median/IQR/range conversion (Wan 2014;
#             Luo 2018), combining groups (Handbook Table 6.5.a), splitting a
#             shared control (Handbook 23.3.4). metafor::conv.fivenum is not
#             used for the conversion itself because its default SD method for
#             the full five-number summary is Shi et al. 2020, not Wan 2014;
#             the plain-R code is cross-checked against conv.fivenum in
#             reference_r/tests/check_reference.R for the cases where they
#             should agree.
#   Methods details and citations: reference_r/METHODS.txt.

# ---------------------------------------------------------------------------
# Library path: prefer the kit-local user library.
# ---------------------------------------------------------------------------
get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
SCRIPT_DIR <- get_script_dir()
local_lib <- file.path(SCRIPT_DIR, "lib")
env_lib <- Sys.getenv("R_LIBS_USER")
.libPaths(unique(c(if (dir.exists(local_lib)) local_lib, if (nzchar(env_lib) && dir.exists(env_lib)) env_lib, .libPaths())))

suppressPackageStartupMessages({
  library(metafor)
  library(jsonlite)
})

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
num <- function(x) suppressWarnings(as.numeric(x))
ok <- function(x) length(x) == 1L && !is.na(x)
ok_vec <- function(x) !is.na(num(x))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
first_ok <- function(...) { for (v in list(...)) if (ok(v)) return(v); NA_real_ }

STD_COLS <- c("study_id", "outcome", "timepoint", "design", "arm", "is_control",
              "n", "mean", "sd", "se", "median", "q1", "q3", "min", "max",
              "mean_change", "sd_change", "baseline_mean", "baseline_sd",
              "corr_prepost", "corr_within", "paired_md", "paired_se",
              "dose_mg_day", "duration_weeks", "population", "risk_of_bias", "notes")
# Optional extension columns, read when present (se_change: SE of the change from baseline;
# timepoint_weeks: the timepoint in weeks, used instead of parsing the label).
EXT_COLS <- c("se_change", "timepoint_weeks")
NUM_COLS <- c("n", "mean", "sd", "se", "median", "q1", "q3", "min", "max",
              "mean_change", "sd_change", "baseline_mean", "baseline_sd",
              "corr_prepost", "corr_within", "paired_md", "paired_se",
              "dose_mg_day", "duration_weeks", EXT_COLS)
# Arm-level moderators come from the experimental rows only, never from the control row.
ARM_LEVEL_MODS <- c("dose_mg_day")
# tau^2 counts as 0 for hksj_rule = revman below this fraction of the mean within-study
# variance (the same constant as metaproof.pooling.TAU2_ZERO).
TAU2_ZERO <- 1e-8
PI_MIN_K <- 5

POOLED_COLS <- c("estimator", "k", "estimate", "se", "ci_lb", "ci_ub", "pi_lb", "pi_ub",
                 "tau2", "tau2_ci_lb", "tau2_ci_ub", "I2", "H2", "Q", "Q_p")
# stat / pval: the test of the pooled estimate under the configured hksj_rule (z, or t on the
# HKSJ df), the p metaproof writes to pooled_tests.csv and feeds to Benjamini-Hochberg.
POOLED_EXTRA <- c("k_studies", "ci_method", "df", "wald_ci_lb", "wald_ci_ub", "hksj_ci_lb",
                  "hksj_ci_ub", "hksj_trunc_ci_lb", "hksj_trunc_ci_ub", "stat", "pval", "note")

write_out <- function(df, path) write.csv(df, path, row.names = FALSE, na = "NA")

# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
parse_args <- function(a) {
  out <- list(); pos <- character(0); i <- 1L
  while (i <= length(a)) {
    if (a[i] %in% c("--data", "--config", "--out")) {
      out[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L
    } else if (grepl("^--(data|config|out)=", a[i])) {
      kv <- strsplit(sub("^--", "", a[i]), "=", fixed = TRUE)[[1]]
      out[[kv[1]]] <- paste(kv[-1], collapse = "="); i <- i + 1L
    } else { pos <- c(pos, a[i]); i <- i + 1L }
  }
  keys <- c("data", "config", "out")
  for (j in seq_along(pos)) { k <- keys[!keys %in% names(out)][1]; if (!is.na(k)) out[[k]] <- pos[j] }
  if (is.null(out$data) || is.null(out$out))
    stop("usage: Rscript run_meta.R --data extraction.csv [--config config.yaml|.json] --out OUTDIR")
  out
}

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DEFAULTS <- list(
  measure = "MD", estimator = "REML", hksj = TRUE, prediction_interval = TRUE,
  # revman: HKSJ only when tau^2 > 0 and k - p >= 2, else Wald z (Handbook v6.5 10.10.4.4);
  # always: metafor test = "knha"; truncated: test = "adhoc" (scale floored at 1).
  hksj_rule = "revman", min_k_pool = 3,
  pi_type = "riley", basis = "auto",
  corr_prepost = 0.5, corr_prepost_sensitivity = c(0.2, 0.8),
  # corr_within_sensitivity "auto" (or null): rerun at AUTO_CORR_WITHIN_SENSITIVITY when a
  # crossover comparison used the assumed corr_within (Handbook 23.2.7), as metaproof does; a
  # list is always run; [] turns it off.
  corr_within = 0.5, corr_within_sensitivity = "auto",
  # luo2018 = Luo 2018 mean + Wan 2014 SD, the protocol default and the metaproof default.
  median_method = "luo2018", allow_range_only = FALSE,
  multi_arm = "combine_treatments", timepoint = "latest", crossover_smd_j = TRUE,
  rom_geary = "flag",
  report_estimators = c("REML", "DL", "PM", "ML", "FE"),
  subgroups = character(0), subgroup_method = "common", subgroup_min_k = 10,
  subgroup_min_per_level = 2, moderators = character(0), metareg_joint = FALSE,
  # metareg_min_k is studies per moderator (Handbook 10.11.5.1: about ten per characteristic);
  # below 10 it needs metareg_allow_small_k: true, as egger_min_k needs egger_allow_small_k.
  metareg_min_k = 10, metareg_allow_small_k = FALSE,
  level = 95, egger_min_k = 10, egger_allow_small_k = FALSE, rom_bias_correction = FALSE,
  hksj_adhoc = FALSE, smd_variance = "LS")

# Keys the Python engine (metaproof) understands but this engine does not implement. A shared
# config must not change one engine only, so a value that would act in Python stops this
# engine instead of being dropped with a note. bh_alpha and primary only act when bh is
# true, so they are accepted (and ignored) while bh is false.
PY_ONLY_KEYS <- c("combine_arms", "bh_alpha", "bh", "primary", "metareg_permutations")
PY_ONLY_INERT <- list(
  combine_arms = function(v) !length(v),
  bh = function(v) isFALSE(tryCatch(cfg_flag(v, "bh", ""), error = function(e) NA)),
  metareg_permutations = function(v) is.numeric(v) && length(v) == 1L && !is.na(v) && v == 0,
  bh_alpha = function(v) TRUE, primary = function(v) TRUE)

KEY_ALIASES <- c(effect_measure = "measure", tau2_estimator = "estimator", method = "estimator",
                 hartung_knapp = "hksj", knha = "hksj", hk = "hksj",
                 prediction = "prediction_interval", median_conversion = "median_method",
                 median_iqr_method = "median_method", multiarm = "multi_arm",
                 multi_arm_handling = "multi_arm", subgroup = "subgroups",
                 metareg = "moderators", meta_regression = "moderators",
                 estimators = "report_estimators", data_basis = "basis")

AUTO_CORR_WITHIN_SENSITIVITY <- c(0.2, 0.8)

normalise_keys <- function(x) {
  if (!length(x)) return(list())
  nm <- names(x); hit <- nm %in% names(KEY_ALIASES)
  nm[hit] <- KEY_ALIASES[nm[hit]]; names(x) <- nm; x
}

# Strict readers: a value this engine cannot implement stops the run with the key and the
# analysis named, so a typo or an unsupported option cannot silently change one engine only.
cfg_fail <- function(where, key, msg) stop(sprintf("config %s.%s: %s", where, key, msg), call. = FALSE)
# Booleans: true/false, 0/1, or the strings true/false in any case (the vocabulary metaproof
# accepts); anything else, e.g. "yes" in quotes or "maybe", is refused.
cfg_flag <- function(v, key, where) {
  if (length(v) == 1L && !is.na(v)) {
    if (is.logical(v)) return(v)
    if (is.numeric(v) && v %in% c(0, 1)) return(v == 1)
    if (is.character(v) && tolower(trimws(v)) %in% c("true", "false")) return(tolower(trimws(v)) == "true")
  }
  cfg_fail(where, key, sprintf("must be true or false, got %s", deparse(v)))
}
cfg_num <- function(v, key, where, lo = -Inf, hi = Inf, whole = FALSE) {
  if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < lo || v > hi || (whole && v != round(v)))
    cfg_fail(where, key, sprintf("must be a %s in [%s, %s], got %s", if (whole) "whole number" else "number",
                                 format(lo), format(hi), deparse(v)))
  as.numeric(v)
}
cfg_choice <- function(v, key, where, choices, aliases = character(0), prep = tolower) {
  if (!is.character(v) || length(v) != 1L || is.na(v)) cfg_fail(where, key, sprintf("must be one of %s, got %s", paste(choices, collapse = " | "), deparse(v)))
  x <- prep(v)
  if (x %in% names(aliases)) x <- aliases[[x]]
  if (!x %in% choices) cfg_fail(where, key, sprintf("must be one of %s, got %s", paste(choices, collapse = " | "), deparse(v)))
  x
}
EST_ALIASES <- c(EE = "FE", CE = "FE", FIXED = "FE")
ESTIMATORS <- c("REML", "DL", "PM", "ML", "FE")

canon <- function(cfg, where = "defaults", given = names(cfg)) {
  unknown <- setdiff(names(cfg), c(names(DEFAULTS), "outcome", PY_ONLY_KEYS))
  if (length(unknown)) cfg_fail(where, paste(unknown, collapse = ", "), "unknown key (not implemented by run_meta.R or metaproof)")
  for (k in intersect(names(cfg), PY_ONLY_KEYS)) {
    if (!PY_ONLY_INERT[[k]](cfg[[k]]))
      cfg_fail(where, k, sprintf("%s is implemented only in the Python engine; run_meta.R cannot reproduce it. Remove it for the cross-check run.",
                                 deparse(cfg[[k]])))
    cfg[[k]] <- NULL
  }
  cfg$measure <- cfg_choice(cfg$measure, "measure", where, c("MD", "SMD", "ROM"),
                            c(G = "SMD", HEDGESG = "SMD", LOGROM = "ROM", LNRR = "ROM"),
                            prep = function(x) toupper(gsub("[^A-Za-z]", "", x)))
  up <- function(x) toupper(x)
  cfg$estimator <- cfg_choice(cfg$estimator, "estimator", where, ESTIMATORS, EST_ALIASES, prep = up)
  re <- unlist(cfg$report_estimators)
  if (!is.character(re) || !length(re)) cfg_fail(where, "report_estimators", "must be a non-empty list of estimators")
  re <- vapply(re, function(e) cfg_choice(e, "report_estimators", where, ESTIMATORS, EST_ALIASES, prep = up), "", USE.NAMES = FALSE)
  cfg$report_estimators <- unique(c(cfg$estimator, re))
  cfg$median_method <- cfg_choice(cfg$median_method, "median_method", where, c("wan2014", "luo2018"),
                                  c(wan = "wan2014", luo = "luo2018", luo2016 = "luo2018", luowan = "luo2018",
                                    luo2018wan2014 = "luo2018"),
                                  prep = function(x) tolower(gsub("[^A-Za-z0-9]", "", x)))
  cfg$multi_arm <- cfg_choice(cfg$multi_arm, "multi_arm", where, c("split_control", "combine_treatments"),
                              c(split = "split_control", combine = "combine_treatments"),
                              prep = function(x) gsub("-", "_", tolower(x)))
  cfg$hksj_rule <- cfg_choice(cfg$hksj_rule, "hksj_rule", where, c("revman", "always", "truncated"),
                              c(knha = "always", adhoc = "truncated", trunc = "truncated"))
  cfg$subgroup_method <- cfg_choice(cfg$subgroup_method, "subgroup_method", where, c("common", "separate"),
                                    c(common_tau2 = "common", separate_tau2 = "separate"))
  cfg$rom_geary <- cfg_choice(cfg$rom_geary, "rom_geary", where, c("flag", "exclude"))
  cfg$basis <- cfg_choice(cfg$basis, "basis", where, c("auto", "final", "change"))
  cfg$pi_type <- cfg_choice(cfg$pi_type, "pi_type", where, c("riley", "metafor", "revman"))
  cfg$smd_variance <- cfg_choice(cfg$smd_variance, "smd_variance", where, c("LS", "LS2", "UB"), prep = up)
  if (cfg$measure == "ROM" && cfg$basis == "change") cfg_fail(where, "basis", "ROM cannot be computed on change scores")
  for (k in c("hksj", "prediction_interval", "metareg_joint", "rom_bias_correction", "hksj_adhoc",
              "allow_range_only", "crossover_smd_j", "egger_allow_small_k", "metareg_allow_small_k"))
    cfg[[k]] <- cfg_flag(cfg[[k]], k, where)
  if (cfg$hksj_adhoc) {
    # R spelling of hksj_rule: truncated. An explicit different rule is a contradiction, refused
    # as metaproof refuses it, so the two engines cannot resolve it differently.
    if ("hksj_rule" %in% given && cfg$hksj_rule != "truncated")
      cfg_fail(where, "hksj_adhoc", sprintf("true conflicts with hksj_rule: %s; remove one of them", cfg$hksj_rule))
    cfg$hksj_rule <- "truncated"
  }
  for (k in c("corr_prepost", "corr_within")) cfg[[k]] <- cfg_num(cfg[[k]], k, where, -1, 1)
  for (k in c("corr_prepost_sensitivity", "corr_within_sensitivity")) {
    v <- unlist(cfg[[k]])
    if (k == "corr_within_sensitivity" && (is.null(cfg[[k]]) || identical(tolower(as.character(v)), "auto"))) {
      cfg[[k]] <- "auto"; next }
    cfg[[k]] <- vapply(v %||% numeric(0), function(x) cfg_num(x, k, where, -1, 1), 0.0, USE.NAMES = FALSE)
  }
  for (k in c("egger_min_k", "subgroup_min_k", "subgroup_min_per_level", "metareg_min_k"))
    cfg[[k]] <- cfg_num(cfg[[k]], k, where, 0, Inf, whole = TRUE)
  cfg$min_k_pool <- cfg_num(cfg$min_k_pool, "min_k_pool", where, 1, Inf, whole = TRUE)
  # level may be given as 95 (metafor style) or 0.95 (metaproof style); 1 and 100 are refused.
  lv <- cfg_num(cfg$level, "level", where, 0, 100)
  if (lv <= 1) lv <- 100 * lv
  if (!(lv > 0 && lv < 100)) cfg_fail(where, "level", sprintf("must be in (0, 1) or (0, 100), got %s", deparse(cfg$level)))
  cfg$level <- lv
  if (cfg$egger_min_k < 10 && !cfg$egger_allow_small_k)
    cfg_fail(where, "egger_min_k", "< 10 needs egger_allow_small_k: true (tests for funnel asymmetry need at least 10 studies; Sterne 2011)")
  if (cfg$metareg_min_k < 10 && !cfg$metareg_allow_small_k)
    cfg_fail(where, "metareg_min_k", "< 10 needs metareg_allow_small_k: true (about ten studies per characteristic modelled; Handbook 10.11.5.1)")
  for (k in c("subgroups", "moderators")) {
    v <- unlist(cfg[[k]])
    if (length(v) && !is.character(v)) cfg_fail(where, k, "must be a list of column names")
    cfg[[k]] <- as.character(v) %||% character(0)
  }
  # timepoint: last | first | one label or number | a list of labels in priority order. A
  # one-element YAML/JSON list is still a list (labels only, as in metaproof), not a single rule.
  tp <- cfg$timepoint
  is_list <- isTRUE(attr(tp, "cfg_seq")) || length(tp) > 1L
  tp <- unlist(tp)
  if (!length(tp) || !(is.character(tp) || is.numeric(tp)) || anyNA(tp))
    cfg_fail(where, "timepoint", "must be last, first, a label, a number of weeks or a non-empty list of labels")
  tp <- as.character(tp)
  if (!is_list) {
    tl <- tolower(tp)
    if (tl %in% c("last", "latest")) tp <- "latest"
    if (tl %in% c("first", "earliest")) tp <- "earliest"
  }
  attr(tp, "tp_list") <- is_list
  cfg$timepoint <- tp
  if (!is.null(cfg$outcome) && !(is.character(cfg$outcome) && length(cfg$outcome) == 1L))
    cfg_fail(where, "outcome", "must be one outcome name")
  cfg
}

# Mark YAML/JSON sequences so a one-element list stays distinguishable from a scalar.
mark_seq <- function(x) {
  if (!length(x)) return(structure(list(), cfg_seq = TRUE))
  v <- unlist(x)
  if (all(vapply(x, function(e) is.atomic(e) && length(e) == 1L, TRUE)) && length(v) == length(x)) {
    attr(v, "cfg_seq") <- TRUE; v
  } else x
}
simplify_json <- function(x) {
  if (!is.list(x)) return(x)
  if (is.null(names(x))) return(mark_seq(lapply(x, simplify_json)))
  lapply(x, simplify_json)
}

read_config <- function(path, outcomes_in_data) {
  raw <- list()
  if (!is.null(path)) {
    if (grepl("\\.ya?ml$", path, ignore.case = TRUE)) {
      suppressPackageStartupMessages(library(yaml)); raw <- yaml::read_yaml(path, handlers = list(seq = mark_seq))
    } else raw <- simplify_json(jsonlite::fromJSON(path, simplifyVector = FALSE))
    raw <- raw %||% list()
  }
  top <- setdiff(names(raw), c("defaults", "outcomes"))
  if (length(top)) stop("unknown top-level config keys: ", paste(top, collapse = ", "), call. = FALSE)
  raw_def <- normalise_keys(raw$defaults %||% list())
  if ("outcome" %in% names(raw_def))
    stop("config defaults.outcome: not allowed; outcome selects the rows of one analysis, so set it under outcomes.<name>", call. = FALSE)
  defaults <- modifyList(DEFAULTS, raw_def)
  invisible(canon(defaults, "defaults", names(raw_def)))  # the defaults block must be valid on its own (as in metaproof)
  oc <- raw$outcomes
  if (is.null(oc) || !length(oc)) {
    oc <- setNames(rep(list(list()), length(outcomes_in_data)), outcomes_in_data)
  } else if (is.null(names(oc))) {  # list of dicts each with outcome/name
    nm <- vapply(oc, function(o) as.character(o$outcome %||% o$name), "")
    oc <- setNames(lapply(oc, function(o) { o$outcome <- NULL; o$name <- NULL; o }), nm)
  }
  res <- lapply(names(oc), function(nm) {
    own <- normalise_keys(oc[[nm]] %||% list())
    for (k in intersect(c("bh", "bh_alpha"), names(own)))
      cfg_fail(paste0("outcomes.", nm), k, sprintf("%s applies across analyses; set it under defaults", k))
    canon(modifyList(defaults, own), paste0("outcomes.", nm), union(names(raw_def), names(own)))
  })
  names(res) <- names(oc)
  # An analysis key may point at rows of another outcome with `outcome:` (as in metaproof), so
  # one outcome can be analysed twice, e.g. at two timepoints.
  for (nm in names(res)) res[[nm]]$outcome <- trimws(as.character(res[[nm]]$outcome %||% nm))
  res
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
# Values that cannot be used as entered stop the run with the row and column named, the same
# cells that stop metaproof (templates/check_extraction.py rules the Python pipeline treats as
# fatal): a blank study_id/outcome/design/arm/is_control/n, design not parallel/crossover,
# is_control not 0/1, risk_of_bias not low/some concerns/high, a number that does not parse or
# is not finite, n not a whole number >= 1, a negative SD/SE/dose, duration_weeks or paired_se
# <= 0, a correlation outside [-1, 1], and min <= q1 <= median <= q3 <= max violated. A typo
# read as missing would otherwise change which data a study contributes.
NONNEG_COLS <- c("sd", "se", "sd_change", "se_change", "baseline_sd", "dose_mg_day", "timepoint_weeks")
POS_COLS <- c("paired_se", "duration_weeks")
CORR_COLS <- c("corr_prepost", "corr_within")
ROB_LEVELS <- c("low", "some concerns", "high")
check_extraction <- function(raw, line = seq_len(nrow(raw)) + 1L) {
  # line: spreadsheet row of each data row (header is row 1), as the Python checker numbers them
  errs <- character(0)
  add <- function(i, msg) if (length(i)) errs <<- c(errs, sprintf("row %d: %s", line[i], msg))
  blank <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))
  for (col in c("study_id", "outcome", "design", "arm", "is_control", "n")) {
    b <- which(blank(raw[[col]])); for (i in b) add(i, sprintf("%s is required", col))
  }
  ds <- tolower(trimws(raw$design))
  for (i in which(!blank(raw$design) & !ds %in% c("parallel", "crossover")))
    add(i, sprintf("design must be one of ['parallel', 'crossover'], got '%s'", raw$design[i]))
  rb <- tolower(trimws(raw$risk_of_bias))
  for (i in which(!blank(raw$risk_of_bias) & !rb %in% ROB_LEVELS))
    add(i, sprintf("risk_of_bias must be one of ['low', 'some concerns', 'high'], got '%s'", raw$risk_of_bias[i]))
  ic <- num(raw$is_control)
  for (i in which(!blank(raw$is_control) & !(ic %in% c(0, 1))))
    add(i, sprintf("is_control must be 0 or 1, got '%s'", raw$is_control[i]))
  vals <- list()
  for (col in NUM_COLS) {
    x <- raw[[col]]; v <- num(x); filled <- !blank(x)
    for (i in which(filled & is.na(v))) add(i, sprintf("%s must be a number, got '%s'", col, x[i]))
    for (i in which(filled & !is.na(v) & !is.finite(v))) add(i, sprintf("%s must be finite, got '%s'", col, x[i]))
    v[!is.finite(v)] <- NA; vals[[col]] <- v
  }
  n <- vals$n
  for (i in which(!is.na(n) & n != round(n))) add(i, sprintf("n must be a whole number, got %s", format(n[i])))
  for (i in which(!is.na(n) & n < 1)) add(i, sprintf("n must be >= 1, got %s", format(n[i])))
  for (col in NONNEG_COLS) for (i in which(!is.na(vals[[col]]) & vals[[col]] < 0))
    add(i, sprintf("%s must be >= 0, got %s", col, format(vals[[col]][i])))
  for (col in POS_COLS) for (i in which(!is.na(vals[[col]]) & vals[[col]] <= 0))
    add(i, sprintf("%s must be > 0, got %s", col, format(vals[[col]][i])))
  for (col in CORR_COLS) for (i in which(!is.na(vals[[col]]) & abs(vals[[col]]) > 1))
    add(i, sprintf("%s must be in [-1, 1], got %s", col, format(vals[[col]][i])))
  five <- do.call(cbind, vals[c("min", "q1", "median", "q3", "max")])
  for (i in seq_len(nrow(raw))) {
    f <- five[i, ]; f <- f[!is.na(f)]
    if (length(f) > 1L && any(diff(f) < 0)) add(i, "min <= q1 <= median <= q3 <= max is violated")
  }
  errs
}

read_extraction <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character",
                na.strings = c("", "NA", "NaN", "na", "N/A", "n/a"), strip.white = TRUE)
  filled <- rowSums(!is.na(d) & nzchar(trimws(as.matrix(d)))) > 0L  # rows of empty cells are dropped
  line <- which(filled) + 1L
  d <- d[filled, , drop = FALSE]
  miss <- setdiff(STD_COLS, names(d))
  for (m in miss) d[[m]] <- NA_character_
  for (m in setdiff(EXT_COLS, names(d))) d[[m]] <- NA_character_
  errs <- check_extraction(d, line)
  if (length(errs)) {
    message(paste(errs, collapse = "\n"))
    stop(sprintf("%s: %d value(s) cannot be used as entered; first: %s", basename(path), length(errs), errs[1]),
         call. = FALSE)
  }
  for (m in NUM_COLS) d[[m]] <- num(d[[m]])
  d$is_control <- as.integer(num(d$is_control))
  d$design <- tolower(trimws(d$design))
  d$study_id <- as.character(d$study_id); d$outcome <- as.character(d$outcome)
  d$arm <- as.character(d$arm); d$timepoint <- as.character(d$timepoint)
  d$.row <- seq_len(nrow(d))
  attr(d, "missing_cols") <- miss
  d
}

# Timepoint label -> weeks: a number with an optional unit ("wk6", "12 weeks", "day 21",
# "3 months", "12-month", "Wk-12"); a bare number is weeks. A hyphen between the number and
# the unit (or the unit and the number) is a separator, never a minus sign, so a label never
# has negative weeks. Unknown units or no number give NA. The same table as
# metaproof.pipeline.tp_weeks, so 'day 21' (3 weeks) sorts before 'wk6' in both engines.
UNIT_WEEKS <- c(w = 1, wk = 1, wks = 1, week = 1, weeks = 1, d = 1 / 7, day = 1 / 7, days = 1 / 7,
                h = 1 / 168, hr = 1 / 168, hrs = 1 / 168, hour = 1 / 168, hours = 1 / 168,
                m = 365.25 / 84, mo = 365.25 / 84, mon = 365.25 / 84, month = 365.25 / 84,
                months = 365.25 / 84, y = 365.25 / 7, yr = 365.25 / 7, yrs = 365.25 / 7,
                year = 365.25 / 7, years = 365.25 / 7)
TP_RE <- "([0-9]*\\.?[0-9]+)\\s*-?\\s*([a-z]*)|([a-z]+)\\s*[-:]?\\s*([0-9]*\\.?[0-9]+)"
tp_weeks1 <- function(x) {
  s <- tolower(trimws(as.character(x)))
  if (is.na(s)) return(NA_real_)
  m1 <- regmatches(s, regexec(TP_RE, s))[[1]]
  if (!length(m1)) return(NA_real_)
  if (nzchar(m1[2])) { v <- num(m1[2]); u <- m1[3] } else { u <- m1[4]; v <- num(m1[5]) }
  u <- sub("\\.$", "", u)
  if (!nzchar(u)) return(v)
  if (u %in% names(UNIT_WEEKS)) v * UNIT_WEEKS[[u]] else NA_real_
}
tp_weeks <- function(x) vapply(x, tp_weeks1, 0.0, USE.NAMES = FALSE)

# Weeks of each timepoint label of one study (rows i of d): the first filled timepoint_weeks
# of that label's rows wins, else the label is parsed (as metaproof.pipeline._weeks). Used by
# both the latest/earliest rules and a numeric timepoint rule, so they agree with each other.
label_weeks <- function(d, i, labels) {
  vapply(labels, function(t) {
    tw <- if ("timepoint_weeks" %in% names(d)) d$timepoint_weeks[i][!is.na(d$timepoint[i]) & d$timepoint[i] == t] else NA
    tw <- tw[!is.na(tw)]; if (length(tw)) tw[1] else tp_weeks1(t) }, 0.0, USE.NAMES = FALSE)
}

# Labels of one study whose rows give more than one timepoint_weeks: the study cannot be
# placed in time, so it is excluded with this reason (as metaproof.pipeline._weeks).
label_conflicts <- function(d, i, labels) {
  if (!"timepoint_weeks" %in% names(d)) return(character(0))
  out <- character(0)
  for (t in unique(labels)) {
    tw <- unique(d$timepoint_weeks[i][!is.na(d$timepoint[i]) & d$timepoint[i] == t])
    tw <- sort(tw[!is.na(tw)])
    if (length(tw) > 1L) out <- c(out, sprintf("rows labelled '%s' give timepoint_weeks %s", t, paste(format(tw), collapse = ", ")))
  }
  out
}

select_timepoint <- function(d, tp, excl) {
  keep <- logical(nrow(d))
  bad <- function(s, reason) excl[[length(excl) + 1L]] <<- data.frame(study_id = s, reason = reason)
  is_list <- isTRUE(attr(tp, "tp_list")) || length(tp) > 1L
  for (s in unique(d$study_id)) {
    i <- which(d$study_id == s); u <- unique(d$timepoint[i])
    if (is_list) {
      # a list of labels: the first one in list order that the study reports
      pick <- tp[tp %in% u][1]
      if (is.na(pick)) bad(s, paste0("none of the timepoints ", paste(tp, collapse = ", "), " reported"))
      else keep[i[!is.na(d$timepoint[i]) & d$timepoint[i] == pick]] <- TRUE
    } else if (tolower(tp) %in% c("latest", "earliest")) {
      if (length(u) <= 1L) { keep[i] <- TRUE; next }
      cf <- label_conflicts(d, i, u)
      if (length(cf)) { bad(s, paste0(paste(cf, collapse = "; "), "; fix timepoint_weeks")); next }
      wk <- label_weeks(d, i, u)
      if (anyNA(wk)) { bad(s, paste0("timepoint labels ", paste(u[is.na(wk)], collapse = ", "), " have no number of weeks")); next }
      best <- if (tolower(tp) == "latest") max(wk) else min(wk)
      if (sum(wk == best) > 1L) { bad(s, paste0("timepoints ", paste(u[wk == best], collapse = ", "), " are the same number of weeks")); next }
      keep[i[d$timepoint[i] %in% u[wk == best]]] <- TRUE
    } else {
      hit <- !is.na(d$timepoint[i]) & d$timepoint[i] == tp
      if (!any(hit) && grepl("^\\s*[-+]?[0-9]*\\.?[0-9]+\\s*$", tp)) {
        # a bare number is a number of weeks, matched on the same label -> weeks map
        cf <- label_conflicts(d, i, d$timepoint[i])
        if (length(cf)) { bad(s, paste0(paste(cf, collapse = "; "), "; fix timepoint_weeks")); next }
        lw <- label_weeks(d, i, d$timepoint[i])
        hit <- !is.na(lw) & lw == num(tp)
      }
      if (!any(hit)) bad(s, paste0("no rows at timepoint ", tp))
      keep[i[hit]] <- TRUE
    }
  }
  list(d = d[keep, , drop = FALSE], excl = excl)
}

# ---------------------------------------------------------------------------
# Median / quartile / range conversions (plain R)
#   Wan X, Wang W, Liu J, Tong T (2014). BMC Med Res Methodol 14:135.
#   Luo D, Wan X, Liu J, Tong T (2018). Stat Methods Med Res 27(6):1785-1805.
# Scenarios: C1 = {min, median, max}; C3 = {q1, median, q3}; C2 = all five.
# ---------------------------------------------------------------------------
fivenum_case <- function(a, q1, m, q3, b) {
  if (ok(a) && ok(q1) && ok(q3) && ok(b)) "C2" else if (ok(q1) && ok(q3)) "C3" else if (ok(a) && ok(b)) "C1" else NA
}
est_mean_fivenum <- function(a, q1, m, q3, b, n, method) {
  cs <- fivenum_case(a, q1, m, q3, b)
  if (is.na(cs) || !ok(m) || !ok(n)) return(list(v = NA_real_, note = NA))
  if (method == "wan2014") {
    v <- switch(cs, C1 = (a + 2 * m + b) / 4, C3 = (q1 + m + q3) / 3,
                C2 = (a + 2 * q1 + 2 * m + 2 * q3 + b) / 8)
    note <- paste0("mean from median (Wan 2014, ", cs, ")")
  } else {
    v <- switch(cs,
      C1 = { w <- 4 / (4 + n^0.75); w * (a + b) / 2 + (1 - w) * m },
      C3 = { w <- 0.7 + 0.39 / n; w * (q1 + q3) / 2 + (1 - w) * m },
      C2 = { w1 <- 2.2 / (2.2 + n^0.75); w2 <- 0.7 - 0.72 / n^0.55
             w1 * (a + b) / 2 + w2 * (q1 + q3) / 2 + (1 - w1 - w2) * m })
    note <- paste0("mean from median (Luo 2018, ", cs, ")")
  }
  list(v = v, note = note)
}
est_sd_wan <- function(a, q1, q3, b, n) {
  cs <- fivenum_case(a, q1, NA, q3, b)
  if (is.na(cs) || !ok(n)) return(list(v = NA_real_, note = NA))
  xi <- 2 * qnorm((n - 0.375) / (n + 0.25))
  eta <- 2 * qnorm((0.75 * n - 0.125) / (n + 0.25))
  v <- switch(cs, C1 = (b - a) / xi, C3 = (q3 - q1) / eta,
              C2 = (b - a) / (2 * xi) + (q3 - q1) / (2 * eta))
  list(v = v, note = paste0("SD from ", c(C1 = "range", C3 = "IQR", C2 = "range+IQR")[cs], " (Wan 2014, ", cs, ")"))
}

# ---------------------------------------------------------------------------
# Arm-level summaries
# ---------------------------------------------------------------------------
final_stats <- function(r, median_method) {
  n <- r$n; m <- r$mean; s <- r$sd; notes <- character(0)
  if (!ok(s) && ok(r$se) && ok(n)) s_from_se <- TRUE else s_from_se <- FALSE
  need <- (!ok(m) && ok(r$median)) || (!ok(s) && !s_from_se)
  range_only <- need && identical(fivenum_case(r$min, r$q1, r$median, r$q3, r$max), "C1")
  if (!ok(m) && ok(r$median)) {
    e <- est_mean_fivenum(r$min, r$q1, r$median, r$q3, r$max, n, median_method)
    if (ok(e$v)) { m <- e$v; notes <- c(notes, e$note) }
  }
  if (!ok(s) && ok(r$se) && ok(n)) { s <- r$se * sqrt(n); notes <- c(notes, "SD = SE*sqrt(n)") }
  if (!ok(s)) {
    e <- est_sd_wan(r$min, r$q1, r$q3, r$max, n)
    if (ok(e$v)) { s <- e$v; notes <- c(notes, e$note) }
  }
  list(mean = m, sd = s, n = n, notes = notes, range_only = range_only)
}

change_stats <- function(r, fin, rho_default) {
  notes <- character(0); cm <- r$mean_change
  if (!ok(r$sd_change) && ok(r$se_change) && ok(r$n)) {
    r$sd_change <- r$se_change * sqrt(r$n); notes <- c(notes, "change SD from SE of change") }
  if (!ok(cm) && ok(fin$mean) && ok(r$baseline_mean)) { cm <- fin$mean - r$baseline_mean; notes <- c(notes, "change mean = final - baseline") }
  rho_own <- ok(r$corr_prepost); rho <- if (rho_own) r$corr_prepost else rho_default
  imp <- if (ok(r$baseline_sd) && ok(fin$sd)) sqrt(r$baseline_sd^2 + fin$sd^2 - 2 * rho * r$baseline_sd * fin$sd) else NA_real_
  list(mean = cm, sd_rep = r$sd_change, sd_imp = imp, rho = rho, rho_imputed = !rho_own, notes = notes)
}

# Cochrane Handbook Table 6.5.a: combine two groups into one.
combine_groups <- function(n, m, s) {
  N <- n[1]; M <- m[1]; S <- s[1]
  if (length(n) > 1L) for (j in 2:length(n)) {
    n2 <- n[j]; m2 <- m[j]; s2 <- s[j]
    # (M - m2)^2 rather than M^2 + m2^2 - 2 M m2, and a centred mean update, so a common shift
    # of the means (e.g. 1e8) cannot cancel the SD away (as metaproof.conversions.combine_groups).
    S <- sqrt(((N - 1) * S^2 + (n2 - 1) * s2^2 + N * n2 / (N + n2) * (M - m2)^2) / (N + n2 - 1))
    M <- M + n2 * (m2 - M) / (N + n2); N <- N + n2
  }
  list(n = N, mean = M, sd = S)
}

# ---------------------------------------------------------------------------
# One study -> one or more comparisons
# ---------------------------------------------------------------------------
# Moderators of one comparison come from its experimental rows (g). A value that differs
# between combined arms is NA, so the comparison drops out of models on that moderator.
# Study-level columns fall back to any row of the study; arm-level ones (dose) never take
# the control row's value.
moderator_value <- function(var, g, rows) {
  if (!var %in% names(rows)) return(NA)
  u <- unique(g[[var]][!is.na(g[[var]])])
  if (length(u) == 1L) return(u)
  if (length(u) > 1L || var %in% ARM_LEVEL_MODS) return(NA)
  vv <- rows[[var]][!is.na(rows[[var]])]; if (length(vv)) vv[1] else NA
}
hedges_j <- function(m) exp(lgamma(m / 2) - 0.5 * log(m / 2) - lgamma((m - 1) / 2))

# Geary's check for the log ratio of means: mean sqrt(n)/SD >= 3 in both arms (Hedges,
# Gurevitch & Curtis 1999; Lajeunesse 2015). rom_geary = flag writes the note, exclude drops
# the comparison; the message is the one metaproof writes.
GEARY_MIN <- 3
geary_check <- function(tv, cv, cfg) {
  g_min <- min(tv$mean * sqrt(tv$n) / tv$sd, cv$mean * sqrt(cv$n) / cv$sd)
  if (!ok(g_min) || g_min >= GEARY_MIN) return(list(g_min = g_min, exclude = FALSE, msg = NULL))
  list(g_min = g_min, exclude = cfg$rom_geary == "exclude",
       msg = sprintf("Geary's check failed (min mean sqrt(n)/SD = %.2f < 3): the log ratio and its variance are unreliable (Hedges 1999)", g_min))
}

# Exclusion messages shared with metaproof.effects (same wording, so excluded.csv reads alike).
CROSSOVER_MULTI <- paste0(
  "multi-period crossover (%s in the same participants): the comparisons share ",
  "participants and the covariance of their contrasts is not available from summary data, ",
  "so the periods can be neither combined with Table 6.5.a (which is for separate groups) ",
  "nor entered as separate comparisons; splitting N does not supply the covariance ",
  "(Handbook 23.2.7 and 23.3.4). Keep one prespecified NR-vs-placebo pair of rows, or one ",
  "paired estimate the trial reports")
CROSSOVER_SMD_CHANGE <- paste0(
  "a crossover SMD is standardised by the between-person SD of the period measurements, not ",
  "by an SD of change, so it cannot be pooled with change-score SMDs (Handbook 10.5.2 and ",
  "23.2.7.2). Use data_basis: final or auto for SMD analyses with crossover trials")
CROSSOVER_N_DIFFER <- paste0(
  "crossover period n differ (%g vs %g): the number of complete pairs is not known and the ",
  "kit does not guess it. Enter the complete-pair N on both rows with the period summaries of ",
  "those participants, or (MD) the reported paired MD and SE")
RANGE_ONLY_MSG <- paste0(
  "%s: only a range is reported (Wan 2014 scenario S1); a mean or SD estimated from a range ",
  "alone is unstable (Handbook 6.5.2.5-6.5.2.6), so it is excluded unless allow_range_only: ",
  "true (use for a sensitivity analysis)")

# Arm values the formulas would silently accept (a squared negative SD, r outside [-1, 1]).
# read_extraction already stops on these; this guards direct calls (as metaproof.effects._check_row).
row_problem <- function(r) {
  if (!ok(r$is_control) || !r$is_control %in% c(0L, 1L)) return(sprintf("%s: is_control must be 0 or 1", r$arm))
  if (!ok(r$n) || r$n < 1) return(sprintf("%s: n is missing or < 1", r$arm))
  for (col in c("sd", "se", "sd_change", "se_change", "baseline_sd"))
    if (ok(r[[col]]) && r[[col]] < 0) return(sprintf("%s: %s must be >= 0 (got %g)", r$arm, col, r[[col]]))
  if (ok(r$paired_se) && r$paired_se <= 0) return(sprintf("%s: paired_se must be > 0 (got %g)", r$arm, r$paired_se))
  for (col in CORR_COLS)
    if (ok(r[[col]]) && abs(r[[col]]) > 1) return(sprintf("%s: %s must be in [-1, 1] (got %g)", r$arm, col, r[[col]]))
  f <- unlist(r[c("min", "q1", "median", "q3", "max")]); f <- f[!is.na(f)]
  if (length(f) > 1L && any(diff(f) < 0)) return(sprintf("%s: five-number summary out of order; need min <= q1 <= median <= q3 <= max", r$arm))
  NULL
}

process_study <- function(rows, cfg, rho_pp, rho_w, mod_vars) {
  sid <- rows$study_id[1]
  bad <- function(reason) list(comp = NULL, excl = data.frame(study_id = sid, reason = reason))
  design <- unique(rows$design)
  if (length(design) != 1L || is.na(design) || !design %in% c("parallel", "crossover"))
    return(bad("design missing or inconsistent (must be parallel or crossover)"))
  if (anyNA(rows$is_control)) return(bad("is_control missing"))
  if (anyDuplicated(rows$arm)) return(bad(paste0("arm label ", paste(unique(rows$arm[duplicated(rows$arm)]), collapse = ", "),
                                                 " appears more than once (several timepoints or duplicate rows; Handbook 6.2.4)")))
  for (j in seq_len(nrow(rows))) { pr <- row_problem(as.list(rows[j, ])); if (!is.null(pr)) return(bad(pr)) }
  ctl <- rows[rows$is_control == 1L, , drop = FALSE]
  trt <- rows[rows$is_control == 0L, , drop = FALSE]
  # A paired MD/SE belongs on the NR row (data dictionary); a control-row value would be reused
  # for every comparison of the study.
  on_ctl <- ctl$arm[ok_vec(ctl$paired_md) | ok_vec(ctl$paired_se)]
  if (length(on_ctl)) return(bad(sprintf("paired_md/paired_se belong on the NR row, not the control row [%s]", paste(on_ctl, collapse = ", "))))
  if (!nrow(trt)) return(bad("no non-control arm"))
  if (design == "crossover") {
    # One NR period against one placebo period. Several periods in the same people cannot be
    # combined (Table 6.5.a is for separate groups) or split: their covariance is unknown.
    if (!nrow(ctl)) return(bad("need at least one control and one experimental arm (a crossover needs its placebo period row)"))
    if (nrow(trt) > 1L || nrow(ctl) > 1L)
      return(bad(sprintf(CROSSOVER_MULTI, sprintf("%d NR and %d control periods", nrow(trt), nrow(ctl)))))
    if (cfg$measure == "SMD" && cfg$basis == "change") return(bad(CROSSOVER_SMD_CHANGE))
  }

  # Arm summaries
  summ <- lapply(seq_len(nrow(rows)), function(j) {
    r <- as.list(rows[j, ]); fin <- final_stats(r, cfg$median_method)
    list(row = r, fin = fin, chg = change_stats(r, fin, rho_pp))
  })
  names(summ) <- as.character(rows$.row)
  S <- function(r) summ[[as.character(r$.row)]]
  # A mean or SD from a range alone (Wan C1) excludes the study only when that value feeds the
  # effect: final values are used, or (change basis) the change mean or SD was built from final
  # values. An unused median/range next to reported change statistics is ignored (Astra 10).
  range_hit <- function(df, basis) {
    if (cfg$allow_range_only) return(character(0))
    hit <- vapply(seq_len(nrow(df)), function(j) {
      z <- S(as.list(df[j, ]))
      used_final <- basis != "change" || !ok(z$chg$sd_rep) || !ok(z$row$mean_change)
      isTRUE(z$fin$range_only) && used_final }, TRUE)
    df$arm[hit]
  }
  # Several reported between-group MDs against one control share it (Handbook 23.3.4).
  trt_paired <- ok_vec(trt$paired_md) & ok_vec(trt$paired_se)
  if (design == "parallel" && cfg$measure == "MD" && nrow(trt) > 1L && any(trt_paired)) {
    fin_trt <- all(vapply(summ[as.character(trt$.row)], function(z) ok(z$fin$mean) && ok(z$fin$sd), TRUE))
    if (!nrow(ctl) || cfg$multi_arm == "split_control" || !fin_trt)
      return(bad("more than one reported between-group MD (paired_md/paired_se) against one control (Handbook 23.3.4)"))
  }

  # Basis for parallel arm-level data
  basis <- "final"; basis_note <- "final values"
  if (design == "parallel") {
    fin_ok <- all(vapply(summ, function(z) ok(z$fin$mean) && ok(z$fin$sd) && ok(z$fin$n), TRUE))
    chg_rep <- all(vapply(summ, function(z) ok(z$chg$mean) && ok(z$chg$sd_rep) && ok(z$row$n), TRUE))
    chg_imp <- all(vapply(summ, function(z) ok(z$chg$mean) && (ok(z$chg$sd_rep) || ok(z$chg$sd_imp)) && ok(z$row$n), TRUE))
    if (cfg$basis == "change" || (cfg$basis == "auto" && cfg$measure == "MD")) {
      if (chg_rep) { basis <- "change"; basis_note <- "change from baseline (reported SD)"
      } else if (cfg$basis == "auto" && fin_ok) { basis <- "final"
      } else if (chg_imp && cfg$measure != "SMD") { basis <- "change"; basis_note <- "change from baseline (SD imputed, Handbook 6.5.2.8)"
      } else basis <- NA  # SMD never standardises by an imputed change SD (Handbook 23.2.7.2)
    } else basis <- if (fin_ok) "final" else NA
  }

  arm_vals <- function(r) {
    z <- S(r)
    if (basis == "change") {
      sd <- if (ok(z$chg$sd_rep)) z$chg$sd_rep else z$chg$sd_imp
      list(n = z$row$n, mean = z$chg$mean, sd = sd,
           imputed_rho = !ok(z$chg$sd_rep) && z$chg$rho_imputed,
           notes = c(z$fin$notes[!ok(z$chg$sd_rep)], z$chg$notes,
                     if (!ok(z$chg$sd_rep)) sprintf("%s: SD change imputed with r=%.2f%s", z$row$arm, z$chg$rho, if (z$chg$rho_imputed) " (default)" else " (study)")))
    } else list(n = z$row$n, mean = z$fin$mean, sd = z$fin$sd, imputed_rho = FALSE,
                notes = if (length(z$fin$notes)) paste0(z$row$arm, ": ", z$fin$notes) else character(0))
  }
  pool_arms <- function(df) {
    v <- lapply(seq_len(nrow(df)), function(j) arm_vals(as.list(df[j, ])))
    if (length(v) == 1L) return(c(v[[1]], list(label = df$arm[1])))
    cg <- combine_groups(vapply(v, `[[`, 0, "n"), vapply(v, `[[`, 0, "mean"), vapply(v, `[[`, 0, "sd"))
    c(cg, list(imputed_rho = any(vapply(v, `[[`, TRUE, "imputed_rho")),
               notes = c(unlist(lapply(v, `[[`, "notes")), paste0("combined arms ", paste(df$arm, collapse = "+"), " (Handbook Table 6.5.a)")),
               label = paste(df$arm, collapse = "+")))
  }

  # Treatment groups to compare (a crossover has exactly one NR row by now)
  if (cfg$multi_arm == "combine_treatments" && nrow(trt) > 1L && design == "parallel") groups <- list(trt) else
    groups <- lapply(seq_len(nrow(trt)), function(j) trt[j, , drop = FALSE])
  m_comp <- length(groups)
  comps <- list(); excl <- list()

  for (g in groups) {
    trow <- if (nrow(g) == 1L) as.list(g[1, ]) else NULL
    label <- if (m_comp > 1L) paste0(sid, "|", paste(g$arm, collapse = "+")) else sid
    notes <- c(design); yi <- NA_real_; vi <- NA_real_; imp_pp <- FALSE; imp_w <- FALSE
    det <- list(n_t = NA, mean_t = NA, sd_t = NA, n_c = NA, mean_c = NA, sd_c = NA, r_within = NA, basis = NA,
                geary_min = NA, prec_v = NA)

    # Reported paired / adjusted MD, from the NR row only (control rows holding them were refused)
    pmd <- if (!is.null(trow)) trow$paired_md else NA; pse <- if (!is.null(trow)) trow$paired_se else NA
    use_paired <- ok(pmd) && ok(pse) && cfg$measure == "MD"

    if (use_paired) {
      yi <- pmd; vi <- pse^2; det$basis <- "paired_md"
      notes <- c(notes, "reported paired/adjusted MD and SE used directly")
      if (m_comp > 1L) notes <- c(notes, "shared-control correlation not adjusted for reported MD")
    } else if (!nrow(ctl)) {
      excl[[length(excl) + 1L]] <- data.frame(study_id = label, reason = "no control arm and no usable paired_md/paired_se"); next
    } else if (design == "parallel") {
      if (is.na(basis)) { excl[[length(excl) + 1L]] <- data.frame(study_id = label, reason = paste0("arm means/SDs not obtainable for basis=", cfg$basis, ", measure=", cfg$measure)); next }
      rh <- range_hit(rbind(g, ctl), basis)
      if (length(rh)) return(bad(sprintf(RANGE_ONLY_MSG, paste(rh, collapse = ", "))))
      tv <- pool_arms(g); cv <- pool_arms(ctl)
      nc <- cv$n
      if (cfg$multi_arm == "split_control" && m_comp > 1L) {
        nc <- cv$n / m_comp; notes <- c(notes, sprintf("shared control n %.4g split over %d comparisons (Handbook 23.3.4)", cv$n, m_comp))
      }
      notes <- c(notes, basis_note, tv$notes, cv$notes); imp_pp <- isTRUE(tv$imputed_rho) || isTRUE(cv$imputed_rho)
      if (cfg$measure == "ROM" && (tv$mean <= 0 || cv$mean <= 0)) {
        excl[[length(excl) + 1L]] <- data.frame(study_id = label, reason = "ROM needs positive means"); next }
      if (cfg$measure == "ROM") {
        gy <- geary_check(tv, cv, cfg)
        det$geary_min <- gy$g_min
        if (gy$exclude) { excl[[length(excl) + 1L]] <- data.frame(study_id = label, reason = gy$msg); next }
        notes <- c(notes, gy$msg)
      }
      es <- escalc(measure = cfg$measure, m1i = tv$mean, sd1i = tv$sd, n1i = tv$n,
                   m2i = cv$mean, sd2i = cv$sd, n2i = nc,
                   correct = if (cfg$measure == "ROM") cfg$rom_bias_correction else TRUE,
                   vtype = if (cfg$measure == "SMD") cfg$smd_variance else "LS")
      yi <- as.numeric(es$yi); vi <- as.numeric(es$vi)
      notes <- c(notes, switch(cfg$measure, MD = "escalc MD", SMD = paste0("escalc SMD (Hedges g, exact J, vtype ", cfg$smd_variance, ")"),
                                 ROM = paste0("escalc ROM", if (cfg$rom_bias_correction) " (bias-corrected)" else " (no bias correction)")))
      det[c("n_t", "mean_t", "sd_t", "n_c", "mean_c", "sd_c", "basis")] <- list(tv$n, tv$mean, tv$sd, nc, cv$mean, cv$sd, basis)
      det$prec_v <- 1 / tv$n + 1 / nc
    } else {  # crossover, approximate analysis from per-period summaries (Handbook 23.2.7)
      tv <- pool_arms(g); cv <- pool_arms(ctl)
      if (tv$n != cv$n) {
        # min(n_e, n_c) is not the number of people with both periods (Handbook 23.2.7.1 needs
        # complete pairs), so it is not guessed; a reported paired MD and SE (above) needs no N.
        excl[[length(excl) + 1L]] <- data.frame(study_id = label, reason = sprintf(CROSSOVER_N_DIFFER, tv$n, cv$n)); next }
      need_means <- !(cfg$measure == "SMD" && ok(pmd))   # SMD numerator may be the reported paired MD
      if (!ok(tv$sd) || !ok(cv$sd) || (need_means && (!ok(tv$mean) || !ok(cv$mean)))) {
        excl[[length(excl) + 1L]] <- data.frame(study_id = label, reason = "crossover: need paired_md/paired_se (MD) or per-condition mean and SD"); next }
      rh <- range_hit(rbind(g, ctl), "final")
      if (length(rh)) return(bad(sprintf(RANGE_ONLY_MSG, paste(rh, collapse = ", "))))
      n <- tv$n
      rw_row <- first_ok(if (!is.null(trow)) trow$corr_within else NA, ctl$corr_within[1])
      rw <- if (ok(rw_row)) rw_row else rho_w; imp_w <- !ok(rw_row); rw_src <- if (imp_w) " (default)" else " (study)"
      zero_sd <- !ok(rw_row) && ok(pse) && (isTRUE(tv$sd == 0) || isTRUE(cv$sd == 0))
      if (zero_sd && cfg$measure != "MD") {
        excl[[length(excl) + 1L]] <- data.frame(study_id = label, reason = sprintf("crossover %s: a period SD is 0, so the within-person r cannot be derived from the paired SE", cfg$measure)); next }
      if (!ok(rw_row) && ok(pse) && !zero_sd) {
        # A reported paired SE with both period SDs fixes the within-person correlation:
        # SD_diff = SE sqrt(N), r = (SDe^2 + SDc^2 - SD_diff^2) / (2 SDe SDc) (Elbourne et al.
        # 2002; Handbook 23.2.7). Using it beats an assumed r.
        rd <- (tv$sd^2 + cv$sd^2 - pse^2 * n) / (2 * tv$sd * cv$sd)
        if (!ok(rd) || rd < -1 || rd > 1) {
          excl[[length(excl) + 1L]] <- data.frame(study_id = label, reason = sprintf("within-person r derived from paired SE is %.3g, outside [-1, 1]", rd)); next }
        rw <- rd; imp_w <- FALSE; rw_src <- " (derived from paired SE)"
      }
      if (cfg$measure == "ROM") {
        gy <- geary_check(tv, cv, cfg)
        det$geary_min <- gy$g_min
        if (gy$exclude) { excl[[length(excl) + 1L]] <- data.frame(study_id = label, reason = gy$msg); next }
        notes <- c(notes, gy$msg)
      }
      notes <- c(notes, tv$notes, cv$notes, sprintf("within-person r=%.2f%s", rw, rw_src))
      if (cfg$measure == "MD") {
        yi <- tv$mean - cv$mean; vi <- (tv$sd^2 + cv$sd^2 - 2 * rw * tv$sd * cv$sd) / n
        notes <- c(notes, "crossover MD, SE from Handbook 23.2.7 (plain R)")
        if (zero_sd) { vi <- pse^2; notes <- c(notes, "variance from the reported paired SE (a period SD is 0)") }
      } else if (cfg$measure == "SMD") {
        sdp <- sqrt((tv$sd^2 + cv$sd^2) / 2)
        md <- if (ok(pmd)) pmd else tv$mean - cv$mean
        if (ok(pmd)) notes <- c(notes, "SMD numerator = reported paired MD")
        jx <- if (cfg$crossover_smd_j) hedges_j(n - 1) else 1
        yi <- jx * md / sdp; vi <- (1 / n + yi^2 / (2 * n)) * 2 * (1 - rw)
        notes <- c(notes, paste0("crossover SMD = J x MD/SD(measurements), SE from Handbook 23.2.7, ",
                                 if (cfg$crossover_smd_j) "J on N - 1 df" else "no small-sample correction", " (plain R)"))
      } else {
        if (tv$mean <= 0 || cv$mean <= 0) { excl[[length(excl) + 1L]] <- data.frame(study_id = label, reason = "ROM needs positive means"); next }
        yi <- log(tv$mean / cv$mean)
        vi <- tv$sd^2 / (n * tv$mean^2) + cv$sd^2 / (n * cv$mean^2) - 2 * rw * tv$sd * cv$sd / (n * tv$mean * cv$mean)
        notes <- c(notes, "crossover log ROM, delta-method paired variance (plain R)")
        if (cfg$rom_bias_correction) {
          # Lajeunesse 2015 term for each period mean, E[ln Xbar] ~ ln mu - sigma^2/(2 N mu^2);
          # expectation is linear, so r_within does not enter it (as metaproof.effects.crossover_rom).
          yi <- yi + 0.5 * (tv$sd^2 / (n * tv$mean^2) - cv$sd^2 / (n * cv$mean^2))
          notes <- c(notes, "ROM bias-corrected (Lajeunesse 2015, per period mean; independent of r_within)")
        }
      }
      det[c("n_t", "mean_t", "sd_t", "n_c", "mean_c", "sd_c", "r_within", "basis")] <- list(n, tv$mean, tv$sd, n, cv$mean, cv$sd, rw, "final")
      # The part of the crossover SMD variance that does not depend on g, 2(1 - r)/N (Handbook
      # 23.2.7; Elbourne et al. 2002), is the precision term of the modified Egger test; the
      # parallel term 1/n1 + 1/n2 would make a crossover look 1/(1 - r) times less precise.
      det$prec_v <- 2 * (1 - rw) / n
    }
    if (!ok(yi) || !ok(vi) || vi <= 0) {
      excl[[length(excl) + 1L]] <- data.frame(study_id = label, reason = "effect or variance not computable (check inputs)"); next }
    mods <- lapply(mod_vars, function(v) moderator_value(v, g, rows))
    names(mods) <- mod_vars
    comps[[length(comps) + 1L]] <- c(list(study_id = label, source_study = sid, yi = yi, vi = vi,
      method_notes = paste(unique(notes[!is.na(notes) & nzchar(notes)]), collapse = "; "),
      imputed_prepost = imp_pp, imputed_within = imp_w, design = design), det, mods)
  }
  list(comp = comps, excl = if (length(excl)) do.call(rbind, excl) else NULL)
}

build_effects <- function(d, cfg, rho_pp, rho_w, mod_vars) {
  comps <- list(); excl <- list()
  for (s in unique(d$study_id)) {
    res <- process_study(d[d$study_id == s, , drop = FALSE], cfg, rho_pp, rho_w, mod_vars)
    comps <- c(comps, res$comp); if (!is.null(res$excl)) excl[[length(excl) + 1L]] <- res$excl
  }
  eff <- if (length(comps)) do.call(rbind, lapply(comps, function(z) as.data.frame(z, stringsAsFactors = FALSE, check.names = FALSE))) else NULL
  if (!is.null(eff)) eff$sei <- sqrt(eff$vi)
  list(eff = eff, excl = if (length(excl)) do.call(rbind, excl) else NULL)
}

# ---------------------------------------------------------------------------
# Pooling
# ---------------------------------------------------------------------------
# metafor stops iterating when the change in tau^2 is below control$threshold
# (default 1e-5) and confint() solves the Q-profile equations with uniroot at
# tol = .Machine$double.eps^0.25 (about 1.2e-4). Both defaults leave tau^2 and
# its CI limits off the exact solution by up to ~1e-5, which is larger than the
# 1e-6 / 1e-4 tolerances of validation/crosscheck.py. Checked on dat.bcg ~ ablat:
# default threshold gives tau^2 = 0.0763547, threshold 1e-10 gives 0.0763480,
# and direct maximisation of the REML likelihood gives 0.0763480.
# tol also sets the uniroot tolerance metafor uses for the Paule-Mandel estimator (default
# about 1.2e-4, which left PM tau^2 off by ~1e-5 on dat.bcg).
RMA_CONTROL <- list(threshold = 1e-10, maxiter = 1000, tol = 1e-12)
# confint() searches the Q-profile limits only up to tau2.max (default max(100, 10 tau^2))
# and reports tau2.max itself when the root lies beyond it (dat.normand1999, DL: 2054.09
# reported, true limit 2889.45). Q_gen(t) <= sum((y - mean(y))^2) / t, so the bound below is
# always above the upper limit for any level >= 99.98%.
ci_control <- function(f) {
  df <- f$k - f$p
  list(tol = 1e-10, tau2.max = 100 + 10 * sum((f$yi - mean(f$yi))^2) / qchisq(1e-4, df) + 10 * max(f$vi))
}

# test for one fit under cfg$hksj_rule. revman needs tau^2 first, so it fits with z and
# refits with knha when tau^2 > 0 and k - p >= 2 (Handbook v6.5 10.10.4.4).
# The Paule-Mandel root is searched only up to control$tau2.max, whose metafor default
# max(100, 10 mad(yi)^2) one outlying study can exceed ("Value of 'tau2.max' too low").
# Q_gen(t) <= sum((y - mean(y))^2) / t for any design with an intercept, so the PM root is at
# most SS/(k - p) and the bound below always contains it (Paule & Mandel 1982).
rma_control <- function(yi) c(RMA_CONTROL, list(tau2.max = 100 + 10 * sum((yi - mean(yi))^2)))
# The reason the last fit_rma call failed, for the note of an all-NA pooled row.
FIT_STATUS <- new.env()
fit_rma <- function(yi, vi, est, cfg, mods = NULL, slab = NULL, rule = cfg$hksj_rule, intercept = TRUE) {
  hk <- cfg$hksj && est != "FE"
  test <- if (!hk) "z" else switch(rule, always = "knha", truncated = "adhoc", revman = "z")
  ctl <- rma_control(yi)
  FIT_STATUS$error <- NULL
  one <- function(test) {
    args <- list(yi = yi, vi = vi, method = est, test = test, level = cfg$level, control = ctl)
    if (!is.null(mods)) { args$mods <- mods; args$intercept <- intercept }
    if (!is.null(slab)) args$slab <- slab
    f <- tryCatch(do.call(rma, args), error = function(e) NULL)
    if (is.null(f)) {  # retry with damped Fisher scoring for convergence problems
      args$control <- c(ctl, list(stepadj = 0.5))
      f <- tryCatch(do.call(rma, args), error = function(e) {
        FIT_STATUS$error <- conditionMessage(e); message("  rma failed (", est, "): ", conditionMessage(e)); NULL })
    }
    f
  }
  f <- one(test)
  if (hk && rule == "revman" && !is.null(f) && f$k - f$p >= 2 && f$tau2 > TAU2_ZERO * mean(vi)) f <- one("knha")
  f
}

ci_method <- function(f) {
  if (is.null(f) || f$test == "z") return("Wald z")
  paste0(if (f$test == "adhoc") "HKSJ truncated" else "HKSJ", " t(", f$k - f$p, ")")
}

na_pooled <- function(est, k, note = NA_character_) {
  r <- as.data.frame(setNames(as.list(rep(NA_real_, length(POOLED_COLS))), POOLED_COLS))
  r$estimator <- est; r$k <- k
  for (c in POOLED_EXTRA) r[[c]] <- NA
  r$ci_method <- "not pooled"; r$note <- note; r
}

k_studies_of <- function(eff) if (is.null(eff)) 0L else length(unique(eff$source_study))

# Prediction interval written to pooled.csv (and drawn by forest.png, see forest_args).
pred_interval <- function(f, est, cfg) {
  pi <- c(NA_real_, NA_real_)
  if (est != "FE" && cfg$prediction_interval) {
    if (cfg$pi_type == "riley") {
      if (f$k >= 3) { p <- predict(f, predtype = "Riley"); pi <- c(p$pi.lb, p$pi.ub) }
    } else if (cfg$pi_type == "revman") {
      # t on k - 1 df whatever the CI method (the form the review reports for RevMan and
      # Handbook v6.5 10.10.4.3), with the SE used for the CI.
      if (f$k >= 2) {
        h <- qt(1 - (1 - cfg$level / 100) / 2, f$k - 1) * sqrt(f$tau2 + f$se[1]^2)
        pi <- c(f$b[1] - h, f$b[1] + h)
      }
    } else { p <- predict(f); pi <- c(p$pi.lb, p$pi.ub) }
  }
  as.numeric(pi)
}
# metafor predtype that reproduces each pi_type: Riley is t on k - 2 df, "t" is t on k - 1 df
# with the fit's SE (equal to the revman formula above), "default" is metafor's own.
PI_PREDTYPE <- c(riley = "Riley", revman = "t", metafor = "default")

pooled_row <- function(f, est, cfg, k_st = NA, fail = NULL) {
  if (is.null(f)) return(na_pooled(est, NA, if (is.null(fail)) NA_character_ else paste0("fit failed: ", fail)))
  pi <- pred_interval(f, est, cfg)
  tci <- c(NA_real_, NA_real_)
  if (est != "FE" && f$k >= 2) {
    ci <- tryCatch(confint(f, control = ci_control(f)), error = function(e) NULL)
    if (!is.null(ci)) tci <- unname(ci$random["tau^2", c("ci.lb", "ci.ub")])
  }
  # The three intervals hksj_rule chooses between, always written (protocol K6).
  alt <- function(rule, hk = TRUE) {
    g <- fit_rma(as.numeric(f$yi), f$vi, est, modifyList(cfg, list(hksj = hk)), rule = rule)
    if (is.null(g)) c(NA_real_, NA_real_) else c(g$ci.lb[1], g$ci.ub[1])
  }
  wald <- alt("always", hk = FALSE)
  hk <- if (est != "FE" && f$k >= 2) alt("always") else c(NA_real_, NA_real_)
  tr <- if (est != "FE" && f$k >= 2) alt("truncated") else c(NA_real_, NA_real_)
  notes <- character(0)
  if (est != "FE" && ok(k_st) && k_st < PI_MIN_K) {
    if (cfg$prediction_interval && ok(pi[1])) notes <- c(notes, sprintf("prediction interval not interpretable with %d studies (< 5)", k_st))
    notes <- c(notes, sprintf("tau2 CI uninformative with %d studies (< 5)", k_st))
  }
  data.frame(estimator = est, k = f$k, estimate = as.numeric(f$b[1]), se = f$se[1],
             ci_lb = f$ci.lb[1], ci_ub = f$ci.ub[1], pi_lb = pi[1], pi_ub = pi[2],
             tau2 = if (est == "FE") 0 else f$tau2, tau2_ci_lb = tci[1], tau2_ci_ub = tci[2],
             I2 = f$I2, H2 = f$H2, Q = f$QE, Q_p = f$QEp, k_studies = k_st, ci_method = ci_method(f),
             df = if (f$test == "z") NA_real_ else f$k - f$p, wald_ci_lb = wald[1], wald_ci_ub = wald[2],
             hksj_ci_lb = hk[1], hksj_ci_ub = hk[2], hksj_trunc_ci_lb = tr[1], hksj_trunc_ci_ub = tr[2],
             stat = f$zval[1], pval = f$pval[1],
             note = if (length(notes)) paste(notes, collapse = "; ") else NA_character_)
}

# Leave one study out (all comparisons of a split multi-arm trial together), refitting with
# the same estimator and hksj_rule, so the result matches what the pooled analysis would show.
leave1out_df <- function(f, eff, cfg) {
  empty <- data.frame(study_id = character(0), k = integer(0), estimate = numeric(0), se = numeric(0),
                      ci_lb = numeric(0), ci_ub = numeric(0), pval = numeric(0), tau2 = numeric(0),
                      I2 = numeric(0), H2 = numeric(0), Q = numeric(0), Q_p = numeric(0))
  if (is.null(f) || k_studies_of(eff) < 3) return(empty)
  rows <- lapply(unique(eff$source_study), function(s) {
    keep <- eff$source_study != s
    g <- fit_rma(eff$yi[keep], eff$vi[keep], cfg$estimator, cfg)
    if (is.null(g)) return(NULL)
    data.frame(study_id = s, k = g$k, estimate = as.numeric(g$b[1]), se = g$se[1], ci_lb = g$ci.lb[1],
               ci_ub = g$ci.ub[1], pval = g$pval[1], tau2 = if (cfg$estimator == "FE") 0 else g$tau2,
               I2 = g$I2, H2 = g$H2, Q = g$QE, Q_p = g$QEp)
  })
  do.call(rbind, c(list(empty), rows))
}

influence_df <- function(f, eff) {
  empty <- data.frame(study_id = character(0), rstudent = numeric(0), dffits = numeric(0), cooks_d = numeric(0),
                      cov_r = numeric(0), tau2_del = numeric(0), QE_del = numeric(0), hat = numeric(0),
                      weight = numeric(0), dfbetas = numeric(0), influential = logical(0))
  if (is.null(f) || f$k < 3) return(empty)
  inf <- tryCatch(influence(f), error = function(e) NULL)
  if (is.null(inf)) return(empty)
  data.frame(study_id = eff$study_id, rstudent = inf$inf$rstudent, dffits = inf$inf$dffits, cooks_d = inf$inf$cook.d,
             cov_r = inf$inf$cov.r, tau2_del = inf$inf$tau2.del, QE_del = inf$inf$QE.del, hat = inf$inf$hat,
             weight = inf$inf$weight, dfbetas = inf$dfbs[[1]], influential = as.logical(inf$is.infl))
}

subgroups_df <- function(eff, cfg) {
  rows <- list()
  zc <- qnorm(1 - (1 - cfg$level / 100) / 2)
  for (v in cfg$subgroups) {
    if (!v %in% names(eff)) { rows[[length(rows) + 1L]] <- data.frame(moderator = v, level = NA, k = NA, note = "column not in data"); next }
    x <- eff[[v]]; use <- !is.na(x)
    # Levels in code-point (C-locale) order, as Python's sorted(), so both engines pick the same
    # reference level whatever the session collation ('Older adults' sorts before 'adults').
    lv <- sort(unique(as.character(x[use])), method = "radix"); ests <- c(); ses <- c()
    k_st <- length(unique(eff$source_study[use]))
    per <- vapply(lv, function(l) length(unique(eff$source_study[use & as.character(x) == l])), 0L)
    # Protocol I: at least subgroup_min_k studies, and subgroup_min_per_level in every subgroup.
    if (k_st < cfg$subgroup_min_k) {
      rows[[length(rows) + 1L]] <- data.frame(moderator = v, level = NA, k = sum(use),
        note = sprintf("not run: %d studies < subgroup_min_k = %d (protocol I)", k_st, cfg$subgroup_min_k)); next }
    if (length(per) && min(per) < cfg$subgroup_min_per_level) {
      rows[[length(rows) + 1L]] <- data.frame(moderator = v, level = NA, k = sum(use),
        note = sprintf("not run: a subgroup has %d studies < subgroup_min_per_level = %d (protocol I)", min(per), cfg$subgroup_min_per_level)); next }
    if (cfg$subgroup_method == "common") {
      # One model with a dummy per subgroup and no intercept: common tau^2; the test for
      # subgroup differences is the omnibus test of the same model with an intercept.
      fx <- factor(as.character(x[use]), levels = lv)
      X0 <- model.matrix(~ fx - 1)
      f0 <- fit_rma(eff$yi[use], eff$vi[use], cfg$estimator, cfg, mods = X0, intercept = FALSE)
      if (is.null(f0)) { rows[[length(rows) + 1L]] <- data.frame(moderator = v, level = NA, k = sum(use), note = "model did not fit"); next }
      for (j in seq_along(lv))
        rows[[length(rows) + 1L]] <- data.frame(moderator = v, level = lv[j], k = sum(fx == lv[j]),
          estimate = as.numeric(f0$b[j]), se = f0$se[j], ci_lb = f0$ci.lb[j], ci_ub = f0$ci.ub[j],
          pi_lb = NA, pi_ub = NA, tau2 = if (cfg$estimator == "FE") 0 else f0$tau2, I2 = NA, Q = NA, Q_p = NA,
          note = "common tau2 (meta-regression on subgroup dummies)", test = NA, test_stat = NA, test_df1 = NA, test_df2 = NA, test_p = NA)
      if (length(lv) >= 2L && sum(use) > length(lv)) {
        fm <- fit_rma(eff$yi[use], eff$vi[use], cfg$estimator, cfg, mods = model.matrix(~ fx)[, -1, drop = FALSE])
        if (!is.null(fm)) rows[[length(rows) + 1L]] <- data.frame(moderator = v, level = "ALL", k = sum(use), estimate = NA, se = NA,
          ci_lb = NA, ci_ub = NA, pi_lb = NA, pi_ub = NA, tau2 = fm$tau2, I2 = fm$I2, Q = fm$QE, Q_p = fm$QEp,
          note = "test for subgroup differences, common tau2 (omnibus moderator test)",
          test = "QM_common_tau2", test_stat = fm$QM, test_df1 = fm$QMdf[1], test_df2 = fm$QMdf[2], test_p = fm$QMp)
      }
      next
    }
    for (l in lv) {
      e <- eff[use & as.character(x) == l, , drop = FALSE]
      if (nrow(e) == 1L) {
        r <- data.frame(estimate = e$yi, se = e$sei, ci_lb = e$yi - zc * e$sei, ci_ub = e$yi + zc * e$sei,
                        pi_lb = NA, pi_ub = NA, tau2 = NA, I2 = NA, Q = NA, Q_p = NA, note = "single study (z CI)")
      } else {
        f <- fit_rma(e$yi, e$vi, cfg$estimator, cfg); pr <- pooled_row(f, cfg$estimator, cfg)
        r <- data.frame(pr[c("estimate", "se", "ci_lb", "ci_ub", "pi_lb", "pi_ub", "tau2", "I2", "Q", "Q_p")], note = "separate tau2 per subgroup")
        # The test for subgroup differences weights each subgroup estimate by the inverse of
        # its model-based variance 1/sum(w*) (Borenstein et al. 2009, ch. 19), not the
        # HKSJ-scaled variance, because Q_between is referred to a chi-square distribution.
        fz <- fit_rma(e$yi, e$vi, cfg$estimator, modifyList(cfg, list(hksj = FALSE)))
        r$se_model <- if (is.null(fz)) NA_real_ else fz$se[1]
      }
      if (is.null(r$se_model)) r$se_model <- r$se
      ests <- c(ests, r$estimate); ses <- c(ses, r$se_model)
      r$se_model <- NULL
      rows[[length(rows) + 1L]] <- data.frame(moderator = v, level = l, k = nrow(e), r, test = NA, test_stat = NA,
                                              test_df1 = NA, test_df2 = NA, test_p = NA)
    }
    if (length(lv) >= 2L && all(is.finite(ses))) {
      tb <- rma(ests, sei = ses, method = "FE")
      rows[[length(rows) + 1L]] <- data.frame(moderator = v, level = "ALL", k = sum(use), estimate = NA, se = NA, ci_lb = NA,
        ci_ub = NA, pi_lb = NA, pi_ub = NA, tau2 = NA, I2 = NA, Q = NA, Q_p = NA,
        note = "test for subgroup differences: Q between subgroup estimates (separate tau2, model-based SEs)",
        test = "Q_between_separate_tau2", test_stat = tb$QE, test_df1 = length(lv) - 1L, test_df2 = NA, test_p = tb$QEp)
      fx <- factor(as.character(x[use]), levels = lv)
      fm <- fit_rma(eff$yi[use], eff$vi[use], cfg$estimator, cfg, mods = model.matrix(~ fx)[, -1, drop = FALSE])
      if (!is.null(fm)) rows[[length(rows) + 1L]] <- data.frame(moderator = v, level = "ALL", k = sum(use), estimate = NA, se = NA,
        ci_lb = NA, ci_ub = NA, pi_lb = NA, pi_ub = NA, tau2 = fm$tau2, I2 = fm$I2, Q = fm$QE, Q_p = fm$QEp,
        note = "omnibus moderator test, common tau2 (meta-regression on factor)",
        test = "QM_common_tau2", test_stat = fm$QM, test_df1 = fm$QMdf[1], test_df2 = fm$QMdf[2], test_p = fm$QMp)
    }
  }
  cols <- c("moderator", "level", "k", "estimate", "se", "ci_lb", "ci_ub", "pi_lb", "pi_ub", "tau2", "I2", "Q", "Q_p",
            "test", "test_stat", "test_df1", "test_df2", "test_p", "note")
  if (!length(rows)) return(as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols)))
  out <- do.call(rbind, lapply(rows, function(r) { for (c in setdiff(cols, names(r))) r[[c]] <- NA; r[cols] }))
  out
}

metareg_df <- function(eff, cfg) {
  sets <- as.list(cfg$moderators)
  if (cfg$metareg_joint && length(cfg$moderators) > 1L) sets <- c(sets, list(cfg$moderators))
  rows <- list()
  for (s in sets) {
    lab <- paste(s, collapse = "+")
    if (!all(s %in% names(eff))) { rows[[length(rows) + 1L]] <- data.frame(model = lab, note = "column not in data"); next }
    X <- eff[, s, drop = FALSE]
    for (c in s) { nv <- num(X[[c]]); if (all(is.na(nv) == is.na(X[[c]]))) X[[c]] <- nv else
      X[[c]] <- factor(X[[c]], levels = sort(unique(as.character(X[[c]])), method = "radix")) }  # C-locale order
    use <- stats::complete.cases(X)
    e <- eff[use, , drop = FALSE]; Xu <- X[use, , drop = FALSE]
    k_st <- length(unique(e$source_study))
    # metareg_min_k studies per moderator (Handbook 10.11.5.1: about ten per characteristic); a
    # categorical moderator counts once however many levels it has (as metaproof.diagnostics).
    n_mod <- length(s); need <- cfg$metareg_min_k * n_mod
    if (k_st < need) { rows[[length(rows) + 1L]] <- data.frame(model = lab, k = nrow(e),
      note = sprintf("not run: %d studies < %s (Handbook 10.11.5.1: about ten studies per characteristic; protocol I)", k_st,
                     if (n_mod == 1L) sprintf("metareg_min_k = %d", cfg$metareg_min_k)
                     else sprintf("metareg_min_k x %d moderators = %d x %d = %d", n_mod, cfg$metareg_min_k, n_mod, need))); next }
    cats <- s[vapply(s, function(c) is.factor(Xu[[c]]), TRUE)]
    lvk <- lapply(cats, function(c) tapply(e$source_study, droplevels(Xu[[c]]), function(z) length(unique(z))))
    names(lvk) <- cats
    min_level <- if (length(lvk)) min(unlist(lvk)) else NA_real_
    mm <- tryCatch(model.matrix(as.formula(paste("~", paste0("`", s, "`", collapse = "+"))), data = Xu), error = function(e) NULL)
    if (is.null(mm)) { rows[[length(rows) + 1L]] <- data.frame(model = lab, k = nrow(e), note = "moderator has fewer than 2 distinct values"); next }
    if (nrow(e) <= ncol(mm)) { rows[[length(rows) + 1L]] <- data.frame(model = lab, k = nrow(e), note = "too few studies for this model"); next }
    # A constant or collinear moderator: metafor would drop the redundant column and the
    # term labels below would no longer match the coefficients, so stop here as metaproof does.
    if (qr(mm)$rank < ncol(mm)) { rows[[length(rows) + 1L]] <- data.frame(model = lab, k = nrow(e), note = "not estimable: collinear moderators (or a constant moderator)"); next }
    f <- fit_rma(e$yi, e$vi, cfg$estimator, cfg, mods = mm[, -1, drop = FALSE])
    if (is.null(f)) { rows[[length(rows) + 1L]] <- data.frame(model = lab, k = nrow(e), note = "model did not fit"); next }
    notes <- character(0)
    if (k_st < 10 * n_mod) notes <- c(notes, sprintf("%d studies for %d moderator(s), fewer than 10 per characteristic: meta-regression is unreliable (Handbook 10.11.5.1)", k_st, n_mod))
    sparse <- names(lvk)[vapply(lvk, function(z) min(z) < 2, TRUE)]
    if (length(sparse)) notes <- c(notes, sprintf("a category of %s has 1 study: its coefficient rests on that study alone", paste(sparse, collapse = ", ")))
    note <- if (length(notes)) paste(notes, collapse = "; ") else NA
    # Pseudo-R2 from the intercept-only fit on the same studies, with R2 = 0 when that tau2 is 0
    # (below TAU2_ZERO x mean vi, the rule metaproof shares): metafor stops a boundary REML fit at
    # tau2 ~ 1e-12 to 1e-6 instead of 0 and then reports an R2 of 80-100% for nothing explained.
    f0 <- if (cfg$estimator == "FE") NULL else fit_rma(e$yi, e$vi, cfg$estimator, cfg)
    R2 <- if (cfg$estimator == "FE" || is.null(f0)) NA_real_ else if (f0$tau2 <= TAU2_ZERO * mean(e$vi)) 0 else
      max(0, 100 * (f0$tau2 - f$tau2) / f0$tau2)
    terms <- colnames(mm); terms[terms == "(Intercept)"] <- "intrcpt"
    rows[[length(rows) + 1L]] <- data.frame(model = lab, term = terms, k = f$k, estimate = as.numeric(f$b), se = f$se,
      statistic = f$zval, df = if (is.null(f$ddf) || all(is.na(f$ddf))) NA else f$ddf[1], pval = f$pval, ci_lb = f$ci.lb, ci_ub = f$ci.ub,
      tau2 = f$tau2, I2 = f$I2, R2 = R2, QM = f$QM, QM_df1 = f$QMdf[1],
      QM_df2 = f$QMdf[2], QM_p = f$QMp, QE = f$QE, QE_p = f$QEp, n_moderators = n_mod,
      k_per_param = k_st / (ncol(mm) - 1L), min_level_k = min_level, note = note)
  }
  cols <- c("model", "term", "k", "estimate", "se", "statistic", "df", "pval", "ci_lb", "ci_ub", "tau2", "I2", "R2",
            "QM", "QM_df1", "QM_df2", "QM_p", "QE", "QE_p", "n_moderators", "k_per_param", "min_level_k", "note")
  if (!length(rows)) return(as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols)))
  do.call(rbind, lapply(rows, function(r) { for (c in setdiff(cols, names(r))) r[[c]] <- NA; r[cols] }))
}

small_study_df <- function(f, cfg, eff) {
  cols <- c("test", "model", "predictor", "k", "k_studies", "coef_sei", "coef_se", "statistic", "df", "pval",
            "limit_estimate", "limit_ci_lb", "limit_ci_ub", "note")
  k <- if (is.null(f)) 0L else f$k
  k_st <- k_studies_of(eff)
  if (is.null(f) || k_st < cfg$egger_min_k) {
    r <- as.data.frame(setNames(as.list(rep(NA, length(cols))), cols))
    r$test <- "egger"; r$k <- k; r$k_studies <- k_st; r$note <- sprintf("not run: %d studies < %d", k_st, cfg$egger_min_k); return(r)
  }
  # For SMD the classic test is biased (g and its SE are correlated); use the g-free part of
  # the SE as the precision term (Pustejovsky & Rodgers 2019): sqrt(1/n1 + 1/n2) for a
  # parallel comparison, sqrt(2(1 - r)/N) for a crossover (prec_v, set in process_study).
  modified <- cfg$measure == "SMD" && all(ok_vec(eff$prec_v) & num(eff$prec_v) > 0)
  one <- function(model, label) {
    rt <- tryCatch(if (label == "egger" && modified) regtest(eff$yi, vi = eff$prec_v, model = "lm", predictor = "sei", ret.fit = TRUE)
                   else regtest(f, model = model, predictor = "sei", ret.fit = TRUE), error = function(e) NULL)
    if (is.null(rt)) return(data.frame(test = label, model = model, predictor = "sei", k = k, k_studies = k_st, coef_sei = NA, coef_se = NA,
      statistic = NA, df = NA, pval = NA, limit_estimate = NA, limit_ci_lb = NA, limit_ci_ub = NA, note = "regtest failed"))
    if (model == "lm") { cf <- coef(summary(rt$fit)); b <- cf[nrow(cf), 1]; bse <- cf[nrow(cf), 2]
    } else { b <- as.numeric(rt$fit$b[nrow(rt$fit$b)]); bse <- rt$fit$se[length(rt$fit$se)] }
    data.frame(test = label, model = model, predictor = if (label == "egger" && modified) "sqrt(1/n1+1/n2)" else "sei",
               k = k, k_studies = k_st, coef_sei = b, coef_se = bse,
               statistic = rt$zval, df = if (is.null(rt$dfs)) NA else rt$dfs, pval = rt$pval,
               limit_estimate = rt$est, limit_ci_lb = rt$ci.lb, limit_ci_ub = rt$ci.ub,
               note = if (label == "egger" && modified) paste0("modified Egger test for SMD, precision sqrt(1/n1 + 1/n2) (Pustejovsky & Rodgers 2019)",
                                                    if (any(eff$design == "crossover")) "; crossover rows sqrt(2(1 - r)/N) (Handbook 23.2.7)" else "")
                      else if (model == "lm") "classic Egger: weighted regression of yi on SE, multiplicative dispersion (Egger 1997)"
                      else "mixed-effects regression on SE, same estimator/test as primary fit")
  }
  out <- rbind(one("lm", "egger"), one("rma", "egger_rma"))
  if (cfg$measure == "SMD") out$note[out$test == "egger_rma"] <- paste(out$note[out$test == "egger_rma"], "| SMD: SE correlates with g by construction (Pustejovsky & Rodgers 2019)")
  out
}

measure_label <- function(m) switch(m, MD = "Mean difference", SMD = "Standardized mean difference (Hedges' g)", ROM = "Ratio of means")

# forest() arguments: the drawn CI uses cfg$level, and the PI line is the pooled.csv interval
# (same predtype, drawn only when pooled.csv has one), so the figure matches the table
# (Riley et al. 2011; Handbook v6.5 10.10.4.3).
# With fewer than PI_MIN_K studies the PI is not drawn (pooled.csv keeps it with its caveat),
# so the published figure carries the same caveat as the table (Handbook v6.5 10.10.4.3).
forest_pi_note <- function(f, eff, cfg, pooled_ok) {
  if (!pooled_ok || is.null(f)) return(NULL)
  k_st <- if ("source_study" %in% names(eff)) length(unique(eff$source_study)) else nrow(eff)
  if (k_st >= PI_MIN_K || !ok(pred_interval(f, cfg$estimator, cfg)[1])) return(NULL)
  sprintf("prediction interval not shown: %d studies (< %d; Handbook 10.10.4.3); pooled.csv keeps it with that caveat", k_st, PI_MIN_K)
}
forest_args <- function(f, eff, cfg, pooled_ok) {
  header <- c("Study", sprintf("Estimate [%g%% CI]", cfg$level))
  if (pooled_ok && !is.null(f)) {
    has_pi <- ok(pred_interval(f, cfg$estimator, cfg)[1]) && is.null(forest_pi_note(f, eff, cfg, pooled_ok))
    args <- list(x = f, slab = eff$study_id, header = header, level = cfg$level,
                 xlab = measure_label(cfg$measure), addpred = has_pi,
                 mlab = sprintf("%s model (%s, %s)", if (cfg$estimator == "FE") "Common-effect" else "Random-effects",
                                cfg$estimator, ci_method(f)))
    if (has_pi) args$predtype <- PI_PREDTYPE[[cfg$pi_type]]
  } else {
    # Below min_k_pool: the study results only, no pooled diamond (protocol F).
    args <- list(x = eff$yi, vi = eff$vi, slab = eff$study_id, header = header, level = cfg$level,
                 xlab = measure_label(cfg$measure))
  }
  if (cfg$measure == "ROM") { args$atransf <- exp; args$refline <- 0 }
  args
}

make_plots <- function(f, eff, cfg, dir, outcome, pooled_ok) {
  if (is.null(eff)) return(invisible())
  k <- nrow(eff)
  png(file.path(dir, "forest.png"), width = 2000, height = max(900, 300 + 70 * k), res = 200)
  op <- par(mar = c(4, 4, 2, 2))
  tryCatch({
    do.call(forest, forest_args(f, eff, cfg, pooled_ok))
    title(if (pooled_ok) outcome else paste0(outcome, " (not pooled)"))
    pin <- forest_pi_note(f, eff, cfg, pooled_ok)
    if (!is.null(pin)) mtext(pin, side = 1, line = 2.8, cex = 0.6, adj = 0)
  }, error = function(e) message("  forest plot failed: ", conditionMessage(e)))
  par(op); dev.off()
  png(file.path(dir, "funnel.png"), width = 1600, height = 1400, res = 200)
  # Contour-enhanced funnel centred on 0 (Peters et al. 2008): p 0.10, 0.05, 0.01 regions.
  tryCatch({ funnel(eff$yi, eff$vi, level = c(90, 95, 99), shade = c("white", "gray55", "gray75"), refline = 0,
                    legend = TRUE, xlab = measure_label(cfg$measure), atransf = if (cfg$measure == "ROM") exp else NULL); title(outcome) },
           error = function(e) message("  funnel plot failed: ", conditionMessage(e)))
  dev.off()
}

# ---------------------------------------------------------------------------
# Per-outcome driver
# ---------------------------------------------------------------------------
run_outcome <- function(d_all, outcome, cfg, dir) {
  message("outcome: ", outcome, if (cfg$outcome != outcome) paste0(" (rows of ", cfg$outcome, ")") else "")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  d <- d_all[!is.na(d_all$outcome) & d_all$outcome == cfg$outcome, , drop = FALSE]
  excl <- list()
  sel <- select_timepoint(d, cfg$timepoint, excl); d <- sel$d; excl <- sel$excl
  mod_vars <- unique(c(cfg$subgroups, cfg$moderators))
  be <- build_effects(d, cfg, cfg$corr_prepost, cfg$corr_within, mod_vars)
  eff <- be$eff
  excl_df <- do.call(rbind, c(excl, list(be$excl)))
  write_out(if (is.null(excl_df)) data.frame(study_id = character(0), reason = character(0)) else excl_df, file.path(dir, "excluded.csv"))
  writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(dir, "config_resolved.json"))

  se_cols <- c("study_id", "yi", "vi", "sei", "method_notes")
  if (is.null(eff)) {
    write_out(as.data.frame(setNames(replicate(5, character(0), simplify = FALSE), se_cols)), file.path(dir, "study_effects.csv"))
    write_out(na_pooled(cfg$estimator, 0L), file.path(dir, "pooled.csv"))
    message("  no usable comparisons"); return(invisible())
  }
  write_out(eff[, c(se_cols, "source_study")], file.path(dir, "study_effects.csv"))
  write_out(eff[, c("study_id", "source_study", "design", "basis", "n_t", "mean_t", "sd_t", "n_c", "mean_c", "sd_c", "r_within",
                    "geary_min", "imputed_prepost", "imputed_within", mod_vars[mod_vars %in% names(eff)])], file.path(dir, "study_inputs.csv"))

  k <- nrow(eff); k_st <- k_studies_of(eff)
  pooled_ok <- k_st >= cfg$min_k_pool
  fails <- list()
  fits <- lapply(cfg$report_estimators, function(e) {
    if (!pooled_ok) return(NULL)
    f <- fit_rma(eff$yi, eff$vi, e, cfg, slab = eff$study_id)
    if (is.null(f)) fails[[e]] <<- FIT_STATUS$error %||% "rma did not converge"
    f })
  names(fits) <- cfg$report_estimators
  note <- sprintf("not pooled: %d studies < min_k_pool = %d (protocol F)", k_st, cfg$min_k_pool)
  pooled <- do.call(rbind, lapply(cfg$report_estimators, function(e)
    if (pooled_ok) pooled_row(fits[[e]], e, cfg, k_st, fails[[e]]) else { r <- na_pooled(e, k, note); r$k_studies <- k_st; r }))
  write_out(pooled, file.path(dir, "pooled.csv"))
  f <- fits[[cfg$estimator]]

  write_out(leave1out_df(if (pooled_ok) f else NULL, eff, cfg), file.path(dir, "leave1out.csv"))
  write_out(influence_df(if (pooled_ok) f else NULL, eff), file.path(dir, "influence.csv"))
  write_out(subgroups_df(eff, cfg), file.path(dir, "subgroups.csv"))
  write_out(metareg_df(eff, cfg), file.path(dir, "metareg.csv"))
  fe_any <- if (is.null(f)) fit_rma(eff$yi, eff$vi, cfg$estimator, cfg) else f
  write_out(small_study_df(fe_any, cfg, eff), file.path(dir, "small_study.csv"))
  make_plots(f, eff, cfg, dir, outcome, pooled_ok)

  # Sensitivity to imputed correlations (only comparisons that used a default r change)
  sens <- list()
  if (identical(cfg$corr_within_sensitivity, "auto"))
    cfg$corr_within_sensitivity <- if (any(as.logical(eff$imputed_within))) AUTO_CORR_WITHIN_SENSITIVITY else numeric(0)
  grid <- rbind(data.frame(analysis = "primary", corr_prepost = cfg$corr_prepost, corr_within = cfg$corr_within),
                if (length(cfg$corr_prepost_sensitivity)) data.frame(analysis = paste0("corr_prepost=", cfg$corr_prepost_sensitivity),
                  corr_prepost = cfg$corr_prepost_sensitivity, corr_within = cfg$corr_within),
                if (length(cfg$corr_within_sensitivity)) data.frame(analysis = paste0("corr_within=", cfg$corr_within_sensitivity),
                  corr_prepost = cfg$corr_prepost, corr_within = cfg$corr_within_sensitivity))
  if (pooled_ok) for (j in seq_len(nrow(grid))) {
    b <- if (j == 1L) be else build_effects(d, cfg, grid$corr_prepost[j], grid$corr_within[j], mod_vars)
    fj <- if (!is.null(b$eff)) fit_rma(b$eff$yi, b$eff$vi, cfg$estimator, cfg) else NULL
    sens[[j]] <- cbind(grid[j, ], n_imputed_prepost = if (is.null(b$eff)) 0L else sum(b$eff$imputed_prepost),
                       n_imputed_within = if (is.null(b$eff)) 0L else sum(b$eff$imputed_within),
                       pooled_row(fj, cfg$estimator, cfg, k_studies_of(b$eff))[POOLED_COLS])
  }
  write_out(if (length(sens)) do.call(rbind, sens) else data.frame(analysis = character(0)), file.path(dir, "sensitivity.csv"))
  message(sprintf("  k = %d comparisons from %d studies, %s %s estimate = %.4f [%.4f, %.4f]", k, k_st, cfg$measure, cfg$estimator,
                  pooled$estimate[1], pooled$ci_lb[1], pooled$ci_ub[1]))
}

# Output folders. Names are sanitised as metaproof does (runs of other characters become one
# '_', leading and trailing '_' dropped) and refused when the folder would be hidden or outside
# --out, would take a run file's name, or would be shared by two analyses (also when the names
# differ only in case: macOS and Windows file systems ignore case).
RUN_FILES <- c("summary_across_analyses.csv", "bh_across_outcomes.csv", "extraction_check.txt", "run_manifest.json")
MARKER <- ".run_meta_analysis"   # inside every folder run_meta.R writes; only these are ever removed
STAGING <- ".run_meta-staging-"
safe_name <- function(x) { f <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9._-]+", "_", x)); if (nzchar(f)) f else "outcome" }
analysis_folders <- function(nms) {
  f <- vapply(nms, safe_name, "", USE.NAMES = FALSE)
  for (j in seq_along(nms)) {
    if (startsWith(f[j], ".") || tolower(f[j]) %in% RUN_FILES)
      stop(sprintf("analysis name '%s' gives the folder '%s', which %s; rename the analysis", nms[j], f[j],
                   if (startsWith(f[j], ".")) "would be hidden or outside --out" else "is a file the run writes"), call. = FALSE)
    dup <- which(tolower(f[seq_len(j - 1L)]) == tolower(f[j]))
    if (length(dup)) stop(sprintf("analyses '%s' and '%s' would write to the same folder '%s'; rename one of them",
                                  nms[dup[1]], nms[j], f[j]), call. = FALSE)
  }
  setNames(f, nms)
}

main <- function() {
  a <- parse_args(commandArgs(trailingOnly = TRUE))
  d <- read_extraction(a$data)
  if (length(attr(d, "missing_cols"))) message("note: columns absent from input, treated as NA: ", paste(attr(d, "missing_cols"), collapse = ", "))
  outcomes <- unique(d$outcome[!is.na(d$outcome)])
  cfgs <- read_config(a$config, outcomes)
  # A configured outcome that is not in the data, or a subgroup/moderator column that is not in
  # the file, stops the run (as metaproof does), so a typo cannot drop an analysis in one engine.
  gone <- names(cfgs)[!vapply(cfgs, function(cf) cf$outcome %in% outcomes, TRUE)]
  if (length(gone)) {
    msg <- vapply(gone, function(nm) {
      o <- cfgs[[nm]]$outcome; close <- agrep(o, outcomes, max.distance = 0.3, value = TRUE, ignore.case = TRUE)
      sprintf("'%s' selects outcome '%s'%s", nm, o, if (length(close)) sprintf(" (did you mean '%s'?)", close[1]) else "") }, "")
    stop("outcome not in the extraction file: ", paste(msg, collapse = "; "), call. = FALSE)
  }
  for (nm in names(cfgs)) {
    absent <- setdiff(unique(c(cfgs[[nm]]$subgroups, cfgs[[nm]]$moderators)), names(d))
    if (length(absent)) stop(sprintf("analysis '%s': subgroup or moderator columns %s are not in the extraction file",
                                     nm, paste(absent, collapse = ", ")), call. = FALSE)
  }
  folders <- analysis_folders(names(cfgs))
  dir.create(a$out, recursive = TRUE, showWarnings = FALSE)
  for (p in list.files(a$out, all.files = TRUE, full.names = TRUE, no.. = TRUE))
    if (startsWith(basename(p), STAGING)) unlink(p, recursive = TRUE)  # left by an interrupted run
  # Every analysis is written in full to its own staging folder; only when all have finished do
  # they replace the folders of an earlier run, so no plot or diagnostic from an earlier run
  # survives and a failure part-way leaves the earlier run whole (as metaproof).
  stages <- list()
  on.exit(for (st in stages) unlink(st, recursive = TRUE), add = TRUE)
  for (o in names(cfgs)) {
    stages[[o]] <- tempfile(STAGING, tmpdir = a$out)
    run_outcome(d, o, cfgs[[o]], stages[[o]])
    writeLines("written by run_meta.R; the whole folder is replaced on every run", file.path(stages[[o]], MARKER))
  }
  for (o in names(cfgs)) {
    final <- file.path(a$out, folders[[o]])
    if (file.exists(final) && !dir.exists(final)) stop(sprintf("%s is a file, not an analysis folder; move it and rerun", final), call. = FALSE)
    if (dir.exists(final) && length(list.files(final, all.files = TRUE, no.. = TRUE)) && !file.exists(file.path(final, MARKER)))
      stop(sprintf("%s already exists and run_meta.R did not write it; move it or choose another --out", final), call. = FALSE)
    if (dir.exists(final)) unlink(final, recursive = TRUE)
    if (!file.rename(stages[[o]], final)) stop(sprintf("could not move the finished analysis into %s", final), call. = FALSE)
  }
  # Folders of analyses no longer configured, only if run_meta.R wrote them.
  kept <- tolower(unname(folders))
  for (p in list.dirs(a$out, recursive = FALSE))
    if (file.exists(file.path(p, MARKER)) && !tolower(basename(p)) %in% kept) unlink(p, recursive = TRUE)
}

if (sys.nframe() == 0L) main()
