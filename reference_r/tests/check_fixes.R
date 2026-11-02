#!/usr/bin/env Rscript
# check_fixes.R -- regression checks for the review fixes in reference_r/run_meta.R and
# reference_r/tests/check_reference.R. Each block names the finding it covers.
#   CR4-1  Paule-Mandel tau^2 above metafor's default tau2.max; failed fits carry a note.
#   CR4-2  modified Egger test for SMD: crossover precision term 2(1 - r)/N.
#   CR4-3  check_reference.R never deletes the WORKDIR it is given.
#   CR4-4  rom_geary = flag writes the Geary note and geary_min.
#   CR4-5  forest.png draws the pooled.csv prediction interval at the configured level.
#   CR4-6  a paired SE on the NR row is kept when the paired MD is missing; on the control row
#          it is refused (CR1-3), as metaproof refuses it.
#   Integration (handoffs from the Python fixes, so both engines agree):
#     Astra 1  multi-period crossovers are excluded, not split.
#     Astra 2  unequal crossover period n is excluded unless a paired MD and SE are reported.
#     Astra 3  a crossover SMD is refused in a data_basis: change analysis.
#     Astra 10 an unused median/range does not exclude a study.
#     Astra 15 crossover ROM takes the Lajeunesse bias correction.
#     Astra 19 combining groups is unchanged by a common shift of the means.
#     Astra 11 / CR3-02 / CR3-06 cells that cannot be used stop read_extraction with the row.
#     Astra 18 metareg_min_k is per moderator; below 10 it needs metareg_allow_small_k.
#     CR2-2 R2 = 0 when the intercept-only tau2 is 0; CR2-3 levels in C-locale order.
#     CR3-03 / CR3-13 config: outcome typos, missing columns, outcome under defaults,
#              bh under outcomes, hksj_adhoc against an explicit hksj_rule.
#     CR3-04 / CR3-05 / Astra 5 / Astra 23 output folders: staging, collisions, dot names.
#     CR3-12 forest.png hides the PI below 5 studies; CR5-2 pooled.csv has stat and pval.
#     Astra 16 a label whose rows give two timepoint_weeks excludes the study.
#   CR3-01 timepoint labels with a hyphen ("12-month", "Wk-12") parse to positive weeks.
#   Astra 16  a numeric timepoint rule uses timepoint_weeks, as metaproof does.
#   Astra 25  a shared config: unknown keys, Python-only options and invalid values stop
#             run_meta.R; median_method defaults to luo2018 as in metaproof.
# Usage: Rscript reference_r/tests/check_fixes.R [WORKDIR]
#   RUN_META=<path> checks another copy of run_meta.R (used to show the checks fail on the
#   code before the fixes). Exit status 0 = all checks passed.

a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
TESTS <- dirname(normalizePath(f[1])); RREF <- dirname(TESTS)
.libPaths(c(file.path(RREF, "lib"), .libPaths()))
suppressPackageStartupMessages(library(metafor))
RUN_META <- Sys.getenv("RUN_META", file.path(RREF, "run_meta.R"))
RSCRIPT <- file.path(R.home("bin"), "Rscript")
args <- commandArgs(trailingOnly = TRUE)
WORK <- file.path(if (length(args)) args[1] else tempdir(), "nr_meta_fix_check")
unlink(WORK, recursive = TRUE); dir.create(WORK, recursive = TRUE)

n_fail <- 0L; n_pass <- 0L
check <- function(label, got, want, tol = 1e-8) {
  got <- suppressWarnings(as.numeric(got)); want <- as.numeric(want)
  good <- length(got) == length(want) && identical(is.na(got), is.na(want)) &&
    all(abs(got - want)[!is.na(want)] <= tol * pmax(1, abs(want[!is.na(want)])))
  if (!isTRUE(good)) {
    n_fail <<- n_fail + 1L
    cat(sprintf("FAIL %s (got %s; want %s)\n", label, paste(signif(got, 8), collapse = " "), paste(signif(want, 8), collapse = " ")))
  } else n_pass <<- n_pass + 1L
}
check_true <- function(label, cond) {
  if (!isTRUE(cond)) { n_fail <<- n_fail + 1L; cat("FAIL", label, "\n") } else n_pass <<- n_pass + 1L
}
# Error message of expr, or NA when it runs without error.
err_of <- function(expr) tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))

env <- new.env()
sys.source(RUN_META, envir = env)
has <- function(name) exists(name, envir = env, inherits = FALSE)

# Extraction rows in the shared format; unspecified columns are NA.
COLS <- env$STD_COLS
mk <- function(...) {
  rows <- list(...)
  d <- do.call(rbind, lapply(rows, function(r) {
    x <- as.list(setNames(rep(NA, length(COLS)), COLS)); x[names(r)] <- r
    as.data.frame(x, stringsAsFactors = FALSE)
  }))
  p <- tempfile(fileext = ".csv", tmpdir = WORK); write.csv(d, p, row.names = FALSE, na = ""); env$read_extraction(p)
}
write_cfg <- function(txt, ext = ".yaml") { p <- tempfile(fileext = ext, tmpdir = WORK); writeLines(txt, p); p }
cfg_of <- function(txt, outcome = "o", ext = ".yaml") env$read_config(write_cfg(txt, ext), outcome)[[1]]
base_cfg <- function(...) {
  cfg <- env$canon(modifyList(env$DEFAULTS, list(...)))
  cfg$outcome <- "o"; cfg
}

# ---------------------------------------------------------------------------
cat("CR4-1 Paule-Mandel above the default tau2.max\n")
yi <- c(0, .1, -.1, .05, -.05, 50, .02); vi <- rep(.01, 7)
cfg <- base_cfg(hksj = FALSE)
fp <- env$fit_rma(yi, vi, "PM", cfg)
check_true("PM fit succeeds with an outlier (metafor default cap 100 is below the root)", !is.null(fp))
# PM root by hand: generalised Q(tau2) = k - 1 (Paule & Mandel 1982)
qgen <- function(t) { w <- 1 / (vi + t); mu <- sum(w * yi) / sum(w); sum(w * (yi - mu)^2) - (length(yi) - 1) }
t_pm <- uniroot(qgen, c(0, 1e5), tol = 1e-12)$root
check("PM tau2 equals the root of Q_gen = k - 1", if (is.null(fp)) NA else fp$tau2, t_pm, tol = 1e-6)
fr <- env$fit_rma(yi, vi, "REML", cfg)
check("REML tau2 unchanged by the larger tau2.max", fr$tau2, metafor::rma(yi, vi, method = "REML", control = list(threshold = 1e-10, maxiter = 1000))$tau2, tol = 1e-6)
# A fit that cannot run writes its reason into pooled.csv
bad <- suppressMessages(env$fit_rma(c(NA, NA, NA), c(1, 1, 1), "REML", cfg))  # rma: k = 0
pr <- if (has("FIT_STATUS")) env$pooled_row(NULL, "REML", cfg, 3, env$FIT_STATUS$error) else env$pooled_row(NULL, "REML", cfg, 3)
check_true("failed fit leaves a 'fit failed:' note in the pooled row", is.null(bad) && grepl("^fit failed: ", pr$note))

# ---------------------------------------------------------------------------
cat("CR4-2 modified Egger precision for crossover SMD\n")
rows <- list()
for (j in 1:8) rows <- c(rows, list(
  list(study_id = paste0("P", j), outcome = "o", timepoint = "12", design = "parallel", arm = "NR", is_control = 0,
       n = 10 + 3 * j, mean = 52 + j / 3, sd = 8),
  list(study_id = paste0("P", j), outcome = "o", timepoint = "12", design = "parallel", arm = "pl", is_control = 1,
       n = 11 + 2 * j, mean = 50, sd = 8.2)))
for (j in 1:3) rows <- c(rows, list(
  list(study_id = paste0("X", j), outcome = "o", timepoint = "12", design = "crossover", arm = "NR", is_control = 0,
       n = 14 + 4 * j, mean = 51 + j / 2, sd = 6.5, corr_within = 0.4 + 0.15 * j),
  list(study_id = paste0("X", j), outcome = "o", timepoint = "12", design = "crossover", arm = "pl", is_control = 1,
       n = 14 + 4 * j, mean = 49.6, sd = 6.8)))
d <- do.call(mk, rows)
cfg <- base_cfg(measure = "SMD", hksj = FALSE)
be <- env$build_effects(d, cfg, cfg$corr_prepost, cfg$corr_within, character(0))
e <- be$eff
want_prec <- ifelse(e$design == "crossover", 2 * (1 - as.numeric(e$r_within)) / as.numeric(e$n_t),
                    1 / as.numeric(e$n_t) + 1 / as.numeric(e$n_c))
check("precision variance per comparison (parallel 1/n1+1/n2, crossover 2(1-r)/N)",
      if (is.null(e$prec_v)) rep(NA, nrow(e)) else e$prec_v, want_prec)
f <- env$fit_rma(e$yi, e$vi, "REML", cfg)
ss <- env$small_study_df(f, cfg, e)
rt <- regtest(e$yi, vi = want_prec, model = "lm", predictor = "sei")
check("modified Egger p uses the crossover precision term", ss$pval[ss$test == "egger"], rt$pval)

# ---------------------------------------------------------------------------
cat("CR4-3 check_reference.R keeps the WORKDIR it is given\n")
keep_dir <- file.path(WORK, "user_folder"); dir.create(keep_dir)
writeLines("user data", file.path(keep_dir, "precious.txt"))
st <- system2(RSCRIPT, c(file.path(TESTS, "check_reference.R"), keep_dir), stdout = TRUE, stderr = TRUE,
              env = "NR_REF_CHECK_SETUP_ONLY=1", timeout = 120)
check_true("a file in the WORKDIR survives check_reference.R", file.exists(file.path(keep_dir, "precious.txt")))

# ---------------------------------------------------------------------------
cat("CR4-4 Geary flag under the default rom_geary\n")
d <- mk(list(study_id = "G1", outcome = "o", timepoint = "12", design = "parallel", arm = "NR", is_control = 0, n = 10, mean = 2, sd = 3),
        list(study_id = "G1", outcome = "o", timepoint = "12", design = "parallel", arm = "pl", is_control = 1, n = 10, mean = 1.8, sd = 2.5),
        list(study_id = "G2", outcome = "o", timepoint = "12", design = "crossover", arm = "NR", is_control = 0, n = 12, mean = 2, sd = 3),
        list(study_id = "G2", outcome = "o", timepoint = "12", design = "crossover", arm = "pl", is_control = 1, n = 12, mean = 1.8, sd = 2.5))
g_par <- min(2 * sqrt(10) / 3, 1.8 * sqrt(10) / 2.5); g_x <- min(2 * sqrt(12) / 3, 1.8 * sqrt(12) / 2.5)
cfg <- base_cfg(measure = "ROM")
e <- env$build_effects(d, cfg, cfg$corr_prepost, cfg$corr_within, character(0))$eff
check_true("flag: both comparisons kept", identical(nrow(e), 2L))
check_true("flag: parallel note names the Geary statistic",
           grepl(sprintf("Geary's check failed (min mean sqrt(n)/SD = %.2f < 3)", g_par), e$method_notes[e$study_id == "G1"], fixed = TRUE))
check_true("flag: crossover note names the Geary statistic",
           grepl(sprintf("Geary's check failed (min mean sqrt(n)/SD = %.2f < 3)", g_x), e$method_notes[e$study_id == "G2"], fixed = TRUE))
check("flag: geary_min recorded", if (is.null(e$geary_min)) c(NA, NA) else e$geary_min, c(g_par, g_x))
cfg <- base_cfg(measure = "ROM", rom_geary = "exclude")
be <- env$build_effects(d, cfg, cfg$corr_prepost, cfg$corr_within, character(0))
check_true("exclude: both comparisons excluded with the Geary reason",
           is.null(be$eff) && nrow(be$excl) == 2L && all(grepl("Geary's check failed", be$excl$reason)))

# ---------------------------------------------------------------------------
cat("CR4-5 / Astra 9 forest.png prediction interval and level\n")
yi <- c(0.5, 1.2, 0.9, 2.1, 1.6, 0.3); vi <- c(.10, .12, .08, .20, .15, .09)
eff <- data.frame(study_id = paste0("S", 1:6), yi = yi, vi = vi)
drawn_pi <- function(args) {
  got <- NULL
  suppressMessages(trace("predict.rma", where = asNamespace("metafor"), print = FALSE,
                         exit = quote(assign("pp", returnValue(), envir = .GlobalEnv))))
  on.exit(suppressMessages(untrace("predict.rma", where = asNamespace("metafor"))))
  if (exists("pp", envir = .GlobalEnv)) rm("pp", envir = .GlobalEnv)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  do.call(forest, args)
  if (exists("pp", envir = .GlobalEnv)) got <- get("pp", envir = .GlobalEnv)
  if (is.null(got) || is.null(got$pi.lb)) c(NA_real_, NA_real_) else c(got$pi.lb, got$pi.ub)
}
fa <- function(f, cfg, pooled_ok = TRUE) {
  if (has("forest_args")) return(env$forest_args(f, eff, cfg, pooled_ok))
  # code before the fix: the arguments make_plots built inline
  list(x = f, slab = eff$study_id, header = c("Study", "Estimate [95% CI]"), xlab = "x",
       addpred = (cfg$estimator != "FE" && cfg$prediction_interval))
}
for (pt in c("riley", "revman", "metafor")) for (rule in c("always", "revman")) {
  cfg <- base_cfg(pi_type = pt, hksj_rule = rule, level = 90)
  f <- env$fit_rma(yi, vi, "REML", cfg)
  pr <- env$pooled_row(f, "REML", cfg, 6)
  check(sprintf("forest PI = pooled.csv PI (pi_type %s, hksj_rule %s, level 90)", pt, rule),
        drawn_pi(fa(f, cfg)), c(pr$pi_lb, pr$pi_ub), tol = 1e-10)
}
cfg <- base_cfg(level = 90)
f <- env$fit_rma(yi, vi, "REML", cfg)
check_true("forest header shows the configured level", identical(fa(f, cfg)$header[2], "Estimate [90% CI]"))
check_true("unpooled forest passes the configured level", identical(fa(f, cfg, pooled_ok = FALSE)$level, 90))
cfg <- base_cfg(min_k_pool = 1)
f2 <- env$fit_rma(yi[1:2], vi[1:2], "REML", cfg)
check_true("riley with k = 2: no PI line, as pooled.csv has none", isFALSE(fa(f2, cfg)$addpred))

# ---------------------------------------------------------------------------
cat("CR4-6 paired SE on the NR row with no paired MD\n")
one_x <- function(on_nr) {
  d <- mk(list(study_id = "C1", outcome = "o", timepoint = "12", design = "crossover", arm = "NR", is_control = 0,
               n = 16, mean = 50.8, sd = 6.5, paired_se = if (on_nr) 0.9 else NA),
          list(study_id = "C1", outcome = "o", timepoint = "12", design = "crossover", arm = "pl", is_control = 1,
               n = 16, mean = 49.6, sd = 6.8, paired_se = if (on_nr) NA else 0.9))
  cfg <- base_cfg(measure = "SMD")
  env$build_effects(d, cfg, cfg$corr_prepost, cfg$corr_within, character(0))$eff
}
r_want <- (6.5^2 + 6.8^2 - 0.9^2 * 16) / (2 * 6.5 * 6.8)
check("r derived from a paired SE on the NR row", one_x(TRUE)$r_within, r_want)
check_true("a paired SE on the control row is refused (CR1-3: it belongs on the NR row)", is.null(one_x(FALSE)))

# ---------------------------------------------------------------------------
cat("CR3-01 timepoint labels with hyphens\n")
check("tp_weeks", env$tp_weeks(c("12-month", "3-mo", "6-months", "Wk-12", "week-12", "day-21", "12-week",
                                 "wk6", "12 weeks", "day 21", "3 months", "8", "W12", "1.5 months")),
      c(365.25 / 7, 365.25 / 28, 365.25 / 14, 12, 12, 3, 12, 6, 12, 3, 365.25 / 28, 8, 12, 1.5 * 365.25 / 84))
d <- mk(list(study_id = "T1", outcome = "o", timepoint = "6 months", design = "parallel", arm = "NR", is_control = 0, n = 10, mean = 1, sd = 1),
        list(study_id = "T1", outcome = "o", timepoint = "12-month", design = "parallel", arm = "NR", is_control = 0, n = 10, mean = 2, sd = 1))
check_true("latest keeps '12-month' over '6 months'", identical(env$select_timepoint(d, "latest", list())$d$timepoint, "12-month"))
check_true("rule 12 does not match '12-month'", nrow(env$select_timepoint(d, "12", list())$d) == 0L)

# ---------------------------------------------------------------------------
cat("Astra 16 numeric timepoint rule uses timepoint_weeks\n")
d <- mk(list(study_id = "V1", outcome = "o", timepoint = "visit A", design = "parallel", arm = "NR", is_control = 0, n = 10, mean = 1, sd = 1),
        list(study_id = "V1", outcome = "o", timepoint = "visit A", design = "parallel", arm = "pl", is_control = 1, n = 10, mean = 1, sd = 1),
        list(study_id = "V1", outcome = "o", timepoint = "visit B", design = "parallel", arm = "NR", is_control = 0, n = 10, mean = 1, sd = 1),
        list(study_id = "V1", outcome = "o", timepoint = "visit B", design = "parallel", arm = "pl", is_control = 1, n = 10, mean = 1, sd = 1))
d$timepoint_weeks <- c(12, 12, 4, NA)
s12 <- env$select_timepoint(d, "12", list())
check_true("rule 12 keeps the 'visit A' rows (timepoint_weeks 12)", identical(s12$d$timepoint, c("visit A", "visit A")))
check_true("rule 4 keeps 'visit B' (weeks from its first filled row)", identical(env$select_timepoint(d, "4", list())$d$timepoint, c("visit B", "visit B")))
check_true("latest keeps 'visit A'", identical(env$select_timepoint(d, "latest", list())$d$timepoint, c("visit A", "visit A")))

# ---------------------------------------------------------------------------
cat("Astra 25 shared config: R stops on what it does not implement\n")
ok_yaml <- "defaults:\n  estimator: REML\noutcomes:\n  o:\n    measure: MD"
check_true("a valid config reads", is.na(err_of(cfg_of(ok_yaml))))
check_true("median_method defaults to luo2018 (metaproof default luo)", identical(cfg_of(ok_yaml)$median_method, "luo2018"))
refused <- list(
  "unknown key" = "defaults:\n  estimatr: DL\noutcomes:\n  o: {measure: MD}",
  "unknown outcome key" = "outcomes:\n  o: {measure: MD, subgroup_min: 3}",
  "unknown top-level key" = "default:\n  estimator: DL\noutcomes:\n  o: {measure: MD}",
  "bh: true (Python only)" = "defaults:\n  bh: true\noutcomes:\n  o: {measure: MD}",
  "combine_arms (Python only)" = "outcomes:\n  o: {measure: MD, combine_arms: {S1: {NR: [a, b]}}}",
  "metareg_permutations > 0 (Python only)" = "outcomes:\n  o: {measure: MD, metareg_permutations: 1000}",
  "quoted 'yes'" = "defaults:\n  hksj: 'yes'\noutcomes:\n  o: {measure: MD}",
  "hksj: 2" = "defaults:\n  hksj: 2\noutcomes:\n  o: {measure: MD}",
  "non-boolean prediction_interval" = "outcomes:\n  o: {measure: MD, prediction_interval: maybe}",
  "partial multi_arm 's'" = "outcomes:\n  o: {measure: MD, multi_arm: s}",
  "partial pi_type 'r'" = "outcomes:\n  o: {measure: MD, pi_type: r}",
  "measure RR (risk ratio is not implemented)" = "outcomes:\n  o: {measure: RR}",
  "unknown report_estimators" = "outcomes:\n  o: {measure: MD, report_estimators: [REML, XYZ]}",
  "level 1" = "outcomes:\n  o: {measure: MD, level: 1}",
  "corr_within 1.5" = "outcomes:\n  o: {measure: MD, corr_within: 1.5}",
  "corr_prepost_sensitivity outside [-1, 1]" = "outcomes:\n  o: {measure: MD, corr_prepost_sensitivity: [0.2, 2]}",
  "min_k_pool 0" = "outcomes:\n  o: {measure: MD, min_k_pool: 0}",
  "min_k_pool 2.5" = "outcomes:\n  o: {measure: MD, min_k_pool: 2.5}",
  "non-numeric corr_prepost" = "outcomes:\n  o: {measure: MD, corr_prepost: high}",
  "invalid defaults overridden by every outcome" = "defaults:\n  pi_type: wrong\noutcomes:\n  o: {measure: MD, pi_type: riley}")
for (nm in names(refused)) check_true(paste("refused:", nm), !is.na(err_of(cfg_of(refused[[nm]]))))
check_true("quoted 'false' reads as false (metaproof's vocabulary), not true",
           isFALSE(cfg_of("defaults:\n  hksj: 'false'\noutcomes:\n  o: {measure: MD}")$hksj))
check_true("corr_within_sensitivity defaults to auto", identical(cfg_of(ok_yaml)$corr_within_sensitivity, "auto"))
check_true("corr_within_sensitivity [] stays off", identical(cfg_of("outcomes:\n  o: {measure: MD, corr_within_sensitivity: []}")$corr_within_sensitivity, numeric(0)))
check_true("accepted: bh 'false' in quotes (metaproof reads it as false)",
           is.na(err_of(cfg_of("defaults:\n  bh: 'false'\noutcomes:\n  o: {measure: MD}"))))
check_true("accepted: primary false and bh_alpha while bh is false",
           is.na(err_of(cfg_of("defaults:\n  bh: false\n  bh_alpha: 0.1\noutcomes:\n  o: {measure: MD, primary: false}"))))
c1 <- cfg_of("outcomes:\n  o: {measure: MD, timepoint: [wk12]}")
check_true("a one-element timepoint list stays a list (YAML)", isTRUE(attr(c1$timepoint, "tp_list")))
c2 <- cfg_of('{"outcomes": {"o": {"measure": "MD", "timepoint": ["wk12"]}}}', ext = ".json")
check_true("a one-element timepoint list stays a list (JSON)", isTRUE(attr(c2$timepoint, "tp_list")))
c3 <- cfg_of("outcomes:\n  o: {measure: MD, timepoint: 12}")
check_true("a scalar timepoint is a single rule", !isTRUE(attr(c3$timepoint, "tp_list")))
d <- mk(list(study_id = "L1", outcome = "o", timepoint = "visit A", design = "parallel", arm = "NR", is_control = 0, n = 10, mean = 1, sd = 1))
d$timepoint_weeks <- 12
check_true("list rule [12] matches labels only (as metaproof), not weeks",
           nrow(env$select_timepoint(d, cfg_of("outcomes:\n  o: {measure: MD, timepoint: [12]}")$timepoint, list())$d) == 0L)
# The command line exits non-zero on a refused config.
csv <- tempfile(fileext = ".csv", tmpdir = WORK)
write.csv(data.frame(study_id = "A", outcome = "o", timepoint = "12", design = "parallel", arm = c("NR", "pl"),
                     is_control = c(0, 1), n = 10, mean = c(1, 0), sd = 1), csv, row.names = FALSE)
st <- suppressWarnings(system2(RSCRIPT, c(RUN_META, "--data", csv, "--config", write_cfg(refused[["bh: true (Python only)"]]),
                                          "--out", file.path(WORK, "cli_out")), stdout = TRUE, stderr = TRUE))
check_true("run_meta.R exits non-zero on a Python-only option", !is.null(attr(st, "status")) && attr(st, "status") != 0)

# corr_within_sensitivity auto: the crossover rows assumed r, so 0.2 and 0.8 are run
csv2 <- tempfile(fileext = ".csv", tmpdir = WORK)
xs <- do.call(rbind, lapply(1:3, function(j) data.frame(study_id = paste0("X", j), outcome = "o", timepoint = "12",
  design = "crossover", arm = c("NR", "pl"), is_control = c(0, 1), n = 12 + j, mean = c(51 + j, 50), sd = c(6, 6.4))))
write.csv(xs, csv2, row.names = FALSE)
st <- suppressWarnings(system2(RSCRIPT, c(RUN_META, "--data", csv2, "--config", write_cfg(ok_yaml), "--out", file.path(WORK, "auto_out")),
                               stdout = TRUE, stderr = TRUE))
sens <- tryCatch(read.csv(file.path(WORK, "auto_out", "o", "sensitivity.csv")), error = function(e) data.frame(analysis = character(0)))
check_true("auto corr_within sensitivity runs r = 0.2 and 0.8 for assumed crossover r",
           all(c("corr_within=0.2", "corr_within=0.8") %in% sens$analysis))


# ---------------------------------------------------------------------------
cat("Integration: crossover rules shared with metaproof.effects\n")
xrow <- function(sid, arm, ctl, n, ...) c(list(study_id = sid, outcome = "o", timepoint = "12", design = "crossover",
                                               arm = arm, is_control = ctl, n = n), list(...))
eff_of <- function(d, ...) { cfg <- base_cfg(...); env$build_effects(d, cfg, cfg$corr_prepost, cfg$corr_within, character(0)) }
d <- mk(xrow("M3", "NR low", 0, 18, mean = 44, sd = 10, corr_within = 0.5),
        xrow("M3", "NR high", 0, 18, mean = 49, sd = 11, corr_within = 0.5),
        xrow("M3", "pl", 1, 18, mean = 40, sd = 10.5, corr_within = 0.5))
for (ma in c("combine_treatments", "split_control")) {
  b <- eff_of(d, measure = "MD", multi_arm = ma)
  check_true(sprintf("multi-period crossover excluded, not split (multi_arm %s)", ma),
             is.null(b$eff) && any(grepl("^multi-period crossover", b$excl$reason)))
}
d <- mk(xrow("U1", "NR", 0, 14, mean = 48, sd = 11), xrow("U1", "pl", 1, 12, mean = 40, sd = 10))
b <- eff_of(d, measure = "MD")
check_true("unequal crossover n is excluded, not min(n)", is.null(b$eff) && grepl("complete pairs is not known", b$excl$reason[1]))
d <- mk(xrow("U2", "NR", 0, 14, mean = 48, sd = 11, paired_md = 7.5, paired_se = 2.1), xrow("U2", "pl", 1, 12, mean = 40, sd = 10))
check("unequal crossover n with a reported paired MD and SE: used as reported", unlist(eff_of(d, measure = "MD")$eff[c("yi", "vi")]), c(7.5, 2.1^2))
check_true("... but not for SMD (needs the complete-pair N)", is.null(eff_of(d, measure = "SMD")$eff))
d <- mk(xrow("S1", "NR", 0, 14, mean = 48, sd = 11), xrow("S1", "pl", 1, 14, mean = 40, sd = 10))
b <- eff_of(d, measure = "SMD", basis = "change")
check_true("crossover SMD refused under data_basis: change (Handbook 10.5.2)", is.null(b$eff) && grepl("between-person SD", b$excl$reason[1]))
b0 <- eff_of(d, measure = "ROM", rom_bias_correction = FALSE)$eff; b1 <- eff_of(d, measure = "ROM", rom_bias_correction = TRUE)$eff
check("crossover ROM bias correction adds 0.5 (SDe^2/(N Me^2) - SDc^2/(N Mc^2)), variance unchanged",
      c(b1$yi - b0$yi, b1$vi - b0$vi), c(0.5 * (11^2 / (14 * 48^2) - 10^2 / (14 * 40^2)), 0))

cat("Integration: range-only gate, combining groups, input checks\n")
prow <- function(sid, arm, ctl, n, ...) c(list(study_id = sid, outcome = "o", timepoint = "12", design = "parallel",
                                               arm = arm, is_control = ctl, n = n), list(...))
d <- mk(prow("R1", "NR", 0, 20, mean_change = -4, sd_change = 8, median = 70, min = 40, max = 120),
        prow("R1", "pl", 1, 21, mean_change = -1, sd_change = 8.5, median = 32, min = 20, max = 55))
b <- eff_of(d, measure = "MD", basis = "change", allow_range_only = FALSE)
check("reported change stats with an unused median/range are used (allow_range_only false)",
      if (is.null(b$eff)) c(NA, NA) else unlist(b$eff[c("yi", "vi")]), c(-3, 64 / 20 + 8.5^2 / 21))
b <- eff_of(d, measure = "MD", basis = "final", allow_range_only = FALSE)
check_true("a range-only final value that feeds the effect is still excluded", is.null(b$eff) && grepl("only a range", b$excl$reason[1]))
g0 <- env$combine_groups(c(20, 20), c(0, 1), c(1, 1)); g8 <- env$combine_groups(c(20, 20), c(1e8, 1e8 + 1), c(1, 1))
check("combined SD unchanged by shifting both means by 1e8", c(g8$sd, g8$mean - 1e8), c(g0$sd, g0$mean), tol = 1e-12)
csv_of <- function(df) { p <- tempfile(fileext = ".csv", tmpdir = WORK); write.csv(df, p, row.names = FALSE, na = ""); p }
good <- data.frame(study_id = "A", outcome = "o", timepoint = "12", design = "parallel", arm = c("NR", "pl"),
                   is_control = c("0", "1"), n = c("10", "11"), mean = c("1", "0"), sd = c("1", "1"), corr_prepost = NA)
check_true("a clean extraction reads", is.na(err_of(env$read_extraction(csv_of(good)))))
bad_cells <- list("blank n" = list(n = c("", "11")), "n '20 (ITT)'" = list(n = c("20 (ITT)", "11")),
                  "mean '12,5'" = list(mean = c("12,5", "0")), "is_control 'yes'" = list(is_control = c("yes", "1")),
                  "negative SD" = list(sd = c("-1", "1")), "corr_prepost 1.1" = list(corr_prepost = c("1.1", NA)),
                  "n 10.5" = list(n = c("10.5", "11")), "design 'cross-over'" = list(design = c("cross-over", "cross-over")))
for (nm in names(bad_cells)) {
  df <- good; for (c in names(bad_cells[[nm]])) df[[c]] <- bad_cells[[nm]][[c]]
  e <- suppressMessages(err_of(env$read_extraction(csv_of(df))))
  check_true(sprintf("read_extraction stops on %s, naming the row", nm), !is.na(e) && grepl("row 2: ", e))
}

cat("Integration: meta-regression gate, R2 and level order\n")
mr <- function(n_st, mods, ...) {
  set.seed(7)
  e <- data.frame(study_id = paste0("S", 1:n_st), source_study = paste0("S", 1:n_st), yi = rnorm(n_st, 0.3, 0.3),
                  vi = runif(n_st, 0.02, 0.08), dose = rep(c(250, 500, 1000), length.out = n_st),
                  pop = rep(c("adults", "Older adults"), length.out = n_st), stringsAsFactors = FALSE)
  env$metareg_df(e, base_cfg(moderators = mods, metareg_joint = TRUE, ...))
}
m <- mr(12, c("dose", "pop"))
check_true("joint model of 2 moderators with 12 studies is not run (needs 10 x 2)",
           grepl("metareg_min_k x 2 moderators = 10 x 2 = 20", m$note[m$model == "dose+pop"][1]))
check_true("each single-moderator model with 12 studies runs", all(c("dose", "pop") %in% m$model[!is.na(m$term)]))
check_true("categorical levels in C-locale order: 'Older adults' is the reference, as in Python",
           identical(m$term[m$model == "pop"], c("intrcpt", "popadults")))
e0 <- data.frame(study_id = paste0("S", 1:10), source_study = paste0("S", 1:10), yi = rep(c(0.2, 0.21), 5),
                 vi = rep(0.5, 10), dose = c(1:10) * 100, stringsAsFactors = FALSE)
m0 <- env$metareg_df(e0, base_cfg(moderators = "dose"))
check("R2 = 0 when the intercept-only tau2 is 0", m0$R2[1], 0)
check_true("metareg_min_k 0 needs metareg_allow_small_k: true",
           !is.na(err_of(cfg_of("outcomes:\n  o: {measure: MD, metareg_min_k: 0}"))) &&
           is.na(err_of(cfg_of("outcomes:\n  o: {measure: MD, metareg_min_k: 0, metareg_allow_small_k: true}"))))

cat("Integration: config and output folders\n")
for (nm in names(rf <- list(
  "outcome under defaults" = "defaults:\n  outcome: o\noutcomes:\n  o: {measure: MD}",
  "bh_alpha under an outcome" = "outcomes:\n  o: {measure: MD, bh_alpha: 0.1}",
  "hksj_adhoc with an explicit different hksj_rule" = "outcomes:\n  o: {measure: MD, hksj_adhoc: true, hksj_rule: always}")))
  check_true(paste("refused:", nm), !is.na(err_of(cfg_of(rf[[nm]]))))
check_true("hksj_adhoc alone gives truncated", identical(cfg_of("outcomes:\n  o: {measure: MD, hksj_adhoc: true}")$hksj_rule, "truncated"))
cli <- function(csv, yaml, out) suppressWarnings(system2(RSCRIPT, c(RUN_META, "--data", csv, "--config", write_cfg(yaml), "--out", out),
                                                       stdout = TRUE, stderr = TRUE))
status <- function(st) if (is.null(attr(st, "status"))) 0L else attr(st, "status")
three <- do.call(rbind, lapply(1:6, function(j) data.frame(study_id = paste0("P", j), outcome = "o", timepoint = "12",
  design = "parallel", arm = c("NR", "pl"), is_control = c(0, 1), n = 20 + j, mean = c(5 + j / 3, 4), sd = c(2, 2.1), sex = "f")))
csv3 <- csv_of(three)
check_true("an outcome typo stops the run", status(cli(csv3, "outcomes:\n  oo: {measure: MD}", file.path(WORK, "typo"))) != 0L)
check_true("a moderator column not in the file stops the run",
           status(cli(csv3, "outcomes:\n  o: {measure: MD, moderators: [age]}", file.path(WORK, "nocol"))) != 0L)
check_true("two analyses sharing a folder ignoring case are refused",
           status(cli(csv3, "outcomes:\n  SBP: {outcome: o, measure: MD}\n  sbp: {outcome: o, measure: MD}", file.path(WORK, "case"))) != 0L)
check_true("an analysis named '..' is refused",
           status(cli(csv3, "outcomes:\n  '..': {outcome: o, measure: MD}", file.path(WORK, "dots"))) != 0L)
out <- file.path(WORK, "rerun")
st <- cli(csv3, "outcomes:\n  o: {measure: MD}\n  old: {outcome: o, measure: MD}", out)
check_true("first run writes forest.png", file.exists(file.path(out, "o", "forest.png")))
po <- read.csv(file.path(out, "o", "pooled.csv"))
yy <- three$mean[three$is_control == 0] - 4; vv <- 2^2 / (20 + 1:6) + 2.1^2 / (20 + 1:6)
f0 <- rma(yy, vv, method = "REML", control = env$rma_control(yy))  # hksj_rule revman: knha only when tau2 > 0
f <- rma(yy, vv, method = "REML", test = if (f0$tau2 > env$TAU2_ZERO * mean(vv)) "knha" else "z", control = env$rma_control(yy))
check("pooled.csv stat and pval are the configured test (REML, HKSJ)",
      if (all(c("stat", "pval") %in% names(po))) unlist(po[1, c("stat", "pval")]) else c(NA, NA), c(f$zval, f$pval), tol = 1e-6)
writeLines("user", file.path(out, "mine.txt")); dir.create(file.path(out, "user_dir"))
st <- cli(csv3, "outcomes:\n  o: {measure: MD, timepoint: wk99}", out)
check_true("rerun with no usable study leaves no forest.png from the earlier run", !file.exists(file.path(out, "o", "forest.png")))
check_true("rerun removes the folder of an analysis no longer configured", !dir.exists(file.path(out, "old")))
check_true("rerun leaves files and folders run_meta.R did not write", file.exists(file.path(out, "mine.txt")) && dir.exists(file.path(out, "user_dir")))
check_true("no staging folder is left behind", !any(startsWith(list.files(out, all.files = TRUE), ".run_meta-staging-")))

cat("Integration: forest PI below 5 studies; timepoint_weeks conflicts\n")
cfg <- base_cfg(pi_type = "revman")
eff4 <- data.frame(study_id = paste0("S", 1:4), source_study = paste0("S", 1:4), yi = yi[1:4], vi = vi[1:4])
f4 <- env$fit_rma(eff4$yi, eff4$vi, "REML", cfg)
check_true("forest draws no PI with 4 studies (pooled.csv keeps it)",
           isFALSE(env$forest_args(f4, eff4, cfg, TRUE)$addpred) && !is.na(env$pooled_row(f4, "REML", cfg, 4)$pi_lb))
check_true("... and the figure says why", grepl("prediction interval not shown: 4 studies", env$forest_pi_note(f4, eff4, cfg, TRUE)))
d <- mk(list(study_id = "W1", outcome = "o", timepoint = "visit A", design = "parallel", arm = "NR", is_control = 0, n = 10, mean = 1, sd = 1),
        list(study_id = "W1", outcome = "o", timepoint = "visit A", design = "parallel", arm = "pl", is_control = 1, n = 10, mean = 1, sd = 1),
        list(study_id = "W1", outcome = "o", timepoint = "visit B", design = "parallel", arm = "NR", is_control = 0, n = 10, mean = 1, sd = 1),
        list(study_id = "W1", outcome = "o", timepoint = "visit B", design = "parallel", arm = "pl", is_control = 1, n = 10, mean = 1, sd = 1))
d$timepoint_weeks <- c(12, 8, 4, 4)
sw <- env$select_timepoint(d, "latest", list())
check_true("rows of one label with two timepoint_weeks exclude the study", nrow(sw$d) == 0L &&
           grepl("fix timepoint_weeks", sw$excl[[1]]$reason))

cat(sprintf("\n%d checks passed, %d failed\n", n_pass, n_fail))
quit(status = if (n_fail) 1L else 0L)
