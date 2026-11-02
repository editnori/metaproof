#!/usr/bin/env Rscript
# make_reference_data.R -- export metafor reference datasets (real published
# data, NOT nicotinamide riboside data and NOT synthetic) to
# validation/reference_data/, in the shared extraction format where possible,
# together with expected values computed by calling metafor directly.
# This script does not source run_meta.R, so the expected values are an
# independent check of both engines.
#
# Datasets (all shipped with the metadat package that metafor loads):
#   dat.bcg         Colditz et al. (1994) JAMA 271(9):698-702. BCG vaccine and
#                   tuberculosis, 13 trials, 2x2 counts. Exported as log risk
#                   ratios (escalc measure "RR") in a yi-vi CSV, and in the
#                   extraction format by passing yi/sei through paired_md /
#                   paired_se. Moderator: absolute latitude (ablat), as in
#                   Berkey et al. (1995) Stat Med 14(4):395-411.
#   dat.normand1999 Normand (1999) Stat Med 18(3):321-359. Length of hospital
#                   stay (days), specialist vs routine stroke care, 9 studies,
#                   arm-level means/SDs. MD and SMD (Hedges' g).
#   dat.curtis1998  Curtis & Wang (1998) Oecologia 113(3):299-313, as analysed
#                   by Hedges, Gurevitch & Curtis (1999) Ecology 80(4):1150-1156.
#                   Woody plant biomass under elevated vs ambient CO2, 102
#                   effects, arm-level means/SDs. Log ratio of means (ROM) with
#                   meta-regression on duration and subgroups by functional group.
#
# Usage: Rscript reference_r/make_reference_data.R [OUTDIR]
#   default OUTDIR = <kit>/validation/reference_data

get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
SCRIPT_DIR <- get_script_dir()
if (dir.exists(file.path(SCRIPT_DIR, "lib"))) .libPaths(c(file.path(SCRIPT_DIR, "lib"), .libPaths()))
suppressPackageStartupMessages({ library(metafor); library(jsonlite) })

args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args)) args[1] else file.path(dirname(SCRIPT_DIR), "validation", "reference_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
w <- function(df, name) write.csv(df, file.path(OUT, name), row.names = FALSE, na = "NA")
wj <- function(x, name) writeLines(toJSON(x, auto_unbox = TRUE, pretty = TRUE), file.path(OUT, name))

EXTRACTION_COLS <- c("study_id", "outcome", "timepoint", "design", "arm", "is_control", "n", "mean", "sd", "se",
                     "median", "q1", "q3", "min", "max", "mean_change", "sd_change", "baseline_mean", "baseline_sd",
                     "corr_prepost", "corr_within", "paired_md", "paired_se", "dose_mg_day", "duration_weeks",
                     "population", "risk_of_bias", "notes")
blank_extraction <- function(n) as.data.frame(setNames(replicate(length(EXTRACTION_COLS), rep(NA, n), simplify = FALSE), EXTRACTION_COLS))

# Tight convergence so expected values sit on the exact optimum (metafor defaults:
# threshold 1e-5 for tau^2, uniroot tol ~1.2e-4 for the Q-profile CI). With the
# defaults, dat.bcg ~ ablat gives tau^2 = 0.0763547 instead of the REML optimum
# 0.0763480. leave1out() and influence() reuse the fitted model's control.
RMA_CONTROL <- list(threshold = 1e-10, maxiter = 1000, tol = 1e-12)
# tau2.max: the Paule-Mandel root is searched only up to metafor's default max(100,
# 10 mad(yi)^2), which one outlying study can exceed. Q_gen(t) <= sum((y - mean(y))^2) / t,
# so the PM root is at most SS/(k - p) and this bound always contains it.
rma <- function(yi, ...) metafor::rma(yi, ..., control = c(RMA_CONTROL, list(tau2.max = 100 + 10 * sum((yi - mean(yi))^2))))
# tau2.max: confint() otherwise stops at max(100, 10 tau^2) and reports that cap as the
# upper limit (dat.normand1999 DL: 2054.09 instead of 2889.45). Q_gen(t) <=
# sum((y - mean(y))^2) / t keeps this bound above the root. tol in RMA_CONTROL is also the
# uniroot tolerance of the Paule-Mandel estimator.
confint_q <- function(f) confint(f, control = list(tol = 1e-10,
  tau2.max = 100 + 10 * sum((f$yi - mean(f$yi))^2) / qchisq(1e-4, f$k - f$p) + 10 * max(f$vi)))

ESTIMATORS <- c("REML", "DL", "PM", "ML", "FE")
fits_grid <- function() {
  g <- expand.grid(estimator = ESTIMATORS, test = c("z", "knha"), stringsAsFactors = FALSE)
  g[!(g$estimator == "FE" & g$test == "knha"), ]
}

pooled_expected <- function(yi, vi, dataset, outcome, measure) {
  g <- fits_grid()
  do.call(rbind, lapply(seq_len(nrow(g)), function(j) {
    est <- g$estimator[j]; tst <- g$test[j]
    f <- rma(yi, vi, method = est, test = tst)
    re <- est != "FE"
    pr <- if (re) predict(f, predtype = "Riley") else NULL
    pd <- if (re) predict(f) else NULL
    tci <- if (re) confint_q(f)$random["tau^2", c("ci.lb", "ci.ub")] else c(NA, NA)
    data.frame(dataset = dataset, outcome = outcome, measure = measure, estimator = est, test = tst, k = f$k,
               estimate = as.numeric(f$b), se = f$se, ci_lb = f$ci.lb, ci_ub = f$ci.ub,
               pi_lb = if (re) pr$pi.lb else NA, pi_ub = if (re) pr$pi.ub else NA,
               tau2 = if (re) f$tau2 else 0, tau2_ci_lb = unname(tci[1]), tau2_ci_ub = unname(tci[2]),
               I2 = f$I2, H2 = f$H2, Q = f$QE, Q_p = f$QEp,
               pi_lb_metafor = if (re) pd$pi.lb else NA, pi_ub_metafor = if (re) pd$pi.ub else NA,
               pval = f$pval, pi_note = "pi_lb/pi_ub: predtype='Riley' (t, k-2 df); pi_*_metafor: predict() default (t k-1 df under knha, z otherwise)")
  }))
}

metareg_expected <- function(yi, vi, X, dataset, model_label) {
  g <- fits_grid(); g <- g[g$estimator %in% c("REML", "DL", "PM", "ML"), ]
  do.call(rbind, lapply(seq_len(nrow(g)), function(j) {
    f <- rma(yi, vi, mods = X, method = g$estimator[j], test = g$test[j])
    data.frame(dataset = dataset, model = model_label, estimator = g$estimator[j], test = g$test[j],
               term = c("intrcpt", colnames(X)), k = f$k, estimate = as.numeric(f$b), se = f$se, statistic = f$zval,
               df = if (g$test[j] == "knha") f$ddf[1] else NA, pval = f$pval, ci_lb = f$ci.lb, ci_ub = f$ci.ub,
               tau2 = f$tau2, I2 = f$I2, R2 = f$R2, QM = f$QM, QM_df1 = f$QMdf[1], QM_df2 = f$QMdf[2],
               QM_p = f$QMp, QE = f$QE, QE_p = f$QEp)
  }))
}

small_study_expected <- function(yi, vi, dataset) {
  rows <- list()
  f0 <- rma(yi, vi, method = "REML", test = "knha")
  rl <- regtest(f0, model = "lm", ret.fit = TRUE); cf <- coef(summary(rl$fit))
  rows[[1]] <- data.frame(dataset = dataset, test = "egger", model = "lm", estimator = NA, rma_test = NA, k = f0$k,
                          coef_sei = cf[2, 1], coef_se = cf[2, 2], statistic = rl$zval, df = rl$dfs, pval = rl$pval,
                          limit_estimate = rl$est, limit_ci_lb = rl$ci.lb, limit_ci_ub = rl$ci.ub)
  g <- fits_grid(); g <- g[g$estimator %in% c("REML", "DL", "PM"), ]
  for (j in seq_len(nrow(g))) {
    f <- rma(yi, vi, method = g$estimator[j], test = g$test[j])
    r <- regtest(f, model = "rma", ret.fit = TRUE)
    rows[[length(rows) + 1L]] <- data.frame(dataset = dataset, test = "egger_rma", model = "rma", estimator = g$estimator[j],
      rma_test = g$test[j], k = f$k, coef_sei = as.numeric(r$fit$b[2]), coef_se = r$fit$se[2], statistic = r$zval,
      df = if (is.null(r$dfs)) NA else r$dfs, pval = r$pval, limit_estimate = r$est, limit_ci_lb = r$ci.lb, limit_ci_ub = r$ci.ub)
  }
  do.call(rbind, rows)
}

loo_influence_expected <- function(yi, vi, slab, dataset) {
  f <- rma(yi, vi, method = "REML", test = "knha", slab = slab)
  l <- leave1out(f)
  loo <- data.frame(dataset = dataset, estimator = "REML", test = "knha", study_id = slab, k = f$k - 1L,
                    estimate = l$estimate, se = l$se, ci_lb = l$ci.lb, ci_ub = l$ci.ub, pval = l$pval,
                    tau2 = l$tau2, I2 = l$I2, H2 = l$H2, Q = l$Q, Q_p = l$Qp)
  inf <- influence(f)
  infl <- data.frame(dataset = dataset, estimator = "REML", test = "knha", study_id = slab,
                     rstudent = inf$inf$rstudent, dffits = inf$inf$dffits, cooks_d = inf$inf$cook.d,
                     cov_r = inf$inf$cov.r, tau2_del = inf$inf$tau2.del, QE_del = inf$inf$QE.del, hat = inf$inf$hat,
                     weight = inf$inf$weight, dfbetas = inf$dfbs[[1]], influential = as.logical(inf$is.infl))
  list(loo = loo, inf = infl)
}

subgroup_expected <- function(yi, vi, grp, dataset, moderator) {
  rows <- list(); ests <- c(); ses <- c()
  for (l in sort(unique(grp))) {
    i <- grp == l
    f <- rma(yi[i], vi[i], method = "REML", test = "knha")
    p <- if (f$k >= 3) predict(f, predtype = "Riley") else list(pi.lb = NA, pi.ub = NA)
    # Q_between uses the model-based SE of each subgroup estimate (Borenstein et al. 2009,
    # ch. 19): the chi-square test is not defined on HKSJ-scaled variances.
    fz <- rma(yi[i], vi[i], method = "REML", test = "z")
    ests <- c(ests, f$b[1]); ses <- c(ses, fz$se)
    rows[[length(rows) + 1L]] <- data.frame(dataset = dataset, moderator = moderator, level = l, k = f$k,
      estimate = as.numeric(f$b), se = f$se, ci_lb = f$ci.lb, ci_ub = f$ci.ub, pi_lb = p$pi.lb, pi_ub = p$pi.ub,
      tau2 = f$tau2, I2 = f$I2, Q = f$QE, Q_p = f$QEp, test = NA, test_stat = NA, test_df1 = NA, test_df2 = NA, test_p = NA)
  }
  tb <- rma(ests, sei = ses, method = "FE")
  rows[[length(rows) + 1L]] <- data.frame(dataset = dataset, moderator = moderator, level = "ALL", k = length(yi),
    estimate = NA, se = NA, ci_lb = NA, ci_ub = NA, pi_lb = NA, pi_ub = NA, tau2 = NA, I2 = NA, Q = NA, Q_p = NA,
    test = "Q_between_separate_tau2", test_stat = tb$QE, test_df1 = length(ests) - 1L, test_df2 = NA, test_p = tb$QEp)
  fm <- rma(yi, vi, mods = model.matrix(~ factor(grp))[, -1, drop = FALSE], method = "REML", test = "knha")
  rows[[length(rows) + 1L]] <- data.frame(dataset = dataset, moderator = moderator, level = "ALL", k = length(yi),
    estimate = NA, se = NA, ci_lb = NA, ci_ub = NA, pi_lb = NA, pi_ub = NA, tau2 = fm$tau2, I2 = fm$I2, Q = fm$QE, Q_p = fm$QEp,
    test = "QM_common_tau2", test_stat = fm$QM, test_df1 = fm$QMdf[1], test_df2 = fm$QMdf[2], test_p = fm$QMp)
  do.call(rbind, rows)
}

PROVENANCE <- "REAL PUBLISHED DATA from the metafor/metadat package (not NR data, not synthetic); used only to validate the engines"

# ---------------------------------------------------------------------------
# 1. dat.bcg -- log risk ratio
# ---------------------------------------------------------------------------
bcg <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)
bcg_id <- paste0("bcg_", bcg$trial, "_", gsub("[^A-Za-z]+", "", bcg$author), "_", bcg$year)
w(data.frame(study_id = bcg_id, author = bcg$author, year = bcg$year, yi = as.numeric(bcg$yi), vi = as.numeric(bcg$vi),
             sei = sqrt(as.numeric(bcg$vi)), tpos = bcg$tpos, tneg = bcg$tneg, cpos = bcg$cpos, cneg = bcg$cneg,
             ablat = bcg$ablat, alloc = bcg$alloc,
             measure = "log risk ratio (BCG vs control), escalc measure='RR'",
             source = paste("metafor dat.bcg; Colditz et al. 1994 JAMA;", PROVENANCE)),
  "bcg_logrr_yi_vi.csv")
ex <- blank_extraction(nrow(bcg))
ex$study_id <- bcg_id; ex$outcome <- "tb_log_rr"; ex$timepoint <- NA; ex$design <- "parallel"
ex$arm <- "BCG vaccine vs no vaccine (precomputed log RR)"; ex$is_control <- 0
ex$n <- bcg$tpos + bcg$tneg + bcg$cpos + bcg$cneg
ex$paired_md <- as.numeric(bcg$yi); ex$paired_se <- sqrt(as.numeric(bcg$vi))
ex$population <- bcg$alloc
ex$notes <- paste("dat.bcg: log RR from 2x2 counts passed through paired_md/paired_se; run with measure MD (pass-through); n = total participants;",
                  "population column holds treatment allocation method;", PROVENANCE)
ex$ablat <- bcg$ablat
w(ex, "bcg_logrr_extraction.csv")
wj(list(defaults = list(measure = "MD", estimator = "REML", hksj = TRUE, prediction_interval = TRUE, pi_type = "riley", hksj_rule = "always", subgroup_method = "separate",
                        subgroup_min_k = 0, subgroup_min_per_level = 0, metareg_min_k = 0, metareg_allow_small_k = TRUE, min_k_pool = 1),
        outcomes = list(tb_log_rr = list(moderators = list("ablat"), subgroups = list("population")))),
   "bcg_config.json")

bcg_pooled <- pooled_expected(bcg$yi, bcg$vi, "dat.bcg", "tb_log_rr", "log RR")
w(bcg_pooled, "bcg_expected_pooled.csv")
w(metareg_expected(bcg$yi, bcg$vi, cbind(ablat = bcg$ablat), "dat.bcg", "ablat"), "bcg_expected_metareg.csv")
w(small_study_expected(bcg$yi, bcg$vi, "dat.bcg"), "bcg_expected_small_study.csv")
li <- loo_influence_expected(bcg$yi, bcg$vi, bcg_id, "dat.bcg")
w(li$loo, "bcg_expected_leave1out.csv"); w(li$inf, "bcg_expected_influence.csv")
w(subgroup_expected(as.numeric(bcg$yi), as.numeric(bcg$vi), bcg$alloc, "dat.bcg", "population"), "bcg_expected_subgroups.csv")

# ---------------------------------------------------------------------------
# 2. dat.normand1999 -- MD and SMD
# ---------------------------------------------------------------------------
nd <- dat.normand1999
nid <- paste0("normand_", nd$study, "_", gsub("[^A-Za-z]+", "", nd$source))
mk_norm <- function(outcome) {
  e <- blank_extraction(2 * nrow(nd))
  e$study_id <- rep(nid, each = 2); e$outcome <- outcome; e$timepoint <- "discharge"; e$design <- "parallel"
  e$arm <- rep(c("specialist care", "routine care"), nrow(nd)); e$is_control <- rep(c(0, 1), nrow(nd))
  e$n <- c(rbind(nd$n1i, nd$n2i)); e$mean <- c(rbind(nd$m1i, nd$m2i)); e$sd <- c(rbind(nd$sd1i, nd$sd2i))
  e$notes <- paste("dat.normand1999: length of hospital stay (days), Normand 1999 Stat Med;", PROVENANCE)
  e
}
w(rbind(mk_norm("los_days_md"), mk_norm("los_days_smd")), "normand1999_extraction.csv")
wj(list(defaults = list(estimator = "REML", hksj = TRUE, prediction_interval = TRUE, pi_type = "riley", hksj_rule = "always", subgroup_method = "separate",
                        subgroup_min_k = 0, subgroup_min_per_level = 0, metareg_min_k = 0, metareg_allow_small_k = TRUE, min_k_pool = 1),
        outcomes = list(los_days_md = list(measure = "MD"), los_days_smd = list(measure = "SMD"))),
   "normand1999_config.json")
md <- escalc("MD", m1i = m1i, sd1i = sd1i, n1i = n1i, m2i = m2i, sd2i = sd2i, n2i = n2i, data = nd)
smd <- escalc("SMD", m1i = m1i, sd1i = sd1i, n1i = n1i, m2i = m2i, sd2i = sd2i, n2i = n2i, data = nd)
w(rbind(data.frame(dataset = "dat.normand1999", outcome = "los_days_md", measure = "MD", study_id = nid,
                   yi = as.numeric(md$yi), vi = as.numeric(md$vi), sei = sqrt(as.numeric(md$vi))),
        data.frame(dataset = "dat.normand1999", outcome = "los_days_smd", measure = "SMD (Hedges g, exact J, vtype LS)", study_id = nid,
                   yi = as.numeric(smd$yi), vi = as.numeric(smd$vi), sei = sqrt(as.numeric(smd$vi)))),
  "normand1999_expected_study_effects.csv")
w(rbind(pooled_expected(md$yi, md$vi, "dat.normand1999", "los_days_md", "MD"),
        pooled_expected(smd$yi, smd$vi, "dat.normand1999", "los_days_smd", "SMD")),
  "normand1999_expected_pooled.csv")

# ---------------------------------------------------------------------------
# 3. dat.curtis1998 -- log ratio of means, meta-regression and subgroups
# ---------------------------------------------------------------------------
cu <- dat.curtis1998
cid <- paste0("curtis1998_id", cu$id)
e <- blank_extraction(2 * nrow(cu))
e$study_id <- rep(cid, each = 2); e$outcome <- "biomass_rom"; e$timepoint <- "end"; e$design <- "parallel"
e$arm <- rep(c("elevated CO2", "ambient CO2"), nrow(cu)); e$is_control <- rep(c(0, 1), nrow(cu))
e$n <- c(rbind(cu$n1i, cu$n2i)); e$mean <- c(rbind(cu$m1i, cu$m2i)); e$sd <- c(rbind(cu$sd1i, cu$sd2i))
e$duration_weeks <- rep(cu$time / 7, each = 2); e$population <- rep(cu$fungrp, each = 2)
e$notes <- paste0("dat.curtis1998 paper ", rep(cu$paper, each = 2),
                  ": total above- plus below-ground woody plant biomass (g); duration_weeks = maximum CO2 exposure days / 7; population = plant functional group; ",
                  "effects from the same paper are treated as independent, as in the metafor examples; ", PROVENANCE)
e$method <- rep(cu$method, each = 2); e$time_days <- rep(cu$time, each = 2)
w(e, "curtis1998_extraction.csv")
wj(list(defaults = list(measure = "ROM", estimator = "REML", hksj = TRUE, prediction_interval = TRUE, pi_type = "riley", hksj_rule = "always", subgroup_method = "separate",
                        subgroup_min_k = 0, subgroup_min_per_level = 0, metareg_min_k = 0, metareg_allow_small_k = TRUE, min_k_pool = 1,
                        rom_bias_correction = FALSE),
        outcomes = list(biomass_rom = list(moderators = list("duration_weeks"), subgroups = list("population")))),
   "curtis1998_config.json")
rom <- escalc("ROM", m1i = m1i, sd1i = sd1i, n1i = n1i, m2i = m2i, sd2i = sd2i, n2i = n2i, data = cu, correct = FALSE)
w(data.frame(dataset = "dat.curtis1998", outcome = "biomass_rom", measure = "ROM (log ratio of means, no bias correction)",
             study_id = cid, yi = as.numeric(rom$yi), vi = as.numeric(rom$vi), sei = sqrt(as.numeric(rom$vi))),
  "curtis1998_expected_study_effects.csv")
w(pooled_expected(rom$yi, rom$vi, "dat.curtis1998", "biomass_rom", "ROM"), "curtis1998_expected_pooled.csv")
w(metareg_expected(rom$yi, rom$vi, cbind(duration_weeks = cu$time / 7), "dat.curtis1998", "duration_weeks"),
  "curtis1998_expected_metareg.csv")
w(small_study_expected(rom$yi, rom$vi, "dat.curtis1998"), "curtis1998_expected_small_study.csv")
w(subgroup_expected(as.numeric(rom$yi), as.numeric(rom$vi), cu$fungrp, "dat.curtis1998", "population"),
  "curtis1998_expected_subgroups.csv")
li <- loo_influence_expected(rom$yi, rom$vi, cid, "dat.curtis1998")
w(li$loo, "curtis1998_expected_leave1out.csv"); w(li$inf, "curtis1998_expected_influence.csv")

writeLines(c(paste("Generated", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), "by reference_r/make_reference_data.R"),
             paste("R", getRversion(), "| metafor", packageVersion("metafor"), "| metadat", packageVersion("metadat"))),
           file.path(OUT, "GENERATED_BY.txt"))
message("wrote reference data to ", OUT)
