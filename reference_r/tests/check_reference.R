#!/usr/bin/env Rscript
# check_reference.R -- self-checks for reference_r/run_meta.R.
#   1. run_meta.R on the metafor reference datasets reproduces the expected
#      values that make_reference_data.R computed by calling metafor directly.
#   2. The plain-R median conversions agree with metafor::conv.fivenum where the
#      two use the same method (Luo mean; Wan SD for range-only and IQR-only).
#   3. Hand-computed SYNTHETIC cases (Handbook formulas) match run_meta.R.
#   4. Every SYNTHETIC example row is labelled SYNTHETIC.
#   5. The fixes in tests/check_fixes.R (config strictness, timepoints, Geary, forest
#      PI, PM tau2.max, crossover Egger precision, WORKDIR safety).
# Usage: Rscript reference_r/tests/check_reference.R [WORKDIR]
#   Outputs go to WORKDIR/nr_meta_ref_check (default: under tempdir()). Only that
#   subdirectory is ever deleted, never WORKDIR itself.
# Exit status 0 = all checks passed.

a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
TESTS <- dirname(normalizePath(f[1])); RREF <- dirname(TESTS); KIT <- dirname(RREF)
.libPaths(c(file.path(RREF, "lib"), .libPaths()))
suppressPackageStartupMessages(library(metafor))
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(trailingOnly = TRUE)
# The script owns only this subdirectory: a WORKDIR the user passes is never deleted.
WORK <- file.path(if (length(args)) args[1] else tempdir(), "nr_meta_ref_check")
unlink(WORK, recursive = TRUE); dir.create(WORK, recursive = TRUE)
# check_fixes.R uses this to test the WORKDIR handling without running every check.
if (nzchar(Sys.getenv("NR_REF_CHECK_SETUP_ONLY"))) { cat(WORK, "\n"); quit(status = 0L) }
REF <- file.path(KIT, "validation", "reference_data")
RSCRIPT <- file.path(R.home("bin"), "Rscript")
TOL <- 1e-8

n_fail <- 0L; n_pass <- 0L
check <- function(label, got, want, tol = TOL) {
  got <- as.numeric(got); want <- as.numeric(want)
  same_na <- identical(is.na(got), is.na(want))
  d <- if (same_na && any(!is.na(want))) max(abs(got - want)[!is.na(want)] / pmax(1, abs(want[!is.na(want)]))) else 0
  if (length(got) != length(want) || !same_na || d > tol) {
    n_fail <<- n_fail + 1L
    cat(sprintf("FAIL %s (max rel diff %.3g; got %s; want %s)\n", label, d,
                paste(signif(head(got, 4), 8), collapse = " "), paste(signif(head(want, 4), 8), collapse = " ")))
  } else n_pass <<- n_pass + 1L
}
run <- function(data, config, out) {
  st <- system2(RSCRIPT, c(file.path(RREF, "run_meta.R"), "--data", data, "--config", config, "--out", out), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(st, "status")) && attr(st, "status") != 0) { cat(st, sep = "\n"); stop("run_meta.R failed") }
}
rd <- function(...) read.csv(file.path(...), stringsAsFactors = FALSE)
POOLED <- c("k", "estimate", "se", "ci_lb", "ci_ub", "pi_lb", "pi_ub", "tau2", "tau2_ci_lb", "tau2_ci_ub", "I2", "H2", "Q", "Q_p")

# pi = "riley" compares pi_lb/pi_ub with predtype Riley, "metafor" with predict()'s default.
# With yi/vi given, the truncated-HKSJ interval is checked against metafor test = "adhoc".
compare_pooled <- function(name, out_dir, exp_file, outcome, pi = "riley", yi = NULL, vi = NULL) {
  got <- rd(out_dir, "pooled.csv"); exp <- rd(REF, exp_file); exp <- exp[exp$outcome == outcome, ]
  for (est in got$estimator) {
    e <- exp[exp$estimator == est & exp$test == if (est == "FE") "z" else "knha", ]
    if (pi == "metafor") { e$pi_lb <- e$pi_lb_metafor; e$pi_ub <- e$pi_ub_metafor }
    g <- got[got$estimator == est, ]
    for (c in POOLED) check(sprintf("%s pooled %s %s", name, est, c), g[[c]], e[[c]])
    # the three intervals hksj_rule chooses between (pooled.csv extra columns)
    ez <- exp[exp$estimator == est & exp$test == "z", ]
    check(sprintf("%s pooled %s wald_ci", name, est), unlist(g[c("wald_ci_lb", "wald_ci_ub")]), c(ez$ci_lb, ez$ci_ub))
    if (est != "FE") {
      check(sprintf("%s pooled %s hksj_ci", name, est), unlist(g[c("hksj_ci_lb", "hksj_ci_ub")]), c(e$ci_lb, e$ci_ub))
      if (!is.null(yi)) {
        ft <- rma(yi, vi, method = est, test = "adhoc", control = list(threshold = 1e-10, maxiter = 1000, tol = 1e-12,
                                                                     tau2.max = 100 + 10 * sum((yi - mean(yi))^2)))
        check(sprintf("%s pooled %s hksj_trunc_ci", name, est), unlist(g[c("hksj_trunc_ci_lb", "hksj_trunc_ci_ub")]), c(ft$ci.lb, ft$ci.ub))
      }
    }
  }
}
# A config variant of a reference config, written to WORK.
variant <- function(cfg_file, name, defaults = list(), outcome = NULL, outcome_over = list()) {
  cfg <- jsonlite::fromJSON(file.path(REF, cfg_file), simplifyVector = FALSE)
  cfg$defaults <- modifyList(cfg$defaults, defaults)
  if (!is.null(outcome)) cfg$outcomes[[outcome]] <- modifyList(cfg$outcomes[[outcome]], outcome_over)
  p <- file.path(WORK, paste0(name, ".json")); writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE), p); p
}
compare_metareg <- function(name, out_dir, exp_file) {
  got <- rd(out_dir, "metareg.csv"); exp <- rd(REF, exp_file); exp <- exp[exp$estimator == "REML" & exp$test == "knha", ]
  for (c in c("estimate", "se", "statistic", "df", "pval", "ci_lb", "ci_ub", "tau2", "I2", "R2", "QM", "QM_df1", "QM_df2", "QM_p", "QE", "QE_p"))
    check(sprintf("%s metareg %s", name, c), got[[c]][match(exp$term, got$term)], exp[[c]])
}
compare_small <- function(name, out_dir, exp_file) {
  got <- rd(out_dir, "small_study.csv"); exp <- rd(REF, exp_file)
  e1 <- exp[exp$test == "egger", ]; e2 <- exp[exp$test == "egger_rma" & exp$estimator == "REML" & exp$rma_test == "knha", ]
  for (c in c("coef_sei", "coef_se", "statistic", "df", "pval", "limit_estimate", "limit_ci_lb", "limit_ci_ub")) {
    check(sprintf("%s egger %s", name, c), got[got$test == "egger", c], e1[[c]])
    check(sprintf("%s egger_rma %s", name, c), got[got$test == "egger_rma", c], e2[[c]])
  }
}
compare_loo_inf <- function(name, out_dir, prefix) {
  g <- rd(out_dir, "leave1out.csv"); e <- rd(REF, paste0(prefix, "_expected_leave1out.csv"))
  for (c in c("estimate", "se", "ci_lb", "ci_ub", "pval", "tau2", "I2", "H2", "Q", "Q_p")) check(sprintf("%s leave1out %s", name, c), g[[c]][match(e$study_id, g$study_id)], e[[c]])
  g <- rd(out_dir, "influence.csv"); e <- rd(REF, paste0(prefix, "_expected_influence.csv"))
  for (c in c("rstudent", "dffits", "cooks_d", "cov_r", "tau2_del", "QE_del", "hat", "weight", "dfbetas")) check(sprintf("%s influence %s", name, c), g[[c]][match(e$study_id, g$study_id)], e[[c]])
}
compare_subgroups <- function(name, out_dir, exp_file) {
  g <- rd(out_dir, "subgroups.csv"); e <- rd(REF, exp_file)
  key <- function(d) paste(d$level, ifelse(is.na(d$test), "", d$test))
  i <- match(key(e), key(g))
  for (c in c("k", "estimate", "se", "ci_lb", "ci_ub", "pi_lb", "pi_ub", "tau2", "I2", "Q", "Q_p", "test_stat", "test_df1", "test_df2", "test_p"))
    check(sprintf("%s subgroups %s", name, c), g[[c]][i], e[[c]])
}

# 1. Reference datasets ------------------------------------------------------
cat("1. metafor reference datasets\n")
run(file.path(REF, "bcg_logrr_extraction.csv"), file.path(REF, "bcg_config.json"), file.path(WORK, "bcg"))
o <- file.path(WORK, "bcg", "tb_log_rr")
yv <- rd(REF, "bcg_logrr_yi_vi.csv"); se <- rd(o, "study_effects.csv")
check("bcg yi", se$yi, yv$yi); check("bcg vi", se$vi, yv$vi)
compare_pooled("bcg", o, "bcg_expected_pooled.csv", "tb_log_rr", yi = yv$yi, vi = yv$vi)
compare_metareg("bcg", o, "bcg_expected_metareg.csv")
compare_small("bcg", o, "bcg_expected_small_study.csv")
compare_loo_inf("bcg", o, "bcg")
compare_subgroups("bcg", o, "bcg_expected_subgroups.csv")
# pi_type metafor and hksj_rule revman on real data: tau2 > 0 and k = 13, so revman uses
# HKSJ and every pooled column equals the test = "knha" expected row.
run(file.path(REF, "bcg_logrr_extraction.csv"), variant("bcg_config.json", "bcg_pimetafor_revman",
    list(pi_type = "metafor", hksj_rule = "revman")), file.path(WORK, "bcg_v"))
ov <- file.path(WORK, "bcg_v", "tb_log_rr")
compare_pooled("bcg pi_type metafor, hksj_rule revman", ov, "bcg_expected_pooled.csv", "tb_log_rr", pi = "metafor")
check("bcg hksj_rule revman chooses HKSJ", as.numeric(rd(ov, "pooled.csv")$ci_method[1] == "HKSJ t(12)"), 1)

run(file.path(REF, "normand1999_extraction.csv"), file.path(REF, "normand1999_config.json"), file.path(WORK, "normand"))
ese <- rd(REF, "normand1999_expected_study_effects.csv")
for (oc in c("los_days_md", "los_days_smd")) {
  o <- file.path(WORK, "normand", oc); se <- rd(o, "study_effects.csv"); e <- ese[ese$outcome == oc, ]
  check(paste("normand", oc, "yi"), se$yi[match(e$study_id, se$study_id)], e$yi)
  check(paste("normand", oc, "vi"), se$vi[match(e$study_id, se$study_id)], e$vi)
  compare_pooled(paste("normand", oc), o, "normand1999_expected_pooled.csv", oc,
                 yi = e$yi[match(se$study_id, e$study_id)], vi = e$vi[match(se$study_id, e$study_id)])
}
# smd_variance LS2 and UB against escalc directly
for (vt in c("LS2", "UB")) {
  run(file.path(REF, "normand1999_extraction.csv"), variant("normand1999_config.json", paste0("normand_", vt),
      outcome = "los_days_smd", outcome_over = list(smd_variance = vt)), file.path(WORK, paste0("normand_", vt)))
  se <- rd(WORK, paste0("normand_", vt), "los_days_smd", "study_effects.csv")
  es <- escalc("SMD", m1i = m1i, sd1i = sd1i, n1i = n1i, m2i = m2i, sd2i = sd2i, n2i = n2i, data = dat.normand1999, vtype = vt)
  nid <- paste0("normand_", dat.normand1999$study, "_", gsub("[^A-Za-z]+", "", dat.normand1999$source))
  check(paste("normand SMD vtype", vt, "yi"), se$yi[match(nid, se$study_id)], es$yi)
  check(paste("normand SMD vtype", vt, "vi"), se$vi[match(nid, se$study_id)], es$vi)
}

run(file.path(REF, "curtis1998_extraction.csv"), file.path(REF, "curtis1998_config.json"), file.path(WORK, "curtis"))
o <- file.path(WORK, "curtis", "biomass_rom"); se <- rd(o, "study_effects.csv"); e <- rd(REF, "curtis1998_expected_study_effects.csv")
check("curtis yi", se$yi[match(e$study_id, se$study_id)], e$yi); check("curtis vi", se$vi[match(e$study_id, se$study_id)], e$vi)
compare_pooled("curtis", o, "curtis1998_expected_pooled.csv", "biomass_rom", yi = e$yi[match(se$study_id, e$study_id)],
               vi = e$vi[match(se$study_id, e$study_id)])
compare_metareg("curtis", o, "curtis1998_expected_metareg.csv")
compare_small("curtis", o, "curtis1998_expected_small_study.csv")
compare_loo_inf("curtis", o, "curtis1998")
compare_subgroups("curtis", o, "curtis1998_expected_subgroups.csv")
for (p in c("forest.png", "funnel.png")) if (!file.exists(file.path(o, p))) { n_fail <- n_fail + 1L; cat("FAIL missing", p, "\n") }
# rom_bias_correction: escalc(correct = TRUE), the Lajeunesse 2015 bias-corrected log ROM
run(file.path(REF, "curtis1998_extraction.csv"), variant("curtis1998_config.json", "curtis_rombc", list(rom_bias_correction = TRUE)),
    file.path(WORK, "curtis_bc"))
se <- rd(WORK, "curtis_bc", "biomass_rom", "study_effects.csv")
es <- escalc("ROM", m1i = m1i, sd1i = sd1i, n1i = n1i, m2i = m2i, sd2i = sd2i, n2i = n2i, data = dat.curtis1998, correct = TRUE)
cid <- paste0("curtis1998_id", dat.curtis1998$id)
check("curtis ROM bias-corrected yi", se$yi[match(cid, se$study_id)], es$yi)
check("curtis ROM bias-corrected vi", se$vi[match(cid, se$study_id)], es$vi)

# 2. Median conversion vs conv.fivenum ---------------------------------------
cat("2. median conversion vs metafor::conv.fivenum\n")
env <- new.env(); sys.source(file.path(RREF, "run_meta.R"), envir = env)
ns <- c(10, 25, 60, 200)
for (n in ns) {
  cf <- conv.fivenum(min = 30, median = 50, max = 75, n = n, test = FALSE)
  check(sprintf("C1 Luo mean n=%d", n), env$est_mean_fivenum(30, NA, 50, NA, 75, n, "luo2018")$v, cf$mean)
  check(sprintf("C1 Wan SD n=%d", n), env$est_sd_wan(30, NA, NA, 75, n)$v, cf$sd)
  cf <- conv.fivenum(q1 = 44, median = 50, q3 = 58, n = n, test = FALSE)
  check(sprintf("C3 Luo mean n=%d", n), env$est_mean_fivenum(NA, 44, 50, 58, NA, n, "luo2018")$v, cf$mean)
  check(sprintf("C3 Wan SD n=%d", n), env$est_sd_wan(NA, 44, 58, NA, n)$v, cf$sd)
  cf <- conv.fivenum(min = 30, q1 = 44, median = 50, q3 = 58, max = 75, n = n, test = FALSE)
  check(sprintf("C2 Luo mean n=%d", n), env$est_mean_fivenum(30, 44, 50, 58, 75, n, "luo2018")$v, cf$mean)
}

# 3. SYNTHETIC hand checks ----------------------------------------------------
cat("3. SYNTHETIC hand checks (Handbook formulas)\n")
ex <- file.path(RREF, "examples")
run(file.path(ex, "SYNTHETIC_example_extraction.csv"), file.path(ex, "SYNTHETIC_example_config.yaml"), file.path(WORK, "synth"))
s <- rd(WORK, "synth", "SYNTHETIC_outcome_MD", "study_effects.csv"); g <- function(id, c) s[[c]][s$study_id == id]
check("S01 latest timepoint MD", g("SYNTH-S01", "yi"), 53.1 - 50.2)
check("S02 SD from SE", g("SYNTH-S02", "vi"), 1.6^2 + 1.7^2)
xi <- function(n) 2 * qnorm((n - 0.375) / (n + 0.25)); eta <- function(n) 2 * qnorm((0.75 * n - 0.125) / (n + 0.25))
check("S03 Wan C3", c(g("SYNTH-S03", "yi"), g("SYNTH-S03", "vi")),
      c((46 + 51 + 57) / 3 - (44.5 + 49 + 55) / 3, ((57 - 46) / eta(30))^2 / 30 + ((55 - 44.5) / eta(31))^2 / 31))
check("S04 Wan C1", c(g("SYNTH-S04", "yi"), g("SYNTH-S04", "vi")),
      c((35 + 100 + 68) / 4 - (33 + 96 + 66) / 4, ((68 - 35) / xi(18))^2 / 18 + ((66 - 33) / xi(18))^2 / 18))
sd5 <- function(a, q1, q3, b, n) (b - a) / (2 * xi(n)) + (q3 - q1) / (2 * eta(n))
check("S05 Wan C2", c(g("SYNTH-S05", "yi"), g("SYNTH-S05", "vi")),
      c((38 + 94 + 104 + 116 + 70) / 8 - (36 + 90 + 100 + 112 + 67) / 8,
        sd5(38, 47, 58, 70, 22)^2 / 22 + sd5(36, 45, 56, 67, 21)^2 / 21))
check("S06 reported change", c(g("SYNTH-S06", "yi"), g("SYNTH-S06", "vi")), c(2.1, 25 / 45 + 5.2^2 / 44))
imp <- function(b, f, r) b^2 + f^2 - 2 * r * b * f
check("S07 imputed change SD r=0.5", g("SYNTH-S07", "vi"), imp(7.5, 7.8, 0.5) / 28 + imp(7.7, 7.6, 0.5) / 27)
check("S12 imputed change SD study r", g("SYNTH-S12", "vi"), imp(8.1, 8.3, 0.65) / 35 + imp(8.0, 8.0, 0.65) / 35)
check("S08 split control", c(g("SYNTH-S08|NR 500 mg", "vi"), g("SYNTH-S08|NR 1000 mg", "vi")),
      c(7.9^2 / 30 + 8^2 / 16, 8.1^2 / 31 + 8^2 / 16))
check("S09 paired MD", c(g("SYNTH-S09", "yi"), g("SYNTH-S09", "vi")), c(1.4, 0.81))
check("S10 crossover MD r=0.5", g("SYNTH-S10", "vi"), (6.5^2 + 6.8^2 - 2 * 0.5 * 6.5 * 6.8) / 16)
check("S11 crossover MD r=0.7", g("SYNTH-S11", "vi"), (7.4^2 + 7.1^2 - 2 * 0.7 * 7.4 * 7.1) / 24)
exc <- rd(WORK, "synth", "SYNTHETIC_outcome_MD", "excluded.csv")
if (!identical(exc$study_id, "SYNTH-S13")) { n_fail <- n_fail + 1L; cat("FAIL S13 should be the only exclusion\n") }
sens <- rd(WORK, "synth", "SYNTHETIC_outcome_MD", "sensitivity.csv")
check("sensitivity r=0.2 differs from primary", abs(sens$estimate[2] - sens$estimate[1]) > 0, 1)
s2 <- rd(WORK, "synth", "SYNTHETIC_outcome_SMD", "study_effects.csv")
# Crossover SMD carries Hedges' J on N - 1 df (crossover_smd_j default true), N = 16
sdp <- sqrt((6.5^2 + 6.8^2) / 2); J15 <- exp(lgamma(15 / 2) - 0.5 * log(15 / 2) - lgamma(14 / 2))
y10 <- J15 * 1.2 / sdp
check("S10 crossover SMD", c(s2$yi[s2$study_id == "SYNTH-S10"], s2$vi[s2$study_id == "SYNTH-S10"]), c(y10, (1 / 16 + y10^2 / 32) * 2 * 0.5))
s3 <- rd(WORK, "synth", "SYNTHETIC_outcome_ROM", "study_effects.csv")
check("S10 crossover ROM", s3$vi[s3$study_id == "SYNTH-S10"],
      6.5^2 / (16 * 50.8^2) + 6.8^2 / (16 * 49.6^2) - 2 * 0.5 * 6.5 * 6.8 / (16 * 50.8 * 49.6))
# Pooled MD from the study effects with metafor directly (REML + Riley PI), with the
# hksj_rule revman choice made here independently: knha only when tau2 > 0 and k >= 3.
fz <- rma(s$yi, s$vi, method = "REML", test = "z")
use_knha <- fz$k >= 3 && fz$tau2 > 1e-8 * mean(s$vi)
fm <- if (use_knha) rma(s$yi, s$vi, method = "REML", test = "knha") else fz
pm <- predict(fm, predtype = "Riley")
check("SYNTH pooled ci_method", as.numeric(identical(rd(WORK, "synth", "SYNTHETIC_outcome_MD", "pooled.csv")$ci_method[1],
      if (use_knha) sprintf("HKSJ t(%d)", fz$k - 1) else "Wald z")), 1)
pp <- rd(WORK, "synth", "SYNTHETIC_outcome_MD", "pooled.csv")
check("SYNTH pooled REML", unlist(pp[1, c("estimate", "se", "ci_lb", "ci_ub", "pi_lb", "pi_ub", "tau2")]),
      c(fm$b, fm$se, fm$ci.lb, fm$ci.ub, pm$pi.lb, pm$pi.ub, fm$tau2))

# 4. SYNTHETIC labelling -------------------------------------------------------
cat("4. SYNTHETIC labelling\n")
d <- read.csv(file.path(ex, "SYNTHETIC_example_extraction.csv"), stringsAsFactors = FALSE)
lab_ok <- all(grepl("SYNTHETIC", d$notes)) && all(grepl("^SYNTH-", d$study_id)) && all(grepl("^SYNTHETIC_", d$outcome))
if (!lab_ok) { n_fail <- n_fail + 1L; cat("FAIL some SYNTHETIC rows are not labelled\n") } else n_pass <- n_pass + 1L

# 5. Review fixes ---------------------------------------------------------------
cat("5. review fixes (tests/check_fixes.R)\n")
fx <- system2(RSCRIPT, c(file.path(TESTS, "check_fixes.R"), file.path(WORK, "fixes")), stdout = TRUE, stderr = TRUE)
fx_status <- attr(fx, "status") %||% 0L
cat(grep("^(FAIL|[0-9]+ checks)", fx, value = TRUE), sep = "\n")
if (fx_status != 0L) { n_fail <- n_fail + 1L; cat("FAIL check_fixes.R exited with status", fx_status, "\n") } else n_pass <- n_pass + 1L

cat(sprintf("\n%d checks passed, %d failed (outputs in %s)\n", n_pass, n_fail, WORK))
quit(status = if (n_fail) 1L else 0L)
