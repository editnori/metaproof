#!/usr/bin/env Rscript
# make_SYNTHETIC_example.R -- writes SYNTHETIC_example_extraction.csv.
#
# SYNTHETIC DATA. Every number below was invented to exercise the code paths
# of the meta-analysis engines (SE->SD, median/IQR/range conversion, change
# scores with reported and imputed SDs, multi-arm splitting, crossover with
# and without a reported paired MD, timepoint selection, exclusion). These are
# NOT nicotinamide riboside trial results and must never be presented as such.

args <- commandArgs(trailingOnly = TRUE)
a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
here <- if (length(f)) dirname(normalizePath(f[1])) else getwd()
out <- if (length(args)) args[1] else file.path(here, "SYNTHETIC_example_extraction.csv")

COLS <- c("study_id", "outcome", "timepoint", "design", "arm", "is_control", "n", "mean", "sd", "se",
          "median", "q1", "q3", "min", "max", "mean_change", "sd_change", "baseline_mean", "baseline_sd",
          "corr_prepost", "corr_within", "paired_md", "paired_se", "dose_mg_day", "duration_weeks",
          "population", "risk_of_bias", "notes")
TAG <- "SYNTHETIC - invented numbers for software testing only; not real trial data"
rows <- list()
r <- function(study_id, arm, is_control, ..., design = "parallel", timepoint = "12", dose = NA, weeks = 12,
              pop = "SYNTHETIC-older adults", rob = "low", what = "") {
  x <- setNames(as.list(rep(NA, length(COLS))), COLS)
  x$study_id <- study_id; x$timepoint <- timepoint; x$design <- design; x$arm <- arm; x$is_control <- is_control
  x$dose_mg_day <- if (is_control == 1) 0 else dose; x$duration_weeks <- weeks
  x$population <- pop; x$risk_of_bias <- rob
  x$notes <- paste0(TAG, if (nzchar(what)) paste0("; path tested: ", what) else "")
  v <- list(...); for (k in names(v)) x[[k]] <- v[[k]]
  rows[[length(rows) + 1L]] <<- x
}

# S01: final mean/SD at two timepoints (latest = 12 is used)
r("SYNTH-S01", "NR", 0, n = 40, mean = 52.0, sd = 8.0, timepoint = "6", dose = 1000, weeks = 6, what = "earlier timepoint, should be ignored")
r("SYNTH-S01", "placebo", 1, n = 40, mean = 50.5, sd = 8.2, timepoint = "6", weeks = 6, what = "earlier timepoint, should be ignored")
r("SYNTH-S01", "NR", 0, n = 38, mean = 53.1, sd = 8.4, dose = 1000, what = "final mean/SD")
r("SYNTH-S01", "placebo", 1, n = 39, mean = 50.2, sd = 8.1, what = "final mean/SD")
# S02: SE instead of SD
r("SYNTH-S02", "NR", 0, n = 25, mean = 49.0, se = 1.6, dose = 500, weeks = 8, rob = "some concerns", what = "SD = SE*sqrt(n)")
r("SYNTH-S02", "placebo", 1, n = 24, mean = 47.9, se = 1.7, weeks = 8, rob = "some concerns", what = "SD = SE*sqrt(n)")
# S03: median and IQR (Wan C3 / Luo C3)
r("SYNTH-S03", "NR", 0, n = 30, median = 51.0, q1 = 46.0, q3 = 57.0, dose = 1000, pop = "SYNTHETIC-younger adults", what = "median/IQR conversion")
r("SYNTH-S03", "placebo", 1, n = 31, median = 49.0, q1 = 44.5, q3 = 55.0, pop = "SYNTHETIC-younger adults", what = "median/IQR conversion")
# S04: median and range (C1)
r("SYNTH-S04", "NR", 0, n = 18, median = 50.0, min = 35.0, max = 68.0, dose = 2000, weeks = 4, rob = "high", what = "median/range conversion")
r("SYNTH-S04", "placebo", 1, n = 18, median = 48.0, min = 33.0, max = 66.0, weeks = 4, rob = "high", what = "median/range conversion")
# S05: full five-number summary (C2)
r("SYNTH-S05", "NR", 0, n = 22, median = 52.0, q1 = 47.0, q3 = 58.0, min = 38.0, max = 70.0, dose = 500, pop = "SYNTHETIC-younger adults", what = "five-number conversion")
r("SYNTH-S05", "placebo", 1, n = 21, median = 50.0, q1 = 45.0, q3 = 56.0, min = 36.0, max = 67.0, pop = "SYNTHETIC-younger adults", what = "five-number conversion")
# S06: change from baseline with reported SD of change
r("SYNTH-S06", "NR", 0, n = 45, mean_change = 2.4, sd_change = 5.0, baseline_mean = 50.0, baseline_sd = 8.0, dose = 1000, weeks = 24, what = "reported change SD")
r("SYNTH-S06", "placebo", 1, n = 44, mean_change = 0.3, sd_change = 5.2, baseline_mean = 50.4, baseline_sd = 7.9, weeks = 24, what = "reported change SD")
# S07: change mean, baseline SD and final SD only -> SD of change imputed with default r
r("SYNTH-S07", "NR", 0, n = 28, mean_change = 1.9, baseline_mean = 49.5, baseline_sd = 7.5, sd = 7.8, dose = 250, weeks = 12, what = "imputed change SD, default r")
r("SYNTH-S07", "placebo", 1, n = 27, mean_change = 0.6, baseline_mean = 49.9, baseline_sd = 7.7, sd = 7.6, weeks = 12, what = "imputed change SD, default r")
# S08: three arms sharing one control -> control n split in two
r("SYNTH-S08", "NR 500 mg", 0, n = 30, mean = 51.5, sd = 7.9, dose = 500, what = "multi-arm, shared control")
r("SYNTH-S08", "NR 1000 mg", 0, n = 31, mean = 52.8, sd = 8.1, dose = 1000, what = "multi-arm, shared control")
r("SYNTH-S08", "placebo", 1, n = 32, mean = 50.0, sd = 8.0, what = "multi-arm, shared control")
# S09: crossover with a reported paired MD and SE (plus per-condition summaries)
r("SYNTH-S09", "NR", 0, n = 20, mean = 51.2, sd = 7.0, paired_md = 1.4, paired_se = 0.9, design = "crossover", dose = 1000, weeks = 6, what = "crossover, reported paired MD")
r("SYNTH-S09", "placebo", 1, n = 20, mean = 49.9, sd = 7.2, design = "crossover", weeks = 6, what = "crossover, reported paired MD")
# S10: crossover, per-condition means/SDs only, within-person r defaulted
r("SYNTH-S10", "NR", 0, n = 16, mean = 50.8, sd = 6.5, design = "crossover", dose = 2000, weeks = 4, what = "crossover, default within-person r")
r("SYNTH-S10", "placebo", 1, n = 16, mean = 49.6, sd = 6.8, design = "crossover", weeks = 4, what = "crossover, default within-person r")
# S11: crossover with study-reported within-person r
r("SYNTH-S11", "NR", 0, n = 24, mean = 52.4, sd = 7.4, corr_within = 0.7, design = "crossover", dose = 500, weeks = 6, pop = "SYNTHETIC-younger adults", what = "crossover, reported within-person r")
r("SYNTH-S11", "placebo", 1, n = 24, mean = 51.0, sd = 7.1, corr_within = 0.7, design = "crossover", weeks = 6, pop = "SYNTHETIC-younger adults", what = "crossover, reported within-person r")
# S12: change with study-reported pre-post r (not varied in sensitivity)
r("SYNTH-S12", "NR", 0, n = 35, mean_change = 2.0, baseline_mean = 50.0, baseline_sd = 8.1, sd = 8.3, corr_prepost = 0.65, dose = 1000, weeks = 16, rob = "some concerns", what = "imputed change SD, study r")
r("SYNTH-S12", "placebo", 1, n = 35, mean_change = 0.1, baseline_mean = 50.2, baseline_sd = 8.0, sd = 8.0, corr_prepost = 0.65, weeks = 16, rob = "some concerns", what = "imputed change SD, study r")
# S13: no dispersion at all -> excluded
r("SYNTH-S13", "NR", 0, n = 15, mean = 53.0, dose = 1000, rob = "high", what = "no SD/SE: must be excluded")
r("SYNTH-S13", "placebo", 1, n = 15, mean = 50.0, rob = "high", what = "no SD/SE: must be excluded")

base <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
all <- do.call(rbind, lapply(c("SYNTHETIC_outcome_MD", "SYNTHETIC_outcome_SMD", "SYNTHETIC_outcome_ROM"),
                             function(o) { b <- base; b$outcome <- o; b[COLS] }))
write.csv(all, out, row.names = FALSE, na = "")
message("wrote ", out, " (SYNTHETIC)")
