validation/reference_data: metafor reference datasets and expected values
==========================================================================

These are REAL PUBLISHED datasets that ship with the R packages metafor and
metadat. They are not nicotinamide riboside (NR) data and are not synthetic.
They are here only to check that the kit's engines reproduce metafor.

Built by reference_r/make_reference_data.R, which calls metafor directly and
does not use reference_r/run_meta.R. GENERATED_BY.txt records the R, metafor
and metadat versions. Rebuild with:
  Rscript reference_r/make_reference_data.R

Expected-value conventions
  estimator  REML, DL, PM, ML, FE
  test       z (no adjustment) or knha (Hartung-Knapp-Sidik-Jonkman). FE uses z only.
  pi_lb/pi_ub          metafor predict(predtype = "Riley"): t with k-2 df
  pi_*_metafor         metafor predict() default: t with k-1 df under knha, z otherwise
  tau2_ci_*            Q-profile CI from confint(), for every random-effects estimator
  I2                   percent (0-100)
  leave1out, influence, subgroups: REML with knha
  small_study: "egger" = regtest(model = "lm") (classic Egger 1997);
               "egger_rma" = regtest(model = "rma") for each estimator and test
  convergence  rma(control = list(threshold = 1e-10, maxiter = 1000, tol = 1e-12)) and
               confint(control = list(tol = 1e-10, tau2.max = <above the upper limit>)),
               so values sit on the exact optimum and Q-profile roots. metafor's defaults
               (threshold 1e-5, uniroot tol ~1.2e-4, tau2.max = max(100, 10 tau2)) give
               tau2 values up to ~3e-5 off and can report the tau2.max cap as the upper
               limit (dat.normand1999, DL). For the BCG latitude model (REML, z) the default
               threshold gives tau2 0.0763547 and QM 16.3571; the optimum is tau2 0.0763480.
  subgroups    Q_between_separate_tau2 weights each subgroup estimate by its model-based
               (non-HKSJ) variance (Borenstein et al. 2009, ch. 19).


1. dat.bcg: BCG vaccine and tuberculosis, 13 trials
   Colditz GA, Brewer TF, Berkey CS, et al. (1994). JAMA 271(9):698-702.
   Moderator: absolute latitude (ablat), as in Berkey CS, Hoaglin DC,
   Mosteller F, Colditz GA (1995). Stat Med 14(4):395-411.
   bcg_logrr_yi_vi.csv           log risk ratio (escalc "RR") as yi, vi, sei, plus counts,
                                 ablat and alloc
   bcg_logrr_extraction.csv      the same yi/sei in the extraction format via
                                 paired_md / paired_se (one row per trial,
                                 is_control = 0). Run it with measure MD, which passes
                                 yi through. population = allocation method. ablat is an
                                 extra column.
   bcg_config.json               measure MD, REML, HKSJ, moderators [ablat],
                                 subgroups [population]
   bcg_expected_pooled.csv, _metareg.csv (~ ablat), _small_study.csv,
   _leave1out.csv, _influence.csv, _subgroups.csv (by allocation)

2. dat.normand1999: length of hospital stay (days), specialist vs routine
   stroke care, 9 studies, arm-level means and SDs.
   Normand S-LT (1999). Stat Med 18(3):321-359.
   normand1999_extraction.csv    two rows per study; the same rows appear under
                                 outcome los_days_md and outcome los_days_smd
   normand1999_config.json       los_days_md: MD; los_days_smd: SMD (Hedges' g)
   normand1999_expected_study_effects.csv, normand1999_expected_pooled.csv
   Note: several SMDs exceed 2 in absolute value; metafor prints a warning
   about this. That is a property of the data.

3. dat.curtis1998: total woody plant biomass under elevated vs ambient CO2,
   102 effects, arm-level means and SDs.
   Curtis PS, Wang X (1998). Oecologia 113(3):299-313; analysed as response
   ratios by Hedges LV, Gurevitch J, Curtis PS (1999). Ecology 80(4):1150-1156.
   curtis1998_extraction.csv     two rows per effect; duration_weeks = maximum CO2
                                 exposure in days / 7; population = plant
                                 functional group (ANGIO, GYMNO, N2FIX); extra
                                 columns method and time_days. Effects from one
                                 paper are treated as independent, as in the
                                 metafor examples.
   curtis1998_config.json        measure ROM (log ratio of means, no bias
                                 correction), REML, HKSJ, moderators
                                 [duration_weeks], subgroups [population]
   curtis1998_expected_study_effects.csv, _pooled.csv, _metareg.csv
   (~ duration_weeks), _small_study.csv, _subgroups.csv, _leave1out.csv,
   _influence.csv

Data paths these sets do not cover (crossover, change scores, median
conversion, multi-arm splitting) are covered by the SYNTHETIC example in
reference_r/examples/ and the hand checks in reference_r/tests/check_reference.R.

Engine configs (bcg_config.json, normand1999_config.json, curtis1998_config.json)
  These pin the settings the expected values were generated with, which since the
  2026-11-02 review differ from the kit defaults: hksj_rule: always (metafor test =
  "knha" for every random-effects fit), subgroup_method: separate, subgroup_min_k,
  subgroup_min_per_level and metareg_min_k: 0, min_k_pool: 1. They exist to compare
  formulas with metafor, not to model the protocol.
