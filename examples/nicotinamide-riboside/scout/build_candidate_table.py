"""Build examples/nicotinamide-riboside/scout/nr_candidate_trials.csv from hand-screened records, verifying each
PMID with NCBI esummary and each NCT ID with the ClinicalTrials.gov API v2.

Usage:
    .venv/bin/python examples/nicotinamide-riboside/scout/build_candidate_table.py            # live verification
    .venv/bin/python examples/nicotinamide-riboside/scout/build_candidate_table.py --offline  # skip network, mark unchecked

What is verified and what is not
--------------------------------
* first_author, year, title, journal and doi come from esummary at build time
  (column pmid_check says "esummary ok <date>"). A PMID that esummary cannot
  resolve, or whose first author does not match the screened record, is an error.
* ctgov_* columns come from the ClinicalTrials.gov API v2 at build time. Fetching a
  registry record checks that the NCT ID exists; it does not check that the
  publication reports that trial. registration_link says how the link was made:
  "stated" (in the abstract or the PubMed DataBank field), "probable (unverified)"
  (matched by the screener on title, design, n or population; registration_id starts
  with UNVERIFIED), "registry record" (registry-only rows) or blank (none found).
  ctgov_check repeats the caveat for probable links.
* Every other field was read from the PubMed abstract (or the registry record
  where the column says so) during title/abstract screening on 2026-11-01/02.
  Full texts were not read. A value that could not be confirmed from those
  sources starts with "UNVERIFIED" and the field name is listed in
  unverified_fields. Treat the whole table as a screening aid: confirm every
  design, n, dose and outcome against the full text before extraction.

record_id follows the order of ROWS, then REGISTRY_ONLY. examples/nicotinamide-riboside/scout/pooling_issues.txt
cites these IDs, so add new screened records at the end of ROWS and rerun; then
update the R-numbers in pooling_issues.txt for the registry-only rows, which shift.

Outcome family codes (see examples/nicotinamide-riboside/scout/pooling_issues.txt):
NAD = blood NAD+ metabolome; BP = blood pressure and arterial stiffness;
GLU = glucose and insulin sensitivity; LIP = lipids; INF = inflammatory markers;
BODY = body composition; MUSC = muscle and physical function; COG = cognition;
OTHER = anything else (named in outcomes_reported).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "nr_candidate_trials.csv"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
CTGOV = "https://clinicaltrials.gov/api/v2/studies/"

U = "UNVERIFIED"
NOT_IN_ABS = f"{U}: not stated in PubMed abstract"

# Suggested screening status values (a suggestion for the two screeners, not a decision):
#   RCT_NR_vs_control      randomized; an NR-alone vs placebo/no-NR contrast exists
#   RCT_NR_combination     randomized; NR given only as part of a mixture, effect not separable
#   RCT_check_design       randomized per abstract but design details need the full text
#   secondary_report       another report of a trial listed elsewhere (link by trial_key)
#   protocol_only          protocol; no results yet
#   nonrandomized          single-arm, open-label pre/post, PK or cohort; likely excluded from RCT synthesis
#   no_eligible_comparator randomized, but no placebo or no-treatment arm (e.g. doses of NR only); PRISMA reason: no eligible comparator
#   registry_only          completed or terminated registered trial with no PubMed report linked

ROWS: list[dict] = [
    # ---- NR alone vs placebo (or a factorial/multi-arm trial with a clean NR contrast) ----
    dict(trial_key="NOPARK", report_role="primary", pmid="42848371", author_check="Brakedal",
         randomized="yes", design="parallel, 2 arms, double-blind, phase 3, 11 centres",
         population="early Parkinson disease (diagnosis within 2 years, Hoehn and Yahr < 3), Norway",
         population_category="disease_parkinson", n_randomized="410 (NR 206, placebo 204); 393 in primary analysis",
         nr_dose_mg_day="1000 (500 mg twice daily)", duration="52 weeks", comparator="matching placebo",
         intervention_type="NR alone",
         outcomes_reported="MDS-UPDRS I-III total change (primary, MMRM); MDS-UPDRS parts I, II, III; DAT-SPECT; Non-Motor Symptoms Scale",
         outcome_families="OTHER (neurological rating scales, imaging)",
         registration_id="NCT03568968", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_vs_control", unverified_fields="", notes="Largest NR trial found. Change-from-baseline primary outcome."),
    dict(trial_key="Martens_aMCI", report_role="primary", pmid="42478598", author_check="Martens",
         randomized="yes", design="parallel, 2 arms, double-blind, phase II pilot",
         population="older adults with amnestic mild cognitive impairment",
         population_category="disease_mci", n_randomized=f"{U}: abstract gives 42 completers (NR 22, placebo 20); registry probable match lists 64 enrolled",
         nr_dose_mg_day=f"{U}: not in abstract; probable registry record gives 4 x 250 mg capsules daily (1000)",
         duration="12 weeks", comparator="placebo", intervention_type="NR alone",
         outcomes_reported="cognitive function (primary); total and regional cerebral blood flow; blood pressure; blood NAD+; adverse events",
         outcome_families="COG; BP; NAD; OTHER (cerebral blood flow)",
         registration_id=f"{U}: probable NCT03482167 (title, population, 12-week primary timepoint match)", registration_basis="registry search; not stated in abstract",
         status="RCT_NR_vs_control", unverified_fields="n_randomized; nr_dose_mg_day; registration_id",
         notes="Registry has results posted. Confirm the n randomized and the dose in the full text."),
    dict(trial_key="Christen_3boosters", report_role="primary", pmid="41540253", author_check="Christen",
         randomized="yes", design="parallel, 4 arms (nicotinamide, NR, NMN, placebo), open-label",
         population="healthy adults", population_category="healthy_adult",
         n_randomized="65 per abstract (registry lists 68 enrolled)",
         nr_dose_mg_day="1000 (per registry; abstract gives no dose)", duration="14 days", comparator="placebo (microcrystalline cellulose, per registry)",
         intervention_type="NR alone (one arm of a multi-precursor trial)",
         outcomes_reported="whole-blood NAD+ and NAD+ metabolome (chronic and acute); ex vivo microbial fermentation",
         outcome_families="NAD; OTHER (microbiome)",
         registration_id="NCT05517122", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_vs_control", unverified_fields="nr_dose_mg_day (registry only)",
         notes="Multi-arm: NR, NMN and nicotinamide share one placebo arm. Open-label. The NMN arm is the only head-to-head NR vs NMN comparison found."),
    dict(trial_key="Thomas_smallfiber", report_role="primary", pmid="41502156", author_check="Thomas",
         randomized="yes", design="parallel, 2 arms, double-blind, phase 2",
         population=f"{U}: experimental capsaicin nerve-injury model; registry eligibility 18-65 y, condition 'small fiber neuropathy'",
         population_category=f"{U}: likely healthy_adult (experimental model)",
         n_randomized=f"{U}: not in abstract; registry lists 45 enrolled",
         nr_dose_mg_day=NOT_IN_ABS, duration=f"{U}: visits to 90 days", comparator="placebo",
         intervention_type="NR alone",
         outcomes_reported="plasma NAD+; skin NAD+; epidermal nerve fibre degeneration and reinnervation (days 60, 90)",
         outcome_families="NAD; OTHER (intraepidermal nerve fibre density)",
         registration_id=f"{U}: probable NCT03912220 (registry title identical to article title)", registration_basis="registry search",
         status="RCT_NR_vs_control", unverified_fields="population; population_category; n_randomized; nr_dose_mg_day; duration; registration_id",
         notes="Registry results report medians for the primary outcomes."),
    dict(trial_key="Wu_longCOVID", report_role="primary", pmid="41357333", author_check="Wu",
         randomized="yes", design="parallel, 2:1, double-blind, placebo-controlled for 10 weeks then delayed start (NR-NR vs placebo-NR)",
         population="community-dwelling adults with long COVID, Boston", population_category="disease_long_covid",
         n_randomized="58 (NR-NR 37, placebo-NR 21)", nr_dose_mg_day="2000", duration="20 weeks NR (randomized placebo comparison for first 10 weeks)",
         comparator="placebo for 10 weeks, then NR", intervention_type="NR alone",
         outcomes_reported="ECog, RBANS, TMT-B (primary cognition); Fatigue Severity Scale; BDI; BAI; PSQI; whole-blood NAD+",
         outcome_families="COG; NAD; OTHER (fatigue, mood, sleep)",
         registration_id="NCT04809974", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_vs_control", unverified_fields="",
         notes="Only the 0-10 week window is a randomized NR vs placebo contrast. Dropout 32.4% NR vs 14.3% placebo at 10 weeks. Registry has results posted."),
    dict(trial_key="Lin_NEET_hypertension", report_role="primary", pmid="40770531", author_check="Lin",
         randomized="yes", design="parallel, 3 arms (NR+exercise, placebo+exercise, NR alone), double-blind",
         population="sedentary adults >= 55 y with daytime SBP >= 130 mmHg", population_category="disease_hypertension",
         n_randomized="54; 49 completed (15, 16, 18)", nr_dose_mg_day="1000", duration="6 weeks",
         comparator="placebo + exercise (primary comparison NR+exercise vs placebo+exercise)", intervention_type="NR with co-intervention (exercise in both compared arms)",
         outcomes_reported="daytime SBP (primary); pulse wave velocity; NAD catabolites; nighttime BP",
         outcome_families="BP; NAD",
         registration_id=f"{U}: probable NCT04112043 (NEET trial; arms and 6-week primary match; registry enrollment 147)", registration_basis="registry search; not stated in abstract",
         status="RCT_NR_vs_control", unverified_fields="registration_id",
         notes="NR-alone arm has no placebo-alone partner, so only NR+exercise vs placebo+exercise gives an NR contrast. Abstract reports change +/- SD."),
    dict(trial_key="Shoji_Werner", report_role="primary", pmid="40459998", author_check="Shoji",
         randomized="yes", design="crossover, 2 periods of 26 weeks, double-blind",
         population="patients with Werner syndrome", population_category="disease_werner_syndrome",
         n_randomized=NOT_IN_ABS, nr_dose_mg_day="1000", duration="26 weeks per period", comparator="placebo",
         intervention_type="NR alone",
         outcomes_reported="safety (primary); plasma NAD+; skin ulcer number and area; blood tests; sarcopenia; heel pad thickness; CAVI; ankle-brachial index; plasma metabolome (creatinine)",
         outcome_families="NAD; BP (arterial stiffness, CAVI); MUSC; OTHER (skin ulcers, renal)",
         registration_id=NOT_IN_ABS, registration_basis="",
         status="RCT_NR_vs_control", unverified_fields="n_randomized; registration_id",
         notes="Washout length not given in abstract; check carry-over handling."),
    dict(trial_key="Wu_SCD_MCI_crossover", report_role="primary", pmid="39817194", author_check="Wu",
         randomized="yes", design="crossover, 2 sequences (NR-placebo, placebo-NR), double-blind",
         population="adults > 55 y with subjective cognitive decline or MCI", population_category="disease_mci",
         n_randomized="46; 41 completed baseline; 37 completed", nr_dose_mg_day="1000", duration="8 weeks per active period",
         comparator="placebo", intervention_type="NR alone",
         outcomes_reported="RBANS (primary); plasma pTau217, GFAP, NfL; Lumosity scores; step counts",
         outcome_families="COG; MUSC (step counts); OTHER (plasma AD biomarkers)",
         registration_id=f"{U}: probable NCT04078178 (crossover SCD/MCI Niagen trial; registry enrollment 61)", registration_basis="registry search; not stated in abstract",
         status="RCT_NR_vs_control", unverified_fields="registration_id",
         notes="Reports both MMRM between-group and paired within-person comparisons. pTau217 reported as percent change. Registry has results posted."),
    dict(trial_key="Norheim_COPD", report_role="primary", pmid="39548320", author_check="Norheim",
         randomized="yes", design="parallel, 2 arms, double-blind",
         population="stable COPD", population_category="disease_copd",
         n_randomized="40", nr_dose_mg_day="2000 (1 g twice daily, per registry; abstract gives no dose)", duration="6 weeks, follow-up 12 weeks after end",
         comparator="placebo", intervention_type="NR alone",
         outcomes_reported="sputum IL-8 change (primary, reported as % difference); whole-blood NAD+; plasma IL-6; airway gene expression; epigenetic age",
         outcome_families="INF; NAD; OTHER (epigenetic age)",
         registration_id="NCT04990869", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_vs_control", unverified_fields="nr_dose_mg_day (registry only)",
         notes="Primary effect is a ratio (log-scale) estimate. Registry enrollment (60) differs from the abstract n (40); check the CONSORT diagram."),
    dict(trial_key="NICE_PAD", report_role="primary", pmid="38871717", author_check="McDermott",
         randomized="yes", design="parallel, 3 arms (NR+resveratrol, NR+placebo, placebo+placebo), double-blind",
         population="lower extremity peripheral artery disease", population_category="disease_pad",
         n_randomized="90", nr_dose_mg_day="1000 (per registry); resveratrol 125 mg/day in one arm", duration="6 months",
         comparator="placebo + placebo", intervention_type="NR alone (NR+placebo arm); NR+resveratrol arm is a combination",
         outcomes_reported="6-minute walk distance (primary); registry lists 34 outcome measures",
         outcome_families="MUSC; OTHER (registry outcomes not checked)",
         registration_id="NCT03743636", registration_basis="registry record title names the NICE trial; enrollment 90 matches; not in abstract",
         status="RCT_NR_vs_control", unverified_fields="nr_dose_mg_day (registry only)",
         notes="Publisher correction PMID 39134546. Primary CI is one-sided 90%. Registry has results posted."),
    dict(trial_key="NR-SAFE", report_role="primary", pmid="38016950", author_check="Berven",
         randomized="yes", design="parallel, 1:1, double-blind, phase I",
         population="Parkinson disease, Bergen, Norway", population_category="disease_parkinson",
         n_randomized="20 (NR 10, placebo 10)", nr_dose_mg_day="3000 (1500 mg twice daily)", duration="4 weeks",
         comparator="placebo", intervention_type="NR alone",
         outcomes_reported="moderate and severe adverse events (primary); mild adverse events; whole-blood and urine NAD metabolome; MDS-UPDRS; serum homocysteine",
         outcome_families="NAD; OTHER (adverse events, MDS-UPDRS, homocysteine)",
         registration_id="NCT05344404", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_vs_control", unverified_fields="", notes="Highest dose found."),
    dict(trial_key="Orr_MCI", report_role="primary", pmid="37994989", author_check="Orr",
         randomized="yes", design="parallel, 2 arms, placebo-controlled pilot",
         population="older adults with mild cognitive impairment", population_category="disease_mci",
         n_randomized="20", nr_dose_mg_day="dose escalation to 1000", duration="10 weeks",
         comparator="placebo", intervention_type="NR alone",
         outcomes_reported="MoCA change (primary); cerebral blood flow; blood NAD+; neurocognitive and psychometric tests; physical performance (walking speed); PBMC DNA methylation",
         outcome_families="COG; NAD; MUSC; OTHER (cerebral blood flow, epigenetic age)",
         registration_id="NCT02942888", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_vs_control", unverified_fields="",
         notes="Registry lists a FACTORIAL model and 46 enrolled; abstract says 20 randomized. Check which participants this report covers."),
    dict(trial_key="Ahmadi_CKD", report_role="primary", pmid="37159264", author_check="Ahmadi",
         randomized="yes", design="crossover, 3 periods (NR, CoQ10, placebo), 6 sequences, double-blind",
         population="CKD with eGFR < 60 mL/min/1.73 m2", population_category="disease_ckd",
         n_randomized="25 (registry lists 26 enrolled)", nr_dose_mg_day="1000", duration="6 weeks per period",
         comparator="placebo period (CoQ10 period is a second active comparator)", intervention_type="NR alone",
         outcomes_reported="VO2 peak and work efficiency (primary); submaximal VO2; eGFR; plasma metabolomics and lipidomics (triglycerides, ceramides)",
         outcome_families="MUSC; LIP; OTHER (eGFR, metabolomics)",
         registration_id="NCT03579693", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_vs_control", unverified_fields="",
         notes="Second report PMID 39847432 (inflammation, oxidative stress); preprint PMID 39228730. Registry has results posted."),
    dict(trial_key="Wang_HFrEF", report_role="primary", pmid="36644285", author_check="Wang",
         randomized=f"{U}: registry says randomized with NR and placebo groups; abstract does not say", design=f"{U}: NR vs placebo; allocation ratio not in abstract",
         population="clinically stable heart failure with reduced ejection fraction", population_category="disease_heart_failure",
         n_randomized="30 (registry 30 enrolled)", nr_dose_mg_day="2000 (1000 mg twice daily)", duration=f"{U}: registry primary timeframe up to 12 weeks",
         comparator="matching placebo (registry)", intervention_type="NR alone",
         outcomes_reported="safety and tolerability; whole-blood NAD+; PBMC basal and maximal respiration; NLRP3 expression",
         outcome_families="NAD; INF; OTHER (PBMC respiration, adverse events)",
         registration_id="NCT03423342", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_check_design", unverified_fields="randomized; design; duration",
         notes="Registry has results posted. Registry allocation reads RANDOMIZED but model SINGLE_GROUP despite two groups, and the abstract never says randomized: settle randomization and design from the full text before inclusion."),
    dict(trial_key="NADPARK", report_role="primary", pmid="35235774", author_check="Brakedal",
         randomized="yes", design="parallel, 2 arms, double-blind, phase I",
         population="newly diagnosed, treatment-naive Parkinson disease", population_category="disease_parkinson",
         n_randomized="30", nr_dose_mg_day="1000", duration="30 days", comparator="placebo", intervention_type="NR alone",
         outcomes_reported="cerebral NAD (31P-MRS); CSF NAD metabolites; FDG-PET; MDS-UPDRS; blood NAD metabolome; blood and muscle transcriptome; serum and CSF inflammatory cytokines",
         outcome_families="NAD; INF; OTHER (brain imaging, MDS-UPDRS)",
         registration_id="NCT03816020", registration_basis="PubMed DataBank (not in abstract text)",
         status="RCT_NR_vs_control", unverified_fields="", notes="Second report PMID 36936793 (DNA methylation)."),
    dict(trial_key="Remie_NRMetabolic", report_role="primary", pmid="32320006", author_check="Remie",
         randomized="yes", design="crossover, 2 periods, double-blind",
         population="healthy overweight or obese men and women", population_category="obesity_insulin_resistance",
         n_randomized="13 (registry lists 15 enrolled)", nr_dose_mg_day="1000", duration="6 weeks per period",
         comparator="placebo", intervention_type="NR alone",
         outcomes_reported="insulin sensitivity (clamp); mitochondrial function; fat-free mass %; sleeping metabolic rate; muscle acetylcarnitine; hepatic and intramyocellular lipid; cardiac energy status and ejection fraction; ambulatory BP; plasma inflammatory markers",
         outcome_families="GLU; BODY; MUSC; BP; INF; OTHER (energy metabolism, cardiac)",
         registration_id="NCT02835664", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_vs_control", unverified_fields="",
         notes="Abstract reports NR vs placebo period values with change +/- SE. BAT report PMID 33524145 shares this registration."),
    dict(trial_key="Nascimento_BAT", report_role=f"{U}: primary or substudy of NCT02835664", pmid="33524145", author_check="Nascimento",
         randomized="yes (human part)", design="crossover, double-blind (human part); also mouse and in vitro work",
         population="volunteers 45-65 y, BMI 27-35", population_category="obesity_insulin_resistance",
         n_randomized=NOT_IN_ABS, nr_dose_mg_day="1000", duration="6 weeks", comparator="placebo", intervention_type="NR alone",
         outcomes_reported="BAT [18F]FDG uptake; cold-induced thermogenesis",
         outcome_families="OTHER (brown adipose tissue activity, thermogenesis)",
         registration_id="NCT02835664; NCT03111719", registration_basis="PubMed DataBank",
         status="RCT_check_design", unverified_fields="report_role; n_randomized",
         notes="Shares NCT02835664 with Remie 2020. Check participant overlap before counting as a separate trial (Handbook 5.2.1)."),
    dict(trial_key="Stocks_exercise", report_role="primary", pmid="33492681", author_check="Stocks",
         randomized=f"{U}: abstract does not say", design=f"{U}: NR vs placebo for 1 week; crossover or parallel not stated in abstract",
         population="young trained men (age 23 +/- 4 y)", population_category="healthy_adult",
         n_randomized="8", nr_dose_mg_day="1000", duration="7 days", comparator="cellulose placebo",
         intervention_type="NR alone",
         outcomes_reported="substrate utilisation at rest and during exercise; muscle mitochondrial respiration; protein acetylation and PARylation; muscle NAD+ metabolome",
         outcome_families="NAD; MUSC; OTHER (substrate metabolism)",
         registration_id=NOT_IN_ABS, registration_basis="",
         status="RCT_check_design", unverified_fields="randomized; design; registration_id", notes=""),
    dict(trial_key="Dollerup_obese_men", report_role="primary", pmid="29992272", author_check="Dollerup",
         randomized="yes", design="parallel, 2 arms, double-blind",
         population="sedentary obese insulin-resistant men, BMI > 30, 40-70 y", population_category="obesity_insulin_resistance",
         n_randomized="40", nr_dose_mg_day="2000 (1000 mg twice daily)", duration="12 weeks", comparator="placebo", intervention_type="NR alone",
         outcomes_reported="insulin sensitivity (clamp); endogenous glucose production; glucose disposal and oxidation; resting energy expenditure; lipolysis and lipid oxidation; body composition (DXA, MRI); intrahepatic lipid; safety",
         outcome_families="GLU; BODY; OTHER (energy metabolism, liver fat)",
         registration_id="NCT02303483", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_vs_control", unverified_fields="",
         notes="Further reports PMID 31390002 (OGTT, beta-cell, incretins) and 31710095 (muscle mitochondria)."),
    dict(trial_key="Elhassan_aged_muscle", report_role="primary", pmid="31412242", author_check="Elhassan",
         randomized="yes", design="crossover, double-blind",
         population="aged men", population_category="healthy_older",
         n_randomized="12", nr_dose_mg_day="1000", duration="21 days", comparator="placebo", intervention_type="NR alone",
         outcomes_reported="muscle NAD+ metabolome; muscle RNA-seq; mitochondrial bioenergetics; circulating inflammatory cytokines",
         outcome_families="NAD; INF; MUSC",
         registration_id=f"{U}: probable NCT02950441 (crossover, 12 enrolled, NR and skeletal muscle)", registration_basis="registry search; not stated in abstract",
         status="RCT_NR_vs_control", unverified_fields="registration_id", notes=""),
    dict(trial_key="Conze_overweight", report_role="primary", pmid="31278280", author_check="Conze",
         randomized="yes", design="parallel, 4 arms (NR 100, 300, 1000 mg/day, placebo), double-blind",
         population="healthy overweight men and women", population_category="healthy_adult",
         n_randomized=f"{U}: not in abstract; probable registry record lists 140 enrolled",
         nr_dose_mg_day="100; 300; 1000 (three arms)", duration="8 weeks", comparator="placebo", intervention_type="NR alone (three doses)",
         outcomes_reported="whole-blood NAD+ (+22%, +51%, +142%) and NAD+ metabolites; adverse events; LDL cholesterol; one-carbon metabolism",
         outcome_families="NAD; LIP; OTHER (adverse events, one-carbon metabolism)",
         registration_id=f"{U}: probable NCT02712593 (Niagen 100/300/1000/placebo in healthy adults)", registration_basis="registry search; not stated in abstract",
         status="RCT_NR_vs_control", unverified_fields="n_randomized; registration_id",
         notes="Multiple doses share one placebo arm. NAD+ effects given as percent change (ratio outcome)."),
    dict(trial_key="Dolopikou_acute", report_role="primary", pmid="30725213", author_check="Dolopikou",
         randomized=f"{U}: abstract says double-blind crossover; randomization not stated", design="crossover, single acute dose, measured 2 h after dose",
         population="12 young and 12 old men", population_category="healthy_older (plus young comparison group)",
         n_randomized="24 (12 young, 12 old)", nr_dose_mg_day=NOT_IN_ABS, duration="single dose (2 hours)", comparator="placebo",
         intervention_type="NR alone",
         outcomes_reported="erythrocyte NADH and NADPH; urine F2-isoprostanes; erythrocyte glutathione; VO2max; concentric and isometric peak torque; fatigue index",
         outcome_families="NAD; MUSC; OTHER (oxidative stress)",
         registration_id=NOT_IN_ABS, registration_basis="",
         status="RCT_check_design", unverified_fields="randomized; nr_dose_mg_day; registration_id",
         notes="Single-dose study; may fall outside a minimum-duration criterion."),
    dict(trial_key="Martens_2018_crossover", report_role="primary", pmid="29599478", author_check="Martens",
         randomized="yes", design="crossover, 2 x 6 weeks, double-blind",
         population="healthy middle-aged and older adults", population_category="healthy_older",
         n_randomized=f"{U}: not in abstract; registry NCT02921659 lists 30 enrolled",
         nr_dose_mg_day="1000 (500 mg twice daily, per the abstract of secondary report PMID 36515353)", duration="6 weeks per period",
         comparator="placebo", intervention_type="NR alone",
         outcomes_reported="NAD+ metabolism; safety and tolerability; blood pressure and arterial stiffness (exploratory); other physiological function",
         outcome_families="NAD; BP; OTHER",
         registration_id="NCT02921659 (stated in the abstract of secondary report PMID 36515353)", registration_basis="via secondary report",
         status="RCT_NR_vs_control", unverified_fields="n_randomized",
         notes="Secondary report PMID 36515353 (neuronal extracellular vesicles, 22 participants)."),
    dict(trial_key="Han_psoriasis", report_role="primary", pmid="42048163", author_check="Han",
         randomized="yes", design="parallel, 2 arms, placebo-controlled",
         population="mild-to-moderate psoriasis", population_category="disease_psoriasis",
         n_randomized=f"{U}: not in abstract; registry NCT04271735 lists 29 enrolled",
         nr_dose_mg_day="1000 (500 mg twice daily)", duration="4 weeks", comparator="matching placebo", intervention_type="NR alone",
         outcomes_reported="Th17 responsiveness (IL-17 secretion change, registry primary); CD4+ T-cell RNA-seq; circulating SLIT2",
         outcome_families="INF; OTHER (T-cell transcriptome)",
         registration_id="NCT04271735 (also lists NCT01143454, NCT01778569, NCT00001846)", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_vs_control", unverified_fields="n_randomized", notes="Registry has results posted."),
    dict(trial_key="NIH_NR_immunity_crossover", report_role=f"{U}: report using participants of NCT02812238", pmid="37586364", author_check="Han",
         randomized=f"{U}: registry NCT02812238 is a randomized crossover; abstract describes 'oral NR supplementation in healthy volunteers'", design=f"{U}: crossover per registry (7 days NR, 7 days placebo)",
         population="healthy volunteers (in vivo part); psoriasis cells ex vivo", population_category="healthy_adult",
         n_randomized=f"{U}: registry lists 38 enrolled", nr_dose_mg_day="1000 (per registry)", duration="7 days per period (per registry)",
         comparator="placebo (per registry)", intervention_type="NR alone",
         outcomes_reported="serum arginine; sqstm1 and antioxidant gene expression; Th17 responsiveness",
         outcome_families="INF; OTHER (gene expression, amino acids)",
         registration_id="NCT02812238 (PubMed DataBank also lists NCT01934660, NCT01143454, NCT00001846)", registration_basis="PubMed DataBank",
         status="RCT_check_design", unverified_fields="report_role; randomized; design; n_randomized",
         notes="PMID 35025762 (Wu 2022, monocytes) also uses NCT02812238 participants. Registry has results posted (IL-1 beta release). Find the primary report or use registry results."),
    dict(trial_key="NIH_NR_immunity_crossover", report_role="secondary_report", pmid="35025762", author_check="Wu",
         randomized=f"{U}: 'in vivo placebo or NR administration' in healthy volunteers", design=f"{U}: see NCT02812238",
         population="healthy volunteers; SLE monocytes ex vivo", population_category="healthy_adult",
         n_randomized=NOT_IN_ABS, nr_dose_mg_day=NOT_IN_ABS, duration=NOT_IN_ABS, comparator="placebo",
         intervention_type="NR alone",
         outcomes_reported="monocyte RNA-seq; LPS-induced IFN-beta; autophagy; inosine",
         outcome_families="INF; OTHER",
         registration_id="NCT02812238; NCT00001846; NCT00001372", registration_basis="stated in abstract",
         status="secondary_report", unverified_fields="randomized; design; n_randomized; nr_dose_mg_day; duration", notes=""),

    # ---- Randomized, NR given as part of a combination or factorial with exercise ----
    dict(trial_key="Lin_Friedreich", report_role="primary", pmid="42009009", author_check="Lin",
         randomized="yes", design="2 x 2 factorial (placebo, NR, exercise+placebo, exercise+NR), phase 2, single centre",
         population="Friedreich ataxia, age 10-40 y, ejection fraction >= 45%", population_category="disease_friedreich_ataxia",
         n_randomized=f"{U}: not in the abstract text retrieved; registry lists 75 enrolled including non-randomized patients",
         nr_dose_mg_day="300, 600 or 900 by body weight", duration="12 weeks",
         comparator="placebo (with attention control) and placebo + exercise", intervention_type="NR in factorial with exercise",
         outcomes_reported="change in peak VO2 (L/min, primary); safety",
         outcome_families="MUSC; OTHER (safety)",
         registration_id="NCT04192136", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_vs_control", unverified_fields="n_randomized",
         notes="Two NR contrasts (with and without exercise). Registry results report median within-participant change. Weight-based dose."),
    dict(trial_key="Bhandari_CCS", report_role="primary", pmid="39387327", author_check="Bhandari",
         randomized="yes", design="parallel, 1:1, telehealth exercise with vs without NR",
         population="adult survivors of childhood cancer with prediabetes", population_category="disease_prediabetes_cancer_survivor",
         n_randomized="20", nr_dose_mg_day=NOT_IN_ABS, duration="6 weeks",
         comparator=f"{U}: exercise alone; whether a placebo was used is not stated", intervention_type="NR with co-intervention (exercise in both arms)",
         outcomes_reported="feasibility (primary); glucose homeostasis biomarkers; myostatin; intramuscular adipose tissue; muscle cross-sectional area",
         outcome_families="GLU; BODY; MUSC; OTHER (feasibility)",
         registration_id=f"{U}: probable NCT05023993 (exercise and NR in survivors; 20 enrolled)", registration_basis="registry search; not stated in abstract",
         status="RCT_NR_vs_control", unverified_fields="nr_dose_mg_day; comparator; registration_id", notes="Feasibility pilot."),
    dict(trial_key="NRPT_Dellinger_2017", report_role="primary", pmid="29184669", author_check="Dellinger",
         randomized="yes", design="parallel, 3 arms (placebo, NRPT 1X, NRPT 2X), double-blind",
         population="healthy adults 60-80 y", population_category="healthy_older",
         n_randomized="120", nr_dose_mg_day=f"{U}: not in abstract; registry arm labels 'Basis 250' and 'Basis 500' suggest 250 and 500 mg NR with pterostilbene",
         duration="8 weeks", comparator="placebo", intervention_type="combination product (NR + pterostilbene)",
         outcomes_reported="whole-blood NAD+ (about +40% and +90%); safety",
         outcome_families="NAD; OTHER (safety)",
         registration_id=f"{U}: probable NCT02678611 (Basis in elderly, 120 enrolled)", registration_basis="registry search; not stated in abstract",
         status="RCT_NR_combination", unverified_fields="nr_dose_mg_day; registration_id",
         notes="Author correction PMID 30155270. Pterostilbene in active arms only."),
    dict(trial_key="NRPT_NAFLD", report_role="primary", pmid="36082508", author_check="Dellinger",
         randomized="yes", design="parallel, 3 arms (placebo, NRPT 1x, NRPT 2x), double-blind",
         population="adults with NAFLD", population_category="disease_nafld",
         n_randomized="111", nr_dose_mg_day="250 (1x) and 500 (2x) with pterostilbene 50 and 100 mg (per registry capsule contents)",
         duration="6 months", comparator="placebo", intervention_type="combination product (NR + pterostilbene)",
         outcomes_reported="hepatic fat fraction (primary); ALT; GGT; ceramide 14:0",
         outcome_families="LIP (ceramides); OTHER (liver fat, liver enzymes)",
         registration_id="NCT03513523", registration_basis="PubMed DataBank",
         status="RCT_NR_combination", unverified_fields="nr_dose_mg_day (registry only)", notes=""),
    dict(trial_key="NRPT_muscle_injury", report_role="primary", pmid="35998039", author_check="Jensen",
         randomized="yes", design="parallel, 2 arms, placebo-controlled",
         population="adults 55-80 y with experimentally induced muscle injury", population_category="healthy_older",
         n_randomized="32", nr_dose_mg_day="1000 NR + 200 mg pterostilbene", duration="2 weeks before injury, biopsies to 30 days after",
         comparator="matched placebo", intervention_type="combination product (NR + pterostilbene)",
         outcomes_reported="muscle stem cell content and proliferation; fibre area; central nuclei; embryonic myosin; myoglobin, creatine kinase, soreness, strength",
         outcome_families="MUSC; OTHER (muscle regeneration histology)",
         registration_id="NCT03754842", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_combination", unverified_fields="", notes="Later report PMID 40473693 uses biopsies from this trial (not an NR comparison)."),
    dict(trial_key="NRPT_AKI", report_role="primary", pmid="32791973", author_check="Simic",
         randomized="yes", design="parallel, stepwise dose escalation, 5 NRPT : 1 placebo per step, double-blind",
         population="hospitalized patients with acute kidney injury", population_category="disease_aki",
         n_randomized="24 (NRPT 20, placebo 4)", nr_dose_mg_day="500, 1000, 1500, 2000 by step (250-1000 mg twice daily) + pterostilbene",
         duration="2 days", comparator="placebo", intervention_type="combination product (NR + pterostilbene)",
         outcomes_reported="whole-blood NAD+ at 48 h; safety laboratory tests (creatinine, eGFR, electrolytes, liver tests, blood counts); adverse events",
         outcome_families="NAD; OTHER (renal function, safety)",
         registration_id="NCT03176628", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_combination", unverified_fields="", notes="Four dose steps share four placebo patients."),
    dict(trial_key="EH301_ALS", report_role="primary", pmid="30668199", author_check="de la Rubia",
         randomized="yes", design="parallel, 2 arms, double-blind pilot",
         population="amyotrophic lateral sclerosis", population_category="disease_als",
         n_randomized="32", nr_dose_mg_day=NOT_IN_ABS, duration="4 months", comparator="placebo",
         intervention_type="combination product (EH301 = NR + pterostilbene)",
         outcomes_reported="ALSFRS-R; forced vital capacity; MRC muscle strength; skeletal muscle/fat weight ratio",
         outcome_families="MUSC; BODY; OTHER (ALSFRS-R, lung function)",
         registration_id=f"{U}: probable NCT03489200 (EH301 for ALS, 32 enrolled); PMID 37810917 cites the same NCT for a 4-group, 40-patient study",
         registration_basis="registry search; not stated in abstract",
         status="RCT_NR_combination", unverified_fields="nr_dose_mg_day; registration_id",
         notes="Not retrieved by NR text words (abstract says only 'EH301'); found by product-name search."),
    dict(trial_key="ALS_MeDi_NRPT", report_role="primary", pmid="37810917", author_check="Carrera-Juli",
         randomized="yes (abstract: 'prospective, mixed, randomized')", design="4 groups (Mediterranean diet control, + NR/pterostilbene, + coconut oil, + both)",
         population="amyotrophic lateral sclerosis", population_category="disease_als",
         n_randomized="40 (4 groups)", nr_dose_mg_day=NOT_IN_ABS, duration=NOT_IN_ABS,
         comparator="Mediterranean diet without supplements", intervention_type="combination product with diet",
         outcomes_reported="weight; fat and muscle mass %; skinfolds; perimeters; BMI; waist-hip and waist-height ratios",
         outcome_families="BODY",
         registration_id="NCT03489200", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_combination", unverified_fields="nr_dose_mg_day; duration",
         notes="Same NCT as EH301_ALS (PMID 30668199) but different group structure and n; resolve whether these are one trial. Further report PMID 39086545 (ketone bodies)."),
    dict(trial_key="CMA_AD_PD", report_role="primary", pmid="36703196", author_check="Yulug",
         randomized="yes", design="parallel, 2 arms, double-blind, phase II",
         population="Alzheimer disease", population_category="disease_alzheimer",
         n_randomized=NOT_IN_ABS, nr_dose_mg_day="1000 per CMA dose; once daily days 0-28, twice daily days 28-84 (so 1000 then 2000)",
         duration="84 days", comparator="placebo (sorbitol, per registry)", intervention_type="combination product (CMA: L-serine, NAC, NR, L-carnitine)",
         outcomes_reported="ADAS-Cog and daily living scores (primary); safety; hippocampal volume and cortical thickness; plasma metabolome and proteome",
         outcome_families="COG; OTHER (imaging, omics)",
         registration_id="NCT04044131", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_combination", unverified_fields="n_randomized",
         notes="Registry NCT04044131 covers both the AD and the PD trial (120 enrolled in total)."),
    dict(trial_key="CMA_AD_PD", report_role="primary (PD arm of registration)", pmid="39816194", author_check="Yulug",
         randomized="yes", design="parallel, 2 arms, double-blind, phase II",
         population="Parkinson disease", population_category="disease_parkinson",
         n_randomized=NOT_IN_ABS, nr_dose_mg_day="1000 per CMA dose; once daily 28 days then twice daily 56 days",
         duration="84 days", comparator="placebo", intervention_type="combination product (CMA)",
         outcomes_reported="UPDRS (primary); MoCA; functional MRI; plasma metabolomics and proteomics; safety",
         outcome_families="COG; OTHER (UPDRS, imaging, omics)",
         registration_id="NCT04044131", registration_basis="stated in abstract and PubMed DataBank",
         status="RCT_NR_combination", unverified_fields="n_randomized", notes=""),
    dict(trial_key="CMA_COVID", report_role="primary", pmid="34180141", author_check="Altay",
         randomized="yes (phase 3 double-blind; phase 2 placebo-controlled open-label)", design="two trials: phase 2 open-label and phase 3 double-blind, placebo-controlled",
         population="ambulatory mild-to-moderate COVID-19", population_category="disease_covid19_acute",
         n_randomized=NOT_IN_ABS, nr_dose_mg_day=NOT_IN_ABS, duration=NOT_IN_ABS,
         comparator="placebo", intervention_type="combination product (CMA)",
         outcomes_reported="time to symptom-free recovery; plasma metabolome and proteome",
         outcome_families="OTHER (recovery time, omics)",
         registration_id=NOT_IN_ABS, registration_basis="",
         status="RCT_NR_combination", unverified_fields="n_randomized; nr_dose_mg_day; duration; registration_id",
         notes="Two trials in one report; treat as two studies if both are included."),
    dict(trial_key="CMA_NAFLD", report_role="primary", pmid="34694070", author_check="Zeybel",
         randomized="yes (PubMed publication type RCT)", design="placebo-controlled", population="NAFLD", population_category="disease_nafld",
         n_randomized=f"{U}: registry NCT04330326 lists 32 enrolled", nr_dose_mg_day=NOT_IN_ABS, duration="10 weeks",
         comparator="placebo (sorbitol, per registry)", intervention_type="combination product (CMA)",
         outcomes_reported="hepatic steatosis; AST; ALT; uric acid; creatinine; metabolomics; proteomics; oral and gut metagenomics",
         outcome_families="OTHER (liver fat, liver enzymes, omics)",
         registration_id="NCT04330326 (registry record lists this PMID as derived)", registration_basis="registry record",
         status="RCT_NR_combination", unverified_fields="n_randomized; nr_dose_mg_day", notes=""),
    dict(trial_key="CMA_acute_precursors", report_role="primary", pmid="37271226", author_check="Li",
         randomized=f"{U}: PubMed type 'Controlled Clinical Trial'; randomization not stated", design="one-day, double-blind, placebo-controlled; six CMA formulations with different NAD+ precursors",
         population=f"{U}: not stated in abstract (likely healthy volunteers)", population_category=f"{U}",
         n_randomized=NOT_IN_ABS, nr_dose_mg_day="1000 (single dose, NR-containing formulation)", duration="1 day",
         comparator="placebo; other precursor formulations", intervention_type="combination product (CMA)",
         outcomes_reported="plasma metabolomics (NAD+ products); flushing; safety",
         outcome_families="NAD; OTHER",
         registration_id=NOT_IN_ABS, registration_basis="",
         status="RCT_NR_combination", unverified_fields="randomized; population; population_category; n_randomized; registration_id", notes="Single-day study."),
    dict(trial_key="BREAK-AD", report_role="primary", pmid="42344884", author_check="Roy",
         randomized="yes", design="parallel, 2 arms", population="mild cognitive impairment", population_category="disease_mci",
         n_randomized="30 (placebo 15, active 15)", nr_dose_mg_day=NOT_IN_ABS, duration="6 months",
         comparator="placebo", intervention_type="combination product (beta-hydroxybutyrate salts + NR)",
         outcomes_reported="brain ketone and glucose uptake (PET); white-matter myelin density (MRI); processing speed and cognition",
         outcome_families="COG; OTHER (brain imaging)",
         registration_id=NOT_IN_ABS, registration_basis="",
         status="RCT_NR_combination", unverified_fields="nr_dose_mg_day; registration_id", notes=""),

    # ---- Further reports of trials above ----
    dict(trial_key="Ahmadi_CKD", report_role="secondary_report", pmid="39847432", author_check="Ahmadi",
         randomized="yes", design="crossover, 3 periods", population="CKD eGFR < 60", population_category="disease_ckd",
         n_randomized="25; bioenergetics subsample 14", nr_dose_mg_day="1000", duration="6 weeks per period", comparator="placebo",
         intervention_type="NR alone",
         outcomes_reported="whole-blood transcriptome; plasma IL-2, IL-13, CRP; F2-isoprostanes; monocyte and lymphocyte bioenergetics",
         outcome_families="INF; OTHER (oxidative stress, transcriptome)",
         registration_id="NCT03579693", registration_basis="stated in abstract", status="secondary_report", unverified_fields="", notes="Same trial as PMID 37159264."),
    dict(trial_key="Ahmadi_CKD", report_role="preprint of secondary report", pmid="39228730", author_check="Ahmadi",
         randomized="yes", design="crossover", population="CKD", population_category="disease_ckd",
         n_randomized="25", nr_dose_mg_day="1000", duration="6 weeks per period", comparator="placebo", intervention_type="NR alone",
         outcomes_reported="as PMID 39847432", outcome_families="INF; OTHER",
         registration_id=f"{U}: not in abstract; same trial as NCT03579693", registration_basis="",
         status="secondary_report", unverified_fields="registration_id", notes="medRxiv preprint; do not double count."),
    dict(trial_key="Dollerup_obese_men", report_role="secondary_report", pmid="31390002", author_check="Dollerup",
         randomized="yes", design="parallel", population="obese insulin-resistant men", population_category="obesity_insulin_resistance",
         n_randomized="40 (20, 20)", nr_dose_mg_day="2000", duration="12 weeks", comparator="placebo", intervention_type="NR alone",
         outcomes_reported="OGTT glucose, insulin, C-peptide, glucagon, GLP-1, GIP; beta-cell indices; bile acids; adipsin",
         outcome_families="GLU; OTHER (incretins, bile acids)",
         registration_id="NCT02303483", registration_basis="PubMed DataBank", status="secondary_report", unverified_fields="", notes=""),
    dict(trial_key="Dollerup_obese_men", report_role="secondary_report", pmid="31710095", author_check="Dollerup",
         randomized="yes", design="parallel", population="obese insulin-resistant men", population_category="obesity_insulin_resistance",
         n_randomized="40", nr_dose_mg_day="2000", duration="12 weeks", comparator="placebo", intervention_type="NR alone",
         outcomes_reported="muscle mitochondrial respiration, content and morphology; NAMPT protein; muscle NAD+",
         outcome_families="MUSC; NAD", registration_id="NCT02303483", registration_basis="PubMed DataBank",
         status="secondary_report", unverified_fields="", notes=""),
    dict(trial_key="NADPARK", report_role="secondary_report", pmid="36936793", author_check="Gaare",
         randomized="yes", design="parallel", population="Parkinson disease", population_category="disease_parkinson",
         n_randomized="29 analysed", nr_dose_mg_day=f"{U}: 'high-dose NR' for up to 30 days", duration="up to 30 days", comparator="placebo",
         intervention_type="NR alone", outcomes_reported="DNA methylation in blood and muscle",
         outcome_families="OTHER (DNA methylation)", registration_id=f"{U}: not in abstract; probably NADPARK (NCT03816020)",
         registration_basis="", status="secondary_report", unverified_fields="nr_dose_mg_day; registration_id", notes=""),
    dict(trial_key="Martens_2018_crossover", report_role="secondary_report", pmid="36515353", author_check="Vreones",
         randomized="yes", design="crossover", population="healthy older adults", population_category="healthy_older",
         n_randomized="22 analysed", nr_dose_mg_day="1000 (500 mg twice daily)", duration="6 weeks per period", comparator="placebo",
         intervention_type="NR alone", outcomes_reported="NAD(H) and Abeta42, pJNK, pERK1/2 in neuron-enriched extracellular vesicles",
         outcome_families="NAD; OTHER (neuronal EV biomarkers)",
         registration_id="NCT02921659", registration_basis="stated in abstract", status="secondary_report", unverified_fields="", notes=""),
    dict(trial_key="NICE_PAD", report_role="erratum", pmid="39134546", author_check="McDermott",
         randomized="", design="", population="", population_category="", n_randomized="", nr_dose_mg_day="", duration="",
         comparator="", intervention_type="", outcomes_reported="publisher correction", outcome_families="",
         registration_id="", registration_basis="", status="secondary_report", unverified_fields="", notes="Read with PMID 38871717."),
    dict(trial_key="NRPT_Dellinger_2017", report_role="erratum", pmid="30155270", author_check="Dellinger",
         randomized="", design="", population="", population_category="", n_randomized="", nr_dose_mg_day="", duration="",
         comparator="", intervention_type="", outcomes_reported="author correction", outcome_families="",
         registration_id="", registration_basis="", status="secondary_report", unverified_fields="", notes="Read with PMID 29184669."),
    dict(trial_key="NRPT_muscle_injury", report_role="secondary_report", pmid="40473693", author_check="Brorson",
         randomized="", design="", population="older adults, muscle injury model", population_category="healthy_older",
         n_randomized="spatial RNA-seq n = 3; flow cytometry n = 9", nr_dose_mg_day="", duration="", comparator="",
         intervention_type="", outcomes_reported="muscle regeneration cell dynamics (no NR comparison reported in abstract)",
         outcome_families="OTHER", registration_id="NCT03754842", registration_basis="stated in abstract",
         status="secondary_report", unverified_fields="", notes="Probably not usable for an NR effect."),
    dict(trial_key="ALS_MeDi_NRPT", report_role="secondary_report", pmid="39086545", author_check="Carrera-Juli",
         randomized="", design="", population="ALS", population_category="disease_als", n_randomized="", nr_dose_mg_day="",
         duration="", comparator="", intervention_type="", outcomes_reported="circulating ketone bodies with coconut oil",
         outcome_families="OTHER", registration_id=NOT_IN_ABS, registration_basis="", status="secondary_report",
         unverified_fields="registration_id", notes="Report from the same pilot as PMID 37810917."),

    # ---- Protocols (no results yet) ----
    dict(trial_key="Freeberg_SBP", report_role="protocol", pmid="35620522", author_check="Freeberg",
         randomized="yes (planned)", design="parallel, 2 arms, double-blind (planned)", population="adults >= 50 y with SBP 120-159 mmHg",
         population_category="healthy_older", n_randomized="94 planned (registry estimate 118)", nr_dose_mg_day="1000 (500 mg twice daily)",
         duration="3 months", comparator="placebo", intervention_type="NR alone",
         outcomes_reported="casual SBP (primary); 24-h SBP; aortic stiffness; DBP; NAD+; oxidative stress and inflammation markers",
         outcome_families="BP; NAD; INF", registration_id="NCT03821623", registration_basis="stated in abstract",
         status="protocol_only", unverified_fields="", notes="Registry status active, not recruiting; watch for results."),
    dict(trial_key="Glaucoma_NR", report_role="protocol", pmid="35039056", author_check="",
         randomized="yes (planned)", design="parallel, multicentre, double-blind (planned)", population="primary open-angle glaucoma",
         population_category="disease_glaucoma", n_randomized="125 planned", nr_dose_mg_day="300", duration="24 months",
         comparator="placebo", intervention_type="NR alone", outcomes_reported="RNFL thinning rate (primary); visual field progression",
         outcome_families="OTHER (ophthalmic)", registration_id="ChiCTR1900021998", registration_basis="stated in abstract",
         status="protocol_only", unverified_fields="", notes="Correction PMID 35144657. Not on ClinicalTrials.gov; WHO ICTRP covers ChiCTR."),
    dict(trial_key="MOONRISE", report_role="protocol", pmid="42843881", author_check="Khanal",
         randomized="yes (planned)", design="parallel pilot feasibility (planned)", population="postpartum mothers of NICU infants",
         population_category="disease_lactation_insufficiency_risk", n_randomized="40 planned", nr_dose_mg_day=NOT_IN_ABS,
         duration="19 days", comparator="placebo", intervention_type="NR alone",
         outcomes_reported="feasibility (primary); milk volume and composition; serum and urine NAD+ metabolites",
         outcome_families="NAD; OTHER (lactation)", registration_id=f"{U}: probable NCT04614714 (registry lists this PMID)",
         registration_basis="registry record", status="protocol_only", unverified_fields="nr_dose_mg_day; registration_id", notes=""),
    dict(trial_key="HCT_AYA", report_role="protocol", pmid="35854224", author_check="Song",
         randomized="yes (planned)", design="2 x 2 factorial (exercise + NR, exercise, NR, control) (planned)",
         population="adolescent and young adult HCT survivors", population_category="disease_hct_survivor",
         n_randomized="80 planned", nr_dose_mg_day=NOT_IN_ABS, duration="16 weeks", comparator="control / placebo",
         intervention_type="NR in factorial with exercise",
         outcomes_reported="muscle strength (primary); cardiopulmonary fitness; DXA muscle mass; muscle OXPHOS; NAD+ metabolomics",
         outcome_families="MUSC; BODY; NAD", registration_id="NCT05194397", registration_basis="stated in abstract",
         status="protocol_only", unverified_fields="nr_dose_mg_day", notes="Registry status recruiting."),

    # ---- Non-randomized human NR studies (likely excluded from an RCT synthesis) ----
    dict(trial_key="Lapatto_twins", report_role="primary", pmid="36638183", author_check="Lapatto",
         randomized="no", design="open-label, single-arm in BMI-discordant monozygotic twin pairs", population="20 BMI-discordant MZ twin pairs",
         population_category="obesity_insulin_resistance", n_randomized="not randomized; 20 twin pairs", nr_dose_mg_day="250 escalating to 1000",
         duration="5 months", comparator="none (cotwin comparison)", intervention_type="NR alone",
         outcomes_reported="systemic NAD+ metabolism; muscle mitochondrial number; myoblast differentiation; gut microbiota; DNA methylation; adiposity; metabolic health",
         outcome_families="NAD; BODY; GLU; MUSC; OTHER", registration_id=f"{U}: probable NCT03951285 (NR in BMI-discordant twins)",
         registration_basis="registry search", status="nonrandomized", unverified_fields="registration_id",
         notes="Not found by the RCT filter lines alone; found by the 'supplementation' line. PMID 42481916 re-analyses this and other studies."),
    dict(trial_key="Airhart_PK", report_role="primary", pmid="29211728", author_check="Airhart",
         randomized="no", design="open-label, non-randomized PK", population="healthy volunteers", population_category="healthy_adult",
         n_randomized="not randomized; 8", nr_dose_mg_day="250 escalating to 2000 (1000 mg twice daily)", duration="9 days",
         comparator="none", intervention_type="NR alone", outcomes_reported="NR pharmacokinetics; whole-blood NAD+; clinical chemistry",
         outcome_families="NAD", registration_id="NCT02689882 (registry lists this PMID as result)", registration_basis="registry record",
         status="nonrandomized", unverified_fields="", notes=""),
    dict(trial_key="Trammell_PK", report_role="primary", pmid="27721479", author_check="Trammell",
         randomized=f"{U}: abstract does not say; probable registry NCT02191462 is a randomized crossover of three doses without placebo",
         design="single-dose PK: 100, 300, 1000 mg (plus one-person pilot and mouse work)", population="healthy volunteers",
         population_category="healthy_adult", n_randomized=f"{U}: not in abstract (registry 12)", nr_dose_mg_day="100; 300; 1000 single doses",
         duration="single doses", comparator="none (dose comparison)", intervention_type="NR alone",
         outcomes_reported="blood NAD+ metabolome including NAAD", outcome_families="NAD",
         registration_id=f"{U}: probable NCT02191462", registration_basis="registry search", status="no_eligible_comparator",
         unverified_fields="randomized; n_randomized; registration_id",
         notes="No placebo arm: the probable registry record is a randomized crossover of three NR doses, so the PRISMA exclusion reason is no eligible comparator, not non-randomized design."),
    dict(trial_key="AT_openlabel_4mo", report_role="primary", pmid="34515380", author_check="",
         randomized="no", design="open-label proof-of-concept, single-arm", population="ataxia telangiectasia", population_category="disease_ataxia_telangiectasia",
         n_randomized="not randomized; 24", nr_dose_mg_day=NOT_IN_ABS, duration="4 months", comparator="none (withdrawal phase)",
         intervention_type="NR alone", outcomes_reported="SARA and ICARS ataxia scores; dysarthria; quality of life; serum IgG; plasma metabolomics",
         outcome_families="OTHER (ataxia scores, immunoglobulins)", registration_id=f"{U}: probable NCT03962114 (24 enrolled)",
         registration_basis="registry search", status="nonrandomized", unverified_fields="nr_dose_mg_day; registration_id", notes=""),
    dict(trial_key="AT_openlabel_2y", report_role="primary", pmid="37899683", author_check="Presterud",
         randomized="no", design="single-arm, open-label, 2 years, historical comparison", population="ataxia telangiectasia",
         population_category="disease_ataxia_telangiectasia", n_randomized="not randomized; n not in abstract", nr_dose_mg_day=NOT_IN_ABS,
         duration="2 years", comparator="historical data", intervention_type="NR alone",
         outcomes_reported="blood NAD+; safety; neuromotor test panels (coordination, eye movements)", outcome_families="NAD; OTHER",
         registration_id=NOT_IN_ABS, registration_basis="", status="nonrandomized", unverified_fields="nr_dose_mg_day; registration_id",
         notes="Further report PMID 42449929 (PBMC transcriptome)."),
    dict(trial_key="AT_openlabel_2y", report_role="secondary_report", pmid="42449929", author_check="Suaste",
         randomized="no", design="longitudinal transcriptomics within the 24-month open-label trial", population="ataxia telangiectasia",
         population_category="disease_ataxia_telangiectasia", n_randomized="not randomized", nr_dose_mg_day="", duration="24 months",
         comparator="none", intervention_type="NR alone", outcomes_reported="PBMC interferon gene signatures", outcome_families="INF; OTHER",
         registration_id=NOT_IN_ABS, registration_basis="", status="nonrandomized", unverified_fields="registration_id", notes=""),
    dict(trial_key="NRPT_menopause", report_role="primary", pmid="42211736", author_check="Holmes",
         randomized="no", design="open-label, single-arm pilot", population="healthy women > 35 y (32 with menopause symptoms, 8 without)",
         population_category="healthy_adult", n_randomized="not randomized; 40", nr_dose_mg_day="250 NR + 50 mg pterostilbene",
         duration="7 days", comparator="none (baseline)", intervention_type="combination product (NR + pterostilbene)",
         outcomes_reported="menopause symptom survey; urinary estrone, estradiol and E2/E1 ratio; urinary vitamin B3 catabolome",
         outcome_families="NAD (urinary catabolome); OTHER (symptoms, hormones)", registration_id="NCT04841499", registration_basis="stated in abstract",
         status="nonrandomized", unverified_fields="", notes=""),
    dict(trial_key="PAD_openlabel", report_role="primary", pmid="40479886", author_check="Szarvas",
         randomized="no", design="open-label, single-arm pilot", population="older adults with PAD", population_category="disease_pad",
         n_randomized="not randomized; 8", nr_dose_mg_day=NOT_IN_ABS, duration="4 weeks", comparator="none",
         intervention_type="NR alone", outcomes_reported="peripheral endothelial function; cerebrovascular response (fNIRS); cognitive performance; ex vivo endothelial assays with patient sera",
         outcome_families="BP (endothelial function); COG; OTHER", registration_id="NCT06534944", registration_basis="stated in abstract",
         status="nonrandomized", unverified_fields="nr_dose_mg_day", notes=""),
    dict(trial_key="NAD-brain_PK", report_role="primary", pmid="41858901", author_check="Berven",
         randomized=f"{U}: abstract does not say; probable registry NCT05698771 is NON_RANDOMIZED", design="phase I PK; NR or NMN",
         population="healthy individuals (6) and Parkinson disease (6)", population_category="healthy_adult; disease_parkinson",
         n_randomized="12 (not randomized per probable registry)", nr_dose_mg_day="1200 (NR or NMN)", duration=f"{U}: at least 4 weeks plus washout",
         comparator="none", intervention_type="NR alone or NMN alone", outcomes_reported="blood NAD and metabolite kinetics; cerebral NAD",
         outcome_families="NAD", registration_id=f"{U}: probable NCT05698771", registration_basis="registry search",
         status="nonrandomized", unverified_fields="randomized; duration; registration_id", notes=""),
    dict(trial_key="NR_cerebral_MRS", report_role="primary", pmid="39044608", author_check="",
         randomized="no", design="single-arm acute dose, test-retest baseline scans", population="healthy volunteers",
         population_category="healthy_adult", n_randomized="not randomized; 10", nr_dose_mg_day="900 single dose", duration="single dose (4 h)",
         comparator="baseline scans", intervention_type="NR alone", outcomes_reported="cerebral NAD+ by downfield 1H-MRS at 7 T",
         outcome_families="NAD (brain)", registration_id=NOT_IN_ABS, registration_basis="", status="nonrandomized",
         unverified_fields="registration_id", notes=""),
    dict(trial_key="Zhou_HF_LVAD", report_role="primary", pmid="32790648", author_check="Zhou",
         randomized="no", design="open-label NR before LVAD; plus case-control and in vitro work", population="stage D heart failure",
         population_category="disease_heart_failure", n_randomized="not randomized; 4 treated (registry 5 enrolled)", nr_dose_mg_day=NOT_IN_ABS,
         duration="5-9 days", comparator="none (pre/post)", intervention_type="NR alone",
         outcomes_reported="PBMC respiration; proinflammatory cytokine gene expression", outcome_families="INF; OTHER",
         registration_id="NCT03727646", registration_basis="stated in abstract", status="nonrandomized",
         unverified_fields="nr_dose_mg_day", notes="Registry has results posted."),
    dict(trial_key="HF-AF_ENERGY", report_role="design paper", pmid="36227441", author_check="Pool",
         randomized="no", design="prospective pre/post intervention (4-month observation, 4-month NR)", population="ischaemic heart disease with AF and ICD",
         population_category="disease_heart_failure", n_randomized="not randomized; n not in abstract", nr_dose_mg_day=NOT_IN_ABS,
         duration="4 months", comparator="pre-treatment period", intervention_type="NR alone",
         outcomes_reported="AF burden (ICD monitoring); echocardiography; mitochondrial function markers; NAD metabolome",
         outcome_families="NAD; OTHER (AF burden)", registration_id=NOT_IN_ABS, registration_basis="",
         status="nonrandomized", unverified_fields="nr_dose_mg_day; registration_id", notes=""),
    dict(trial_key="Traba_fasting_NR", report_role="primary", pmid="26529255", author_check="Traba",
         randomized=f"{U}: not stated; 'additional group of 8 healthy volunteers' given NR", design=f"{U}: NR substudy within fasting/refeeding study",
         population="healthy volunteers", population_category="healthy_adult", n_randomized="8 in NR substudy", nr_dose_mg_day=NOT_IN_ABS,
         duration=NOT_IN_ABS, comparator=f"{U}", intervention_type="NR alone",
         outcomes_reported="NLRP3 inflammasome activation in leukocytes", outcome_families="INF",
         registration_id="NCT00442195; NCT02122575", registration_basis="stated in abstract", status="nonrandomized",
         unverified_fields="randomized; design; nr_dose_mg_day; duration; comparator", notes=""),
    dict(trial_key="CMA_kinetics", report_role="primary", pmid="32337855", author_check="",
         randomized=f"{U}: not stated", design="human calibration and control study of combined cofactors; rat toxicology",
         population="healthy subjects", population_category="healthy_adult", n_randomized=NOT_IN_ABS, nr_dose_mg_day=NOT_IN_ABS,
         duration="acute", comparator="control study without supplementation", intervention_type="combination product (CMA)",
         outcomes_reported="plasma metabolite kinetics; inflammatory proteins; clinical parameters", outcome_families="NAD; INF; OTHER",
         registration_id=f"{U}: probable NCT03838822 (kinetics of metabolic cofactors)", registration_basis="registry search",
         status="nonrandomized", unverified_fields="randomized; n_randomized; nr_dose_mg_day; registration_id", notes=""),
    dict(trial_key="Epigenetic_reanalysis", report_role="secondary analysis", pmid="42481916", author_check="",
         randomized="mixed (three source studies)", design="re-analysis of NR (5-month) and HIIT studies", population="human skeletal muscle across three studies",
         population_category="mixed", n_randomized=NOT_IN_ABS, nr_dose_mg_day=NOT_IN_ABS, duration="5 months NR",
         comparator="pre/post", intervention_type="NR alone", outcomes_reported="muscle epigenetic age acceleration (seven clocks)",
         outcome_families="OTHER (epigenetic age)", registration_id=NOT_IN_ABS, registration_basis="", status="nonrandomized",
         unverified_fields="n_randomized; nr_dose_mg_day; registration_id", notes="The 5-month NR source is probably the twin study (PMID 36638183); confirm."),
    dict(trial_key="Blood_NAD_cohorts", report_role="secondary analysis", pmid="42135539", author_check="",
         randomized=f"{U}", design="assay validation across seven cohorts, some with NR supplementation", population="mixed human cohorts",
         population_category="mixed", n_randomized=NOT_IN_ABS, nr_dose_mg_day=NOT_IN_ABS, duration=NOT_IN_ABS, comparator=f"{U}",
         intervention_type=f"{U}", outcomes_reported="whole-blood NAD+ (UHPLC-HRMS)", outcome_families="NAD", registration_id=NOT_IN_ABS,
         registration_basis="", status="nonrandomized", unverified_fields="randomized; n_randomized; nr_dose_mg_day; comparator",
         notes="Useful for assay comparability of whole-blood NAD+ across trials; check whether it re-reports any RCT."),
]

# Completed or terminated randomized registrations with actual enrollment and no PubMed
# report linked by NCT ID or by the NR search (checked 2026-11-02 UTC). Fields come from
# the registry record only.
REGISTRY_ONLY = [
    ("NCT04818216", "NIRVANA: NR in SARS-CoV-2 patients for renal protection", "disease_covid19_acute"),
    ("NCT03685253", "NR for diabetic neuropathy (terminated)", "disease_diabetic_neuropathy"),
    ("NCT04110028", "NR in hospitalized patients (sequential assignment)", "disease_hospitalized"),
    ("NCT04342975", "Protection from AKI with Basis (NR + pterostilbene)", "disease_aki"),
    ("NCT03818802", "NR effects on bone, skeletal muscle and metabolism", f"{U}"),
    ("NCT03501433", "NR effects on metabolism and vascular function (crossover)", f"{U}"),
    ("NCT03562468", "TRU NIAGEN 300 and 1000 mg/day on cognition, mood and sleep (crossover)", f"{U}"),
    ("NCT02721537", "31P-MRS brain NAD+ in current and former collegiate athletes", "healthy_adult"),
    ("NCT04907110", "NR supplementation and exercise", f"{U}"),
    ("NCT04166396", "Exercise intolerance in cystic fibrosis: resveratrol, NR, placebo (crossover)", "disease_cystic_fibrosis"),
    ("NCT05589766", "N-DOSE: NR dose optimisation in Parkinson disease", "disease_parkinson"),
    ("NCT05617508", "N-DOSE AD: NR dose optimisation in Alzheimer disease", "disease_alzheimer"),
    ("NCT05593939", "Slow Age: exercise, time-restricted feeding, NR", f"{U}"),
    ("NCT05703074", "Mental intervention and NR in long COVID (factorial)", "disease_long_covid"),
    ("NCT05561738", "NR in ulcerative colitis (terminated)", "disease_ulcerative_colitis"),
    ("NCT06018467", "Senolytics (dasatinib, quercetin) and NR to improve osteoporosis therapy", "disease_osteoporosis"),
    ("NCT07284225", "Single-dose NR on sleep", f"{U}"),
    ("NCT04573153", "CMA + hydroxychloroquine in COVID-19 (terminated)", "disease_covid19_acute"),
    ("NCT02300740", "Pharmacokinetic analysis of NR (crossover)", "healthy_adult"),
    ("NCT06382688", "IV Niagen vs NAD+ (status unknown)", f"{U}"),
    ("NCT07802548", "Sirtuin activators (NR named in record) + weight-loss diet vs placebo + diet", f"{U}"),
    ("NCT05411835", "Oral ketone salts + NR chloride vs placebo in long-chain FAOD (crossover; combination)",
     "disease_faod"),
]

COLUMNS = [
    "record_id", "trial_key", "report_role", "suggested_status", "pmid", "pmid_check",
    "first_author", "year", "title", "journal", "doi",
    "randomized", "design", "population", "population_category", "n_randomized",
    "intervention_type", "nr_dose_mg_day", "duration", "comparator",
    "outcomes_reported", "outcome_families",
    "registration_id", "registration_basis", "registration_link",
    "ctgov_overall_status", "ctgov_enrollment", "ctgov_allocation", "ctgov_model", "ctgov_has_results", "ctgov_check",
    "unverified_fields", "notes",
]


def _get_json(url: str, params: dict | None = None, tries: int = 3):
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "metaproof/scout"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode())
        except Exception as exc:  # network errors are retried, then re-raised
            last = exc
            time.sleep(1.0 + i)
    raise RuntimeError(f"request failed: {url}: {last}")


def esummary(pmids: list[str]) -> dict:
    out = {}
    for i in range(0, len(pmids), 150):
        chunk = pmids[i:i + 150]
        d = _get_json(EUTILS, {"db": "pubmed", "id": ",".join(chunk), "retmode": "json",
                                "tool": "metaproof"})["result"]
        for u in d.get("uids", []):
            out[u] = d[u]
        time.sleep(0.4)
    return out


def ctgov(nct: str) -> dict:
    s = _get_json(CTGOV + nct, {"format": "json"})
    p = s["protocolSection"]
    de = p.get("designModule", {})
    di = de.get("designInfo", {})
    return {
        "ctgov_overall_status": p.get("statusModule", {}).get("overallStatus", ""),
        "ctgov_enrollment": f"{de.get('enrollmentInfo', {}).get('count', '')} "
                            f"{de.get('enrollmentInfo', {}).get('type', '')}".strip(),
        "ctgov_allocation": di.get("allocation", ""),
        "ctgov_model": di.get("interventionModel", ""),
        "ctgov_has_results": "yes" if s.get("hasResults") else "no",
        "brief_title": p.get("identificationModule", {}).get("briefTitle", ""),
    }


def registration_link(registration_id: str, registry_row: bool = False) -> str:
    """How a record is tied to its registration: stated, probable (unverified), registry
    record, or blank when no NCT ID is given."""
    if registry_row:
        return "registry record"
    if not first_nct(registration_id):
        return ""
    return "probable (unverified)" if registration_id.strip().startswith(U) else "stated"


def first_nct(text: str) -> str | None:
    import re
    m = re.search(r"NCT\d{8}", text or "")
    return m.group(0) if m else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--offline", action="store_true", help="skip esummary and registry checks")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)
    today = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d")

    rows = [dict(r) for r in ROWS]
    pmids = [r["pmid"] for r in rows]
    summ = {} if args.offline else esummary(pmids)
    errors = []
    for i, r in enumerate(rows, 1):
        r["record_id"] = f"R{i:03d}"
        r["suggested_status"] = r.pop("status")
        s = summ.get(r["pmid"])
        if args.offline:
            r["pmid_check"] = "UNVERIFIED: offline build"
        elif not s or "error" in s:
            r["pmid_check"] = "FAILED: esummary returned no record"
            errors.append(r["pmid"])
        else:
            fa = s.get("sortfirstauthor", "")
            want = r.get("author_check", "")
            if want and not fa.lower().startswith(want.lower()):
                r["pmid_check"] = f"MISMATCH: esummary first author {fa!r}, screened as {want!r}"
                errors.append(r["pmid"])
            else:
                r["pmid_check"] = f"esummary ok {today}"
            r["first_author"] = fa
            r["year"] = (s.get("pubdate") or "")[:4]
            r["title"] = s.get("title", "")
            r["journal"] = s.get("source", "")
            r["doi"] = next((a["value"] for a in s.get("articleids", []) if a.get("idtype") == "doi"), "")
        r.pop("author_check", None)

    for nct, title, popcat in REGISTRY_ONLY:
        rows.append(dict(record_id=f"R{len(rows) + 1:03d}", trial_key=nct, report_role="registry record",
                         suggested_status="registry_only", pmid="", pmid_check="no PubMed report linked (searched 2026-11-02)",
                         title=title, population_category=popcat, registration_id=nct,
                         registration_basis="ClinicalTrials.gov API v2 search",
                         unverified_fields="all design, n, dose and outcome fields (registry only; see record)",
                         notes="Check registry results and contact investigators (Handbook chapter 13, missing results)."))

    for r in rows:
        r["registration_link"] = registration_link(r.get("registration_id", ""),
                                                   r.get("suggested_status") == "registry_only")

    if not args.offline:
        cache: dict[str, dict] = {}
        for r in rows:
            nct = first_nct(r.get("registration_id", ""))
            if not nct:
                continue
            if nct not in cache:
                try:
                    cache[nct] = ctgov(nct)
                except RuntimeError as exc:
                    cache[nct] = {"ctgov_check": f"FAILED: {exc}"}
                time.sleep(0.3)
            info = dict(cache[nct])
            info.pop("brief_title", None)
            r.update(info)
            if r["registration_link"] == "probable (unverified)":
                r.setdefault("ctgov_check", f"registry record fetched {today} ({nct}); that this "
                                            "report belongs to this trial is UNVERIFIED")
            else:
                r.setdefault("ctgov_check", f"registry ok {today} ({nct})")

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in COLUMNS})
    print(f"wrote {len(rows)} rows to {args.out}")
    if errors:
        print("PMID problems:", ", ".join(errors))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
