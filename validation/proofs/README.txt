validation/proofs -- formal checks on the meta-analysis formulas

What this folder proves
  The formulas the kit uses must have certain properties no matter what the
  data are: weights are positive, pooled estimates sit between the study
  estimates, tau^2 and I^2 respect their bounds, and so on. This folder
  proves those properties and keeps wrong versions (mutants) that must fail,
  so the checks are shown to tell a right formula from a wrong one.
  RECEIPTS.txt states the exact domain of every claim (real or integer,
  general k or k = 2..8, which estimator); quote it rather than this summary.

  sympy_proofs.py   real-valued proofs (sympy 1.14). Ten properties P01-P10,
                    15 mutants M01-M15. Writes receipt_sympy.txt.
  engine_conformance.py
                    ties the proofs to the code: each formula is built in sympy,
                    evaluated on random SYNTHETIC inputs and compared with the
                    metaproof function (E01-E06), plus engine properties no formula
                    proof covers (E07-E08): control participants counted once
                    on the arm-level paths the engine accepts (split, combine,
                    several control rows, explicit combine_arms; a
                    multi-period crossover must be excluded or count the
                    placebo period once), each comparison's yi and vi equal to
                    an independent formula at the split n, each row's dose
                    from its own experimental arm(s), one timepoint per study.
                    These are sampled SYNTHETIC checks, not proofs.
                    Writes receipt_engine.txt.
  test_proofs.py    pytest wrapper over the proofs, mutants and conformance
                    checks, plus tests that each check fails on a broken engine
                    or formula and tests of the Bend runner's rules.
  bend/             Bend2 v2.0.20 laws over Nat (non-negative integers):
                    integer models of the truncation step, the integer split,
                    and the sum identity behind the pooled lower bound. They
                    are not proofs of the real-valued formulas or of the
                    engine. Checked on a Linux machine by bend/run_bend.py, which
                    writes bend/receipt_bend.txt.
  RECEIPTS.txt      what each proof shows, and its domain.
  receipt.txt       the combined table.

How to run
  .venv/bin/python validation/proofs/sympy_proofs.py
  .venv/bin/python validation/proofs/engine_conformance.py
  .venv/bin/python -m pytest validation/proofs -q
  validation/proofs/bend/run_bend_on_vm.sh      (runs bend/run_bend.py)

  Bend2 needs Linux with glibc, so the runner copies only the .bend files
  (formulas and synthetic counts) over ssh to a Linux machine and runs
  `./bend.sh FILE.bend --check-only` in a folder there. Set BEND_HOST (the ssh
  host), BEND_DIR (the folder holding an executable bend.sh that runs Bend2
  v2.0.20, https://github.com/bendlang/bend) and BEND_HOME. Exit status 0
  is a pass, 1 a failed law or mutant, 2 an infrastructure failure.

What a PASS means
  The sympy and Bend receipts are formula proofs: they show that the formulas
  as written in this folder have the stated properties. They say nothing about
  the engine by themselves. engine_conformance.py is the separate test that the
  engine computes those formulas. P07 proves the approximate Hedges J while the
  engine uses the exact gamma form (E03), and P04 proves the Q-based I^2 while
  random-effects rows report tau2/(tau2 + s^2) (E02). P10 shows only that a
  split keeps the total n; split comparisons stay correlated (Handbook v6.5
  23.3.4 and 6.2.9), which is why the engine default is to combine arms.

What counts as proven
  sympy: a step passes only when sympy decides it (an identity simplifies to
  0, or the assumption system returns True for the sign). "Cannot decide"
  is a failure. Steps marked NUMERIC are scans, not proofs: the t and normal
  quantile order (P06) and the exact gamma form of Hedges' J (P07).
  Steps marked k=2..8 are proven for each number of studies from 2 to 8.
  Steps marked "general k" use sympy Sum over indexed symbols with sign
  assumptions and hold for every k.

  Bend2: each law is a type and its def is a checked proof term over all
  natural numbers and all list lengths, and nothing wider. Every law is over
  Nat: Nat.sub is truncated subtraction, Q and df are integers and nothing
  is divided. pooled_not_below_min proves only the identity
  sum w(m + d) = m sum w + sum w d over naturals. The split laws model an
  integer floor/ceil and round-robin split that the engine does not run (it
  divides the shared control as n/k: P10 and E07). The real-valued
  statements are P02, P03, P04 and P10.

Mutants
  sympy: a mutant is refuted when it cannot be proven and a SYNTHETIC
  counterexample is found. Every counterexample value is synthetic and is
  not trial data.

  Bend: a rejected copy of a proof is not evidence, because Bend also
  rejects a wrong proof of a true statement. A Bend mutant is refuted only
  when (1) meta_laws.bend checks a law refute_mutant_<X> whose hypothesis is
  the mutant file's own law (run_bend.py compares the text, and compares
  every function the mutant defines with meta_laws.bend), and (2) Bend
  reports a type error at the mutant's law. A parse error, an unknown name,
  a missing file, a type error in another def, or an ssh or wrapper failure
  is an ERROR and fails the run; it is never counted as a rejection.

Edges the proofs surfaced (the kit should guard these)
  - Hedges J = 1 - 3/(4 df - 1) is 0 at df = 1; J is in (0,1) only for df > 1.
  - Change-score SD is 0 when r = 1 and SD_baseline = SD_final; the formula is
    a valid SD only for r in [-1, 1].
  - DL C > 0 needs k >= 2 studies.
  - The prediction interval uses t with k - 2 df (pi_type riley), so it needs
    k >= 3; pi_type revman uses k - 1. Below 5 studies the engine flags any PI
    as not interpretable (Handbook v6.5 10.10.4.3).
  - PI-contains-CI (P06) holds when the PI and the CI use the same SE and
    the PI's critical value is at least the CI's. The engine's PI uses the
    SE of its CI (HK-scaled under HKSJ) and a critical value no smaller
    (riley t_(k-2) or revman t_(k-1), against the CI's t_(k-1) or z; metafor
    uses the CI's own), so the engine PI always contains its CI. E06 checks
    this for DL and REML, with HKSJ off and under every hksj_rule.
  - "Q <= df gives tau^2 = 0" is a DerSimonian-Laird property (P03). REML
    can give tau^2 > 0 with Q < df: SYNTHETIC y = (0, 1, 2),
    v = (0.5, 0.1, 10) gives Q = 1.80 < 2, REML tau^2 = 0.149 and the
    reported I^2 = tau2/(tau2 + s^2) = 20.7%.

Sources
  Cochrane Handbook for Systematic Reviews of Interventions v6: Chapter 10
  (section 10.3 inverse-variance; section 10.10 heterogeneity, I^2,
  random effects, prediction intervals); Chapter 6 section 6.5.2.8 (SD of
  change from baseline); Chapter 23 section 23.2 (cross-over trials) and
  section 23.3.4 (multiple groups, splitting a shared group).
  DerSimonian R, Laird N. Controlled Clinical Trials 1986;7:177-188.
  Higgins JPT, Thompson SG. Statistics in Medicine 2002;21:1539-1558.
  Higgins JPT, Thompson SG, Spiegelhalter DJ. J R Stat Soc A 2009;172:137-159.
  Hedges LV. J Educational Statistics 1981;6:107-128.
  Elbourne DR et al. Int J Epidemiol 2002;31:140-149.
