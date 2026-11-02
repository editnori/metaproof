# meta_laws.bend -- discrete and structural laws for the NR meta-analysis kit.
#
# Checked by Bend2 v2.0.20. "All terms check." means every law below is proven.
# Only formulas live here: no trial data, no study values. Literals in the
# CHECK section are synthetic counts chosen to exercise the definitions.
#
# Model: every law here is over Nat (non-negative integers). Nat.sub is
# truncated subtraction, Nat.sub(a, b) = max(0, a - b), the max(0, .) step of DL
# tau^2 and of I^2 = max(0, (Q - df) / Q), but with integer Q and df and with no
# division by C or by Q. These laws are integer models, not proofs about the
# real-valued formulas or the engine: the real-valued versions are P02, P03,
# P04 and P10 in ../sympy_proofs.py. ../RECEIPTS.txt states the exact domain
# of each law.

import Base

# ---------------------------------------------------------------------------
# 1. Arithmetic support lemmas
# ---------------------------------------------------------------------------

law add_zero_r:
  for +a: Nat
  {Nat.add(a, 0n) == a : Nat}

def add_zero_r(a):
  match a:
    case 0n:
      {==}
    case 1n+p:
      %add_zero_r(p) : {1n+Nat.add(p, 0n) == 1n+_ : Nat}
      {==}

law add_succ_r:
  for +a: Nat
  for +b: Nat
  {Nat.add(a, 1n+b) == 1n+Nat.add(a, b) : Nat}

def add_succ_r(a, b):
  match a:
    case 0n:
      {==}
    case 1n+p:
      %add_succ_r(p, b) : {1n+Nat.add(p, 1n+b) == 1n+_ : Nat}
      {==}

# the same fact read right to left, for rewriting a goal that contains a + (1 + b)
law succ_add_r:
  for +a: Nat
  for +b: Nat
  {1n+Nat.add(a, b) == Nat.add(a, 1n+b) : Nat}

def succ_add_r(a, b):
  Equal.sym(Nat, Nat.add(a, 1n+b), 1n+Nat.add(a, b), add_succ_r(a, b))

# a + (x + s) == x + (a + s)
law add_left_comm:
  for +a: Nat
  for +x: Nat
  for +s: Nat
  {Nat.add(a, Nat.add(x, s)) == Nat.add(x, Nat.add(a, s)) : Nat}

def add_left_comm(a, x, s):
  match a:
    case 0n:
      {==}
    case 1n+p:
      %succ_add_r(x, Nat.add(p, s)) : {1n+Nat.add(p, Nat.add(x, s)) == _ : Nat}
      %add_left_comm(p, x, s) : {1n+Nat.add(p, Nat.add(x, s)) == 1n+_ : Nat}
      {==}

# ---------------------------------------------------------------------------
# 2. Splitting a shared control group keeps the total n: INTEGER split model
#    (Cochrane Handbook v6, Chapter 23, section 23.3.4: a shared group may be
#    split into two or more groups with smaller sample sizes.)
#    The engine does not run this split: it divides the shared control as a
#    fraction, n / k (metaproof/effects.py), whose total is P10 in sympy_proofs.py
#    and whose use in each comparison is checked by E07 in engine_conformance.py.
#    These laws would apply if an integer split were ever added.
# ---------------------------------------------------------------------------

# two comparisons share one control arm of n people: floor(n/2) and ceil(n/2)
def half(n: Nat) -> Nat:
  match n:
    case 0n:
      0n
    case 1n+p:
      match p:
        case 0n:
          0n
        case 1n+q:
          1n+half(q)

def halfup(n: Nat) -> Nat:
  match n:
    case 0n:
      0n
    case 1n+p:
      match p:
        case 0n:
          1n
        case 1n+q:
          1n+halfup(q)

law split2_keeps_n:
  for +n: Nat
  {Nat.add(half(n), halfup(n)) == n : Nat}

def split2_keeps_n(n):
  match n:
    case 0n:
      {==}
    case 1n+p:
      match p:
        case 0n:
          {==}
        case 1n+q:
          %succ_add_r(half(q), halfup(q)) : {1n+_ == 2n+q : Nat}
          %split2_keeps_n(q) : {2n+Nat.add(half(q), halfup(q)) == 2n+_ : Nat}
          {==}

# the two halves differ by at most one person: halfup(n) = half(n) + (n mod 2)
def parity(n: Nat) -> Nat:
  match n:
    case 0n:
      0n
    case 1n+p:
      match p:
        case 0n:
          1n
        case 1n+q:
          parity(q)

law split2_balanced:
  for +n: Nat
  {halfup(n) == Nat.add(half(n), parity(n)) : Nat}

def split2_balanced(n):
  match n:
    case 0n:
      {==}
    case 1n+p:
      match p:
        case 0n:
          {==}
        case 1n+q:
          %split2_balanced(q) : {1n+halfup(q) == 1n+_ : Nat}
          {==}

# k comparisons share one control arm: deal its n people one at a time to the
# k sub-arms in rotation. h is the first sub-arm's count, t the others'.
def snoc(xs: List<&2, Nat>, x: Nat) -> List<&2, Nat>:
  match xs:
    case Nil{}:
      Con{x, Nil{}}
    case Con{h, t}:
      Con{h, snoc(t, x)}

def total(xs: List<&2, Nat>) -> Nat:
  match xs:
    case Nil{}:
      0n
    case Con{h, t}:
      Nat.add(h, total(t))

def zeros(k: Nat) -> List<&2, Nat>:
  match k:
    case 0n:
      Nil{}
    case 1n+p:
      Con{0n, zeros(p)}

def deal(n: Nat, h: Nat, t: List<&2, Nat>) -> List<&2, Nat>:
  match n:
    case 0n:
      Con{h, t}
    case 1n+p:
      match t:
        case Nil{}:
          deal(p, 1n+h, Nil{})
        case Con{h2, t2}:
          deal(p, h2, snoc(t2, 1n+h))

law total_snoc:
  for +t: List<&2, Nat>
  for +x: Nat
  {total(snoc(t, x)) == Nat.add(x, total(t)) : Nat}

def total_snoc(t, x):
  match t:
    case Nil{}:
      {==}
    case Con{a, r}:
      %add_left_comm(a, x, total(r)) : {Nat.add(a, total(snoc(r, x))) == _ : Nat}
      %total_snoc(r, x) : {Nat.add(a, total(snoc(r, x))) == Nat.add(a, _) : Nat}
      {==}

law deal_keeps_n:
  for +n: Nat
  for +h: Nat
  for +t: List<&2, Nat>
  {total(deal(n, h, t)) == Nat.add(n, Nat.add(h, total(t))) : Nat}

def deal_keeps_n(n, h, t):
  match n:
    case 0n:
      {==}
    case 1n+p:
      match t:
        case Nil{}:
          %add_succ_r(p, Nat.add(h, 0n)) : {total(deal(p, 1n+h, Nil{})) == _ : Nat}
          deal_keeps_n(p, 1n+h, Nil{})
        case Con{h2, t2}:
          %add_succ_r(p, Nat.add(h, Nat.add(h2, total(t2)))) : {total(deal(p, h2, snoc(t2, 1n+h))) == _ : Nat}
          %add_left_comm(h2, 1n+h, total(t2)) : {total(deal(p, h2, snoc(t2, 1n+h))) == Nat.add(p, _) : Nat}
          %total_snoc(t2, 1n+h) : {total(deal(p, h2, snoc(t2, 1n+h))) == Nat.add(p, Nat.add(h2, _)) : Nat}
          deal_keeps_n(p, h2, snoc(t2, 1n+h))

law total_zeros:
  for +k: Nat
  {total(zeros(k)) == 0n : Nat}

def total_zeros(k):
  match k:
    case 0n:
      {==}
    case 1n+p:
      total_zeros(p)

# starting from 1 + k empty sub-arms, the sub-arms add back to the shared n
law deal_from_empty_keeps_n:
  for +n: Nat
  for +k: Nat
  {total(deal(n, 0n, zeros(k))) == n : Nat}

def deal_from_empty_keeps_n(n, k):
  %add_zero_r(n) : {total(deal(n, 0n, zeros(k))) == _ : Nat}
  %total_zeros(k) : {total(deal(n, 0n, zeros(k))) == Nat.add(n, _) : Nat}
  deal_keeps_n(n, 0n, zeros(k))

# with two sub-arms that start equal, the deal gives floor(n/2) and ceil(n/2)
law deal2_floor_ceil:
  for +n: Nat
  for +a: Nat
  {deal(n, a, Con{a, Nil{}}) == Con{Nat.add(half(n), a), Con{Nat.add(halfup(n), a), Nil{}}} : List<&2, Nat>}

def deal2_floor_ceil(n, a):
  match n:
    case 0n:
      {==}
    case 1n+p:
      match p:
        case 0n:
          {==}
        case 1n+q:
          %add_succ_r(half(q), a) : {deal(2n+q, a, Con{a, Nil{}}) == Con{_, Con{1n+Nat.add(halfup(q), a), Nil{}}} : List<&2, Nat>}
          %add_succ_r(halfup(q), a) : {deal(2n+q, a, Con{a, Nil{}}) == Con{Nat.add(half(q), 1n+a), Con{_, Nil{}}} : List<&2, Nat>}
          deal2_floor_ceil(q, 1n+a)

# so a two-way deal from empty sub-arms is (floor(n/2), ceil(n/2)); with
# split2_balanced the two sub-arms differ by at most one person. A balance law
# for three or more sub-arms is NOT proven here.
law deal2_from_empty:
  for +n: Nat
  {deal(n, 0n, zeros(1n)) == Con{half(n), Con{halfup(n), Nil{}}} : List<&2, Nat>}

def deal2_from_empty(n):
  %add_zero_r(half(n)) : {deal(n, 0n, zeros(1n)) == Con{_, Con{halfup(n), Nil{}}} : List<&2, Nat>}
  %add_zero_r(halfup(n)) : {deal(n, 0n, zeros(1n)) == Con{Nat.add(half(n), 0n), Con{_, Nil{}}} : List<&2, Nat>}
  deal2_floor_ceil(n, 0n)

# ---------------------------------------------------------------------------
# 3. Truncation in DerSimonian-Laird tau^2 and in I^2
#    tau2_DL = max(0, (Q - df) / C)   (DerSimonian & Laird 1986)
#    I2      = max(0, (Q - df) / Q)   (Higgins & Thompson 2002)
#    Both numerators are Nat.sub(Q, df) here.
# ---------------------------------------------------------------------------

# Q <= df (written df = Q + s) gives a zero numerator: for the DL estimator and the
# Q-based I^2 this means tau^2 = 0 and I^2 = 0. Integer Q and df only; REML and
# other iterative estimators can give tau^2 > 0 when Q < df (RECEIPTS.txt).
law trunc_zero_when_q_le_df:
  for +q: Nat
  for +s: Nat
  {Nat.sub(q, Nat.add(q, s)) == 0n : Nat}

def trunc_zero_when_q_le_df(q, s):
  match q:
    case 0n:
      match s:
        case 0n:
          {==}
        case 1n+sp:
          {==}
    case 1n+p:
      trunc_zero_when_q_le_df(p, s)

# Q > df (written Q = df + s) gives exactly Q - df: truncation changes nothing
law trunc_exact_above_df:
  for +d: Nat
  for +s: Nat
  {Nat.sub(Nat.add(d, s), d) == s : Nat}

def trunc_exact_above_df(d, s):
  match d:
    case 0n:
      match s:
        case 0n:
          {==}
        case 1n+sp:
          {==}
    case 1n+p:
      trunc_exact_above_df(p, s)

# max(0, Q - df) + min(Q, df) == Q, so the I^2 numerator never exceeds Q
law i2_slack:
  for +q: Nat
  for +d: Nat
  {Nat.add(Nat.sub(q, d), Nat.min(q, d)) == q : Nat}

def i2_slack(q, d):
  match q d:
    case 0n 0n:
      {==}
    case 0n 1n+dp:
      {==}
    case 1n+qp 0n:
      add_zero_r(1n+qp)
    case 1n+qp 1n+dp:
      %succ_add_r(Nat.sub(qp, dp), Nat.min(qp, dp)) : {_ == 1n+qp : Nat}
      %i2_slack(qp, dp) : {1n+Nat.add(Nat.sub(qp, dp), Nat.min(qp, dp)) == 1n+_ : Nat}
      {==}

law i2_strict_slack:
  for +q: Nat
  for +d: Nat
  {Nat.add(Nat.sub(q, d), 1n+Nat.min(q, d)) == 1n+q : Nat}

def i2_strict_slack(q, d):
  %succ_add_r(Nat.sub(q, d), Nat.min(q, d)) : {_ == 1n+q : Nat}
  %i2_slack(q, d) : {1n+Nat.add(Nat.sub(q, d), Nat.min(q, d)) == 1n+_ : Nat}
  {==}

# for INTEGER Q >= 1 and df >= 1 the numerator max(0, Q - df) is at most Q - 1,
# so (Q - df)/Q < 1. Over the reals the bound needs df >= 1 as well (Q = 1.5,
# df = 0.2 gives 1.3 > 0.5);
# the real-valued I^2 < 1 is P04 in sympy_proofs.py.
law i2_below_one:
  for +q: Nat
  for +d: Nat
  exs s: Nat
  {Nat.add(Nat.sub(1n+q, 1n+d), 1n+s) == 1n+q : Nat}

def i2_below_one(q, d):
  (Nat.min(q, d), i2_strict_slack(q, d))

# ---------------------------------------------------------------------------
# 4. The identity behind the pooled-estimate lower bound, for any number of
#    studies k (structural induction over the study list), over Nat only.
#    Each study has an integer weight w and an estimate y = m + d with integer
#    d >= 0. The law proves only
#        sum(w * y) == m * sum(w) + sum(w * d)
#    It has no division, does not make m the smallest estimate, has no upper
#    bound, and sum(w * d) >= 0 holds because the type is Nat, not by proof.
#    The real-valued bound min y <= pooled <= max y for any k, with FE and RE
#    weights, is P02 in sympy_proofs.py. The law keeps its historical name.
# ---------------------------------------------------------------------------

law add_assoc:
  for +a: Nat
  for +b: Nat
  for +c: Nat
  {Nat.add(Nat.add(a, b), c) == Nat.add(a, Nat.add(b, c)) : Nat}

def add_assoc(a, b, c):
  match a:
    case 0n:
      {==}
    case 1n+p:
      %add_assoc(p, b, c) : {1n+Nat.add(Nat.add(p, b), c) == 1n+_ : Nat}
      {==}

law add_assoc_r:
  for +a: Nat
  for +b: Nat
  for +c: Nat
  {Nat.add(a, Nat.add(b, c)) == Nat.add(Nat.add(a, b), c) : Nat}

def add_assoc_r(a, b, c):
  match a:
    case 0n:
      {==}
    case 1n+p:
      %add_assoc_r(p, b, c) : {1n+Nat.add(p, Nat.add(b, c)) == 1n+_ : Nat}
      {==}

# (a + b) + (c + d) == (a + c) + (b + d)
law add_swap4:
  for +a: Nat
  for +b: Nat
  for +c: Nat
  for +d: Nat
  {Nat.add(Nat.add(a, b), Nat.add(c, d)) == Nat.add(Nat.add(a, c), Nat.add(b, d)) : Nat}

def add_swap4(a, b, c, d):
  %add_assoc_r(a, c, Nat.add(b, d)) : {Nat.add(Nat.add(a, b), Nat.add(c, d)) == _ : Nat}
  %add_left_comm(b, c, d) : {Nat.add(Nat.add(a, b), Nat.add(c, d)) == Nat.add(a, _) : Nat}
  %add_assoc(a, b, Nat.add(c, d)) : {Nat.add(Nat.add(a, b), Nat.add(c, d)) == _ : Nat}
  {==}

# w * (m + d) == w*m + w*d
law mul_add_l:
  for +w: Nat
  for +m: Nat
  for +d: Nat
  {Nat.mul(w, Nat.add(m, d)) == Nat.add(Nat.mul(w, m), Nat.mul(w, d)) : Nat}

def mul_add_l(w, m, d):
  match w:
    case 0n:
      {==}
    case 1n+p:
      %add_swap4(m, d, Nat.mul(p, m), Nat.mul(p, d)) : {Nat.add(Nat.add(m, d), Nat.mul(p, Nat.add(m, d))) == _ : Nat}
      %mul_add_l(p, m, d) : {Nat.add(Nat.add(m, d), Nat.mul(p, Nat.add(m, d))) == Nat.add(Nat.add(m, d), _) : Nat}
      {==}

# a*m + b*m == (a + b) * m
law mul_add_r:
  for +a: Nat
  for +b: Nat
  for +m: Nat
  {Nat.add(Nat.mul(a, m), Nat.mul(b, m)) == Nat.mul(Nat.add(a, b), m) : Nat}

def mul_add_r(a, b, m):
  match a:
    case 0n:
      {==}
    case 1n+p:
      %mul_add_r(p, b, m) : {Nat.add(Nat.add(m, Nat.mul(p, m)), Nat.mul(b, m)) == Nat.add(m, _) : Nat}
      %add_assoc(m, Nat.mul(p, m), Nat.mul(b, m)) : {Nat.add(Nat.add(m, Nat.mul(p, m)), Nat.mul(b, m)) == _ : Nat}
      {==}

type Study is Data:
  Study{w: Nat, d: Nat}

def sum_w(xs: List<&2, Study>) -> Nat:
  match xs:
    case Nil{}:
      0n
    case Con{Study{w, d}, t}:
      Nat.add(w, sum_w(t))

def sum_wd(xs: List<&2, Study>) -> Nat:
  match xs:
    case Nil{}:
      0n
    case Con{Study{w, d}, t}:
      Nat.add(Nat.mul(w, d), sum_wd(t))

def sum_wy(+m: Nat, xs: List<&2, Study>) -> Nat:
  match xs:
    case Nil{}:
      0n
    case Con{Study{w, d}, t}:
      Nat.add(Nat.mul(w, Nat.add(m, d)), sum_wy(m, t))

law pooled_not_below_min:
  for +m: Nat
  for +xs: List<&2, Study>
  {sum_wy(m, xs) == Nat.add(Nat.mul(sum_w(xs), m), sum_wd(xs)) : Nat}

def pooled_not_below_min(m, xs):
  match xs:
    case Nil{}:
      {==}
    case Con{Study{w, d}, t}:
      %mul_add_r(w, sum_w(t), m) : {sum_wy(m, xs) == Nat.add(_, Nat.add(Nat.mul(w, d), sum_wd(t))) : Nat}
      %add_swap4(Nat.mul(w, m), Nat.mul(w, d), Nat.mul(sum_w(t), m), sum_wd(t)) : {sum_wy(m, xs) == _ : Nat}
      %mul_add_l(w, m, d) : {sum_wy(m, xs) == Nat.add(_, Nat.add(Nat.mul(sum_w(t), m), sum_wd(t))) : Nat}
      %pooled_not_below_min(m, t) : {sum_wy(m, xs) == Nat.add(Nat.mul(w, Nat.add(m, d)), _) : Nat}
      {==}

# ---------------------------------------------------------------------------
# 5. Refutations of the mutant formulas (see mutant_*.bend)
#    A rejected copy of a proof is not evidence that a statement is false: Bend
#    also rejects a wrong proof of a true statement. The evidence is here.
#    refute_<case> shows the mutant's equation false at one synthetic count;
#    refute_mutant_<X> takes the mutant's law, quantified exactly as
#    mutant_<X>.bend states it, and derives Empty, so that law is false.
#    run_bend.py checks that each hypothesis matches the mutant file's law.
# ---------------------------------------------------------------------------

def disc(n: Nat) -> Type:
  match n:
    case 0n:
      Empty
    case 1n+p:
      Unit

law succ_ne_zero:
  for -p: Nat
  for e: {1n+p == 0n : Nat}
  Empty

def succ_ne_zero(p, e):
  %e : disc(_)
  Unit{}

# mutant A: each of two comparisons gets floor(n/2) of the shared control.
# synthetic n = 3: floor + floor = 2, so the claim "== 3" is contradictory.
law refute_floor_twice:
  for e: {Nat.add(half(3n), half(3n)) == 3n : Nat}
  Empty

def refute_floor_twice(e):
  succ_ne_zero(0n, Equal.cong(Nat, Nat, x => Nat.sub(3n, x), Nat.add(half(3n), half(3n)), 3n, e))

# mutant B: the shared control is counted in full in both comparisons.
# synthetic n = 1: 1 + 1 == 1 is contradictory.
law refute_full_n_twice:
  for e: {Nat.add(1n, 1n) == 1n : Nat}
  Empty

def refute_full_n_twice(e):
  succ_ne_zero(0n, Equal.cong(Nat, Nat, x => Nat.sub(x, 1n), Nat.add(1n, 1n), 1n, e))

# mutant C: numerator max(0, df - Q) instead of max(0, Q - df).
# synthetic Q = 5, df = 2: the correct numerator is 3, the mutant's is 0.
law refute_swapped_numerator:
  for e: {Nat.sub(2n, 5n) == Nat.sub(5n, 2n) : Nat}
  Empty

def refute_swapped_numerator(e):
  succ_ne_zero(2n, Equal.sym(Nat, Nat.sub(2n, 5n), Nat.sub(5n, 2n), e))

# mutant D: I^2 slack written as max(0, Q - df) + df == Q (min dropped).
# synthetic Q = 1, df = 2: 0 + 2 == 1 is contradictory.
law refute_i2_no_min:
  for e: {Nat.add(Nat.sub(1n, 2n), 2n) == 1n : Nat}
  Empty

def refute_i2_no_min(e):
  succ_ne_zero(0n, Equal.cong(Nat, Nat, x => Nat.sub(x, 1n), Nat.add(Nat.sub(1n, 2n), 2n), 1n, e))

# mutant E: a k-way deal that drops the person on rotation (the 1n+ is lost
# when the dealt-to sub-arm moves to the back). Same body as mutant_E's deal.
def deal_drop(n: Nat, h: Nat, t: List<&2, Nat>) -> List<&2, Nat>:
  match n:
    case 0n:
      Con{h, t}
    case 1n+p:
      match t:
        case Nil{}:
          deal_drop(p, 1n+h, Nil{})
        case Con{h2, t2}:
          deal_drop(p, h2, snoc(t2, h))

# synthetic n = 1 split over two sub-arms: the mutant deals [0, 0], total 0.
law refute_deal_drops:
  for e: {total(deal_drop(1n, 0n, zeros(1n))) == 1n : Nat}
  Empty

def refute_deal_drops(e):
  succ_ne_zero(0n, Equal.sym(Nat, total(deal_drop(1n, 0n, zeros(1n))), 1n, e))

law refute_mutant_A:
  for claim: @n:Nat -> {Nat.add(half(n), half(n)) == n : Nat}
  Empty

def refute_mutant_A(claim):
  refute_floor_twice(claim(3n))

law refute_mutant_B:
  for claim: @n:Nat -> {Nat.add(n, n) == n : Nat}
  Empty

def refute_mutant_B(claim):
  refute_full_n_twice(claim(1n))

law refute_mutant_C:
  for claim: @d:Nat -> @s:Nat -> {Nat.sub(d, Nat.add(d, s)) == s : Nat}
  Empty

def refute_mutant_C(claim):
  refute_swapped_numerator(claim(2n, 3n))

law refute_mutant_D:
  for claim: @q:Nat -> @d:Nat -> {Nat.add(Nat.sub(q, d), d) == q : Nat}
  Empty

def refute_mutant_D(claim):
  refute_i2_no_min(claim(1n, 2n))

law refute_mutant_E:
  for claim: @n:Nat -> @h:Nat -> @t:List<&2, Nat> -> {total(deal_drop(n, h, t)) == Nat.add(n, Nat.add(h, total(t))) : Nat}
  Empty

def refute_mutant_E(claim):
  refute_deal_drops(claim(1n, 0n, zeros(1n)))

# ---------------------------------------------------------------------------
# 6. Evaluation checks on synthetic counts (SYNTHETIC values, not trial data)
# ---------------------------------------------------------------------------

law synthetic_split2_n7:
  {Nat.add(half(7n), halfup(7n)) == 7n : Nat}

def synthetic_split2_n7():
  {==}

law synthetic_halves_n7:
  {half(7n) == 3n : Nat}

def synthetic_halves_n7():
  {==}

law synthetic_deal_n7_k3:
  {total(deal(7n, 0n, zeros(2n))) == 7n : Nat}

def synthetic_deal_n7_k3():
  {==}

law synthetic_pooled_lower:
  {sum_wy(2n, Con{Study{3n, 0n}, Con{Study{1n, 4n}, Nil{}}}) == 12n : Nat}

def synthetic_pooled_lower():
  {==}
