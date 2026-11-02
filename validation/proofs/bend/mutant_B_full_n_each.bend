# MUTANT B (must FAIL): the shared control arm is entered in full in each of
# two comparisons, so the meta-analysis counts n + n people. Claim: n + n == n.
import Base

law succ_add_r:
  for +a: Nat
  for +b: Nat
  {1n+Nat.add(a, b) == Nat.add(a, 1n+b) : Nat}

def succ_add_r(a, b):
  match a:
    case 0n:
      {==}
    case 1n+p:
      %succ_add_r(p, b) : {2n+Nat.add(p, b) == 1n+_ : Nat}
      {==}

law full_n_each_keeps_n:
  for +n: Nat
  {Nat.add(n, n) == n : Nat}

def full_n_each_keeps_n(n):
  match n:
    case 0n:
      {==}
    case 1n+p:
      %full_n_each_keeps_n(p) : {1n+Nat.add(p, 1n+p) == 1n+_ : Nat}
      {==}
