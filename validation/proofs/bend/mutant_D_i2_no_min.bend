# MUTANT D (must FAIL): I^2 slack written as max(0, Q - df) + df == Q.
# Correct form is max(0, Q - df) + min(Q, df) == Q.
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

law i2_slack_no_min:
  for +q: Nat
  for +d: Nat
  {Nat.add(Nat.sub(q, d), d) == q : Nat}

def i2_slack_no_min(q, d):
  match q d:
    case 0n 0n:
      {==}
    case 0n 1n+dp:
      {==}
    case 1n+qp 0n:
      {==}
    case 1n+qp 1n+dp:
      %succ_add_r(Nat.sub(qp, dp), dp) : {_ == 1n+qp : Nat}
      %i2_slack_no_min(qp, dp) : {1n+Nat.add(Nat.sub(qp, dp), dp) == 1n+_ : Nat}
      {==}
