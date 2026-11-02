# MUTANT C (must FAIL): DL tau^2 / I^2 numerator written as max(0, df - Q).
# Claim checked: when Q = df + s the numerator equals s (true for the correct
# max(0, Q - df), false for the mutant).
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

law swapped_exact_above_df:
  for +d: Nat
  for +s: Nat
  {Nat.sub(d, Nat.add(d, s)) == s : Nat}

def swapped_exact_above_df(d, s):
  match d:
    case 0n:
      match s:
        case 0n:
          {==}
        case 1n+sp:
          {==}
    case 1n+p:
      swapped_exact_above_df(p, s)
