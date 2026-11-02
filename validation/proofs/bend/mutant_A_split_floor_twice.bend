# MUTANT A (must FAIL): two comparisons share a control arm of n people and
# each gets floor(n/2). Odd n loses a person. The proof below copies the
# correct split2 proof; Bend must reject it.
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

law split_floor_twice_keeps_n:
  for +n: Nat
  {Nat.add(half(n), half(n)) == n : Nat}

def split_floor_twice_keeps_n(n):
  match n:
    case 0n:
      {==}
    case 1n+p:
      match p:
        case 0n:
          {==}
        case 1n+q:
          %succ_add_r(half(q), half(q)) : {1n+_ == 2n+q : Nat}
          %split_floor_twice_keeps_n(q) : {2n+Nat.add(half(q), half(q)) == 2n+_ : Nat}
          {==}
