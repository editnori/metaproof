# MUTANT E (must FAIL): k-way split of a shared control in which the sub-arm
# that just received a person goes to the back of the rotation without that
# person (snoc(t2, h) instead of snoc(t2, 1n+h)). Every split with two or more
# sub-arms loses people: n = 1 over two sub-arms gives [0, 0].
# The proof below copies the correct deal_keeps_n proof; Bend must reject it at
# deal_keeps_n. The evidence that the law is false is refute_mutant_E in
# meta_laws.bend, which derives Empty from this file's law.
import Base

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

law succ_add_r:
  for +a: Nat
  for +b: Nat
  {1n+Nat.add(a, b) == Nat.add(a, 1n+b) : Nat}

def succ_add_r(a, b):
  Equal.sym(Nat, Nat.add(a, 1n+b), 1n+Nat.add(a, b), add_succ_r(a, b))

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

def deal(n: Nat, h: Nat, t: List<&2, Nat>) -> List<&2, Nat>:
  match n:
    case 0n:
      Con{h, t}
    case 1n+p:
      match t:
        case Nil{}:
          deal(p, 1n+h, Nil{})
        case Con{h2, t2}:
          deal(p, h2, snoc(t2, h))

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
          %add_succ_r(p, Nat.add(h, Nat.add(h2, total(t2)))) : {total(deal(p, h2, snoc(t2, h))) == _ : Nat}
          %add_left_comm(h2, 1n+h, total(t2)) : {total(deal(p, h2, snoc(t2, h))) == Nat.add(p, _) : Nat}
          %total_snoc(t2, 1n+h) : {total(deal(p, h2, snoc(t2, h))) == Nat.add(p, Nat.add(h2, _)) : Nat}
          deal_keeps_n(p, h2, snoc(t2, h))
