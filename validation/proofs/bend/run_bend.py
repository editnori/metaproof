"""Check the Bend2 laws on a Linux machine over ssh and write receipt_bend.txt beside this file.

Run from the repo root (or through run_bend_on_vm.sh):
    .venv/bin/python validation/proofs/bend/run_bend.py

What passes
-----------
meta_laws.bend must print "All terms check." with exit status 0. That checks every law in
it, including the refute_mutant_<X> laws.

A mutant counts as rejected only when all of these hold:
  1. static tie (checked here, no VM): meta_laws.bend has law refute_mutant_<X> whose
     hypothesis is the mutant file's law, quantified the same way, with the mutant's
     functions renamed as MUTANTS says; and every function the mutant file defines has the
     same body in meta_laws.bend (under that rename). So the checked refutation is about
     the mutant's own statement, not a look-alike.
  2. Bend exits 1 with a type error ("- expected :" / "- observed :") whose Location is the
     mutant's law. A parse error, an unknown name, a bad import, a missing file, an ssh or
     wrapper failure, or a type error in any other def is NOT a rejection: it is reported as
     ERROR and the run fails.
The rejected copy of a proof only illustrates the failure. Bend also rejects a wrong proof
of a true statement, so the evidence that a mutant is false is its refute_mutant_<X> law.

Exit status: 0 all pass, 1 a law or mutant check failed, 2 infrastructure failure (scp,
ssh, wrapper, Bend not run). Only formulas and synthetic counts are sent to the VM.
"""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

HERE = Path(__file__).resolve().parent
# Where Bend2 runs. Bend2 needs Linux with glibc, so the laws are checked over ssh on a
# Linux machine whose BEND_DIR holds an executable wrapper bend.sh that runs `bend`.
VM = os.environ.get("BEND_HOST", "bend-host")
VMDIR = os.environ.get("BEND_DIR", "~/metaproof-bend")
VM_HOME = os.environ.get("BEND_HOME", "~")
LAWS = "meta_laws.bend"
PASS_LINE = "All terms check."
EXIT_MARK = "__BEND_EXIT__="


@dataclass(frozen=True)
class MutantSpec:
    law: str                     # the mutant's law; Bend must reject its proof here
    refute: str                  # law in meta_laws.bend that derives Empty from it
    renames: dict = field(default_factory=dict)   # mutant def name -> meta_laws def name


MUTANTS = {
    "mutant_A_split_floor_twice.bend": MutantSpec("split_floor_twice_keeps_n", "refute_mutant_A"),
    "mutant_B_full_n_each.bend": MutantSpec("full_n_each_keeps_n", "refute_mutant_B"),
    "mutant_C_swapped_numerator.bend": MutantSpec("swapped_exact_above_df", "refute_mutant_C"),
    "mutant_D_i2_no_min.bend": MutantSpec("i2_slack_no_min", "refute_mutant_D"),
    "mutant_E_deal_drops_on_rotation.bend": MutantSpec("deal_keeps_n", "refute_mutant_E",
                                                       {"deal": "deal_drop"}),
}

# "- expected :" texts that are not type errors in a proof
NOT_TYPE_ERRORS = ("a defined name", "an import", "end of input", "'")


# --- parsing the .bend sources -------------------------------------------------------------

def _blocks(text: str) -> dict[tuple[str, str], list[str]]:
    """Top-level 'law NAME:' and 'def NAME(...)' blocks: (kind, name) -> lines."""
    out, key = {}, None
    for line in text.splitlines():
        m = re.match(r"(law|def|type) (\w+)", line)
        if m:
            key = (m.group(1), m.group(2))
            out[key] = [line]
        elif line and not line[0].isspace():
            key = None
        elif key is not None:
            out[key].append(line)
    return {k: [x for x in v if x.strip() and not x.strip().startswith("#")] for k, v in out.items()}


def _rename(text: str, renames: dict) -> str:
    for a, b in renames.items():
        text = re.sub(rf"\b{re.escape(a)}\b", b, text)
    return text


def law_as_hypothesis(law_lines: list[str]) -> str:
    """'law f:\\n  for +d: Nat\\n  for +s: Nat\\n  {P}' -> '@d:Nat -> @s:Nat -> {P}'."""
    binders, claim = [], None
    for line in law_lines[1:]:
        m = re.match(r"\s*for [+-]?(\w+): (.+)$", line)
        if m:
            binders.append(f"@{m.group(1)}:{m.group(2).strip()}")
        elif re.match(r"\s*exs ", line):
            raise ValueError("a law with a witness (exs) cannot be a mutant hypothesis")
        else:
            claim = line.strip()
    if claim is None:
        raise ValueError(f"{law_lines[0]}: no claim line")
    return " -> ".join(binders + [claim])


def static_ties(here: Path = HERE) -> list[str]:
    """Problems that make a mutant's rejection worthless as evidence; [] when all tie."""
    problems = []
    laws = _blocks((here / LAWS).read_text())
    files = sorted(p.name for p in here.glob("mutant_*.bend"))
    for f in files:
        if f not in MUTANTS:
            problems.append(f"{f}: not listed in run_bend.MUTANTS (no refute law named)")
    for f, spec in MUTANTS.items():
        if f not in files:
            problems.append(f"{f}: listed in MUTANTS but missing")
            continue
        mut = _blocks((here / f).read_text())
        if ("law", spec.law) not in mut:
            problems.append(f"{f}: has no law {spec.law}")
            continue
        refute = laws.get(("law", spec.refute))
        if refute is None:
            problems.append(f"{f}: {LAWS} has no law {spec.refute}")
            continue
        want = _rename(law_as_hypothesis(mut[("law", spec.law)]), spec.renames)
        hyp = [m.group(1).strip() for line in refute[1:]
               if (m := re.match(r"\s*for \w+: (.+)$", line))]
        if refute[-1].strip() != "Empty" or hyp != [want]:
            problems.append(f"{f}: {spec.refute} must be 'for claim: {want}' -> Empty; "
                            f"found {hyp} -> {refute[-1].strip()}")
        for (kind, name), lines in mut.items():
            if kind != "def" or "->" not in lines[0]:
                continue                      # proof defs (no return type) are the copies
            target = spec.renames.get(name, name)
            theirs = laws.get(("def", target))
            if theirs is None or _rename("\n".join(lines), spec.renames) != "\n".join(theirs):
                problems.append(f"{f}: def {name} differs from {LAWS} def {target}, so "
                                f"{spec.refute} would refute a different function")
    return problems


# --- classifying checker output -------------------------------------------------------------

def classify(out: str, status: int | None, law: str | None) -> tuple[str, str]:
    """(verdict, detail). law None means the file must check.

    verdicts: CHECKS, REJECTED (type error at law), ACCEPTED (mutant checked), ERROR."""
    lines = [x.rstrip() for x in out.strip().splitlines()]
    if status is None:
        return "ERROR", "Bend did not run (ssh, wrapper or marker missing)"
    if status == 0 and lines and lines[-1] == PASS_LINE:
        return ("CHECKS", PASS_LINE) if law is None else ("ACCEPTED", "mutant law checked")
    if status != 1 or not lines or lines[0] != "Error:":
        return "ERROR", f"exit {status}, not a checker verdict: {(lines or [''])[0][:80]}"
    if any(x.startswith("- message") for x in lines):
        return "ERROR", "checker message (not a type error): " + \
            next(x for x in lines if x.startswith("- message"))
    exp = next((x for x in lines if x.startswith("- expected :")), None)
    obs = next((x for x in lines if x.startswith("- observed :")), None)
    loc = next((x for x in lines if x.startswith("Location:")), None)
    if exp is None or obs is None or loc is None:
        return "ERROR", "checker error without expected/observed/Location"
    if exp.split(":", 1)[1].strip().startswith(NOT_TYPE_ERRORS):
        return "ERROR", f"not a type error: {exp}"
    where = loc.split(":", 1)[1].strip()
    if law is None:
        return "ERROR", f"type error at {where}"
    if where != law:
        return "ERROR", f"type error at {where}, not at the mutant law {law}"
    return "REJECTED", f"type error at {law}"


# --- running on the VM -----------------------------------------------------------------------

Runner = Callable[[list[str]], tuple[int, str]]


SSH_OPTS = ["-o", "ConnectTimeout=20", "-o", "BatchMode=yes"]


def _sh(cmd: list[str]) -> tuple[int, str]:
    if cmd[0] in ("ssh", "scp"):
        cmd = [cmd[0], *SSH_OPTS, *cmd[1:]]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    except subprocess.TimeoutExpired:
        return 124, f"timed out: {' '.join(cmd[:2])}"
    return p.returncode, p.stdout + p.stderr


def check_file(name: str, run: Runner, vm: str, vmdir: str) -> tuple[int | None, str]:
    code, out = run(["ssh", vm, f"cd {vmdir} && HOME={VM_HOME} ./bend.sh {name} --check-only "
                     f"2>&1; echo {EXIT_MARK}$?"])
    m = re.search(rf"^{EXIT_MARK}(\d+)$", out, re.M)
    if code != 0 or m is None:
        return None, out
    return int(m.group(1)), out[:m.start()]


def main(run: Runner = _sh, vm: str = VM, vmdir: str = VMDIR, here: Path = HERE) -> int:
    files = [LAWS] + sorted(p.name for p in here.glob("mutant_*.bend"))
    rows, logs, infra = [], [], False
    ties = static_ties(here)
    code, _ = run(["ssh", vm, f"cd {vmdir} && test -x bend.sh && rm -f meta_laws.bend mutant_*.bend"])
    if code == 0:
        code, _ = run(["scp", "-q", *[str(here / f) for f in files], f"{vm}:{vmdir}/"])
    if code != 0:
        infra = True
        version = "not reached"
    else:
        code, version = run(["ssh", vm, f"cd {vmdir} && HOME={VM_HOME} ./bend.sh version 2>&1"])
        version = version.strip() if code == 0 else "not reached"
        infra = code != 0
    for f in files:
        sha = hashlib.sha256((here / f).read_bytes()).hexdigest()[:16]
        law = None if f == LAWS else MUTANTS.get(f, MutantSpec("?", "?")).law
        if infra:
            verdict, detail, out = "ERROR", "not run: VM unreachable", ""
        else:
            status, out = check_file(f, run, vm, vmdir)
            verdict, detail = classify(out, status, law)
            infra = infra or status is None
        expected = "CHECKS" if law is None else f"REJECTED at {law}"
        tie_ok = law is None or not [p for p in ties if p.startswith(f + ":")]
        result = "PASS" if verdict == expected.split(" ")[0] and tie_ok else "FAIL"
        refute = "" if law is None else f" | refuted by {MUTANTS[f].refute}" if f in MUTANTS else ""
        rows.append(f"{f} | {sha} | {expected} | {verdict}: {detail}{refute} | {result}")
        logs += ["", f"--- {f} output (first 12 lines)", *out.strip().splitlines()[:12]]
    ok = all(r.endswith("| PASS") for r in rows) and not ties
    status = 0 if ok else (2 if infra else 1)
    text = "\n".join([
        "Bend2 receipt for validation/proofs/bend (written by run_bend.py)",
        f"date_utc: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"host: {vm}:{vmdir}", f"bend: {version}", "",
        "A mutant passes only if Bend reports a type error at the mutant's own law AND",
        "meta_laws.bend checks a refute_mutant_<X> law whose hypothesis is that law.",
        "", "file | sha256 | expected | observed | result", *rows, "",
        "static ties (mutant law == refute hypothesis; mutant defs == meta_laws defs):",
        *(f"  PROBLEM {p}" for p in ties), *([] if ties else ["  all tie"]), "",
        "laws proven in meta_laws.bend:",
        *(f"  {k[1]}" for k in _blocks((here / LAWS).read_text()) if k[0] == "law"),
        "", "checker output:", *logs, "",
        f"overall: {'PASS' if ok else 'FAIL (infrastructure)' if infra else 'FAIL'}",
    ]) + "\n"
    (here / "receipt_bend.txt").write_text(text)
    print(text)
    return status


if __name__ == "__main__":
    import os
    raise SystemExit(main(vm=VM, vmdir=VMDIR))
