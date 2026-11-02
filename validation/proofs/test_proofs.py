"""pytest wrapper: every sympy proof must pass, every mutant must be refuted, and the metaproof
engine must compute the proven formulas (engine_conformance.py).

Run: .venv/bin/python -m pytest validation/proofs -q
The Bend2 laws run on a Linux machine over ssh through bend/run_bend_on_vm.sh. The rules that decide what
counts as a Bend rejection (bend/run_bend.py) are tested here without the VM.
"""

import shutil
import sys
from pathlib import Path

import pytest
import sympy as sp

import sympy_proofs as sm

sys.path.insert(0, str(Path(__file__).resolve().parent / "bend"))
import run_bend as rb  # noqa: E402


@pytest.mark.parametrize("proof", sm.PROOFS, ids=lambda f: f.__name__)
def test_proof_holds(proof):
    r = proof()
    failed = [text for text, ok in r.steps if not ok]
    assert r.ok, f"{r.pid} failed steps: {failed}"


@pytest.mark.parametrize("make", sm.MUTANTS, ids=lambda f: f.__name__)
def test_mutant_refuted(make):
    m = make()
    assert not m.proven, f"{m.mid}: sympy proved the wrong formula ({m.wrong})"
    assert m.counterexample, f"{m.mid}: no counterexample found for {m.wrong}"


import engine_conformance as ec  # noqa: E402


@pytest.mark.parametrize("check", ec.CHECKS, ids=lambda f: f.__name__)
def test_engine_conforms(check):
    c = check()
    assert c.ok, f"{c.cid} {c.prop}: max scaled diff {c.max_rel:.3g}; {c.detail}"


# --- the proofs must fail when the formula they claim to cover is wrong ----------------------

def test_p04_fails_for_wrong_i2(monkeypatch):
    """CR6-7: P04's 'I2 < 1' was the constant 0 < 1, so P04 passed with M06's I2."""
    monkeypatch.setattr(sm, "i2", lambda Q, df: sp.Max(0, (Q - df) / df))
    assert not sm.p04_i2().ok


def test_general_k_between_fails_without_division():
    """The general-k P02 steps must reject pooled = sum(w y) (M03's formula)."""
    out = sm.general_k_between(pool=lambda w, y, i, k: sp.Sum(w * y, (i, 1, k)))
    assert not any(lo or hi for _, lo, hi in out)
    assert all(lo and hi for _, lo, hi in sm.general_k_between())


# --- engine conformance must fail when the engine is wrong ------------------------------------

def test_e07_catches_wrong_dose(monkeypatch):
    """CR6-5: E07 passed when every dose was replaced by an invented positive number."""
    from metaproof import effects as ef
    orig = ef.moderator_values

    def wrong_dose(rows, members, carried):
        out = orig(rows, members, carried)
        out["dose_min"] = out["dose_max"] = 12345.0
        return out
    monkeypatch.setattr(ef, "moderator_values", wrong_dose)
    assert not ec.e07_control_counted_once().ok


def test_e07_catches_variance_from_unsplit_control(monkeypatch):
    """CR6-5: E07 passed when vi used the full shared control while n_c said split."""
    from metaproof import effects as ef
    orig = ef._parallel_effect

    def full_control(e, c, s, notes, basis):
        c2 = ef.Arm(**{**c.__dict__, "n": c.n * 4})
        return orig(e, c2, s, notes, basis)
    monkeypatch.setattr(ef, "_parallel_effect", full_control)
    assert not ec.e07_control_counted_once().ok


def test_e06_catches_metafor_pi_using_z_under_hksj(monkeypatch):
    """CR6-8: E06 never checked pi_type metafor or fits without HKSJ."""
    from scipy import stats
    from metaproof import pooling as pl
    orig = pl.prediction_interval

    def z_always(f, level=0.95, pi_type="riley"):
        if pi_type != "metafor" or f.method == "FE":
            return orig(f, level, pi_type)
        half = stats.norm.ppf(1 - (1 - level) / 2) * (f.tau2 + f.se[0] ** 2) ** 0.5
        return f.estimate - half, f.estimate + half
    monkeypatch.setattr(pl, "prediction_interval", z_always)
    assert not ec.e06_prediction_interval().ok


# --- Bend runner: only a type error at the mutant's own law is a rejection -------------------

TYPE_ERROR = """Error:
- expected : 0n
- observed : 1n
Location: split_floor_twice_keeps_n
40 |         case 0n:"""
LAW = "split_floor_twice_keeps_n"


@pytest.mark.parametrize("out,status,verdict", [
    (TYPE_ERROR, 1, "REJECTED"),
    (TYPE_ERROR.replace(LAW, "succ_add_r"), 1, "ERROR"),                 # wrong def failed
    ("Error:\n- message  : no such file: mutant_A.bend", 1, "ERROR"),    # missing file
    ("Error:\n- expected : '}'\n- observed : end of input\nLocation:\n3 |", 1, "ERROR"),  # parse
    ("Error:\n- expected : a defined name\n- observed : bar\nLocation: " + LAW, 1, "ERROR"),
    ("Error:\n- expected : an import ('import Base')\n- observed : 'import Nope'\nLocation:", 1,
     "ERROR"),
    ("ssh: connect to host bend-host port 22: Operation timed out", None, "ERROR"),
    ("bwrap: execvp bend: No such file or directory", 127, "ERROR"),
    (TYPE_ERROR, 0, "ERROR"),                                             # exit status disagrees
    ("All terms check.", 0, "ACCEPTED"),
])
def test_classify_mutant_output(out, status, verdict):
    assert rb.classify(out, status, LAW)[0] == verdict


def test_classify_laws_file():
    assert rb.classify("All terms check.", 0, None)[0] == "CHECKS"
    assert rb.classify(TYPE_ERROR, 1, None)[0] == "ERROR"


def test_static_ties_hold_for_shipped_files():
    assert rb.static_ties() == []
    assert set(rb.MUTANTS) == {p.name for p in rb.HERE.glob("mutant_*.bend")}


def _copy_bend(tmp_path):
    for p in rb.HERE.glob("*.bend"):
        shutil.copy(p, tmp_path / p.name)
    return tmp_path


def test_static_ties_catch_missing_refutation(tmp_path):
    """CR6-4: a mutant whose falsity no checked law shows must not count as refuted."""
    d = _copy_bend(tmp_path)
    laws = (d / rb.LAWS).read_text().replace("law refute_mutant_E:", "law refute_other_E:")
    (d / rb.LAWS).write_text(laws.replace("def refute_mutant_E(", "def refute_other_E("))
    assert any("refute_mutant_E" in p for p in rb.static_ties(d))


def test_static_ties_catch_refuting_a_different_statement(tmp_path):
    d = _copy_bend(tmp_path)
    laws = (d / rb.LAWS).read_text().replace(
        "for claim: @n:Nat -> {Nat.add(half(n), half(n)) == n : Nat}",
        "for claim: @n:Nat -> {Nat.add(half(n), halfup(n)) == n : Nat}")
    (d / rb.LAWS).write_text(laws)
    assert any("refute_mutant_A" in p for p in rb.static_ties(d))


def test_static_ties_catch_diverging_mutant_function(tmp_path):
    d = _copy_bend(tmp_path)
    f = d / "mutant_E_deal_drops_on_rotation.bend"
    f.write_text(f.read_text().replace("deal(p, h2, snoc(t2, h))\n\nlaw total_snoc",
                                       "deal(p, h2, snoc(t2, 1n+h))\n\nlaw total_snoc"))
    assert any("def deal" in p for p in rb.static_ties(d))


def test_static_ties_catch_unlisted_mutant(tmp_path):
    d = _copy_bend(tmp_path)
    shutil.copy(d / "mutant_B_full_n_each.bend", d / "mutant_F_new.bend")
    assert any("mutant_F_new.bend" in p for p in rb.static_ties(d))


def _fake(responses):
    """A runner that answers each ssh/scp call from a per-file script."""
    def run(cmd):
        text = " ".join(cmd)
        for key, (code, out) in responses.items():
            if key in text:
                return code, out
        return 0, ""
    return run


def _good(name, law):
    if law is None:
        return 0, f"All terms check.\n{rb.EXIT_MARK}0\n"
    return 0, TYPE_ERROR.replace(LAW, law) + f"\n{rb.EXIT_MARK}1\n"


def test_runner_passes_only_real_rejections(tmp_path):
    d = _copy_bend(tmp_path)
    resp = {"bend.sh version": (0, "bend 2.0.20")}
    resp[f"bend.sh {rb.LAWS} "] = _good(rb.LAWS, None)
    for f, spec in rb.MUTANTS.items():
        resp[f"bend.sh {f} "] = _good(f, spec.law)
    assert rb.main(run=_fake(resp), here=d) == 0


def test_runner_ssh_failure_on_one_mutant_is_not_a_rejection(tmp_path):
    """Astra 21: an ssh failure after the positive control used to be recorded as PASS."""
    d = _copy_bend(tmp_path)
    resp = {"bend.sh version": (0, "bend 2.0.20")}
    resp[f"bend.sh {rb.LAWS} "] = _good(rb.LAWS, None)
    for f, spec in rb.MUTANTS.items():
        resp[f"bend.sh {f} "] = _good(f, spec.law)
    resp["bend.sh mutant_C_swapped_numerator.bend "] = (255, "ssh: Connection reset by peer")
    assert rb.main(run=_fake(resp), here=d) == 2
    assert "mutant_C_swapped_numerator.bend" in next(
        x for x in (d / "receipt_bend.txt").read_text().splitlines() if x.endswith("| FAIL"))


def test_runner_unreachable_vm_fails_everything(tmp_path):
    d = _copy_bend(tmp_path)
    assert rb.main(run=lambda cmd: (255, "ssh: Could not resolve hostname"), here=d) == 2
    assert "| PASS" not in (d / "receipt_bend.txt").read_text()
