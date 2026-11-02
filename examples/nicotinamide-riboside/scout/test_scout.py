"""Offline checks of the scout tables: trial-identity links and the registry-only list."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import build_candidate_table as bct  # noqa: E402
import query_ctgov as qc  # noqa: E402


def _read(name):
    with open(HERE / name, newline="") as fh:
        return list(csv.DictReader(fh))


def test_probable_registry_links_keep_their_uncertainty(tmp_path):
    p = tmp_path / "cand.csv"
    p.write_text("record_id,pmid,registration_id\n"
                 "R001,111,NCT00000001\n"
                 "R002,222,UNVERIFIED: probable NCT00000002 (title match)\n"
                 "R003,,NCT00000003\n")
    links = qc.candidate_links(p)
    assert links["NCT00000001"] == ["R001 PMID 111"]
    assert links["NCT00000002"] == ["R002 PMID 222 [probable, unverified]"]
    assert qc.link_status(links["NCT00000001"]) == "stated report"
    assert qc.link_status(links["NCT00000002"]) == "probable report only"
    assert qc.link_status(links["NCT00000003"]) == "registry-only row"
    assert qc.link_status([]) == "none"


def test_registration_link_kind():
    assert bct.registration_link("NCT03568968") == "stated"
    assert bct.registration_link("UNVERIFIED: probable NCT03482167 (title match)") == \
        "probable (unverified)"
    assert bct.registration_link(bct.NOT_IN_ABS) == ""
    assert bct.registration_link("NCT05411835", registry_row=True) == "registry record"


def test_committed_tables_carry_link_status():
    cand = _read("nr_candidate_trials.csv")
    for r in cand:
        want = bct.registration_link(r["registration_id"], r["suggested_status"] == "registry_only")
        assert r["registration_link"] == want, r["record_id"]
        if want == "probable (unverified)":
            assert "UNVERIFIED" in r["ctgov_check"], r["record_id"]
    links = qc.candidate_links()
    for r in _read("ctgov_nr_registrations.csv"):
        assert r["linked_records"] == "; ".join(links.get(r["nct_id"], [])), r["nct_id"]
        assert r["link_status"] == qc.link_status(links.get(r["nct_id"], [])), r["nct_id"]


def test_registry_only_list_covers_every_unlinked_completed_randomized_trial():
    # CR7-3: the list's own definition is completed or terminated randomized registrations
    # with actual enrollment and no PubMed report linked.
    linked = set()
    for r in _read("nr_candidate_trials.csv"):
        linked |= set(re.findall(r"NCT\d{8}", r["registration_id"]))
    missing = [r["nct_id"] for r in _read("ctgov_nr_registrations.csv")
               if r["allocation"] == "RANDOMIZED" and r["enrollment_type"] == "ACTUAL"
               and r["overall_status"] in ("COMPLETED", "TERMINATED") and r["nct_id"] not in linked]
    assert missing == []


def test_search_strategy_counts_match_the_registry_table():
    # CR7-8b and the second review: counts quoted in search_strategy.txt come from the table.
    reg = _read("ctgov_nr_registrations.csv")
    rand = [r for r in reg if r["allocation"] == "RANDOMIZED"]
    stated = [r for r in rand if r["link_status"] == "stated report"]
    probable = [r for r in rand if r["link_status"] == "probable report only"]
    text = " ".join((HERE / "search_strategy.txt").read_text().split())
    assert f"{len(rand)} of the {len(reg)} registrations are randomized" in text
    assert f"{len(stated)} of those {len(rand)} have a stated link" in text
    assert f"{len(probable)} more have only a probable link" in text
    n = len(_read("nr_candidate_trials.csv"))
    first = n - len(bct.REGISTRY_ONLY) + 1
    assert f"R{first:03d}-R{n:03d}" in text
    assert f"R001-R{n:03d}" in (HERE / "pooling_issues.txt").read_text()
