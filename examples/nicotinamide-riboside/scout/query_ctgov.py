"""Query ClinicalTrials.gov (API v2) for nicotinamide riboside studies.

Writes:
    examples/nicotinamide-riboside/scout/ctgov_nr_registrations.csv   every study the query returns
    examples/nicotinamide-riboside/scout/ctgov_results_posted.csv     the subset with results posted

Each row is linked to examples/nicotinamide-riboside/scout/nr_candidate_trials.csv records that carry the same
NCT ID (column linked_records), so unpublished trials stand out. A link the screener
only judged probable (registration_id starts with UNVERIFIED) is labelled
"[probable, unverified]", and link_status says whether the registration has a stated
link to a PubMed report, only probable ones, only a registry-only row, or none. A
probable link is not evidence that the trial has been published.

Usage:
    .venv/bin/python examples/nicotinamide-riboside/scout/query_ctgov.py            # query the registry and relink
    .venv/bin/python examples/nicotinamide-riboside/scout/query_ctgov.py --relink   # relink the saved CSVs, no network
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
API = "https://clinicaltrials.gov/api/v2/studies"
# Two queries, unioned. The intervention-field query catches product names; the all-fields
# query catches trials that name NR only in the description or arms (e.g. NCT07802548).
QUERIES = [
    {"query.intr": 'nicotinamide riboside OR Niagen OR "NR chloride" OR EH301 OR NRPT OR "Tru Niagen"'},
    {"query.term": '"nicotinamide riboside"'},
]

COLUMNS = ["nct_id", "brief_title", "overall_status", "study_type", "phases", "allocation",
           "intervention_model", "masking", "enrollment", "enrollment_type", "start_date",
           "completion_date", "has_results", "results_first_posted", "conditions",
           "interventions", "primary_outcomes", "registry_pmids", "linked_records", "link_status",
           "retrieved_utc"]
PROBABLE = " [probable, unverified]"


def fetch(query: dict) -> list[dict]:
    studies, token = [], None
    while True:
        params = {**query, "pageSize": 200, "format": "json"}
        if token:
            params["pageToken"] = token
        url = API + "?" + urllib.parse.urlencode(params)
        req = urllib.request.Request(url, headers={"User-Agent": "metaproof/scout"})
        with urllib.request.urlopen(req, timeout=60) as r:
            d = json.loads(r.read().decode())
        studies.extend(d.get("studies", []))
        token = d.get("nextPageToken")
        if not token:
            return studies


def fetch_all() -> list[dict]:
    seen: dict[str, dict] = {}
    for q in QUERIES:
        for s in fetch(q):
            seen.setdefault(s["protocolSection"]["identificationModule"]["nctId"], s)
    return list(seen.values())


def candidate_links(path: Path | None = None) -> dict[str, list[str]]:
    path = path or HERE / "nr_candidate_trials.csv"
    links: dict[str, list[str]] = {}
    if not path.exists():
        return links
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            reg = row.get("registration_id", "")
            probable = reg.strip().startswith("UNVERIFIED")
            for nct in sorted(set(re.findall(r"NCT\d{8}", reg))):
                label = row["record_id"] + (f" PMID {row['pmid']}" if row.get("pmid") else " registry")
                links.setdefault(nct, []).append(label + (PROBABLE if probable else ""))
    return links


def link_status(labels: list[str]) -> str:
    """stated report | probable report only | registry-only row | none."""
    reports = [x for x in labels if " PMID " in x]
    if any(not x.endswith(PROBABLE) for x in reports):
        return "stated report"
    if reports:
        return "probable report only"
    return "registry-only row" if labels else "none"


def flatten(s: dict, links: dict, now: str) -> dict:
    p = s["protocolSection"]
    idm, st, de = p["identificationModule"], p.get("statusModule", {}), p.get("designModule", {})
    di = de.get("designInfo", {})
    arms = p.get("armsInterventionsModule", {}).get("interventions", [])
    outs = p.get("outcomesModule", {}).get("primaryOutcomes", [])
    refs = p.get("referencesModule", {}).get("references", [])
    nct = idm["nctId"]
    return {
        "nct_id": nct,
        "brief_title": idm.get("briefTitle", ""),
        "overall_status": st.get("overallStatus", ""),
        "study_type": de.get("studyType", ""),
        "phases": ";".join(de.get("phases", []) or []),
        "allocation": di.get("allocation", ""),
        "intervention_model": di.get("interventionModel", ""),
        "masking": di.get("maskingInfo", {}).get("masking", ""),
        "enrollment": de.get("enrollmentInfo", {}).get("count", ""),
        "enrollment_type": de.get("enrollmentInfo", {}).get("type", ""),
        "start_date": st.get("startDateStruct", {}).get("date", ""),
        "completion_date": st.get("completionDateStruct", {}).get("date", ""),
        "has_results": "yes" if s.get("hasResults") else "no",
        "results_first_posted": st.get("resultsFirstPostDateStruct", {}).get("date", ""),
        "conditions": "; ".join(p.get("conditionsModule", {}).get("conditions", []) or []),
        "interventions": " | ".join(f"{i.get('type', '')}: {i.get('name', '')}" for i in arms),
        "primary_outcomes": " | ".join(o.get("measure", "") for o in outs),
        "registry_pmids": "; ".join(f"{r['pmid']} ({r.get('type', '')})" for r in refs if r.get("pmid")),
        "linked_records": "; ".join(links.get(nct, [])),
        "link_status": link_status(links.get(nct, [])),
        "retrieved_utc": now,
    }


def relink(rows: list[dict], links: dict) -> list[dict]:
    """Recompute linked_records and link_status on saved rows (no registry query)."""
    for r in rows:
        r["linked_records"] = "; ".join(links.get(r["nct_id"], []))
        r["link_status"] = link_status(links.get(r["nct_id"], []))
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--relink", action="store_true",
                    help="relink the saved registrations to nr_candidate_trials.csv; no network")
    args = ap.parse_args(argv)
    links = candidate_links()
    if args.relink:
        with open(HERE / "ctgov_nr_registrations.csv", newline="") as fh:
            rows = relink(list(csv.DictReader(fh)), links)
    else:
        now = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%MZ")
        rows = sorted((flatten(s, links, now) for s in fetch_all()), key=lambda r: r["nct_id"])
    for name, subset in (("ctgov_nr_registrations.csv", rows),
                         ("ctgov_results_posted.csv", [r for r in rows if r["has_results"] == "yes"])):
        with open(HERE / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS)
            w.writeheader()
            w.writerows(subset)
        print(f"wrote {len(subset)} rows to {HERE / name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
