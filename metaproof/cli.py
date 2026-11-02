"""Command line: python -m metaproof run <extraction.csv> --config <config.yaml> --out out/"""

from __future__ import annotations

import argparse
import sys

from .config import UserError
from .pipeline import run


def _one_line(exc: BaseException) -> str:
    return " ".join(str(exc).split())


def main(argv=None) -> int:
    """Exit 0 when the analyses ran, 2 for a problem in the inputs (printed as one line,
    without a traceback), 1 when the extraction file had no outcome to analyse."""
    ap = argparse.ArgumentParser(prog="metaproof", description="Meta-analysis of continuous outcomes")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run every configured analysis")
    r.add_argument("extraction", help="extraction CSV, one row per arm x outcome x timepoint")
    r.add_argument("--config", help="YAML or JSON config (defaults used when omitted)")
    r.add_argument("--out", default="out", help="output directory (default: out)")
    args = ap.parse_args(argv)
    try:
        res = run(args.extraction, args.config, args.out)
    except UserError as exc:
        print(f"metaproof: error: {_one_line(exc)}", file=sys.stderr)
        return 2
    except OSError as exc:            # unreadable file, full disk, no permission
        where = f"{exc.filename}: " if exc.filename else ""
        print(f"metaproof: error: {where}{exc.strerror or _one_line(exc)}", file=sys.stderr)
        return 2
    for line in res.attrs.get("notices", []):
        print(f"metaproof: {line}", file=sys.stderr)
    if not len(res):
        print("metaproof: no analyses: the extraction file has no outcomes", file=sys.stderr)
        return 1
    print(res.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
