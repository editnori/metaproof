"""Draw a PRISMA 2020 flow diagram (new review) from a small counts CSV.

Usage:
    .venv/bin/python templates/prisma_flow.py counts.csv --out prisma_flow.png
    .venv/bin/python templates/prisma_flow.py counts.csv --out prisma_flow.pdf --title "..."

Counts CSV: columns key,label,n. One row per box value. Keys:

  databases, registers                    records identified
  duplicates_removed, automation_ineligible, other_removed
                                          records removed before screening
  records_screened, records_excluded
  reports_sought, reports_not_retrieved
  reports_assessed
  reports_excluded_reason                 one row per reason; label = reason
  studies_included, reports_included
  other_source                            optional, one row per source
                                          (label = Websites, Organisations,
                                          Citation searching, ...)
  other_reports_sought, other_reports_not_retrieved, other_reports_assessed
  other_reports_excluded_reason           one row per reason

A blank n prints as "(n = )", so the template draws an empty diagram. If any
other_reports_* key has a value, or any other_source or
other_reports_excluded_reason row is present (even with a blank n), the "other
methods" column of the 2020 template is drawn. The previous-studies column
(updated reviews) is not drawn.

Layout follows the PRISMA 2020 flow diagram for new systematic reviews which
included searches of databases, registers and other sources (Page MJ et al.
BMJ 2021;372:n71, figure 1). Arithmetic is checked and any mismatch is printed
as a warning; the diagram is still drawn with the numbers given.
"""

from __future__ import annotations

import argparse
import csv
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

SINGLE_KEYS = {
    "databases", "registers", "duplicates_removed", "automation_ineligible",
    "other_removed", "records_screened", "records_excluded", "reports_sought",
    "reports_not_retrieved", "reports_assessed", "studies_included",
    "reports_included", "other_reports_sought", "other_reports_not_retrieved",
    "other_reports_assessed",
}
LIST_KEYS = {"reports_excluded_reason", "other_source", "other_reports_excluded_reason"}

FONT = 8.0
LINE_H = 0.165      # inches per text line at FONT
PAD = 0.09          # inches inside a box
ROW_GAP = 0.38      # inches between rows
MIN_ROW_H = 0.75    # inches; keeps the rotated 'Included' label inside its bar
PHASE_FILL = "#c6dbef"
HEADER_FILL = "#f5c242"
GREY_FILL = "#e6e6e6"


@dataclass
class Counts:
    single: dict[str, int | None] = field(default_factory=dict)
    lists: dict[str, list[tuple[str, int | None]]] = field(default_factory=dict)

    def get(self, key):
        return self.single.get(key)

    def items(self, key):
        return self.lists.get(key, [])

    def has_other(self) -> bool:
        if any(self.single.get(k) is not None for k in SINGLE_KEYS if k.startswith("other_reports")):
            return True
        return bool(self.items("other_source") or self.items("other_reports_excluded_reason"))


def _parse_n(raw: str, where: str) -> int | None:
    raw = (raw or "").strip().replace(",", "")
    if raw == "":
        return None
    try:
        v = float(raw)
    except ValueError:
        raise ValueError(f"{where}: n must be a whole number, got {raw!r}")
    if v < 0 or v != int(v):
        raise ValueError(f"{where}: n must be a whole number >= 0, got {raw!r}")
    return int(v)


def read_counts(path: Path) -> Counts:
    c = Counts()
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        missing = {"key", "label", "n"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        for lineno, row in enumerate(reader, start=2):
            key = (row["key"] or "").strip()
            if not key or key.startswith("#"):
                continue
            where = f"{path.name} line {lineno}"
            n = _parse_n(row["n"], where)
            if key in SINGLE_KEYS:
                if key in c.single:
                    raise ValueError(f"{where}: {key} given twice")
                c.single[key] = n
            elif key in LIST_KEYS:
                label = (row["label"] or "").strip()
                if not label:
                    raise ValueError(f"{where}: {key} needs a label")
                c.lists.setdefault(key, []).append((label, n))
            else:
                raise ValueError(f"{where}: unknown key {key!r}")
    return c


def _sum(values):
    vals = list(values)
    if not vals or any(v is None for v in vals):
        return None
    return sum(vals)


def check_counts(c: Counts) -> list[str]:
    """Return warnings for arithmetic that does not add up. A check is skipped when a value
    it needs is blank, except that a blank databases or registers count is read as 0 when the
    other is filled, and blank removal counts are read as 0."""
    w = []

    def eq(label, lhs, rhs):
        if lhs is not None and rhs is not None and lhs != rhs:
            w.append(f"{label}: expected {rhs}, got {lhs}")

    if c.get("databases") is not None or c.get("registers") is not None:
        identified = (c.get("databases") or 0) + (c.get("registers") or 0)
        removed = sum(c.get(k) or 0 for k in ("duplicates_removed", "automation_ineligible",
                                                "other_removed"))
        eq("records_screened = identified - removed before screening",
           c.get("records_screened"), identified - removed)
    if c.get("records_screened") is not None and c.get("records_excluded") is not None:
        eq("reports_sought = records_screened - records_excluded",
           c.get("reports_sought"), c.get("records_screened") - c.get("records_excluded"))
    if c.get("reports_sought") is not None and c.get("reports_not_retrieved") is not None:
        eq("reports_assessed = reports_sought - reports_not_retrieved",
           c.get("reports_assessed"), c.get("reports_sought") - c.get("reports_not_retrieved"))
    if c.get("other_reports_sought") is not None and c.get("other_reports_not_retrieved") is not None:
        eq("other_reports_assessed = other_reports_sought - other_reports_not_retrieved",
           c.get("other_reports_assessed"),
           c.get("other_reports_sought") - c.get("other_reports_not_retrieved"))

    excl = _sum(n for _, n in c.items("reports_excluded_reason"))
    other_excl = _sum(n for _, n in c.items("other_reports_excluded_reason")) \
        if c.items("other_reports_excluded_reason") else 0
    assessed = c.get("reports_assessed")
    # With an other-methods column, a blank other_reports_assessed means the sum cannot be
    # checked (it is not 0); without one, the other-methods terms are 0.
    other_assessed = c.get("other_reports_assessed") if c.has_other() else 0
    if assessed is not None and excl is not None and other_excl is not None \
            and other_assessed is not None:
        eq("reports_included = (reports_assessed - excluded) + (other_reports_assessed - other excluded)",
           c.get("reports_included"), assessed - excl + other_assessed - other_excl)
    si, ri = c.get("studies_included"), c.get("reports_included")
    if si is not None and ri is not None and si > ri:
        w.append(f"studies_included ({si}) is larger than reports_included ({ri})")
    return w


_KEEP = "\ue000"  # stands in for spaces inside "(n = x)" so wrapping never splits it


def _n(v) -> str:
    text = "(n = )" if v is None else f"(n = {v:,})"
    return text.replace(" ", _KEEP)


def _wrap(lines, width_in):
    chars = max(12, int((width_in - 2 * PAD) / (FONT * 0.52 / 72)))
    out = []
    for line in lines:
        indent = len(line) - len(line.lstrip(" "))
        out.extend(textwrap.wrap(line, chars, subsequent_indent=" " * (indent + 2)) or [""])
    return [o.replace(_KEEP, " ") for o in out]


@dataclass
class Box:
    x: float
    w: float
    lines: list[str]
    fill: str = "white"
    y: float = 0.0   # top, set by layout
    h: float = 0.0

    def wrap(self):
        self.lines = _wrap(self.lines, self.w)
        self.h = len(self.lines) * LINE_H + 2 * PAD


def build_boxes(c: Counts, other: bool):
    if other:
        cols = dict(main=(0.65, 2.75), excl=(3.75, 2.6), oth=(6.85, 2.75), oexcl=(9.95, 2.6))
        width = 12.75
    else:
        cols = dict(main=(0.65, 3.6), excl=(4.85, 3.6))
        width = 8.65

    def box(col, lines, fill="white"):
        x, w = cols[col]
        return Box(x, w, lines, fill)

    reasons = c.items("reports_excluded_reason") or [("Reason 1", None), ("Reason 2", None)]
    rows = [
        {"main": box("main", ["Records identified from:",
                              f"  Databases {_n(c.get('databases'))}",
                              f"  Registers {_n(c.get('registers'))}"]),
         "excl": box("excl", ["Records removed before screening:",
                              f"  Duplicate records removed {_n(c.get('duplicates_removed'))}",
                              "  Records marked as ineligible by automation tools "
                              f"{_n(c.get('automation_ineligible'))}",
                              f"  Records removed for other reasons {_n(c.get('other_removed'))}"])},
        {"main": box("main", [f"Records screened {_n(c.get('records_screened'))}"]),
         "excl": box("excl", [f"Records excluded {_n(c.get('records_excluded'))}"])},
        {"main": box("main", [f"Reports sought for retrieval {_n(c.get('reports_sought'))}"]),
         "excl": box("excl", [f"Reports not retrieved {_n(c.get('reports_not_retrieved'))}"])},
        {"main": box("main", [f"Reports assessed for eligibility {_n(c.get('reports_assessed'))}"]),
         "excl": box("excl", ["Reports excluded:"] + [f"  {lab} {_n(n)}" for lab, n in reasons])},
        {"main": box("main", [f"Studies included in review {_n(c.get('studies_included'))}",
                              f"Reports of included studies {_n(c.get('reports_included'))}"])},
    ]
    if other:
        sources = c.items("other_source") or [("Websites", None), ("Organisations", None),
                                              ("Citation searching", None)]
        oreasons = c.items("other_reports_excluded_reason") or [("Reason 1", None), ("Reason 2", None)]
        rows[0]["oth"] = box("oth", ["Records identified from:"] + [f"  {lab} {_n(n)}" for lab, n in sources])
        rows[2]["oth"] = box("oth", [f"Reports sought for retrieval {_n(c.get('other_reports_sought'))}"])
        rows[2]["oexcl"] = box("oexcl", [f"Reports not retrieved {_n(c.get('other_reports_not_retrieved'))}"])
        rows[3]["oth"] = box("oth", ["Reports assessed for eligibility "
                                     f"{_n(c.get('other_reports_assessed'))}"])
        rows[3]["oexcl"] = box("oexcl", ["Reports excluded:"] + [f"  {lab} {_n(n)}" for lab, n in oreasons])
    return rows, cols, width


def draw(c: Counts, out: Path, title: str | None = None) -> Path:
    other = c.has_other()
    rows, cols, width = build_boxes(c, other)

    header_h = 0.42
    title_h = 0.45 if title else 0.0
    y = 0.15 + title_h + header_h + ROW_GAP * 0.6
    for row in rows:
        for b in row.values():
            b.wrap()
        h = max(b.h for b in row.values())
        if row is rows[-1]:
            h = max(h, MIN_ROW_H)
        for b in row.values():
            b.y, b.h = y, h
        y += h + ROW_GAP
    height = y - ROW_GAP + 0.2

    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.axis("off")

    if title:
        ax.text(width / 2, 0.15 + title_h / 2, title, ha="center", va="center",
                fontsize=FONT + 2, fontweight="bold")
    hy = 0.15 + title_h
    mx, _ = cols["main"]
    ex, ew = cols["excl"]
    ax.add_patch(FancyBboxPatch((mx, hy), ex + ew - mx, header_h, boxstyle="round,pad=0,rounding_size=0.08",
                                fc=HEADER_FILL, ec="black", lw=0.8))
    ax.text((mx + ex + ew) / 2, hy + header_h / 2, "Identification of studies via databases and registers",
            ha="center", va="center", fontsize=FONT + 0.5, fontweight="bold")
    if other:
        ox, _ = cols["oth"]
        oex, oew = cols["oexcl"]
        ax.add_patch(FancyBboxPatch((ox, hy), oex + oew - ox, header_h,
                                    boxstyle="round,pad=0,rounding_size=0.08",
                                    fc=GREY_FILL, ec="black", lw=0.8))
        ax.text((ox + oex + oew) / 2, hy + header_h / 2, "Identification of studies via other methods",
                ha="center", va="center", fontsize=FONT + 0.5, fontweight="bold")

    for row in rows:
        for b in row.values():
            ax.add_patch(Rectangle((b.x, b.y), b.w, b.h, fc=b.fill, ec="black", lw=0.8))
            ax.text(b.x + PAD, b.y + b.h / 2, "\n".join(b.lines), ha="left", va="center",
                    fontsize=FONT, linespacing=1.25)

    phases = [("Identification", rows[0], rows[0]), ("Screening", rows[1], rows[3]),
              ("Included", rows[4], rows[4])]
    for label, first, last in phases:
        top = first["main"].y
        bottom = last["main"].y + last["main"].h
        ax.add_patch(FancyBboxPatch((0.12, top), 0.36, bottom - top,
                                    boxstyle="round,pad=0,rounding_size=0.08",
                                    fc=PHASE_FILL, ec="black", lw=0.8))
        ax.text(0.30, (top + bottom) / 2, label, rotation=90, ha="center", va="center",
                fontsize=FONT + 0.5, fontweight="bold")

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", lw=0.9, color="black",
                                    shrinkA=0, shrinkB=0, mutation_scale=9))

    def down(a: Box, b: Box):
        xm = a.x + a.w / 2
        arrow(xm, a.y + a.h, xm, b.y)

    def right(a: Box, b: Box):
        ym = a.y + a.h / 2
        arrow(a.x + a.w, ym, b.x, ym)

    for i in range(4):
        down(rows[i]["main"], rows[i + 1]["main"])
    for i in range(4):
        right(rows[i]["main"], rows[i]["excl"])
    if other:
        down(rows[0]["oth"], rows[2]["oth"])
        down(rows[2]["oth"], rows[3]["oth"])
        right(rows[2]["oth"], rows[2]["oexcl"])
        right(rows[3]["oth"], rows[3]["oexcl"])
        src, dst = rows[3]["oth"], rows[4]["main"]
        xm = src.x + src.w / 2
        ym = dst.y + dst.h / 2
        ax.plot([xm, xm], [src.y + src.h, ym], color="black", lw=0.9)
        arrow(xm, ym, dst.x + dst.w, ym)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("counts", type=Path)
    ap.add_argument("--out", type=Path, default=Path("prisma_flow.png"))
    ap.add_argument("--title", default=None)
    args = ap.parse_args(argv)
    try:
        counts = read_counts(args.counts)
    except ValueError as e:
        print(f"ERROR {e}", file=sys.stderr)
        return 1
    for w in check_counts(counts):
        print(f"WARNING {w}")
    print(f"wrote {draw(counts, args.out, args.title)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
