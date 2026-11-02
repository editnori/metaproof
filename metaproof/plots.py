"""Forest and funnel plots (matplotlib, Agg backend)."""

from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import transforms  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#d9d8d4"
ACCENT = "#2a78d6"
SURFACE = "#fcfcfb"
RATIO_TICKS = [0.1, 0.2, 0.25, 0.5, 0.67, 0.8, 0.9, 1, 1.1, 1.25, 1.5, 2, 3, 4, 5, 10]

MEASURE_LABEL = {"MD": "Mean difference", "SMD": "Standardised mean difference (Hedges g)",
                 "ROM": "Ratio of means"}


def _fmt(x, ratio):
    return f"{math.exp(x):.2f}" if ratio else f"{x:.2f}"


# Fixed margins in inches, so the footer never runs into the x-axis label at any height.
TOP_IN = 0.45
FOOT_Y_IN = 0.06       # footer baseline above the bottom edge
FOOT_LINE_IN = 0.17    # one footer line at 8.5 pt
AXIS_BELOW_IN = 0.62   # tick labels plus the x-axis label, below the axes


def forest(path, eff, pooled, measure, level=0.95, title="", pi_note=None):
    """Write the forest plot to path. See forest_figure."""
    fig = forest_figure(eff, pooled, measure, level, title, pi_note)
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def forest_figure(eff, pooled, measure, level=0.95, title="", pi_note=None):
    """eff: DataFrame with study_id, yi, vi; pooled: pooled.csv row (dict) of the primary
    model, or None when the analysis is below min_k_pool (no diamond, no weights).
    pi_note: when given, the prediction interval is not drawn and the note is printed under
    the heterogeneity line (the pipeline passes it with fewer than 5 studies)."""
    ratio = measure == "ROM"
    k = len(eff)
    z = stats.norm.ppf(1 - (1 - level) / 2)
    y, se = eff["yi"].to_numpy(float), np.sqrt(eff["vi"].to_numpy(float))
    lb, ub = y - z * se, y + z * se
    show_pool = pooled is not None and np.isfinite(pooled.get("estimate", np.nan))
    tau2 = pooled["tau2"] if show_pool and np.isfinite(pooled["tau2"]) else 0.0
    w = 1 / (eff["vi"].to_numpy(float) + tau2)
    wpct = 100 * w / w.sum()
    has_pi = show_pool and np.isfinite(pooled["pi_lb"]) and not pi_note

    rows = k + (1 if has_pi else 0) + (2 if show_pool else 1)
    n_foot = 2 if pi_note and show_pool else 1
    height = 0.38 * rows + 1.6 + FOOT_LINE_IN * (n_foot - 1)
    fig, ax = plt.subplots(figsize=(9, height), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    base = (1 if has_pi else 0) + (0.5 if show_pool else -0.5)
    ypos = np.arange(k, 0, -1) + base
    ax.hlines(ypos, lb, ub, color=INK, lw=1.2)
    sizes = 20 + 80 * wpct / wpct.max() if show_pool else np.full(k, 40.0)
    ax.scatter(y, ypos, s=sizes, marker="s", color=INK, zorder=3, linewidths=0)

    labels, ticks = list(eff["study_id"]), list(ypos)
    xs_lo, xs_hi = [lb], [ub]
    if show_pool:
        est, clb, cub = pooled["estimate"], pooled["ci_lb"], pooled["ci_ub"]
        yd = 1.0 + (1 if has_pi else 0) - 0.3
        ax.fill([clb, est, cub, est], [yd, yd + 0.3, yd, yd - 0.3], color=ACCENT, zorder=3)
        ax.axvline(est, color=ACCENT, lw=0.8, ls=":")
        labels.append(f"Pooled ({pooled['estimator']}, {pooled.get('ci_method', '')})")
        ticks.append(yd)
        xs_lo.append([clb])
        xs_hi.append([cub])
        if has_pi:
            ax.hlines(0.7, pooled["pi_lb"], pooled["pi_ub"], color=ACCENT, lw=2)
            ax.vlines([pooled["pi_lb"], pooled["pi_ub"]], 0.55, 0.85, color=ACCENT, lw=2)
            labels.append(f"{int(level * 100)}% prediction interval")
            ticks.append(0.7)
            xs_lo.append([pooled["pi_lb"]])
            xs_hi.append([pooled["pi_ub"]])
    ax.axvline(0, color=INK_2, lw=0.8, ls="-")
    ax.set_yticks(ticks, labels, color=INK, fontsize=9)
    ax.set_ylim(0.2 if show_pool else 0.0, ypos.max() + 0.8)
    xmin = np.nanmin(np.concatenate([np.atleast_1d(a) for a in xs_lo]))
    xmax = np.nanmax(np.concatenate([np.atleast_1d(a) for a in xs_hi]))
    pad = 0.05 * (xmax - xmin or 1)
    ax.set_xlim(min(xmin, 0) - pad, max(xmax, 0) + pad)
    if ratio:
        lo, hi = ax.get_xlim()
        nice = [r for r in RATIO_TICKS if lo <= math.log(r) <= hi]
        ax.set_xticks([math.log(r) for r in nice], [f"{r:g}" for r in nice])
    ax.set_xlabel(MEASURE_LABEL[measure] + (" (log scale)" if ratio else ""), color=INK)

    # Text columns: x in figure coordinates, y in data coordinates.
    trans = transforms.blended_transform_factory(fig.transFigure, ax.transData)
    x_eff, x_w = 0.67, 0.92
    ax.text(x_eff, ypos.max() + 0.6, "Effect [CI]", transform=trans, fontsize=8.5,
            color=INK_2, va="bottom")
    if show_pool:
        ax.text(x_w, ypos.max() + 0.6, "Weight", transform=trans, fontsize=8.5, color=INK_2,
                va="bottom", ha="right")
    for yy, a, b, c, wp in zip(ypos, y, lb, ub, wpct):
        ax.text(x_eff, yy, f"{_fmt(a, ratio)} [{_fmt(b, ratio)}, {_fmt(c, ratio)}]",
                transform=trans, fontsize=8.5, color=INK, va="center")
        if show_pool:
            ax.text(x_w, yy, f"{wp:.1f}%", transform=trans, fontsize=8.5, color=INK,
                    va="center", ha="right")
    if show_pool:
        ax.text(x_eff, yd, f"{_fmt(est, ratio)} [{_fmt(clb, ratio)}, {_fmt(cub, ratio)}]",
                transform=trans, fontsize=8.5, color=INK, va="center", fontweight="bold")
        if has_pi:
            ax.text(x_eff, 0.7, f"[{_fmt(pooled['pi_lb'], ratio)}, "
                    f"{_fmt(pooled['pi_ub'], ratio)}]", transform=trans, fontsize=8.5,
                    color=INK, va="center")
        het = (f"k = {k}; tau² = {pooled['tau2']:.3g}; I² = {pooled['I2']:.1f}%; "
               f"Q = {pooled['Q']:.2f} (p = {pooled['Q_p']:.3g})")
        if pi_note:
            het += "\n" + pi_note
    else:
        het = f"k = {k}: no pooled estimate (below min_k_pool); study results only"
    ax.set_title(title, loc="left", color=INK, fontsize=11)
    fig.text(0.01, FOOT_Y_IN / height, het, color=INK_2, fontsize=8.5, va="bottom",
             linespacing=1.3)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(INK_2)
    ax.tick_params(colors=INK_2, length=0, axis="y")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    bottom = FOOT_Y_IN + FOOT_LINE_IN * n_foot + AXIS_BELOW_IN
    fig.subplots_adjust(left=0.26, right=0.64, bottom=bottom / height, top=1 - TOP_IN / height)
    return fig


CONTOURS = ((0.10, "#e9e8e4"), (0.05, "#d5d4cf"), (0.01, "#bdbcb6"))


def funnel(path, eff, center, measure, level=0.95, title=""):
    """Contour-enhanced funnel plot (Peters JL, Sutton AJ, Jones DR, Abrams KR, Rushton L.
    J Clin Epidemiol 2008;61(10):991-996): effect against SE (inverted), with shaded regions
    where a study's own two-sided p is between 0.10 and 0.05, 0.05 and 0.01, and below 0.01
    (centred on 0), plus pseudo-confidence limits about the pooled estimate when there is one.
    """
    y, se = eff["yi"].to_numpy(float), np.sqrt(eff["vi"].to_numpy(float))
    fig, ax = plt.subplots(figsize=(6, 5), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    smax = se.max() * 1.1
    s_grid = np.linspace(0, smax, 50)
    span = max(np.abs(y).max() + stats.norm.ppf(0.995) * smax, 1e-9) * 1.1
    for p_lev, col in CONTOURS:
        zc = stats.norm.ppf(1 - p_lev / 2)
        ax.fill_betweenx(s_grid, zc * s_grid, span, color=col, lw=0, zorder=0)
        ax.fill_betweenx(s_grid, -span, -zc * s_grid, color=col, lw=0, zorder=0)
    ax.axvline(0, color=INK_2, lw=0.8)
    if center is not None and np.isfinite(center):
        z = stats.norm.ppf(1 - (1 - level) / 2)
        ax.plot(center - z * s_grid, s_grid, color=INK_2, lw=0.8, ls="--")
        ax.plot(center + z * s_grid, s_grid, color=INK_2, lw=0.8, ls="--")
        ax.axvline(center, color=ACCENT, lw=1.2)
    ax.scatter(y, se, s=40, color=INK, zorder=3, edgecolors=SURFACE, linewidths=1.5)
    ax.set_xlim(-span, span)
    ax.set_ylim(smax, 0)
    ax.set_ylabel("Standard error", color=INK)
    ax.set_xlabel("Log ratio of means" if measure == "ROM" else MEASURE_LABEL[measure], color=INK)
    ax.set_title(title, loc="left", color=INK, fontsize=11)
    ax.text(0.99, 0.02, "shading: p 0.10-0.05, 0.05-0.01, < 0.01", transform=ax.transAxes,
            ha="right", fontsize=7.5, color=INK_2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=INK_2)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
