"""
Plot Figure 3 (TFIM lattice) from precomputed JSON results.
"""
import json
import os

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterSciNotation

ROOT = os.path.dirname(__file__)
RESULTS_PATH = os.path.join(ROOT, "results", "data", "results_tfim.json")
OUTPUT_PATH = os.path.join(ROOT, "results", "figures", "fig3.png")
REFERENCE_M_VALUES = [5, 8, 12, 15, 19]
TITLE_FONTSIZE = 18
LABEL_FONTSIZE = 16
TICK_FONTSIZE = 13
LEGEND_FONTSIZE = 11
PANEL_LABEL_FONTSIZE = 20
Y_AXIS_LABEL = r"Trotter steps $N$ and gate complexity $g$ (log scale)"


def _extract_series(rows, method_key, value_key):
    return [row["methods"][method_key][value_key] for row in rows]


def _plot_method_subplot(ax, rows, method_key, title):
    m_values = [row["M"] for row in rows]

    n_analytic = _extract_series(rows, method_key, "N_analytic")
    n_min = _extract_series(rows, method_key, "N_min")
    g_analytic = _extract_series(rows, method_key, "g_analytic")
    g_min = _extract_series(rows, method_key, "g_min")

    ax.plot(
        m_values,
        n_analytic,
        linestyle="--",
        marker="s",
        linewidth=1.2,
        markersize=5,
        color="#1b9e77",
        label=r"$N^{analytic}_{1,\mathrm{det}}$" if method_key == "det1"
        else r"$N^{analytic}_{1,\mathrm{ran}}$" if method_key == "ran1"
        else r"$N^{analytic}_{2,\mathrm{det}}$" if method_key == "det2"
        else r"$N^{analytic}_{2,\mathrm{ran}}$",
    )
    ax.plot(
        m_values,
        n_min,
        linestyle="--",
        marker="x",
        linewidth=1.2,
        markersize=5,
        color="#66a61e",
        label=r"$N^{min}_{1,\mathrm{det}}$" if method_key == "det1"
        else r"$N^{min}_{1,\mathrm{ran}}$" if method_key == "ran1"
        else r"$N^{min}_{2,\mathrm{det}}$" if method_key == "det2"
        else r"$N^{min}_{2,\mathrm{ran}}$",
    )
    ax.plot(
        m_values,
        g_analytic,
        linestyle=":",
        marker="o",
        linewidth=1.2,
        markersize=4,
        color="#6a3d9a",
        label=r"$g^{analytic}_{1,\mathrm{det}}$" if method_key == "det1"
        else r"$g^{analytic}_{1,\mathrm{ran}}$" if method_key == "ran1"
        else r"$g^{analytic}_{2,\mathrm{det}}$" if method_key == "det2"
        else r"$g^{analytic}_{2,\mathrm{ran}}$",
    )
    ax.plot(
        m_values,
        g_min,
        linestyle=":",
        marker="D",
        linewidth=1.2,
        markersize=4,
        color="#7570b3",
        label=r"$g^{min}_{1,\mathrm{det}}$" if method_key == "det1"
        else r"$g^{min}_{1,\mathrm{ran}}$" if method_key == "ran1"
        else r"$g^{min}_{2,\mathrm{det}}$" if method_key == "det2"
        else r"$g^{min}_{2,\mathrm{ran}}$",
    )

    ax.set_title(title, fontsize=TITLE_FONTSIZE)
    ax.set_xlabel("Number of Liouvillian terms, $M$", fontsize=LABEL_FONTSIZE)
    ax.set_xticks(m_values)
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE)

    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=8))
    ax.yaxis.set_major_formatter(LogFormatterSciNotation(base=10))

    ax.legend(loc="upper left", fontsize=LEGEND_FONTSIZE, frameon=True)
    ax.grid(alpha=0.3, which="both")


def main():
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    with open(RESULTS_PATH, "r", encoding="utf-8") as f:
        payload = json.load(f)

    rows = [row for row in payload["rows"] if row["M"] in REFERENCE_M_VALUES]

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    _plot_method_subplot(
        axes[0, 0], rows, "det1", "First Order Deterministic"
    )
    _plot_method_subplot(
        axes[0, 1], rows, "ran1", "First Order Randomised"
    )
    _plot_method_subplot(
        axes[1, 0], rows, "det2", "Second Order Deterministic"
    )
    _plot_method_subplot(
        axes[1, 1], rows, "ran2", "Second Order Randomised"
    )
    fig.supylabel(Y_AXIS_LABEL, fontsize=LABEL_FONTSIZE)

    axes[0, 0].text(
        0.5, -0.18, "(a)", transform=axes[0, 0].transAxes,
        ha="center", va="top", fontsize=PANEL_LABEL_FONTSIZE, family="serif", clip_on=False
    )
    axes[0, 1].text(
        0.5, -0.18, "(b)", transform=axes[0, 1].transAxes,
        ha="center", va="top", fontsize=PANEL_LABEL_FONTSIZE, family="serif", clip_on=False
    )
    axes[1, 0].text(
        0.5, -0.20, "(c)", transform=axes[1, 0].transAxes,
        ha="center", va="top", fontsize=PANEL_LABEL_FONTSIZE, family="serif", clip_on=False
    )
    axes[1, 1].text(
        0.5, -0.20, "(d)", transform=axes[1, 1].transAxes,
        ha="center", va="top", fontsize=PANEL_LABEL_FONTSIZE, family="serif", clip_on=False
    )

    fig.subplots_adjust(left=0.08, right=0.98, top=0.93, bottom=0.10, hspace=0.42, wspace=0.25)
    fig.savefig(OUTPUT_PATH, dpi=300)

    print(f"Wrote {OUTPUT_PATH}")
    print(
        "The number of Trotter steps, N, and gate complexities, g, required by each method "
        "for the transverse-field Ising model with varying M, the number of terms in the "
        "model's Liouvillian generator. Parameters: lambda = 4.00, epsilon = 10^-5, t = 5."
    )


if __name__ == "__main__":
    main()
