"""
Plot the TFIM ZZ-coupling architectures used in the lattice model.
"""
import math
import os

import matplotlib.pyplot as plt

from models.tfim_lattice import _LATTICE_PAIRS


ROOT = os.path.dirname(__file__)
OUTPUT_PATH = os.path.join(ROOT, "results", "figures", "Fig1.png")
SITE_LABEL_FONTSIZE = 16
TITLE_FONTSIZE = 15
PANEL_LABEL_FONTSIZE = 18
SUPTITLE_FONTSIZE = 20


def _ring_positions(n_spins):
    return {
        site: (
            math.cos(math.pi / 2 + 2 * math.pi * site / n_spins),
            math.sin(math.pi / 2 + 2 * math.pi * site / n_spins),
        )
        for site in range(n_spins)
    }


def _positions_for(n_spins):
    positions = {
        2: {0: (0.0, 0.0), 1: (1.0, 0.0)},
        3: {0: (0.0, 0.0), 1: (1.0, 0.0), 2: (2.0, 0.0)},
        4: {0: (0.0, 1.0), 1: (1.0, 1.0), 2: (0.0, 0.0), 3: (1.0, 0.0)},
        6: {
            0: (0.0, 1.0),
            1: (1.0, 1.0),
            2: (2.0, 1.0),
            3: (0.0, 0.0),
            4: (1.0, 0.0),
            5: (2.0, 0.0),
        },
    }
    if n_spins == 5:
        return _ring_positions(n_spins)
    return positions[n_spins]


def _draw_architecture(ax, n_spins, pairs):
    positions = _positions_for(n_spins)

    for i, j in pairs:
        x_values = [positions[i][0], positions[j][0]]
        y_values = [positions[i][1], positions[j][1]]
        ax.plot(x_values, y_values, color="#444444", linewidth=2.4, zorder=1)

    xs = [positions[site][0] for site in range(n_spins)]
    ys = [positions[site][1] for site in range(n_spins)]
    ax.scatter(xs, ys, s=520, color="#f5f5f5", edgecolor="#222222", linewidth=1.6, zorder=2)

    for site, (x_coord, y_coord) in positions.items():
        ax.text(
            x_coord,
            y_coord,
            rf"$q_{{{site}}}$",
            ha="center",
            va="center",
            fontsize=SITE_LABEL_FONTSIZE,
            zorder=4,
        )

    n_pairs = len(pairs)
    m_terms = n_pairs + 2 * n_spins
    ax.set_title(rf"$n={n_spins}$, couplings={n_pairs}, $M={m_terms}$", fontsize=TITLE_FONTSIZE)
    ax.set_aspect("equal")
    ax.axis("off")

    x_padding = 0.55
    y_padding = 0.55
    ax.set_xlim(min(xs) - x_padding, max(xs) + x_padding)
    ax.set_ylim(min(ys) - y_padding, max(ys) + y_padding)


def _add_panel_label(ax, label):
    ax.text(
        0.5,
        -0.08,
        label,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=PANEL_LABEL_FONTSIZE,
        family="serif",
        clip_on=False,
    )


def main():
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)

    fig, axes = plt.subplots(1, len(_LATTICE_PAIRS), figsize=(16, 3.8))
    panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    for ax, (n_spins, pairs), label in zip(axes, sorted(_LATTICE_PAIRS.items()), panel_labels):
        _draw_architecture(ax, n_spins, pairs)
        _add_panel_label(ax, label)

    fig.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.18, wspace=0.28)
    fig.savefig(OUTPUT_PATH, dpi=300)
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
