"""Draw Figure 3 (fig_transitions.pdf): population and fertility paths with nu = 0.5.

Usage: python fig_transitions.py owid-co2-data.csv [output.pdf]

Solid lines use the effort equivalent to current measured clean progress; dashed lines, drawn
where that effort falls short, use 1.1 times the required effort m* of equation (36).
Emission-intensity rates are the compound 2005-2022 rates of Section 5.1.
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from owid_checks import COUNTRIES, GY, INVERT, Q_INF, eps_star, simulate

NU = 0.5
GENERATIONS = 12  # 2023 to 2383
STYLE = {
    "Korea": dict(color="#111111", lw=1.6),
    "Italy": dict(color="#555555", lw=1.4),
    "Germany": dict(color="#999999", lw=1.4),
    "United States": dict(color="#2b6ca3", lw=1.6),
}


def main(path, out):
    raw = pd.read_csv(path, usecols=["iso_code", "year", "gdp", "co2"])
    plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7})
    fig, (ax_l, ax_n) = plt.subplots(1, 2, figsize=(7.0, 3.0))
    for name, style in STYLE.items():
        iso, tfr = COUNTRIES[name]
        d = raw[raw.iso_code == iso].set_index("year")
        intensity = d.co2 / d.gdp
        g_int = (intensity[2022] / intensity[2005]) ** (1 / 17) - 1
        g_b = (1 + g_int) ** -30 - 1
        L_star = Q_INF / eps_star(INVERT(tfr / 2.1))
        m_star = (GY - 1) / g_b * L_star ** -NU
        current = simulate(tfr, g_b, m=1.0, nu=NU, generations=GENERATIONS)
        ax_l.plot(current[:, 0], current[:, 1], label=name, **style)
        ax_n.plot(current[:, 0], 2.1 * current[:, 2], **style)
        if m_star > 1:
            raised = simulate(tfr, g_b, m=1.1 * m_star, nu=NU, generations=GENERATIONS)
            ax_l.plot(raised[:, 0], raised[:, 1], ls="--", **style)
            ax_n.plot(raised[:, 0], 2.1 * raised[:, 2], ls="--", **style)
    ax_l.set_yscale("log")
    ax_n.axhline(2.1, color="#bbbbbb", lw=0.8, zorder=0)
    for ax, title in [(ax_l, "(a) Population, 2023 = 1 (log scale)"), (ax_n, "(b) Total fertility rate")]:
        ax.set_title(title, loc="left", fontsize=9)
        ax.set_xticks([2100, 2200, 2300, 2400])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    ax_l.legend(loc="lower left", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "fig_transitions.pdf")
