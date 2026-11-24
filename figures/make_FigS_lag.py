# -*- coding: utf-8 -*-
"""Supplementary figure: lag L between collision and deposit.

(a) schematic (no data): a deposit of age t is compared with the track and the margin at t (L = 0,
    the case of the main text) or at t + L, when the track segment that could have generated the
    magma reached the trench. As the track moves along the margin, the distance and sometimes the
    nearest track change.
(b) number of porphyries (71) and epithermal deposits (24) that change track with L = 1, 2 and 5 Myr.
(c) median change of the distance to the track.
In (b) and (c), point = median of the nine models, bar = range between models
(analysis/outputs/deposit_lag_9models.csv, from deposit_lag_9models.py).
Output: results/FigS_lag.pdf/.png
Usage:  py make_FigS_lag.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CSV = os.path.join(ROOT, "analysis", "outputs", "deposit_lag_9models.csv")
OUT = os.path.join(ROOT, "results")
C_POR, C_EPI = "#b2182b", "#2166ac"


def esquema(ax):
    y = np.linspace(-3, 3, 50)
    ax.plot(0.15 * y ** 2 * 0 + 0, y, color="0.3", lw=1.4)
    ax.text(0.08, 2.75, "trench", fontsize=8, color="0.3")
    ax.fill_betweenx(y, 0, 4, color="#f2efe6", lw=0, zorder=0)
    ax.text(3.9, 2.75, "South America", fontsize=8, ha="right", color="0.35")
    for yc, ls, c, lab in ((0.0, "-", "k", r"track at $t$"), (-1.2, "--", "#7a7a7a", r"track at $t+L$")):
        ax.plot([-3.2, 0], [yc + 1.0, yc], color=c, lw=1.6, ls=ls)
        ax.plot([0, 2.6], [yc, yc - 0.55], color=c, lw=1.6, ls=ls, alpha=0.8)
        ax.text(-3.1, yc + 1.15, lab, fontsize=8, color=c)
    xd, yd = 1.9, 0.35
    ax.scatter([xd], [yd], s=55, color=C_POR, zorder=5)
    ax.text(xd + 0.12, yd + 0.15, r"deposit, age $t$", fontsize=8, color=C_POR)
    ax.annotate("", xy=(xd, -0.55 * xd / 2.6), xytext=(xd, yd), arrowprops=dict(arrowstyle="<->", color="k", lw=0.8))
    ax.annotate("", xy=(xd + 0.35, -1.2 - 0.55 * (xd + 0.35) / 2.6), xytext=(xd + 0.35, yd),
                arrowprops=dict(arrowstyle="<->", color="#7a7a7a", lw=0.8, ls="--"))
    ax.text(xd - 0.1, -0.15, r"$d_0$", fontsize=8.5, ha="right")
    ax.text(xd + 0.45, -0.9, r"$d_L$", fontsize=8.5, color="#7a7a7a")
    ax.text(-3.1, -2.85, "Nazca plate", fontsize=8, color="0.35")
    ax.set_xlim(-3.3, 4); ax.set_ylim(-3, 3)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("(a) Schematic", fontsize=10)


def main():
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9.5})
    T = pd.read_csv(CSV)
    fig, axs = plt.subplots(1, 3, figsize=(11.0, 3.4), gridspec_kw={"width_ratios": [1.3, 1, 1]})
    esquema(axs[0])
    for k, (col, ylab, tit) in enumerate((("cambian", "Deposits that change track", "(b) Change of track"),
                                          ("dd_mediana_km", "Median distance change (km)",
                                           "(c) Change of distance to the track"))):
        ax = axs[k + 1]
        for dx, tipo, c, nm in ((-0.12, "porphyries", C_POR, "porphyries (71)"),
                                (0.12, "epithermal", C_EPI, "epithermal deposits (24)")):
            g = T[T.tipo == tipo].groupby("lag")[col].agg(["median", "min", "max"])
            x = g.index.values + dx
            ax.errorbar(x, g["median"], yerr=[g["median"] - g["min"], g["max"] - g["median"]],
                        fmt="o", color=c, ms=5, capsize=3, lw=1.1, label=nm)
        ax.set_xticks([1, 2, 5]); ax.set_xlim(0, 6)
        ax.set_xlabel(r"Lag $L$ (Myr)"); ax.set_ylabel(ylab); ax.set_title(tit, fontsize=10)
        ax.set_ylim(bottom=0)
        ax.grid(True, color="#e5e5e5", lw=0.5)
    axs[1].legend(loc="upper left", frameon=False, fontsize=7.6)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"FigS_lag.{ext}"), dpi=300)
    print(T.groupby(["tipo", "lag"])[["cambian", "dd_mediana_km"]].agg(["min", "max"]).round(0))


if __name__ == "__main__":
    main()
