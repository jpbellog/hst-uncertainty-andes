# -*- coding: utf-8 -*-
"""Supplementary figure: frame of Wagenaar et al. (2025) compared with the nine models of the
main text.

One panel per chain: envelope of the 68% bands of the nine models (union p16-p84, f_PAC
0.44, valid ages of each model), nominal collision of Quiero + WK08-A and nominal collision with
the WVH25 frame (Wagenaar et al. 2025 transferred to the Pacific with the Quiero circuit). Reads
the table written by analysis/sens9_wagenaar.py (outputs/sens9/wagenaar_9models.csv).
Output: results/FigS_wagenaar.pdf/.png
Usage:  py make_FigS_wagenaar.py
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
A = os.path.join(ROOT, "analysis", "outputs", "sens9")
OUT = os.path.join(ROOT, "results")
CAD = [("Easter", "Easter track", "#d95f02"), ("TT", "San Félix track", "#7570b3"),
       ("JFR", "Juan Fernández track", "#1b9e77"), ("Foundation", "Foundation track", "#8c510a")]


def main():
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9.5})
    T = pd.read_csv(os.path.join(A, "wagenaar_9models.csv"))
    fig, axs = plt.subplots(1, 4, figsize=(11.0, 3.4))
    for ax, (r, lab, col) in zip(axs, CAD):
        t = T[T.ridge == r].sort_values("age")
        ax.fill_between(t.age, t.p16, t.p84, color="#bdbdbd", alpha=0.6, lw=0,
                        label="68% bands, envelope of 9 models")
        ax.plot(t.age, t.lat_wk08, color="k", lw=1.3, label="Q22 + WK08-A")
        ax.plot(t.age, t.lat, color=col, lw=1.6, ls="--", label="Q22 + WVH25")
        ax.text(0.03, 0.04, f"inside envelope at {100 * t.dentro68.mean():.0f}% of ages",
                transform=ax.transAxes, fontsize=7.6)
        ax.set_title(lab, color=col, fontweight="bold", fontsize=10)
        ax.set_xlim(t.age.max(), 0)
        ax.set_xlabel("Age (Ma)")
        ax.grid(True, color="#e5e5e5", lw=0.5)
    axs[0].set_ylabel("Collision latitude (°)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"FigS_wagenaar.{ext}"), dpi=300)
    print(T.groupby("ridge").dentro68.mean().round(2).to_string())


if __name__ == "__main__":
    main()
