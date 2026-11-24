# -*- coding: utf-8 -*-
"""Supplementary figure: restoration of Andean shortening with the rigid blocks of
Schepers et al. (2017) / Boschman (2021) (Utrecht model) compared with the map-view restoration
of Arriagada et al. (2008) used in the main text.

(a) porphyries: position restored with Arriagada (database v7, Latitude_ROT) and with the Utrecht
    blocks (outputs/andean_shortening/porphyries_utrecht.csv), joined by a segment.
(b) trench restored at 20 and 40 Ma with the two restorations (Margin_SAM_0.1.csv and
    Margin_SAM_0.1_utrecht.csv) and present-day trench.
Both tables are written by analysis/andean_shortening_utrecht.py.
Output: results/FigS_shortening.pdf/.png
Usage:  py make_FigS_shortening.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
U = os.path.join(ROOT, "analysis", "outputs", "andean_shortening")
DEP = os.path.join(ROOT, "inputs", "deposits")
MSAM = os.path.join(ROOT, "inputs", "margin", "Margin_SAM_0.1.csv")
OUT = os.path.join(ROOT, "results")
C_A, C_U = "#2166ac", "#b2182b"


def fosa(df, t):
    j = int(np.argmin(np.abs(df.age.values - t)))
    k = sorted(int(c[4:]) for c in df.columns if c.startswith("Lat_"))
    return (np.array([df.loc[j, f"Lon_{i}"] for i in k], float),
            np.array([df.loc[j, f"Lat_{i}"] for i in k], float))


def main():
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9.5})
    a = pd.read_csv(os.path.join(DEP, "Porphyry_Plate_ID_geodyn_v7.csv"))
    u = pd.read_csv(os.path.join(U, "porphyries_utrecht.csv"))
    assert len(a) == len(u) and (a.Cod.values == u.Cod.values).all()
    la_a, lo_a = a.Latitude_ROT.astype(float).values, a.Longitude_ROT.astype(float).values
    la_u, lo_u = u.Latitude_ROT.astype(float).values, u.Longitude_ROT.astype(float).values
    ok = np.isfinite(la_a + lo_a + la_u + lo_u)
    dk = np.hypot(la_u - la_a, (lo_u - lo_a) * np.cos(np.radians(la_a))) * 111.2
    M0, MU = pd.read_csv(MSAM), pd.read_csv(os.path.join(U, "Margin_SAM_0.1_utrecht.csv"))

    fig, axs = plt.subplots(1, 2, figsize=(6.4, 5.6), sharey=True)
    ax = axs[0]
    x0, y0 = fosa(M0, 0)
    for ax_ in axs:
        ax_.plot(x0, y0, color="0.55", lw=1.0, label="Present-day trench")
    for i in np.where(ok)[0]:
        ax.plot([lo_a[i], lo_u[i]], [la_a[i], la_u[i]], color="0.6", lw=0.6, zorder=1)
    ax.scatter(lo_a[ok], la_a[ok], s=14, facecolor="white", edgecolor=C_A, lw=0.9, zorder=2,
               label="Arriagada et al. (2008)")
    ax.scatter(lo_u[ok], la_u[ok], s=12, color=C_U, zorder=3, label="Rigid blocks (Utrecht model)")
    ax.text(0.03, 0.03, f"{ok.sum()} porphyries\nshift: median {np.nanmedian(dk[ok]):.0f} km, "
            f"maximum {np.nanmax(dk[ok]):.0f} km", transform=ax.transAxes, fontsize=7.8)
    ax.set_title("(a) Restored porphyry positions", fontsize=10)
    ax.legend(loc="upper right", frameon=False, fontsize=7.6)
    ax = axs[1]
    for t, ls in ((20, "-"), (40, "--")):
        for df, c, nm in ((M0, C_A, "Arriagada et al. (2008)"), (MU, C_U, "rigid blocks")):
            x, y = fosa(df, t)
            ax.plot(x, y, color=c, lw=1.2, ls=ls, label=f"{nm}, {t} Ma")
    ax.set_title("(b) Restored trench", fontsize=10)
    ax.legend(loc="upper right", frameon=False, fontsize=7.6)
    ylo, yhi = np.nanmin(la_a[ok]) - 4, np.nanmax(la_a[ok]) + 4
    for ax_ in axs:
        ax_.set_ylim(ylo, yhi)
        ax_.set_xlim(-82, -62)
        ax_.set_aspect(1 / np.cos(np.radians(0.5 * (ylo + yhi))))
        ax_.set_xlabel("Longitude (°)")
        ax_.grid(True, color="#e5e5e5", lw=0.5)
    axs[0].set_ylabel("Latitude (°)")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"FigS_shortening.{ext}"), dpi=300)
    print(f"n={ok.sum()} median={np.nanmedian(dk[ok]):.1f} max={np.nanmax(dk[ok]):.1f} km")


if __name__ == "__main__":
    main()
