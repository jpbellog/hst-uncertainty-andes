# -*- coding: utf-8 -*-
"""Supplementary figure: distance to the track by chain, nine models of the main text.

One panel per chain (Easter, San Felix, Juan Fernandez, Foundation): CDF of the distance of the
porphyries and of the background igneous samples assigned to that chain (nearest of the four
tracks, step 07), line = median of the CDFs of the nine models, fill = range between models.
Background with the same filter as step 12 (method, 0-63.5 Ma). In each panel: n of porphyries
(range between models) and number of models with one-sided Mann-Whitney p < 0.05 (from
models/Q22_WK08-A/deposit_track_test_by_model.csv).
Output: results/FigS_deposit_cdf_chain.pdf/.png
Usage:  PAPER_FIG_TAG=_hsb_cc_g13 py make_FigS_deposit_cdf_chain.py
"""
import importlib
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "models", "Q22_WK08-A"))
import paper_models as MP
from make_F4b_deposit_summary import NOMBRE_TEST

OUT = os.path.join(ROOT, "results")
CAD = [(14, "Easter", "Easter", "#d95f02"), (15, "San Felix", "San Félix", "#7570b3"),
       (13, "Juan Fernandez", "Juan Fernández", "#1b9e77"),
       (16, "Foundation", "Foundation", "#8c510a")]
C_POR, C_BG = "#b2182b", "#6b6b6b"
GRID = np.linspace(0, 900, 451)


def cargar(d, exc):
    """As cargar() of step 12: porphyries and background with chain (HS4) and distance (Dist4_km)."""
    por = pd.read_csv(os.path.join(ROOT, "models", d, "deposit_tracks_porphyries.csv"))
    g = pd.read_csv(os.path.join(ROOT, "models", d, "deposit_tracks_samples.csv"), low_memory=False)
    age = pd.to_numeric(g["Age (Ma)"], errors="coerce")
    meth = g["Method"].astype(str).str.strip().str.lower()
    bad = meth.isin(exc) | g["Method"].isna() | (meth == "") | (meth == "nan")
    bg = g[((~bad) | (age < 1)) & (age >= 0) & (age <= 63.5)]
    return [t[t.HS4.notna()] for t in (por, bg)]


def cdf(x):
    x = np.sort(np.asarray(x, float))
    return np.searchsorted(x, GRID, side="right") / len(x) if len(x) else np.full(GRID.shape, np.nan)


def main():
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9.5})
    os.chdir(os.path.join(ROOT, "models", "Q22_WK08-A"))
    exc = importlib.import_module("07_deposit_tracks").EXCLUDED
    D = {c: cargar(c, exc) for c in MP.CARPETAS}
    T = pd.read_csv(os.path.join(ROOT, "models", "Q22_WK08-A", "deposit_track_test_by_model.csv"))
    nombres = [NOMBRE_TEST[(ci, ma)] for _, ci, ma in MP.MODELOS]
    T = T[(T.tipo == "porphyries") & T.modelo.isin(nombres)]
    fig, axs = plt.subplots(1, 4, figsize=(11.0, 3.3), sharey=True)
    for ax, (pid, cad, lab, col) in zip(axs, CAD):
        for k, c, nm in ((1, C_BG, "Background igneous rocks"), (0, C_POR, "Porphyries")):
            M = np.array([cdf(D[m][k].Dist4_km[D[m][k].HS4 == pid]) for m in MP.CARPETAS])
            ax.fill_between(GRID, np.nanmin(M, 0), np.nanmax(M, 0), color=c, alpha=0.22, lw=0)
            ax.plot(GRID, np.nanmedian(M, 0), color=c, lw=1.8, label=nm)
        t = T[T.cadena == cad]
        nmin, nmax = int(t.n_dep.min()), int(t.n_dep.max())
        ntxt = f"{nmin}" if nmin == nmax else f"{nmin}–{nmax}"
        bmin, bmax = int(t.n_bg.min()), int(t.n_bg.max())
        ax.text(0.97, 0.05, f"porphyries n = {ntxt}\nbackground n = {bmin}–{bmax}\n"
                f"$p<0.05$ in {int((t.mw_p_menor < 0.05).sum())}/{len(t)} models",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7.6)
        ax.set_title(lab, color=col, fontweight="bold", fontsize=10)
        ax.set_xlim(0, 900); ax.set_ylim(0, 1)
        ax.set_xlabel("Distance to the track (km)")
        ax.grid(True, color="#e5e5e5", lw=0.5)
    axs[0].set_ylabel("Cumulative fraction")
    axs[0].legend(loc="upper left", frameon=False, fontsize=7.8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"FigS_deposit_cdf_chain.{ext}"), dpi=300)
    print(T.groupby("cadena")[["n_dep", "n_bg"]].agg(["min", "max"]))


if __name__ == "__main__":
    main()
