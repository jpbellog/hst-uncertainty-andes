# -*- coding: utf-8 -*-
"""Simplified deposit figure (one summary per chain instead of one line per model and all
deposits against age).

(A) Restored distance to the track: porphyries against the background igneous samples, all
    chains. Line = median of the CDFs of the nine models; fill = range between models.
    Grey band = median 68% half-width of the collision band at the porphyries (median over
    the models). Epithermal deposits are not shown (no association in most models; they are
    given in the text).
(B) By chain: median distance of the porphyries and of the background assigned to it (dot =
    median over the models, bar = range between models), with the number of porphyries and
    the number of models in which the one-sided Mann-Whitney test gives p < 0.05.

Data: g13 outputs of step 07 (deposit_tracks_*.csv) of each folder of paper_models.TEXTO,
background with the method filter of make_F4_deposit_distance (0-63.5 Ma) and
models/Q22_WK08-A/deposit_track_test_by_model.csv (step 12).
Output: figures/F4b_deposit_summary.pdf/.png
Usage:  PAPER_FIG_TAG=_hsb_cc_g13 py make_F4b_deposit_summary.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.environ.setdefault("PAPER_FIG_TAG", "_hsb_cc_g13")
import make_F4_deposit_distance as F4
import paper_models as MP

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TAG = os.environ["PAPER_FIG_TAG"]
C_POR, C_BG = "#b2182b", "#6b6b6b"
# model name in deposit_track_test_by_model.csv
NOMBRE_TEST = {("Quiero", "WK08-A"): "WK08", ("Quiero", "GMHRF-B"): "GMHRF-B",
               ("Quiero", "GWG20"): "GWG20", ("WD25", "WK08-A"): "WD25",
               ("WD25", "GMHRF-B"): "WD25+GMHRF-B", ("WD25", "GWG20"): "WD25+GWG20",
               ("SG12", "WK08-A"): "SG12", ("SG12", "GMHRF-B"): "SG12+GMHRF-B",
               ("SG12", "GWG20"): "SG12+GWG20"}
CADENAS = [("Juan Fernandez", "Juan Fernández"), ("San Felix", "San Félix"),
           ("Foundation", "Foundation"), ("Easter", "Easter")]
GRID = np.linspace(0, 900, 451)


def distancias(carpeta):
    P = os.path.join(ROOT, "models", carpeta)
    bb = F4.load_bands(os.path.join(P, f"mc_collision_bands_total{TAG}.csv"))
    por = F4.prep(pd.read_csv(os.path.join(P, "deposit_tracks_porphyries.csv")), "Age MY")
    g = pd.read_csv(os.path.join(P, "deposit_tracks_samples.csv"), low_memory=False)
    meth = g["Method"].astype(str).str.strip().str.lower()
    age = pd.to_numeric(g["Age (Ma)"], errors="coerce")
    bad = meth.isin(F4.EXCLUDED) | g["Method"].isna() | (meth == "") | (meth == "nan")
    bg = F4.prep(g[(~bad) | (age < 1)], "Age (Ma)")
    bg = bg[bg.age <= 63.5]
    h68 = np.array([np.interp(a, bb[r]["age"], bb[r]["h68"]) for r, a in zip(por.ridge, por.age)])
    return por.dist.values, bg.dist.values, float(np.median(h68))


def cdf_en_grilla(x):
    x = np.sort(x)
    return np.searchsorted(x, GRID, side="right") / len(x)


def main():
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9.5,
                         "axes.linewidth": 0.8})
    CP, CB, H68 = [], [], []
    for c in MP.CARPETAS:
        dp, db, h = distancias(c)
        CP.append(cdf_en_grilla(dp)); CB.append(cdf_en_grilla(db)); H68.append(h)
    CP, CB = np.array(CP), np.array(CB)
    h68 = float(np.median(H68))

    T = pd.read_csv(os.path.join(ROOT, "models", "Q22_WK08-A", "deposit_track_test_by_model.csv"))
    nombres = [NOMBRE_TEST[(ci, ma)] for _, ci, ma in MP.MODELOS]
    T = T[(T.tipo == "porphyries") & T.modelo.isin(nombres)]
    assert T.modelo.nunique() == len(nombres)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.3, 3.4),
                                   gridspec_kw={"width_ratios": [1.0, 1.1]})
    # ---- A ----
    axA.axvspan(0, h68, color="#d9d9d9", zorder=0, lw=0)
    axA.text(h68 / 2, 0.97, f"68% half-width of the track ({h68:.0f} km)", rotation=90,
             fontsize=7.0, color="0.3", ha="center", va="top")
    for M, col, lab in ((CP, C_POR, "Porphyry deposits"), (CB, C_BG, "Background igneous rocks")):
        axA.fill_between(GRID, M.min(0), M.max(0), color=col, alpha=0.22, lw=0, zorder=2)
        axA.plot(GRID, np.median(M, 0), color=col, lw=2.0, zorder=3, label=lab)
    axA.set_xlim(0, 900); axA.set_ylim(0, 1)
    axA.set_xlabel("Distance to the nearest track (km)")
    axA.set_ylabel("Cumulative fraction")
    axA.grid(True, color="#e5e5e5", lw=0.5); axA.set_axisbelow(True)
    axA.legend(loc="lower right", frameon=False, fontsize=7.8,
               title=f"line: median of {len(MP.CARPETAS)} models\nshading: range",
               title_fontsize=7.2)
    axA.set_title("A", loc="left", fontweight="bold", fontsize=10)

    # ---- B ----
    for i, (cad, lab) in enumerate(CADENAS):
        t = T[T.cadena == cad]
        y = len(CADENAS) - 1 - i
        for col, dcol, dy, mk in ((C_BG, "d_med_bg", -0.12, "o"), (C_POR, "d_med_dep", 0.12, "*")):
            v = t[dcol].values
            axB.hlines(y + dy, v.min(), v.max(), color=col, lw=2.2, alpha=0.55, zorder=2)
            axB.plot(np.median(v), y + dy, mk, color=col, ms=9 if mk == "*" else 6,
                     mec="white", mew=0.6, zorder=3)
        nmin, nmax = int(t.n_dep.min()), int(t.n_dep.max())
        nsig = int((t.mw_p_menor < 0.05).sum())
        ntxt = f"{nmin}" if nmin == nmax else f"{nmin}–{nmax}"
        axB.text(505, y, f"n = {ntxt}\n$p<0.05$: {nsig}/{len(t)}", fontsize=7.4, va="center",
                 ha="left", color="0.15")
    axB.set_yticks(range(len(CADENAS)))
    axB.set_yticklabels([lab for _, lab in CADENAS][::-1])
    axB.set_xlim(0, 660); axB.set_ylim(-0.6, len(CADENAS) - 0.4)
    axB.set_xticks([0, 100, 200, 300, 400, 500])
    axB.axvline(500, color="0.8", lw=0.6)
    axB.set_xlabel("Median distance to the track (km)")
    axB.grid(True, axis="x", color="#e5e5e5", lw=0.5); axB.set_axisbelow(True)
    axB.plot([], [], "*", color=C_POR, ms=9, label="Porphyries")
    axB.plot([], [], "o", color=C_BG, ms=6, label="Background")
    axB.legend(loc="lower left", bbox_to_anchor=(0.12, 1.0), ncol=2, fontsize=7.6, frameon=False,
               handletextpad=0.3, borderaxespad=0.2)
    axB.set_title("B", loc="left", fontweight="bold", fontsize=10)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fn = os.path.join(HERE, f"F4b_deposit_summary.{ext}")
        fig.savefig(fn, dpi=300 if ext == "png" else None)
        print("written:", fn)
    print(f"median 68% half-width {h68:.0f} km")
    print(T.groupby("cadena").agg(n_min=("n_dep", "min"), n_max=("n_dep", "max"),
                                  dp_min=("d_med_dep", "min"), dp_max=("d_med_dep", "max"),
                                  db_min=("d_med_bg", "min"), db_max=("d_med_bg", "max"),
                                  sig=("mw_p_menor", lambda p: int((p < 0.05).sum()))))


if __name__ == "__main__":
    main()
