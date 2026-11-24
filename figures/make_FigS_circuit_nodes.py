# -*- coding: utf-8 -*-
"""Supplementary figure: nodes of the plate circuit, their source and their covariance.
Uses the current .rot files (GTS2020, Granot et al. 2013, Cande et al. 2000 at 43.5 Ma) and
includes the SG12 row (South America-Nazca of Somoza and Ghidella 2012). The branch-change
ages (active interval of each link) are those of the .rot files.
Age on the horizontal axis (0-84 Ma), links on the vertical axis in circuit order. Each node
of the .rot is a square in the colour of its source (comment of the .rot); the fill shows the
covariance used by the Monte Carlo (PAPER_CIRCUIT_COV=b, circuit_covariance.aplicar).

Usage:  py make_FigS_circuit_nodes.py
Outputs: results/FigS_circuit_nodes.{png,pdf} and
         analysis/outputs/circuit_nodes_g13.csv
"""
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

AQUI = os.path.dirname(os.path.abspath(__file__))
P = (__import__("os").environ["PAPER_ROOT"])
sys.path.insert(0, P + "/models/Q22_WK08-A")
import importlib
eng = importlib.import_module("02_rot_engine")
sys.path.insert(0, P + "/analysis")
import circuit_covariance as cc

ROT = {"base": P + "/models/Q22_WK08-A/rot/Global_Model_R_BG_2025_hsb.rot",
       "WD25": P + "/models/WD25_WK08-A/rot/Global_Model_R_BG_2025_WD25_hsb.rot",
       "SG12": P + "/models/SG12_WK08-A/rot/Global_Model_R_BG_2025_SG12_hsb.rot"}
TMAX = 84.0
TOPE = 71.07   # upper age limit of the circuit with covariance

# (label, moving plate, fixed plate, .rot model, active interval in the NZ -> SAM chain)
ESLABONES = [  # branch changes read from the .rot files
    ("Nazca–Pacific (Q22, SG12)", 911, 901, "base", (0, 21.806)),
    ("Nazca–Pacific (WD25)", 911, 901, "WD25", (0, 33.726)),
    ("Farallon–Pacific (Q22, SG12)", 902, 901, "base", (21.806, TMAX)),
    ("Farallon–Pacific (WD25)", 902, 901, "WD25", (33.726, TMAX)),
    ("Pacific–West Antarctica (Q22, SG12)", 901, 804, "base", (0, 44.846)),
    ("Pacific–West Antarctica (WD25)", 901, 804, "WD25", (0, 44.846)),
    ("West–East Antarctica", 804, 802, "base", (0, 44.846)),
    ("Pacific–Lord Howe Rise", 901, 833, "base", (44.846, 83.65)),
    ("Lord Howe Rise–Australia", 833, 801, "base", (44.846, TMAX)),
    ("Australia–East Antarctica", 801, 802, "base", (44.846, TMAX)),
    ("East Antarctica–Nubia", 802, 701, "base", (0, 28.278)),
    ("East Antarctica–Somalia", 802, 709, "base", (28.278, TMAX)),
    ("Somalia–Nubia", 709, 701, "base", (28.278, TMAX)),
    ("South America–Nubia (Q22, WD25)", 201, 701, "base", (0, TMAX)),
    ("South America–Nazca (SG12)", 201, 911, "SG12", (0, 71.567)),
    ("South America–Nubia (SG12)", 201, 701, "SG12", (71.567, TMAX)),
]

# source from the text of the .rot; order matters (technical nodes first)
FUENTES = [
    (r"(?i)interpolated|interpolada|cruce|No relative movements|crossover|NZ=FAR|FAR=NZ|^!\s*$", "Technical node (crossover, interpolated, no motion)"),
    (r"Wilson & DeMets", "Wilson and DeMets 2025"),
    (r"Somoza", "Somoza and Ghidella 2012"),
    (r"@REF \?\?\?|@REF: \?\?\?", "Quiero et al. 2022"),
    (r"Rowan", "Rowan and Rowley 2014"),
    (r"Croon", "Croon et al. 2008"),
    (r"Sutherland|LHR locked to PAC", "Sutherland 1995 (Lord Howe Rise fixed to Pacific)"),
    (r"Gaina", "Gaina et al. 1998"),
    (r"Cande & Stock", "Cande and Stock 2004"),
    (r"Tikku", "Tikku and Cande 2000"),
    (r"Granot_\+\+_2013", "Granot et al. 2013"),
    (r"Granot", "Granot and Dyment 2018"),
    (r"Cande et al, Nature 2000", "Cande et al. 2000"),
    (r"Cande et al\. \(2010\)", "Cande et al. 2010"),
    (r"Bernnard", "Bernard et al. 2005"),
    (r"DeMets & Merkouriev", "DeMets and Merkouriev 2019"),
    (r"Horner-Johnson", "Horner-Johnson et al. 2007"),
    (r"Lemaux", "Lemaux et al. 2002"),
    (r"Mueller", "Müller et al. 1999"),
]
TECNICO = FUENTES[0][1]
GRIS = "#8a8984"
# colours: 8 validated slots of the reference palette + 8 dark steps of
# the same ramps; identity does not depend on colour alone (each run carries
# the name of its source written above the row)
COLOR = {"Quiero et al. 2022": "#2a78d6", "Wilson and DeMets 2025": "#eb6834",
         "Rowan and Rowley 2014": "#1baf7a", "Croon et al. 2008": "#eda100",
         "Granot and Dyment 2018": "#e87ba4", "Cande et al. 2000": "#008300",
         "Sutherland 1995 (Lord Howe Rise fixed to Pacific)": "#6250d6", "Gaina et al. 1998": "#e34948",
         "Cande and Stock 2004": "#5f6b2a", "Tikku and Cande 2000": "#8a5d00",
         "DeMets and Merkouriev 2019": "#184f95", "Cande et al. 2010": "#0f6e4c",
         "Bernard et al. 2005": "#a3355f", "Horner-Johnson et al. 2007": "#3b2e8a",
         "Lemaux et al. 2002": "#b02f2e", "Müller et al. 1999": "#a8431a",
         "Granot et al. 2013": "#9c27b0", "Somoza and Ghidella 2012": "#00838f"}
PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
          "#6250d6", "#e34948", "#184f95", "#a8431a", "#0f6e4c", "#8a5d00",
          "#a3355f", "#3b2e8a", "#5f6b2a", "#b02f2e"]

COV_CLASES = ["covariance stored in the rotation file", "published covariance added here",
              "covariance derived or fitted here", "no covariance"]


def fuente(texto):
    for pat, nom in FUENTES:
        if re.search(pat, texto):
            return nom
    return "?"


def filas_rot(fn):
    out = []
    for l in open(fn, errors="replace"):
        s = l.split()
        if len(s) < 6 or not s[0].isdigit():
            continue
        try:
            age = float(s[1])
        except ValueError:
            continue
        fx = re.match(r"\d+", s[5])
        if not fx:
            continue
        com = l[l.find("!"):] if "!" in l else ""
        out.append((int(s[0]), age, int(fx.group()), com.strip()))
    return out


def clases_cov(fn):
    """{(plate, age, fixed): class} from circuit_covariance.aplicar (variant b)."""
    M0 = eng.RotModel(fn)
    M = eng.RotModel(fn)
    ent = cc.aplicar(M, hellinger="b")
    # the ages in ent are on the timescale of each source; _find gives the GTS2020 node of the .rot
    agregada = {(p, round(cc._find(M, p, a, f)[1][0], 3), f):
                ("derivada o ajustada" if ("deriv" in s + n or "interpol" in s + n
                                           or "Hellinger" in s + n or "propagada" in s + n) else "publicada")
                for (p, a, f, C, s, n, m) in ent}
    out = {}
    for p, rows in M.table.items():
        for (age, la, lo, ang, fx, Cv) in rows:
            k = (p, round(age, 3), fx)
            orig = [r for r in M0.table[p] if r[0] == age and r[4] == fx]
            if orig and orig[0][5] is not None:
                out[k] = COV_CLASES[0]
            elif k in agregada:
                out[k] = COV_CLASES[1] if agregada[k] == "publicada" else COV_CLASES[2]
            else:
                out[k] = COV_CLASES[3]
    return out


def main():
    filas = {m: filas_rot(fn) for m, fn in ROT.items()}
    clases = {m: clases_cov(fn) for m, fn in ROT.items()}
    reg = []
    for i, (lab, p, fx, mod, act) in enumerate(ESLABONES):
        for (pp, age, ff, com) in filas[mod]:
            if pp != p or ff != fx or age > TMAX or age == 0.0:
                continue
            reg.append(dict(fila=i, eslabon=lab, movil=p, fija=fx, modelo=mod, edad=age,
                            fuente=fuente(com),
                            cov=clases[mod].get((p, round(age, 3), fx), COV_CLASES[3]),
                            activo=act[0] <= age <= act[1], comentario_rot=com[:120]))
    D = pd.DataFrame(reg)
    D.to_csv(os.path.join(P, "analysis", "outputs", "circuit_nodes_g13.csv"), index=False)
    if (D.fuente == "?").any():
        raise RuntimeError("unclassified source:\n" + D[D.fuente == "?"].to_string())

    orden = [f for f in D.drop_duplicates("fuente").fuente if f != TECNICO]
    color = dict(COLOR, **{TECNICO: GRIS})
    if set(orden) - set(COLOR):
        raise RuntimeError(f"source without colour: {set(orden) - set(COLOR)}")

    n = len(ESLABONES)
    fig, ax = plt.subplots(figsize=(12.5, 0.32 * n + 1.9))   # landscape (horizontal page of the supplement)
    ink, muted = "#2b2b2a", "#6b6a64"
    for i, (lab, p, fx, mod, (a0, a1)) in enumerate(ESLABONES):
        y = n - 1 - i
        ax.add_patch(Rectangle((a0, y - 0.07), a1 - a0, 0.14, color="#d9d8d2", lw=0, zorder=1))
        g = D[D.fila == i].sort_values("edad")
        for r in g.itertuples():
            c = color[r.fuente]
            alpha = 1.0 if r.activo else 0.35
            kw = dict(marker="s", ms=6.5, zorder=3, alpha=alpha, ls="none")
            if r.cov == COV_CLASES[0]:
                ax.plot(r.edad, y, mfc=c, mec=c, mew=1.0, **kw)
            elif r.cov == COV_CLASES[1]:
                ax.plot(r.edad, y, mfc=c, mec=ink, mew=1.6, **kw)
            elif r.cov == COV_CLASES[2]:
                ax.plot(r.edad, y, mfc="white", mec=c, mew=1.6, **kw)
                ax.plot(r.edad, y, marker="x", ms=4.2, mec=c, mew=1.6, zorder=4,
                        alpha=alpha, ls="none")
            else:
                ax.plot(r.edad, y, mfc="white", mec=c, mew=1.6, **kw)
        # source name above each contiguous run (technical nodes not labelled)
        g = g.reset_index(drop=True)
        run = (g.fuente != g.fuente.shift()).cumsum()
        fin_ant = -99.0   # age where the previous label of the row ends
        for _, h in g.groupby(run):
            f = h.fuente.iloc[0]
            if f == TECNICO:
                continue
            # always above the row; if it overlaps the previous label, it moves to the right
            # (approximate width of 0.46 Myr per character at fontsize 7.5 in this figure)
            x = max(h.edad.min() - 0.3, fin_ant + 0.8)
            izq = x < TMAX - 12
            ax.text(min(x, TMAX - 0.3), y + 0.24, f, fontsize=7.5, color=muted,
                    ha="left" if izq else "right", va="bottom", zorder=5)
            fin_ant = x + 0.46 * len(f) if izq else TMAX
        if (p, fx) == (709, 701):
            ax.text(TMAX - 0.3, y + 0.24, "rotation of Lemaux et al. 2002 held constant to the "
                    "250 Ma node (published covariance added) →", fontsize=7.5, color=muted,
                    ha="right", va="bottom", zorder=5)
        ax.text(-0.8, y, lab, ha="right", va="center", fontsize=9.5, color=ink)
    ax.axvline(TOPE, color=ink, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.text(TOPE + 0.3, n - 0.35, f"oldest circuit node\nwith covariance,\n{TOPE} Ma (C32y)", fontsize=8,
            color=ink, va="top")
    ax.set_xlim(0, TMAX)
    ax.set_ylim(-0.7, n - 0.2)
    ax.set_yticks([])
    ax.set_xticks(np.arange(0, TMAX + 1, 5))
    ax.set_xticks(np.arange(0, TMAX + 1, 1), minor=True)
    ax.grid(axis="x", which="major", color="#ecebe6", lw=0.8, zorder=0)
    ax.set_xlabel("Age (Ma)", fontsize=10, color=ink)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(muted)
    ax.tick_params(colors=muted, labelsize=9)

    # legends: sources (colour) and covariance (fill)
    from matplotlib.lines import Line2D
    hf = [Line2D([], [], marker="s", ls="none", ms=8, mfc=color[f], mec=color[f], label=f)
          for f in orden + [TECNICO]]
    hc = [Line2D([], [], marker="s", ls="none", ms=8, mfc="#555", mec="#555", label=COV_CLASES[0]),
          Line2D([], [], marker="s", ls="none", ms=8, mfc="#555", mec=ink, mew=1.6,
                 label=COV_CLASES[1]),
          Line2D([], [], marker="X", ls="none", ms=8, mfc="white", mec="#555", mew=1.2,
                 label=COV_CLASES[2]),
          Line2D([], [], marker="s", ls="none", ms=8, mfc="white", mec="#555", mew=1.6,
                 label=COV_CLASES[3]),
          Rectangle((0, 0), 1, 1, color="#d9d8d2", label="link active in the circuit"),
          Line2D([], [], marker="s", ls="none", ms=8, mfc="#555", mec="#555", alpha=0.35,
                 label="node outside the active interval")]
    l1 = fig.legend(handles=hf, title="Source of the rotation", loc="lower left",
                    bbox_to_anchor=(0.0, 0.005), ncol=4, fontsize=7.4, title_fontsize=8.5,
                    frameon=False)
    fig.legend(handles=hc, title="Covariance", loc="lower right", bbox_to_anchor=(1.0, 0.005),
               ncol=1, fontsize=7.4, title_fontsize=8.5, frameon=False)
    fig.subplots_adjust(left=0.2, right=0.99, top=0.985, bottom=0.25)
    out = os.path.join(P, "results", "FigS_circuit_nodes")
    fig.savefig(out + ".png", dpi=200)
    fig.savefig(out + ".pdf")
    print(out + ".png")
    print(D.groupby(["eslabon", "cov"], sort=False).size().to_string())


if __name__ == "__main__":
    main()
