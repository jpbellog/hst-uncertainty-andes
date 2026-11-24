# -*- coding: utf-8 -*-
"""Test figure and table: arrival of the hotspot tracks beneath the flat-slab segments of
Wilson and DeMets (2025, Table 3) compared with the changes in Nazca/Farallon motion.

Idea: WD25 attribute the changes in Nazca/Farallon velocity to the flattening and steepening of
slab segments; here we test whether the tracks of this work were beneath those segments when the
changes occurred.

Inputs
  - data/wd25_table3_flat_slabs.csv (Table 3 of WD25, transcribed from the PDF).
  - collision bands at the present-day margin, total, f_PAC 0.44, of the nine models of the main text
    (figures/paper_models.py), clipped by the frame_coverage of the frame, in results/.

Outputs (results/)
  - flat_slab_tracks_test.png/.pdf: collision latitude against age; 95% envelope and 68% union
    of the nine models per chain; segments as latitude-time bands (upper half Ramos and
    Folguera 2009, lower half Martinod et al. 2020; light shade = partial or ongoing flattening,
    dark = flat); WD25 motion changes as triangles.
  - flat_slab_tracks_test.csv: per chain and segment, age intervals in which the 68% union and
    the 95% envelope overlap the latitude band of the segment, and in which the medians of
    all valid models fall inside it.
Usage:  py flat_slab_tracks.py
"""
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D

sys.stdout.reconfigure(encoding="utf-8")
A = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(A, ".."))
sys.path.insert(0, os.path.join(P, "figures"))
import paper_models as MP

G13 = os.path.join(P, "results")
TAG, PF, AGE_MAX = "_hsb_cc_g13", "0.44", 50.0
LAT_MIN, LAT_MAX = -40.0, -5.0     # figure up to 5 S
TRK = {"Easter": ("#d95f02", "Easter"), "TT": ("#7570b3", "San Félix"),
       "JFR": ("#1b9e77", "Juan Fernández"), "Foundation": ("#8c510a", "Foundation")}
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9.5})


def nums(txt):
    """Numbers in a cell ('40–32', '14*,8', '>33', '∼8.3', '<3'); [] if N/A or '?'."""
    return [float(x) for x in re.findall(r"\d+(?:\.\d+)?", str(txt))]


def bandas():
    dd = []
    for c in MP.CARPETAS:
        b = pd.read_csv(os.path.join(G13, f"mc_collision_bands_total_present"
                                          f"{MP.sufijo_bandas(c)}{TAG}.csv"), dtype={"pf": str})
        b = MP.recortar(b, c)
        dd.append(b[(b.pf == PF) & (b.age <= AGE_MAX)].assign(modelo=c))
    return pd.concat(dd)


def corridas(edades, mask):
    """Contiguous intervals (1 Myr step) of ages where mask is true, as 'a-b'."""
    e = np.asarray(edades)[np.asarray(mask)]
    if not len(e):
        return ""
    out, ini, prev = [], e[0], e[0]
    for x in e[1:]:
        if x - prev > 1.01:
            out.append((ini, prev)); ini = x
        prev = x
    out.append((ini, prev))
    return "; ".join(f"{a:g}" if a == b else f"{a:g}-{b:g}" for a, b in out)


def intervalos_planos(flat_txt, steep_txt):
    """(start_age, end_age, dark) of what is drawn. Flattening as a range (40-32): light segment
    40-32; steepening as a range (27-18): light segment 27-18; between them, dark. Cell with an
    asterisk: the first number is reduced or inland volcanism (light up to the second
    number or up to the steepening)."""
    f, s = nums(flat_txt), nums(steep_txt)
    if not f:
        return []
    ast = "*" in str(flat_txt)
    s_young = min(s) if s else 0.0
    s_old = max(s) if s else 0.0
    if ast and len(f) == 2:                      # '14*,8': light 14-8, dark 8-steep
        return [(f[0], f[1], False), (f[1], s_old, True)]
    if ast:                                      # '16*' with '5*': all light
        return [(f[0], s_young, False)]
    out = []
    if len(f) == 2:                              # flattening range
        out.append((max(f), min(f), False))
    out.append((min(f), s_old, True))
    if len(s) == 2:                              # steepening range
        out.append((s_old, s_young, False))
    return out


def main():
    T = pd.read_csv(os.path.join(A, "data", "wd25_table3_flat_slabs.csv"), comment="#")
    d = bandas()
    fig, ax = plt.subplots(figsize=(8.6, 6.6))
    # segments (the figure shows only those south of LAT_MAX; the table has all of them)
    for _, s in T[T.lat_n <= LAT_MAX].iterrows():
        lo, hi = s.lat_s, s.lat_n
        mid = (lo + hi) / 2
        ax.axhspan(lo, hi, color="0.93", zorder=0, lw=0)
        for (txt_f, txt_s, y0, y1) in ((s.rf_flat_txt, s.rf_steep_txt, mid, hi),
                                       (s.ma_flat_txt, s.ma_steep_txt, lo, mid)):
            for a0, a1, osc in intervalos_planos(txt_f, txt_s):
                ax.add_patch(Rectangle((a1, y0), a0 - a1, y1 - y0, facecolor="0.55" if osc else "0.80",
                                       edgecolor="none", zorder=0.5))
        ax.text(-0.4, mid, s.segmento, fontsize=7.6, va="center", ha="left", color="0.25",
                clip_on=False)
        for txt, mk, fc in ((s.decel_txt, "v", "black"), (s.accel_txt, "^", "white")):
            v = nums(txt)
            if v:
                ax.plot(v[0], mid, mk, ms=8, mfc=fc, mec="black", mew=1.0, zorder=8)
                if ">" in str(txt):
                    ax.annotate("", (min(v[0] + 3, AGE_MAX), mid), (v[0], mid),
                                arrowprops=dict(arrowstyle="->", lw=0.9), zorder=8)
    # tracks
    rows = []
    for r, (col, lab) in TRK.items():
        x = d[d.ridge == r].dropna(subset=["lat_p50"])
        g = x.groupby("age")
        e95 = pd.DataFrame({"lo": g.lat_p2_5.min(), "hi": g.lat_p97_5.max()})
        e68 = pd.DataFrame({"lo": g.lat_p16.min(), "hi": g.lat_p84.max(),
                            "mlo": g.lat_p50.min(), "mhi": g.lat_p50.max(),
                            "n": g.modelo.nunique()})
        ax.fill_between(e95.index, e95.lo, e95.hi, color=col, alpha=0.18, lw=0, zorder=2)
        ax.fill_between(e68.index, e68.lo, e68.hi, color=col, alpha=0.45, lw=0, zorder=3)
        for _, t in x.groupby("modelo"):
            t = t.sort_values("age")
            ax.plot(t.age, t.lat_p50, color=col, lw=0.6, alpha=0.9, zorder=4)
        for _, s in T.iterrows():
            lo, hi = s.lat_s, s.lat_n
            ov68 = (e68.hi >= lo) & (e68.lo <= hi)
            ov95 = (e95.hi >= lo) & (e95.lo <= hi)
            med = (e68.mlo >= lo) & (e68.mhi <= hi)
            if ov95.any():
                rows.append(dict(cadena=lab, segmento=s.segmento, lat=s.lat_txt,
                                 solape_68=corridas(e68.index, ov68),
                                 solape_95=corridas(e95.index, ov95),
                                 todas_medianas_dentro=corridas(e68.index, med),
                                 modelos_validos_min=int(e68.n[ov95].min()),
                                 rf_plano=f"{s.rf_flat_txt} / {s.rf_steep_txt}",
                                 martinod_plano=f"{s.ma_flat_txt} / {s.ma_steep_txt}",
                                 wd25_decel=s.decel_txt, wd25_accel=s.accel_txt))
    ax.set_xlim(AGE_MAX, 0)
    ax.set_ylim(LAT_MIN, LAT_MAX)
    ax.set_xlabel("Age (Ma)")
    ax.set_ylabel("Collision latitude, present-day margin (°)")
    ax.grid(True, color="#e5e5e5", lw=0.5, zorder=0)
    h = ([Patch(facecolor=c, alpha=0.45, label=f"{l}: 68% envelope") for c, l in TRK.values()]
         + [Patch(facecolor="0.5", alpha=0.18, label="95% envelope (same colours)"),
            Line2D([0], [0], color="0.3", lw=0.6, label=f"median of each of the {len(MP.CARPETAS)} models"),
            Patch(facecolor="0.55", label="flat slab (upper half: Ramos & Folguera 2009;"),
            Patch(facecolor="0.80", label="onset, steepening or partial; lower half: Martinod et al. 2020)"),
            Line2D([0], [0], marker="v", ls="", mfc="black", mec="black", ms=8,
                   label="deceleration (Wilson & DeMets 2025)"),
            Line2D([0], [0], marker="^", ls="", mfc="white", mec="black", ms=8,
                   label="acceleration (Wilson & DeMets 2025)")])
    fig.tight_layout()
    fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=2, frameon=False,
               fontsize=7.6)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(G13, f"flat_slab_tracks_test.{ext}"), dpi=300,
                    bbox_inches="tight")
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(G13, "flat_slab_tracks_test.csv"), index=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 40):
        print(R.to_string(index=False))


if __name__ == "__main__":
    main()
