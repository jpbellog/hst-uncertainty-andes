# -*- coding: utf-8 -*-
"""
Figure S4 (Supplementary) -- flowchart of the complete procedure, from the
inputs (hotspot seed points, nodal covariances of the .rot, geochronology and
coordinates of the reconstructed margin) to the tables and figures of the paper.

The aim is traceability: for each table and each figure, a reader can see which
script produces it and which input it depends on, and reproduce it.

It is drawn by hand with matplotlib (boxes and arrows) so it does not depend on
graphviz or an external binary, and stays reproducible with the rest of the pipeline.

It has two elements with a purpose beyond decoration:
  - a vignette with real data in one corner (the margin, the tracks and their
    68/95% ellipses from step 08), which shows at a glance what the chain
    produces;
  - a glyph of stacked layers next to INPUTS, the classic GIS motif, so that
    the band reads at once as the inputs.

Two orientations:
  py make_FigS4_workflow.py              -> FigS4_workflow.pdf        (portrait)
  PAPER_S4_LAYOUT=landscape py ...       -> FigS4_workflow_landscape.pdf

Usage: py make_FigS4_workflow.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm"})

BASE = os.path.dirname(os.path.abspath(__file__))
PIPE = os.path.join(os.path.dirname(BASE), "models", os.environ.get("PAPER_PIPE", "Q22_WK08-A"))
MARGIN_CSV = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "margin", "Margin_SAM_0.1.csv")  # Margin_SAM (see rotation_models/README.md)
LAYOUT = os.environ.get("PAPER_S4_LAYOUT", "portrait")
LANDSCAPE = LAYOUT == "landscape"

C_IN, C_PROC, C_OUT, C_ALT = "#dbe7f3", "#e9e9e9", "#e2efd9", "#fdebd3"
EC = "#4a4a4a"
TRACK_COL = {"Easter": "#d95f02", "TT": "#7570b3", "JFR": "#1b9e77", "Foundation": "#8c510a"}
# PAPER_FIG_TAG = tag of the step 08 outputs; without it, the untagged outputs are read
# (they do not include Foundation).
TAG = os.environ.get("PAPER_FIG_TAG", "")

# ---------------------------------------------------------------- contents
# Texts describe four chains, the 03b fit, twelve models and steps 07 and 12.
# Order of IN_TXT: chains, .rot, frame covariance, margin, deposits (the last one at one side).
IN_TXT = [
    "Hotspot chains\nseed points and 77\n$^{40}$Ar/$^{39}$Ar ages, 4 chains\n[Supplementary Data S3, S4]",
    "Rotation model .rot\nnodal finite rotations\n+ inline covariances\n(rad$^2$)",
    "Frame covariances\nWK08-A (Wessel & Kroenke\n2008, Table 3), or those\nof GMHRF / GWG20",
    "Reconstructed margin\nMargin_SAM_0.1.csv\n(Arriagada et al. 2008)\nnode coords. per age",
    "Deposits and samples\n71 porphyries, 24 epithermal,\n2355 dated igneous rocks;\nrestored positions",
]
P_TXT = [
    "02_rot_engine.py\nparse .rot and its covariances; plate-circuit tree;\n"
    "stage-quaternion (slerp) interpolation",
    "03b_hotspot_fit.py\n"
    "hotspot position as an unknown (fixed or drifting),\n"
    "fitted to the ages with geological dispersion and\n"
    "outlier rejection; Foundation switch age at the\n"
    "Challenger FZ  $\\rightarrow$  hotspot_fit.json",
    "04_mc_collisions.py   $-$   MONTE CARLO,  N = 10$^3$\n"
    "perturb every covariance-bearing node (Eq. 3);\n"
    "4 sampling variants (total, circuit, frame, hotspot)\n"
    "$\\times$ 3 $f_\\mathrm{PAC}$ models; collision vs age-restored margin",
    "steps 05 / 05b / 09  $-$  convergence velocity\n"
    "analytic propagation of the stage covariances (Eq. 6)\n"
    "and Monte Carlo of the circuit $\\rightarrow$ vc_sigma_*.csv",
    "08_track_ellipses.py,  10_recon_ellipses.py,\n11_recon_model_tracks.py\n"
    "position ellipses and reconstructed tracks;\n"
    "07_deposit_tracks.py,  12_deposit_track_test.py\n"
    "nearest chain of each deposit; deposits vs background",
    "tables/make_present_bands.py\n"
    "map every collision latitude and percentile onto the\n"
    "PRESENT-DAY margin through the node correspondence\n"
    "(margin_map.py) $-$ the convention of the paper",
]
O_TXT = [
    "Tables 1, S1\nhotspot\nfits",
    "Fig. S2\nsource\ndecomposition",
    "Figs. S3, S4\nconvergence velocity\nand its uncertainty",
    "Figs. 1, 4, 6, S7\nposition ellipses,\nreconstruction,\ndeposit distances",
    "Table 2, Data S2\nFigs. 2, 3, 5, S5, S6\ncollision bands\nand comparison",
]
ALT_TXT = ("TWELVE MODELS\n\nthree plate circuits:\nQuiero et al. (2022),\nWilson & DeMets (2025),\n"
           "Somoza & Ghidella (2012)\n\neach with four frames:\nWK08-A, GMHRF-A',\nGMHRF-B, GWG20\n\n"
           "03b re-fitted in each;\nthe whole chain runs\nonce per model")
FOOT = ("Every box is a script of the supplied pipeline and every arrow a file "
        "it writes. PAPER_ROT_FILE,\nPAPER_HS_FIT and PAPER_FRAME_COV switch the "
        "whole chain between the twelve models; PAPER_PIPE\nand PAPER_FIG_TAG "
        "do the same for the table and figure scripts. Reproducing a table "
        "therefore\nrequires only its own script plus the outputs of the boxes "
        "upstream of it.")


# ------------------------------------------------------------------ helpers
def box(ax, x, y, w, h, text, fc, fs=7.6, ec=EC, lw=0.9):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.06,rounding_size=0.12",
                                facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2))
    lines = text.split("\n")
    for i, ln in enumerate(lines):
        ax.text(x + w / 2, y + h - (i + 0.5) * h / len(lines), ln,
                ha="center", va="center", fontsize=fs,
                fontweight="bold" if i == 0 else "normal", zorder=3)
    return (x, y, w, h)


def arrow(ax, p0, p1, rad=0.0, lw=1.1, ls="-", color="#333333"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=11,
                                 linewidth=lw, linestyle=ls, color=color,
                                 connectionstyle=f"arc3,rad={rad}",
                                 shrinkA=0, shrinkB=0, zorder=1))


def bot(b):
    return (b[0] + b[2] / 2, b[1])


def top(b):
    return (b[0] + b[2] / 2, b[1] + b[3])


def lft(b):
    return (b[0], b[1] + b[3] / 2)


def rgt(b):
    return (b[0] + b[2], b[1] + b[3] / 2)


def layer_stack(ax, x, y, w=0.95, h=0.17, n=4, dx=0.10, dy=0.11):
    """Glyph of stacked layers, the classic GIS motif for 'inputs'."""
    for i in range(n):
        ox, oy = x + i * dx, y + i * dy
        ax.add_patch(Polygon([[ox, oy], [ox + w, oy], [ox + w - 0.28, oy + h],
                              [ox - 0.28, oy + h]], closed=True,
                             facecolor=C_IN if i % 2 else "#c9dcee",
                             edgecolor="#2a4b73", linewidth=0.7, zorder=4))


def ellipse_xy(lon, lat, a_km, b_km, az_deg, n=64):
    t = np.linspace(0, 2 * np.pi, n)
    x, y = a_km * np.cos(t), b_km * np.sin(t)
    az = np.radians(az_deg)
    e = x * np.sin(az) + y * np.cos(az)
    nn = x * np.cos(az) - y * np.sin(az)
    return (lon + e / (111.195 * np.cos(np.radians(lat))), lat + nn / 111.195)


def vignette(ax, rect):
    """Mini-map with real data: margin, the tracks and their ellipses.

    `rect` is given in data coordinates of the diagram (x, y, width, height),
    which places it without dealing with figure fractions.
    """
    axv = ax.inset_axes(rect, transform=ax.transData)
    axv.set_facecolor("#f4f8fb")
    try:
        m = pd.read_csv(MARGIN_CSV)
        row = m.iloc[(pd.to_numeric(m["age"], errors="coerce") - 0).abs().idxmin()]
        nn = sum(1 for c in m.columns if c.startswith("Lon_"))
        lo = np.array([row.get(f"Lon_{k}") for k in range(1, nn + 1)], float)
        la = np.array([row.get(f"Lat_{k}") for k in range(1, nn + 1)], float)
        k = np.isfinite(lo) & np.isfinite(la)
        axv.plot(lo[k], la[k], "-", color="#222222", lw=1.2, zorder=4)
    except Exception as exc:                       # the map is decorative
        print(f"  warning: no margin in the vignette ({exc})")
    nom = pd.read_csv(os.path.join(PIPE, f"track_nominal_0.44{TAG}.csv"))
    ell = pd.read_csv(os.path.join(PIPE, f"track_ellipses_0.44{TAG}.csv"))
    for ridge, col in TRACK_COL.items():
        g = nom[nom.ridge == ridge].sort_values("tau")
        axv.plot(g.lon, g.lat, "-", color=col, lw=1.5, zorder=3)
        for _, r in ell[ell.ridge == ridge].iterrows():
            x9, y9 = ellipse_xy(r.lon, r.lat, r.a95_km, r.b95_km, r.az_major_deg)
            axv.fill(x9, y9, color=col, alpha=0.16, lw=0, zorder=2)
            axv.plot(x9, y9, color=col, lw=0.5, alpha=0.8, zorder=2.5)
            x6, y6 = ellipse_xy(r.lon, r.lat, r.a68_km, r.b68_km, r.az_major_deg)
            axv.fill(x6, y6, color=col, alpha=0.34, lw=0, zorder=2.6)
    axv.set_xlim(-110, -66)
    axv.set_ylim(-40, -3)
    axv.set_xticks([]); axv.set_yticks([])
    for s in axv.spines.values():
        s.set_edgecolor("#8aa6bf"); s.set_linewidth(0.8)
    axv.set_title("what the chain produces:\ntracks and their 68 / 95% ellipses",
                  fontsize=6.4, color="#2a4b73", pad=3)
    return axv


# ------------------------------------------------------------------ layouts
def draw_portrait():
    fig, ax = plt.subplots(figsize=(9.2, 12.6))
    ax.set_xlim(0, 12.4); ax.set_ylim(0, 14.8); ax.axis("off")

    layer_stack(ax, 0.15, 14.02)
    ax.text(1.55, 14.32, "INPUTS", fontsize=9, fontweight="bold", color="#2a4b73")
    iw, ih, iy = 2.80, 1.30, 12.55
    ib = [box(ax, x, iy, iw, ih, t, C_IN)
          for x, t in zip((0.10, 3.05, 6.00, 8.95), IN_TXT[:4])]

    px, pw = 3.05, 5.30
    ys = [(10.95, 0.95), (8.95, 1.55), (7.25, 1.35), (5.75, 1.10),
          (3.85, 1.50), (2.25, 1.20)]
    pb = [box(ax, px, y, pw, h, t, C_PROC) for (y, h), t in zip(ys, P_TXT)]
    alt = box(ax, 0.10, 5.20, 2.60, 3.50, ALT_TXT, C_ALT, fs=6.6)
    dep = box(ax, 0.10, 3.65, 2.60, 1.30, IN_TXT[4], C_IN, fs=6.9)

    ax.text(9.55, 10.62, "OUTPUTS", fontsize=9, fontweight="bold", color="#3d6b2e")
    oy = [(9.25, 1.00), (7.45, 1.10), (5.60, 1.30), (3.85, 1.40), (2.00, 1.45)]
    ob = [box(ax, 9.55, y, 2.75, h, t, C_OUT) for (y, h), t in zip(oy, O_TXT)]

    arrow(ax, (ib[0][0] + 0.5, ib[0][1]), (px + 0.45, pb[1][1] + pb[1][3]), rad=0.12)
    arrow(ax, bot(ib[1]), (px + 2.2, pb[0][1] + pb[0][3]))
    arrow(ax, bot(ib[2]), (px + pw - 0.9, pb[0][1] + pb[0][3]), rad=0.06)
    for i in range(5):
        arrow(ax, (px + pw / 2, pb[i][1]), (px + pw / 2, pb[i + 1][1] + pb[i + 1][3]))
    arrow(ax, (alt[0] + alt[2], alt[1] + alt[3] - 0.25),
          (px, pb[1][1] + pb[1][3] * 0.35), rad=-0.15)
    arrow(ax, rgt(dep), (px, pb[4][1] + pb[4][3] * 0.30))
    for pbx, obx in ((pb[1], ob[0]), (pb[2], ob[1]), (pb[3], ob[2]),
                     (pb[4], ob[3]), (pb[5], ob[4])):
        arrow(ax, (px + pw, pbx[1] + pbx[3] / 2), lft(obx))

    # the margin runs down the corridor between the chain and the outputs
    XC = 9.20
    ax.plot([XC, XC], [ib[3][1], pb[5][1] + pb[5][3] * 0.30],
            ls=(0, (4, 2)), lw=1.0, color="#2a4b73", zorder=1)
    for b_, f in ((pb[2], 0.22), (pb[5], 0.30)):
        arrow(ax, (XC, b_[1] + b_[3] * f), (px + pw, b_[1] + b_[3] * f),
              ls=(0, (4, 2)), color="#2a4b73", lw=1.0)

    ax.text(3.05, 1.90, "The margin coordinates enter twice: to restore the "
            "margin to the age of each\ncollision (step 04) and to map every "
            "result back onto the present-day margin.",
            fontsize=6.9, color="#2a4b73", ha="left", va="top", style="italic")
    ax.text(3.05, 1.05, FOOT, fontsize=7.0, ha="left", va="top", color="#333333")
    vignette(ax, [0.15, 0.15, 2.60, 2.05])
    return fig, "FigS4_workflow"


def draw_landscape():
    """Three columns: inputs, chain and outputs, from left to right."""
    fig, ax = plt.subplots(figsize=(14.2, 9.8))
    ax.set_xlim(0, 13.7); ax.set_ylim(0, 11.4); ax.axis("off")

    layer_stack(ax, 0.15, 10.80)
    ax.text(1.55, 11.05, "INPUTS", fontsize=9, fontweight="bold", color="#2a4b73")

    ix, iw = 0.15, 3.15
    iy = [9.42 - i * 1.46 for i in range(5)]
    ib = [box(ax, ix, y, iw, 1.26, t, C_IN, fs=7.2)
          for y, t in zip(iy, [IN_TXT[k] for k in (0, 1, 2, 4, 3)])]
    ib = [ib[k] for k in (0, 1, 2, 4, 3)]   # margin at the bottom: its arrows do not cross boxes

    cx, cw = 3.95, 5.55
    ch = [1.00, 1.40, 1.30, 1.06, 1.40, 1.18]
    cy, yy = [], 10.42
    for h in ch:
        yy -= h
        cy.append(yy)
        yy -= 0.18
    pb = [box(ax, cx, y, cw, h, t, C_PROC, fs=7.3)
          for y, h, t in zip(cy, ch, P_TXT)]

    ax.text(13.55, 10.84, "OUTPUTS", fontsize=9, fontweight="bold",
            color="#3d6b2e", ha="right")
    ox, ow = 9.95, 3.60
    oh = [1.02, 1.10, 1.30, 1.30, 1.34]
    oy, yy = [], 10.42
    for h in oh:
        yy -= h
        oy.append(yy)
        yy -= 0.24
    ob = [box(ax, ox, y, ow, h, t, C_OUT, fs=7.2)
          for y, h, t in zip(oy, oh, O_TXT)]

    alt = box(ax, 0.15, 0.98, 9.35, 1.02,
              "EIGHT MODELS\n"
              "Quiero et al. (2022) circuit with WK08-A, GMHRF-A', GMHRF-B and GWG20; "
              "Wilson & DeMets (2025) with WK08-A, GMHRF-B and GWG20;\n"
              "Somoza & Ghidella (2012) with WK08-A. 03b is re-fitted in each "
              "(SG12 shares the fit of Quiero + WK08-A); the whole chain runs once per model.",
              C_ALT, fs=7.2)

    arrow(ax, rgt(ib[0]), (cx, pb[1][1] + pb[1][3] * 0.60), rad=0.10)
    arrow(ax, rgt(ib[1]), (cx, pb[0][1] + pb[0][3] * 0.65), rad=-0.06)
    arrow(ax, rgt(ib[2]), (cx, pb[0][1] + 0.18), rad=-0.16)
    arrow(ax, rgt(ib[4]), (cx, pb[4][1] + pb[4][3] * 0.40), rad=0.06)
    for i in range(5):
        arrow(ax, (cx + cw / 2, pb[i][1]), (cx + cw / 2, pb[i + 1][1] + pb[i + 1][3]))
    for pbx, obx in ((pb[1], ob[0]), (pb[2], ob[1]), (pb[3], ob[2]),
                     (pb[4], ob[3]), (pb[5], ob[4])):
        arrow(ax, rgt(pbx), lft(obx), rad=-0.08)
    arrow(ax, top(alt), (cx + 1.2, pb[5][1]), rad=-0.12)

    # The margin enters twice, at different points. In landscape both
    # targets lie to the right of the input, so two direct arrows with
    # different curvatures are drawn, without a corridor.
    arrow(ax, rgt(ib[3]), (cx, pb[2][1] + pb[2][3] * 0.28), rad=0.25,
          ls=(0, (4, 2)), color="#2a4b73", lw=1.0)
    arrow(ax, rgt(ib[3]), (cx, pb[5][1] + pb[5][3] * 0.62), rad=-0.05,
          ls=(0, (4, 2)), color="#2a4b73", lw=1.0)
    ax.text(0.15, 3.28, "the margin coordinates are used twice: to restore the "
            "margin at each age\nand to map the results back to the present",
            fontsize=6.6, color="#2a4b73", ha="left", va="top", style="italic")

    ax.text(0.15, 0.80, FOOT.replace("\n", " "), fontsize=6.9, ha="left",
            va="top", color="#333333", wrap=True)
    vignette(ax, [9.95, 0.98, 3.60, 1.70])
    return fig, "FigS4_workflow_landscape"


fig, name = draw_landscape() if LANDSCAPE else draw_portrait()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(os.environ.get("PAPER_FIG_OUT", BASE), f"{name}.{ext}"), dpi=350, bbox_inches="tight")
plt.close(fig)
print(f"{name}: OK")
