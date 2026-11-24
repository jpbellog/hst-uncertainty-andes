# -*- coding: utf-8 -*-
"""
Figure 5 -- Comparison of the collision-latitude bands obtained with the two
versions of the Nazca-Pacific branch of the plate circuit.

Top row, one panel per track: 68% Monte Carlo band and median of the reference
version (Wilder 2003 via Quiero et al. 2022) and of the revised version
(Wilson and DeMets 2025), reference spreading-asymmetry model f_PAC = 0.44.

Bottom row: difference between the two medians at each age, drawn over the
+/- 68% and +/- 95% half-widths of the reference band. The point of the figure
is that the difference curve stays inside the reference uncertainty at almost
every age, i.e. the two kinematics are not distinguishable at the resolution of
the budget.

REAL DATA: tables/mc_collision_bands_total_present.csv       (reference)
           tables/mc_collision_bands_total_present_wd25.csv  (revised)
both from pipeline step 04 (N = 1000 per f_PAC model, joint sampling of circuit
+ WK08-A + geochronology), re-restored to the present-day margin by
tables/make_present_bands.py.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 10.5,
    "axes.linewidth": 0.8,
})

BASE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(os.path.dirname(BASE), "tables")
PF = "0.44"
AGE_MAX = 45 if os.environ.get("PAPER_MODELOS") == "supl" else 66   # GMHRF-A' resolves up to ~45 Ma
XTICKS = [40, 30, 20, 10, 0] if AGE_MAX == 45 else [60, 40, 20, 0]

# Reference keeps the f_PAC = 0.44 blue of Fig. 2, so the two figures read
# together; the revised model gets the orange. The pair is colourblind-safe
# (worst adjacent dE 22.2 protan, 31.7 normal vision).
C_REF, C_NEW = "#2b5f9e", "#d95f02"
C_BAND = "#9a9a9a"

PANELS = [("Easter", "Easter track"),
          ("TT", "San Félix track"),
          ("JFR", "Juan Fernández track"),
          ("Foundation", "Foundation track")]
# PAPER_FIG_TAG = tag of the bands (e.g. _hsb_cc_g13), PAPER_TABLES_DIR = folder of the bands,
# PAPER_FIG_OUT = output folder. Without them, the two-circuit comparison below is drawn.
# Foundation (fourth chain) is included only if it is present in the bands.
TAG = os.environ.get("PAPER_FIG_TAG", "")
TAB = os.environ.get("PAPER_TABLES_DIR", TAB)
FIG_OUT = os.environ.get("PAPER_FIG_OUT", BASE)



def load(sfx):
    d = pd.read_csv(os.path.join(TAB, f"mc_collision_bands_total_present{sfx}{TAG}.csv"),
                    dtype={"pf": str})
    return d[d.pf == PF]


# With PAPER_FIG_TAG the figure compares all models of the set (paper_models.py) equally.
# Top: median of each model (colour = frame, line style = circuit). Bottom: spread of the
# medians between models (total; circuit only = mean over frames of the range between
# circuits; frame only = mean over circuits of the range between frames) against the 68% and
# 95% band widths (median over the models).
FRAME_COL = {"WK08-A": "#2b5f9e", "GMHRF-B": "#d95f02", "GWG20": "#7570b3",
             "GMHRF-A'": "#d95f02"}       # validated palette (CVD dE min 21.4)
C_SCIRC, C_SMARCO = "#1b9e77", "#c51b8a"   # CVD dE 7.8: acceptable with a different line style
CIRC_LS = {"Quiero": "-", "WD25": (0, (5, 2)), "SG12": (0, (1, 1.3))}
CIRC_LAB = {"Quiero": "Quiero et al. (2022)", "WD25": "Wilson & DeMets (2025)",
            "SG12": "Somoza & Ghidella (2012)"}


def make_modelos():
    import paper_models as MP
    # frame_coverage of each model from its frame covariance
    d = pd.concat([MP.recortar(load(MP.sufijo_bandas(c)), c).assign(modelo=c, circ=ci, marco=fr)
                   for c, ci, fr in MP.MODELOS])
    d = d[d.age <= AGE_MAX]
    panels = [p for p in PANELS if p[0] in set(d.ridge)]
    fig, axes = plt.subplots(2, len(panels), figsize=(9.2 / 3 * len(panels), 5.9),
                             sharex=True, gridspec_kw={"height_ratios": [2.1, 1.0]})
    rows = []
    for j, (key, title) in enumerate(panels):
        r = d[d.ridge == key]
        ax = axes[0, j]
        for (c, ci, fr), t in r.groupby(["modelo", "circ", "marco"]):
            t = t.sort_values("age")
            ax.plot(t.age, t.lat_p50, color=FRAME_COL[fr], ls=CIRC_LS[ci], lw=1.1)
        ax.set_title(title, fontsize=10.5, loc="left")
        ax.grid(True, color="#dddddd", lw=0.5)
        ax.set_axisbelow(True)
        if j == 0:
            ax.set_ylabel("Collision latitude (°)")

        ax = axes[1, j]
        p = r.pivot_table(index="age", columns=["circ", "marco"], values="lat_p50")
        # range only where at least two models are valid (with the frame coverage, beyond
        # ~20-40 Ma only the WK08-A models remain and the range between frames is undefined)
        rng = lambda g: (g.max(axis=1) - g.min(axis=1)).where(g.notna().sum(axis=1) >= 2)
        tot = rng(p)
        marcos = sorted(set(r.marco)); circs = sorted(set(r.circ))
        s_circ = (pd.concat([rng(p.xs(m, axis=1, level="marco")) for m in marcos], axis=1)
                  .mean(axis=1) if len(circs) > 1 else None)
        s_marco = (pd.concat([rng(p.xs(c, axis=1, level="circ")) for c in circs], axis=1)
                   .mean(axis=1) if len(marcos) > 1 else None)
        w68 = (r.assign(w=r.lat_p84 - r.lat_p16).groupby("age").w.median())
        w95 = (r.assign(w=r.lat_p97_5 - r.lat_p2_5).groupby("age").w.median())
        ax.fill_between(w68.index, 0, w68.values, color=C_BAND, alpha=0.30, lw=0)
        ax.plot(w95.index, w95.values, color=C_BAND, lw=0.9, ls=(0, (4, 2)))
        ax.plot(tot.index, tot.values, color="#222222", lw=1.5)
        if s_circ is not None:
            ax.plot(s_circ.index, s_circ.values, color=C_SCIRC, lw=1.2, ls=(0, (6, 1.5, 1, 1.5)))
        if s_marco is not None:
            ax.plot(s_marco.index, s_marco.values, color=C_SMARCO, lw=1.2)
        ax.set_ylim(0, None)
        ax.grid(True, color="#dddddd", lw=0.5)
        ax.set_axisbelow(True)
        ax.set_xlabel("Age (Ma)")
        ax.set_xticks(XTICKS)
        ax.set_xlim(AGE_MAX, 0)   # after the ticks: set_xticks widens the axis
        if j == 0:
            ax.set_ylabel("Spread or width (°)")
        k = tot.dropna().index.intersection(w68.dropna().index)
        rows.append(dict(cadena=key, rango_med_deg=float(tot.median()),
                         rango_max_deg=float(tot.max()),
                         edad_max=float(tot.idxmax()),
                         frac_rango_le_w68=float(np.mean(tot.loc[k] <= w68.loc[k])),
                         circ_med_deg=float(s_circ.median()) if s_circ is not None else np.nan,
                         marco_med_deg=float(s_marco.median()) if s_marco is not None else np.nan,
                         w68_med_deg=float(w68.median())))
    ymax = max(ax.get_ylim()[1] for ax in axes[1, :])
    for ax in axes[1, :]:
        ax.set_ylim(0, ymax)
    for ax in axes[1, 1:]:
        ax.tick_params(labelleft=False)
    frames = list(dict.fromkeys(m[2] for m in MP.MODELOS))
    circs = list(dict.fromkeys(m[1] for m in MP.MODELOS))
    h_top = ([Line2D([0], [0], color=FRAME_COL[f], lw=1.6, label=f) for f in frames]
             + [Line2D([0], [0], color="0.3", lw=1.2, ls=CIRC_LS[c], label=CIRC_LAB[c])
                for c in circs])
    h_bot = [Line2D([0], [0], color="#222222", lw=1.5,
                    label=f"range of the {len(MP.MODELOS)} medians")]
    if len(circs) > 1:
        h_bot.append(Line2D([0], [0], color=C_SCIRC, lw=1.2, ls=(0, (6, 1.5, 1, 1.5)),
                            label="range between circuits (mean over frames)"))
    if len(frames) > 1:
        h_bot.append(Line2D([0], [0], color=C_SMARCO, lw=1.2,
                            label="range between frames (mean over circuits)"))
    h_bot += [Patch(facecolor=C_BAND, alpha=0.30, label="68% band width (median of models)"),
              Line2D([0], [0], color=C_BAND, lw=0.9, ls=(0, (4, 2)),
                     label="95% band width (median of models)")]
    fig.tight_layout()
    # two legends, one per row, so the circuit line styles are not confused with those below
    fig.legend(handles=h_top, loc="upper left", bbox_to_anchor=(0.06, 0.0), frameon=False,
               fontsize=7.6, ncol=2, columnspacing=1.4, alignment="left",
               title=f"Top: median collision latitude, $f_\\mathrm{{PAC}}={PF}$\n"
                     "(colour = hotspot frame, line = plate circuit)", title_fontsize=7.8)
    fig.legend(handles=h_bot, loc="upper left", bbox_to_anchor=(0.53, 0.0), frameon=False,
               fontsize=7.6, ncol=1, alignment="left",
               title="Bottom: spread between models and band widths", title_fontsize=7.8)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_OUT, f"F5_model_comparison{MP.SUFIJO_SALIDA}.{ext}"),
                    dpi=300, bbox_inches="tight")
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(FIG_OUT, f"F5_model_comparison{MP.SUFIJO_SALIDA}_numbers.csv"),
             index=False, float_format="%.3f")
    print(R.round(2).to_string(index=False))
    print(f"written: F5_model_comparison{MP.SUFIJO_SALIDA}.pdf/.png")


if TAG:
    make_modelos()
    raise SystemExit

ref, new = load(""), load(os.environ.get("PAPER_F5_NEW_SFX", "_wd25"))
PANELS = [p for p in PANELS if p[0] in set(ref.ridge) & set(new.ridge)]

fig, axes = plt.subplots(2, len(PANELS), figsize=(9.2 / 3 * len(PANELS), 5.6), sharex=True,
                         gridspec_kw={"height_ratios": [2.1, 1.0]})

for j, (key, title) in enumerate(PANELS):
    a = ref[ref.ridge == key].sort_values("age").set_index("age")
    b = new[new.ridge == key].sort_values("age").set_index("age")
    ages = a.index.intersection(b.index).values

    # ---- top: the two bands ----
    ax = axes[0, j]
    ax.fill_between(ages, a.loc[ages, "lat_p16"], a.loc[ages, "lat_p84"],
                    color=C_REF, alpha=0.26, lw=0)
    ax.fill_between(ages, b.loc[ages, "lat_p16"], b.loc[ages, "lat_p84"],
                    color=C_NEW, alpha=0.26, lw=0)
    ax.plot(ages, a.loc[ages, "lat_p50"], color=C_REF, lw=1.5)
    ax.plot(ages, b.loc[ages, "lat_p50"], color=C_NEW, lw=1.5, ls=(0, (5, 2)))
    ax.set_title(title, fontsize=10.5, loc="left")
    ax.grid(True, color="#dddddd", lw=0.5)
    ax.set_axisbelow(True)
    if j == 0:
        ax.set_ylabel("Collision latitude (°)")

    # ---- bottom: difference vs the reference uncertainty ----
    ax = axes[1, j]
    hw68 = (a.loc[ages, "lat_p84"] - a.loc[ages, "lat_p16"]).values / 2.0
    hw95 = (a.loc[ages, "lat_p97_5"] - a.loc[ages, "lat_p2_5"]).values / 2.0
    diff = (b.loc[ages, "lat_p50"] - a.loc[ages, "lat_p50"]).values
    ax.fill_between(ages, -hw68, hw68, color=C_BAND, alpha=0.30, lw=0)
    ax.plot(ages, hw95, color=C_BAND, lw=0.9, ls=(0, (4, 2)))
    ax.plot(ages, -hw95, color=C_BAND, lw=0.9, ls=(0, (4, 2)))
    ax.axhline(0, color="#555555", lw=0.7)
    ax.plot(ages, diff, color="#222222", lw=1.5)
    ax.grid(True, color="#dddddd", lw=0.5)
    ax.set_axisbelow(True)
    ax.set_xlabel("Age (Ma)")
    ax.set_xlim(AGE_MAX, 0)
    ax.set_xticks([60, 40, 20, 0])
    if j == 0:
        ax.set_ylabel("Median difference (°)")

    # fraction of ages inside the 68% band of the reference (printed to the console,
    # it is the number quoted in the text)
    k = np.isfinite(diff) & np.isfinite(hw68)
    print(f"  {key:7s}: |diff| <= half-width 68% in {100*np.mean(np.abs(diff[k]) <= hw68[k]):5.1f}% "
          f"of the ages | <= 95% in {100*np.mean(np.abs(diff[k]) <= hw95[k]):5.1f}% "
          f"| max |diff| {np.nanmax(np.abs(diff)):.2f}°")

ymin = min(ax.get_ylim()[0] for ax in axes[1, :])
ymax = max(ax.get_ylim()[1] for ax in axes[1, :])
for ax in axes[1, :]:
    ax.set_ylim(ymin, ymax)
for ax in axes[1, 1:]:
    ax.tick_params(labelleft=False)

handles_top = [
    Patch(facecolor=C_REF, alpha=0.26, label="68% band, reference circuit"),
    Patch(facecolor=C_NEW, alpha=0.26, label="68% band, revised circuit"),
    Line2D([0], [0], color=C_REF, lw=1.5, label="median, reference"),
    Line2D([0], [0], color=C_NEW, lw=1.5, ls=(0, (5, 2)), label="median, revised"),
]
axes[0, 0].legend(handles=handles_top, loc="lower left", frameon=False,
                  fontsize=7.4, ncol=1, handlelength=1.8, labelspacing=0.28,
                  borderaxespad=0.5)

handles_bot = [
    Line2D([0], [0], color="#222222", lw=1.5, label="revised − reference"),
    Patch(facecolor=C_BAND, alpha=0.30, label="±68% half-width"),
    Line2D([0], [0], color=C_BAND, lw=0.9, ls=(0, (4, 2)), label="±95% half-width"),
]
axes[1, 0].legend(handles=handles_bot, loc="lower left", frameon=False,
                  fontsize=6.6, ncol=1, handlelength=1.6, labelspacing=0.25)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG_OUT, f"F5_model_comparison.{ext}"), dpi=300,
                bbox_inches="tight")
print("written: F5_model_comparison.pdf/.png")
