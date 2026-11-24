# -*- coding: utf-8 -*-
"""
Figure 2 -- Monte Carlo error bands of the ridge-trench collision latitude vs
age, for the three spreading-asymmetry models + 95% envelope.

Reads the pipeline output directly. Output goes to this folder (figures/).

Without PAPER_FIG_TAG, TWO versions are written, one per plate-circuit version; a
single panel with 3 f_PAC models x 2 circuits (6 bands + 2 envelopes per track)
becomes saturated and unreadable.
  F2_collision_bands.pdf       reference circuit
  F2_collision_bands_wd25.pdf  revised circuit (Wilson and DeMets 2025)

REAL DATA: tables/mc_collision_bands_total_present{,_wd25}.csv (pipeline step 04,
N=1000 per f_PAC model, joint sampling of circuit + WK08-A + geochron,
RE-RESTORED to the present-day margin by tables/make_present_bands.py).
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
AGE_MAX = 45 if os.environ.get("PAPER_MODELOS") == "supl" else 66   # GMHRF-A' resolves up to ~45 Ma
MODELS = [("0.5", "#c51b8a", "$f_\\mathrm{PAC}=0.50$"),
          ("0.44", "#2b5f9e", "$f_\\mathrm{PAC}=0.44$"),
          ("0.36", "#31a354", "$f_\\mathrm{PAC}=0.36$")]
ENV_COLOR = "#9a9a9a"
PANELS = [("Easter", "Easter track (Nazca Ridge / Easter Seamount Chain)"),
          ("TT", "San Félix track (Taltal Ridge)"),
          ("JFR", "Juan Fernández track (Juan Fernández Ridge)"),
          ("Foundation", "Foundation track (Iquique Ridge)")]
# PAPER_FIG_TAG = tag of the band files (e.g. _hsb_cc_g13), PAPER_TABLES_DIR = folder of
# the bands, PAPER_FIG_OUT = output folder. Foundation (fourth chain) is plotted only if it
# is present in the bands.
TAG = os.environ.get("PAPER_FIG_TAG", "")
TAB = os.environ.get("PAPER_TABLES_DIR", TAB)
FIG_OUT = os.environ.get("PAPER_FIG_OUT", BASE)

# (input CSV suffix, output file suffix, subtitle)
VERSIONS = [("_Q22_WK08-A", "", "Reference plate circuit"),
            ("_WD25_WK08-A", "_wd25",
             "Revised Nazca–Pacific circuit (Wilson and DeMets 2025)")]


def make(csv_sfx, out_sfx, subtitle):
    src = os.path.join(TAB, f"mc_collision_bands_total_present{csv_sfx}{TAG}.csv")
    if not os.path.exists(src):
        print(f"{src} not found, skipped")
        return
    d = pd.read_csv(src, dtype={"pf": str})
    panels = [p for p in PANELS if p[0] in set(d.ridge)]

    fig, axes = plt.subplots(len(panels), 1, figsize=(6.8, 8.8 / 3 * len(panels)),
                             sharex=True)
    for ax, (key, title) in zip(axes, panels):
        sub = {pf: d[(d.ridge == key) & (d.pf == pf)].sort_values("age")
               for pf, _, _ in MODELS}
        ages = sorted(set.intersection(*[set(s.age) for s in sub.values()]))
        lo = np.min([sub[pf].set_index("age").loc[ages, "lat_p2_5"].values
                     for pf, _, _ in MODELS], axis=0)
        hi = np.max([sub[pf].set_index("age").loc[ages, "lat_p97_5"].values
                     for pf, _, _ in MODELS], axis=0)
        ax.fill_between(ages, lo, hi, color=ENV_COLOR, alpha=0.25, lw=0)
        for pf, col, lab in MODELS:
            t = sub[pf]
            ax.fill_between(t.age, t.lat_p16, t.lat_p84, color=col, alpha=0.28,
                            lw=0)
            ax.plot(t.age, t.lat_p50, color=col, lw=1.5,
                    ls="-" if pf != "0.36" else (0, (5, 2)))
        ax.set_title(title, fontsize=10.5, loc="left")
        ax.set_ylabel("Collision latitude (°)")
        ax.grid(True, color="#dddddd", lw=0.5)
        ax.set_axisbelow(True)
        ax.set_xlim(AGE_MAX, 0)

    handles = ([Patch(facecolor=ENV_COLOR, alpha=0.25,
                      label="95% envelope (all $f_\\mathrm{PAC}$ models)")]
               + [Patch(facecolor=col, alpha=0.28, label=f"68% band, {lab}")
                  for pf, col, lab in MODELS]
               + [Line2D([0], [0], color=col, lw=1.5,
                         ls="-" if pf != "0.36" else (0, (5, 2)),
                         label=f"MC median, {lab}")
                  for pf, col, lab in MODELS])
    axes[0].legend(handles=handles, loc="lower left", frameon=False,
                   fontsize=7.6, ncol=2, columnspacing=1.2)
    axes[-1].set_xlabel("Age (Ma)")
    fig.suptitle(subtitle, fontsize=10.5, y=0.995)

    fig.tight_layout(rect=(0, 0, 1, 0.985))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_OUT, f"F2_collision_bands{out_sfx}.{ext}"),
                    dpi=300)
    plt.close(fig)
    print(f"written: F2_collision_bands{out_sfx}.pdf/.png  [{subtitle}]")


def make_modelos():
    """With PAPER_FIG_TAG the figure shows all models of the set (paper_models.py) equally:
    grey = 95% envelope of all models and f_PAC; colour = union of the models' 68% bands for
    each f_PAC; thin lines = median of each model."""
    import paper_models as MP
    dd = []
    for c in MP.CARPETAS:
        src = os.path.join(TAB, f"mc_collision_bands_total_present{MP.sufijo_bandas(c)}{TAG}.csv")
        # each model is limited to the frame_coverage of its frame covariance
        dd.append(MP.recortar(pd.read_csv(src, dtype={"pf": str}), c).assign(modelo=c))
    d = pd.concat(dd)
    panels = [p for p in PANELS if p[0] in set(d.ridge)]
    fig, axes = plt.subplots(len(panels), 1, figsize=(6.8, 8.8 / 3 * len(panels)),
                             sharex=True)
    for ax, (key, title) in zip(axes, panels):
        r = d[(d.ridge == key) & (d.age <= AGE_MAX)]
        env = r.groupby("age").agg(lo=("lat_p2_5", "min"), hi=("lat_p97_5", "max"))
        ax.fill_between(env.index, env.lo, env.hi, color=ENV_COLOR, alpha=0.25, lw=0)
        for pf, col, lab in MODELS:
            s = r[r.pf == pf]
            u = s.groupby("age").agg(lo=("lat_p16", "min"), hi=("lat_p84", "max"))
            ax.fill_between(u.index, u.lo, u.hi, color=col, alpha=0.22, lw=0)
            for _, t in s.groupby("modelo"):
                t = t.sort_values("age")
                ax.plot(t.age, t.lat_p50, color=col, lw=0.7, alpha=0.85,
                        ls="-" if pf != "0.36" else (0, (5, 2)))
        ax.set_title(title, fontsize=10.5, loc="left")
        ax.set_ylabel("Collision latitude (°)")
        ax.grid(True, color="#dddddd", lw=0.5)
        ax.set_axisbelow(True)
        ax.set_xlim(AGE_MAX, 0)
    n = len(MP.CARPETAS)
    handles = ([Patch(facecolor=ENV_COLOR, alpha=0.25,
                      label=f"95% envelope ({n} models, all $f_\\mathrm{{PAC}}$)")]
               + [Patch(facecolor=col, alpha=0.22, label=f"68% envelope, {lab}")
                  for pf, col, lab in MODELS]
               + [Line2D([0], [0], color=col, lw=0.9,
                         ls="-" if pf != "0.36" else (0, (5, 2)),
                         label=f"median of each model, {lab}")
                  for pf, col, lab in MODELS])
    axes[-1].set_xlabel("Age (Ma)")
    fig.tight_layout()
    # below the figure: inside a panel it covers the medians
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), frameon=False,
               fontsize=7.6, ncol=2, columnspacing=1.2)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_OUT, f"F2_collision_bands{MP.SUFIJO_SALIDA}.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"written: F2_collision_bands{MP.SUFIJO_SALIDA}.pdf/.png  [{n} models]")


if TAG:
    make_modelos()
else:
    for v in VERSIONS:
        make(*v)
