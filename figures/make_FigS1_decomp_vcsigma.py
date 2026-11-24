# -*- coding: utf-8 -*-
"""
Figure S1 (Geology paper, Supplemental) -- two panels:
  (A) Source decomposition of the collision-latitude uncertainty: 95% band width
      contributed by each error source sampled alone (relative plate circuit,
      WK08-A hotspot frame, seamount geochronology; reference model f_PAC=0.44),
      for the three tracks vs age. Colour = source, line style = track.
  (B) 1 sigma of the orthogonal Farallon/Nazca-South America convergence rate as a
      function of latitude and age (hybrid: 0-28 Ma analytic from Quiero et al.
      2022; >28 Ma Monte Carlo of the circuit, Rowan & Rowley 2014 branch only =
      lower bound). Collision latitudes and restored deposits overlaid.

New for this paper; combines the per-source widths of the local
tables/mc_collision_bands_{circuit,wk08,geochron}_present.csv (step 04
RE-RESTORED to the present-day margin by tables/make_present_bands.py) with the
vc-sigma field of make_F5_vc_sigma_field.py. The white collision-latitude lines
of panel B are mapped to the present-day margin with tables/margin_map.py.
Output goes to figures/.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm"})

BASE = os.path.dirname(os.path.abspath(__file__))
TABLES = os.path.join(os.path.dirname(BASE), "tables")
import sys
sys.path.insert(0, TABLES)
import margin_map as MM
# Pipeline with sigma_omega JFR/San Felix = 6.7 km/Myr (recalculated by L.E. Lara).
# PAPER_PIPE selects the kinematic version:
#   py make_FigS1_decomp_vcsigma.py                      -> FigS1_decomp_vcsigma.*
#   PAPER_PIPE=WD25_WK08-A py make_FigS1_...py         -> ..._wd25.*
# In the WD25 version panel B uses the Nazca-South America stage covariances
# published directly by Wilson and DeMets (2025) (step 05b), which reach
# 33.7 Ma instead of 28, so the crossover of the covariance basis shifts.
PIPE_NAME = os.environ.get("PAPER_PIPE", "Q22_WK08-A")
IS_WD25 = PIPE_NAME != "Q22_WK08-A"
OUT_SFX = "_wd25" if IS_WD25 else ""
BANDS_SFX = "_" + PIPE_NAME
PIPE = os.path.join(__import__("os").environ["PAPER_ROOT"], "models", PIPE_NAME)
VC_ANALYTIC = os.path.join(PIPE, "vc_sigma_wd25.csv" if IS_WD25
                           else "vc_sigma_quiero.csv")
DEP = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "deposits")  # deposit database (see rotation_models/README.md)
GEOCHEM = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "geochemistry",
          "Full_Data_v11_tt_with_Clusters_v12_geodyn2.xlsx")
AGE_MAX = 66.0
PF = "0.44"
# Panel B shows sigma_ort, the uncertainty of the orthogonal component (V_perp)
# discussed in the text (eq. 5) and listed in Table S3. sigma_mag (uncertainty of
# the velocity magnitude) differs little in the median (ratio mag/ort 0.94) but
# much at the extreme (max 6.60 vs 11.60 cm/yr in the base model). To plot the
# magnitude instead: SIGMA_COL = "sigma_mag".
SIGMA_COL = "sigma_ort"
# T_SPLIT (age reached by the published stage covariances) is derived from the
# data: 27.9 Ma with Quiero et al. (2022) and 33.7 Ma with Wilson and DeMets (2025).
_halo = [pe.withStroke(linewidth=2.4, foreground="black")]

# ----------------------------------------------------------------- Panel A data
# PAPER_FIG_TAG = tag of the bands (e.g. _hsb_cc_g13), PAPER_TABLES_DIR = folder of the bands,
# PAPER_FIG_OUT = output folder. Without them, the single-model two-panel figure below is drawn.
# Foundation (fourth chain) is included only if it is present in the bands.
TAG = os.environ.get("PAPER_FIG_TAG", "")
TABLES = os.environ.get("PAPER_TABLES_DIR", TABLES)
FIG_OUT = os.environ.get("PAPER_FIG_OUT", BASE)
SRC = [("circuit", "#1b9e77", "Relative plate circuit"),
       ("wk08", "#7570b3", "WK08-A hotspot frame"),
       # with the step 03b fit (option b) the geochronology enters as the "hotspot" variant
       ("hotspot" if TAG else "geochron", "#d95f02", "Hotspot fit (seamount geochronology)")]
TRACKS = [("Easter", "-", "Easter"),
          ("TT", (0, (5, 2)), "San Félix"),
          ("JFR", (0, (1, 1.2)), "Juan Fernández"),
          ("Foundation", (0, (6, 1.5, 1, 1.5)), "Foundation")]
# deposits and drawing helpers (also used by vcsigma_circuito)
por = pd.read_csv(os.path.join(DEP, "Porphyry_Plate_ID_geodyn_v7.csv"))
epi = pd.read_csv(os.path.join(DEP, "Epithermal_Plate_Id_geodyn_v6_thesis.csv"))
for df, ton in ((por, "Tonnage Mt"), (epi, "Total reso")):
    df["Age"] = pd.to_numeric(df["Age MY"], errors="coerce")
    df["Lat"] = pd.to_numeric(df["Latitude_ROT"], errors="coerce")
    df["Ton"] = pd.to_numeric(df[ton], errors="coerce")
por = por.dropna(subset=["Age", "Lat"]); por = por[por.Age <= AGE_MAX]
epi = epi.dropna(subset=["Age", "Lat"]); epi = epi[epi.Age <= AGE_MAX]

def por_area(T):
    T = np.asarray(T, float)
    return 16 + 130 * np.clip(np.log10(np.clip(T, 1, None)) - 2.2, 0, None) ** 1.35

def label_on_line(ax, a, la, lab, anchor):
    j = int(np.argmin(np.abs(a - anchor)))
    j0 = max(0, j - 2); j1 = min(len(a) - 1, j + 2)
    p0 = ax.transData.transform((a[j0], la[j0]))
    p1 = ax.transData.transform((a[j1], la[j1]))
    ang = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
    if ang > 90: ang -= 180
    if ang < -90: ang += 180
    ax.text(a[j], la[j], lab, color="white", fontsize=6.8, fontweight="bold",
            va="bottom", ha="center", rotation=ang, rotation_mode="anchor",
            zorder=6, path_effects=_halo)


def decomp_modelos():
    """Panel A with all models of the set (paper_models.py) weighted equally. One panel per
    chain; for each source, range between models of the 95% width (fill) and its median (line).
    The "frame" source is the frame covariance of each model (WK08-A, GMHRF or GWG20; wk08
    variant of step 04). Panel B (velocity) depends only on the circuit and is drawn
    separately."""
    import paper_models as MP
    src = [("circuit", "#1b9e77", "Relative plate circuit"),
           ("wk08", "#7570b3", "Hotspot frame (WK08-A, GMHRF or GWG20)"),
           ("hotspot", "#d95f02", "Hotspot fit (seamount geochronology)")]
    dd = []
    for c in MP.CARPETAS:
        for key, _, _ in src:
            b = pd.read_csv(os.path.join(TABLES, f"mc_collision_bands_{key}_present"
                                                 f"{MP.sufijo_bandas(c)}{TAG}.csv"),
                            dtype={"pf": str})
            b = MP.recortar(b[(b.pf == PF) & (b.age <= AGE_MAX)], c)   # frame_coverage
            dd.append(b.assign(w95=b.lat_p97_5 - b.lat_p2_5, fuente=key, modelo=c))
    d = pd.concat(dd)
    tr = [t for t in TRACKS if t[0] in set(d.ridge)]
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.6), sharex=True, sharey=True)
    for ax, (rk, _, lab) in zip(axs.flat, tr):
        for key, col, _ in src:
            g = d[(d.ridge == rk) & (d.fuente == key)].groupby("age").w95
            q = pd.DataFrame({"lo": g.min(), "hi": g.max(), "med": g.median()})
            ax.fill_between(q.index, q.lo, q.hi, color=col, alpha=0.22, lw=0)
            ax.plot(q.index, q.med, color=col, lw=1.4)
        ax.set_title(f"{lab} track", fontsize=10, loc="left")
        ax.set_xlim(AGE_MAX, 0); ax.set_ylim(0, None)
        ax.grid(True, color="#dddddd", lw=0.5); ax.set_axisbelow(True)
    for ax in axs[:, 0]:
        ax.set_ylabel("95% band width (°)", fontsize=10)
    for ax in axs[-1, :]:
        ax.set_xlabel("Age (Ma)", fontsize=10)
    from matplotlib.patches import Patch
    h = ([Line2D([0], [0], color=col, lw=1.8, label=lab) for _, col, lab in src]
         + [Line2D([0], [0], color="0.35", lw=1.4, label=f"median of the {len(MP.CARPETAS)} models"),
            Patch(facecolor="0.35", alpha=0.22, label="range between models")])
    fig.suptitle("Source decomposition of the collision-latitude uncertainty "
                 "($f_\\mathrm{PAC}=0.44$)", fontsize=10.5, x=0.02, ha="left")
    fig.tight_layout()
    fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, 0.0), frameon=False,
               fontsize=7.8, ncol=3)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_OUT, f"FigS1_decomp{MP.SUFIJO_SALIDA}.{ext}"), dpi=350,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"FigS1_decomp{MP.SUFIJO_SALIDA}: OK")


def vcsigma_circuito():
    """Panel B (1 sigma of V_perp) for one circuit, PAPER_PIPE = Q22_WK08-A, WD25_WK08-A or
    SG12_WK08-A. Analytic part with the published stage covariances (Quiero et al. 2022,
    step 05; Wilson and DeMets 2025, step 05b); Somoza and Ghidella (2012) do not publish
    stage covariances, so for SG12 the whole range comes from the Monte Carlo of the circuit
    (step 09). Collision lines: the three frames of the main text with that circuit, clipped
    by frame_coverage."""
    import paper_models as MP
    circ = {"Q22_WK08-A": ("Quiero", "", "vc_sigma_quiero.csv"),
            "WD25_WK08-A": ("WD25", "_wd25", "vc_sigma_wd25.csv"),
            "SG12_WK08-A": ("SG12", "_sg12", None)}
    ci, sfx, fa = circ[PIPE_NAME]
    lat_axis = np.arange(-35.75, -5.0, 0.5)
    y_edges = np.concatenate([[lat_axis[0] - 0.25], lat_axis + 0.25])
    _base = plt.cm.jet(np.linspace(0, 1, 256))
    _pastel = 0.55 * _base + 0.45 * np.ones_like(_base)
    _pastel[:, 3] = 1.0
    cmap = LinearSegmentedColormap.from_list("jet_pastel", _pastel)
    fig, ax = plt.subplots(figsize=(7.4, 4.9))
    t_split = 0.0
    if fa:
        q = pd.read_csv(os.path.join(PIPE, fa))
        q["t_young"] = q[["t_ini", "t_fin"]].min(axis=1)
        q["t_old"] = q[["t_ini", "t_fin"]].max(axis=1)
        st = q[["t_young", "t_old"]].drop_duplicates().sort_values("t_young").values
        Q = np.full((len(lat_axis), len(st)), np.nan)
        for j, (t0, t1) in enumerate(st):
            s = q[(q.t_young == t0) & (q.t_old == t1)].sort_values("lat")
            Q[:, j] = np.interp(lat_axis, s.lat, s[SIGMA_COL], left=np.nan, right=np.nan)
        t_split = float(st[-1, 1])
        pcm = ax.pcolormesh(np.concatenate([st[:, 0], [st[-1, 1]]]), y_edges,
                            np.ma.masked_invalid(Q), cmap=cmap, vmin=0, vmax=3.0,
                            shading="flat", zorder=0, rasterized=True)
    dmc = pd.read_csv(os.path.join(PIPE, "vc_sigma_circuit_mc.csv"))
    dmc = dmc[dmc.age > t_split]
    am = np.sort(dmc.age.unique())
    Mz = np.full((len(lat_axis), len(am)), np.nan)
    for j, t in enumerate(am):
        s = dmc[dmc.age == t].sort_values("lat")
        Mz[:, j] = np.interp(lat_axis, s.lat, s[SIGMA_COL], left=np.nan, right=np.nan)
    xm = np.concatenate([[max(t_split, am[0] - 0.5) if fa else am[0] - 1.0], am + 0.5])
    pcm2 = ax.pcolormesh(xm, y_edges, np.ma.masked_invalid(Mz), cmap=cmap, vmin=0, vmax=3.0,
                         shading="flat", zorder=0, rasterized=True)
    if fa:
        ax.axvline(t_split, color="black", lw=1.0, ls=(0, (4, 3)), zorder=4)
    labs = [("Easter", "Easter (Nazca R.)", 13), ("JFR", "Juan Fernández R.", 21),
            ("TT", "San Félix (Taltal R.)", 16), ("Foundation", "Foundation (Iquique R.)", 40)]
    lines = {}
    for c in [c for c, cc, _ in MP.MODELOS if cc == ci]:
        b = MP.recortar(pd.read_csv(os.path.join(
            TABLES, f"mc_collision_bands_total_present{MP.sufijo_bandas(c)}{TAG}.csv"),
            dtype={"pf": str}), c)
        b = b[(b.pf == PF) & (b.age <= AGE_MAX)]
        for k, lab, anc in labs:
            t = b[b.ridge == k].dropna(subset=["lat_nominal"]).sort_values("age")
            if len(t):
                lines.setdefault(lab, []).append((t.age.values, t.lat_nominal.values, anc))
                ax.plot(t.age, t.lat_nominal, "-", color="white", lw=1.0, zorder=3,
                        path_effects=_halo)
    ax.set_xlim(AGE_MAX, 0); ax.set_ylim(-36, -5)
    for lab, ll in lines.items():
        a, la, anc = min(ll, key=lambda x: np.min(np.abs(x[0] - x[2])))
        label_on_line(ax, a, la, lab, anc)
    ax.scatter(por["Age"], por["Lat"], s=por_area(por["Ton"]), marker="*", c="gold",
               edgecolor="black", linewidth=0.5, zorder=8)
    ax.scatter(epi["Age"], epi["Lat"], s=34, marker="P", c="red", edgecolor="black",
               linewidth=0.4, zorder=7)
    ax.set_xlabel("Age (Ma)", fontsize=11)
    ax.set_ylabel("Latitude (°)", fontsize=11)
    nombre = {"Quiero": "Quiero et al. (2022)", "WD25": "Wilson and DeMets (2025)",
              "SG12": "Somoza and Ghidella (2012)"}[ci]
    ax.set_title(f"$1\\sigma$ of orthogonal convergence rate (cm/yr), circuit of {nombre}",
                 fontsize=10.5, fontweight="bold", loc="left")
    cb = fig.colorbar(pcm2, ax=ax, pad=0.015, extend="max", fraction=0.05)
    cb.set_label("$1\\sigma\\;V_\\perp$ (cm/yr)", fontsize=10)
    h = [Line2D([0], [0], color="0.35", lw=2.0, label="Collision latitude, each frame"),
         Line2D([0], [0], marker="*", color="w", markerfacecolor="gold",
                markeredgecolor="k", markersize=12, label="Porphyry Cu"),
         Line2D([0], [0], marker="P", color="w", markerfacecolor="red",
                markeredgecolor="k", markersize=10, label="Epithermal Au–Ag")]
    if fa:
        h.insert(1, Line2D([0], [0], color="black", lw=1.0, ls=(0, (4, 3)),
                           label=f"Covariance basis change ({t_split:.1f} Ma)"))
    ax.legend(handles=h, loc="lower left", fontsize=7.5, framealpha=0.9, edgecolor="0.7",
              ncol=2)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_OUT, f"FigS1b_vcsigma{sfx}.{ext}"), dpi=350)
    plt.close(fig)
    print(f"FigS1b_vcsigma{sfx}: OK  (analytic 0-{t_split:.1f} Ma; MC up to {am.max():.0f} Ma)")


if TAG:
    if os.environ.get("PAPER_S1B", "") == "1":
        vcsigma_circuito()
    else:
        decomp_modelos()
    raise SystemExit

bands = {}
for key, _, _ in SRC:
    b = pd.read_csv(os.path.join(TABLES,
                                 f"mc_collision_bands_{key}_present{BANDS_SFX}{TAG}.csv"),
                    dtype={"pf": str})
    bands[key] = b[b.pf == PF].copy()
    bands[key]["w95"] = bands[key]["lat_p97_5"] - bands[key]["lat_p2_5"]

# ----------------------------------------------------------------- Panel B data
lat_axis = np.arange(-35.75, -5.0, 0.5)
y_edges = np.concatenate([[lat_axis[0] - 0.25], lat_axis + 0.25])
_base = plt.cm.jet(np.linspace(0, 1, 256))
_pastel = 0.55 * _base + 0.45 * np.ones_like(_base)
_pastel[:, 3] = 1.0
CMAP = LinearSegmentedColormap.from_list("jet_pastel", _pastel)

q = pd.read_csv(VC_ANALYTIC)
# The two sources use opposite t_ini/t_fin conventions (Quiero: t_ini is the
# young end; WD25: the old end), so the stages are normalised to young-old
# before building the pcolormesh edges.
q["t_young"] = q[["t_ini", "t_fin"]].min(axis=1)
q["t_old"] = q[["t_ini", "t_fin"]].max(axis=1)
stages = (q[["t_young", "t_old"]].drop_duplicates()
          .sort_values("t_young").values)
xq_edges = np.concatenate([stages[:, 0], [stages[-1, 1]]])
T_SPLIT = float(stages[-1, 1])
Q_abs = np.full((len(lat_axis), len(stages)), np.nan)
for j, (t0, t1) in enumerate(stages):
    s = q[(q.t_young == t0) & (q.t_old == t1)].sort_values("lat")
    Q_abs[:, j] = np.interp(lat_axis, s.lat, s[SIGMA_COL], left=np.nan, right=np.nan)
print(f"analytic coverage: 0-{T_SPLIT:.2f} Ma ({len(stages)} stages)")

dmc = pd.read_csv(os.path.join(PIPE, "vc_sigma_circuit_mc.csv"))
dmc = dmc[dmc.age > T_SPLIT]
ages_mc = np.sort(dmc.age.unique())
M_abs = np.full((len(lat_axis), len(ages_mc)), np.nan)
for j, t in enumerate(ages_mc):
    s = dmc[dmc.age == t].sort_values("lat")
    M_abs[:, j] = np.interp(lat_axis, s.lat, s[SIGMA_COL], left=np.nan, right=np.nan)
xm_edges = np.concatenate([[max(T_SPLIT, ages_mc[0] - 0.5)], ages_mc + 0.5])

_LABS = [("Easter", "Easter (Nazca R.)", 13),
         ("JFR", "Juan Fernández R.", 21),
         ("TT", "San Félix (Taltal R.)", 16),
         ("Foundation", "Foundation (Iquique R.)", 40)]
# The white collision lines are read from the bands of the pipeline itself
# (total variant, f_PAC=0.44, already mapped to the present-day margin), so
# both versions of the figure come from the same source and are comparable.
_tot = pd.read_csv(os.path.join(TABLES,
                                f"mc_collision_bands_total_present{BANDS_SFX}{TAG}.csv"),
                   dtype={"pf": str})
_tot = _tot[_tot.pf == PF]
coll = {}
for t, lab, anc in _LABS:
    c = _tot[(_tot.ridge == t) & (_tot.age <= AGE_MAX)].sort_values("age")
    c = c.dropna(subset=["lat_nominal"])
    if len(c):          # some band tables do not include Foundation
        coll[lab] = (c["age"].values, c["lat_nominal"].values, anc)
TRACKS = [t for t in TRACKS if t[0] in set(_tot.ridge)]

# ----------------------------------------------------------------- figure
fig, (axA, axB) = plt.subplots(2, 1, figsize=(7.4, 8.6),
                               gridspec_kw={"height_ratios": [1.0, 1.25]})
fig.subplots_adjust(left=0.10, right=0.90, top=0.965, bottom=0.09, hspace=0.28)

# Panel A
for key, col, _ in SRC:
    b = bands[key]
    for rk, ls, _ in TRACKS:
        t = b[b.ridge == rk].sort_values("age")
        axA.plot(t.age, t.w95, color=col, ls=ls, lw=1.4)
axA.set_xlim(AGE_MAX, 0)
axA.set_ylim(0, None)
axA.set_xlabel("Age (Ma)", fontsize=11)
axA.set_ylabel("95% band width (°)", fontsize=11)
axA.set_title("A)  Source decomposition of the collision-latitude uncertainty "
              "($f_\\mathrm{PAC}=0.44$)", fontsize=11.5, fontweight="bold", loc="left")
axA.grid(True, color="#dddddd", lw=0.5); axA.set_axisbelow(True)
legA = ([Line2D([0], [0], color=col, lw=2.4, label=lab) for _, col, lab in SRC]
        + [Line2D([0], [0], color="0.35", lw=1.4, ls=ls, label=lab)
           for _, ls, lab in TRACKS])
axA.legend(handles=legA, loc="upper left", frameon=False, fontsize=7.8, ncol=2)

# Panel B
pcm = axB.pcolormesh(xq_edges, y_edges, np.ma.masked_invalid(Q_abs), cmap=CMAP,
                     vmin=0, vmax=3.0, shading="flat", zorder=0, rasterized=True)
axB.pcolormesh(xm_edges, y_edges, np.ma.masked_invalid(M_abs), cmap=CMAP,
               vmin=0, vmax=3.0, shading="flat", zorder=0, rasterized=True)
axB.axvline(T_SPLIT, color="black", lw=1.0, ls=(0, (4, 3)), zorder=4)
for lab, (a, la, anc) in coll.items():
    axB.plot(a, la, "-", color="white", lw=1.6, zorder=3, path_effects=_halo)
    label_on_line(axB, a, la, lab, anc)
axB.scatter(por["Age"], por["Lat"], s=por_area(por["Ton"]), marker="*",
            c="gold", edgecolor="black", linewidth=0.5, zorder=8)
axB.scatter(epi["Age"], epi["Lat"], s=34, marker="P",
            c="red", edgecolor="black", linewidth=0.4, zorder=7)
axB.set_xlim(AGE_MAX, 0); axB.set_ylim(-36, -5)
axB.set_xlabel("Age (Ma)", fontsize=11)
axB.set_ylabel("Latitude (°)", fontsize=11)
axB.set_title("B)  $1\\sigma$ of orthogonal convergence rate (cm/yr)",
              fontsize=11.5, fontweight="bold", loc="left")
cb = fig.colorbar(pcm, ax=axB, pad=0.015, extend="max", fraction=0.05)
cb.set_label("$1\\sigma\\;V_\\perp$ (cm/yr)", fontsize=10)
h_type = [
    Line2D([0], [0], color="0.35", lw=2.0, label="Collision latitude"),
    Line2D([0], [0], color="black", lw=1.0, ls=(0, (4, 3)),
           label=f"Covariance basis change ({T_SPLIT:.1f} Ma)"),
    Line2D([0], [0], marker="*", color="w", markerfacecolor="gold",
           markeredgecolor="k", markersize=12, label="Porphyry Cu"),
    Line2D([0], [0], marker="P", color="w", markerfacecolor="red",
           markeredgecolor="k", markersize=10, label="Epithermal Au–Ag"),
]
axB.legend(handles=h_type, loc="lower left", fontsize=7.5, framealpha=0.9,
           edgecolor="0.7", ncol=2)

for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG_OUT, f"FigS1_decomp_vcsigma{OUT_SFX}.{ext}"), dpi=350)
plt.close(fig)
print(f"FigS1_decomp_vcsigma{OUT_SFX}: OK")
