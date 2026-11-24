# -*- coding: utf-8 -*-
"""
Figure 1 (Geology paper) -- Location map + present-day configuration of the three
modelled hotspot tracks (f_PAC=0.44) with 68/95% confidence ellipses of the
present-day track position at 10-Myr track-age nodes (MC total, step 08).

Reads the model outputs and the input data under PAPER_ROOT. Output goes to this
folder (figures/).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

BASE = os.path.dirname(os.path.abspath(__file__))
# Model folder (sigma_omega JFR/San Felix = 6.7 km/Myr, recalculated by L.E. Lara).
# PAPER_PIPE / PAPER_PIPE_WD select the model folders (name or path) and PAPER_FIG_TAG the
# tag of the step 08 outputs (e.g. _hsb_cc_g13).
PIPE = os.path.join(__import__("os").environ["PAPER_ROOT"], "models", os.environ.get("PAPER_PIPE", "Q22_WK08-A"))
TAG = os.environ.get("PAPER_FIG_TAG", "")
CL_LINE = (__import__("os").environ["PAPER_ROOT"] + "/external/muller2019/Coastlines_USGS_MedRes_AREPS_polyline.shp")
MARGIN = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "margin", "Margin_SAM_0.1.csv")  # copy of Margin_SAM (see rotation_models/README.md)
DEP = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "deposits")  # copy of the deposit database (see rotation_models/README.md)
GEOCHEM = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "geochemistry",
          "Full_Data_v11_tt_with_Clusters_v12_geodyn2.xlsx")

KM = 111.19492664455873
EXT = (-112, -55, -42, 2)
RIDGE_STYLE = {"Easter": ("#d95f02", "Easter track"),
               "TT": ("#7570b3", "San Félix track"),
               "JFR": ("#1b9e77", "Juan Fernández track"),
               # fourth chain; palette validated with the four colours
               "Foundation": ("#8c510a", "Foundation track")}
_wh = [pe.withStroke(linewidth=2.0, foreground="white")]

# Without PAPER_FIG_TAG the figure overlays the TWO circuit versions. Colour identifies
# the track; the model is distinguished by symbology (reference: filled ellipse and
# solid track; revised: dashed outline only).
PIPE_WD = os.path.join(os.path.dirname(BASE), "models", os.environ.get("PAPER_PIPE_WD", "WD25_WK08-A"))

# With PAPER_FIG_TAG the figure shows all models equally: the nominal track of each one
# (thin line) and, at each 10-Myr node, the envelope (union) of their 95% ellipses.
# The full grid has twelve models (3 circuits x 4 frames).
# PAPER_F1_ROOT = folder that contains the model folders
# (default: the project root).
# The model set comes from paper_models.py (PAPER_MODELOS: texto = nine models with
# GMHRF-B, supl = GMHRF-A', todos = twelve)
import paper_models as MP
MODELOS = MP.CARPETAS
OCHO = bool(TAG)
if OCHO:
    ROOT = os.environ.get("PAPER_F1_ROOT", os.path.dirname(BASE))
    ell_m = {d: pd.read_csv(os.path.join(ROOT, "models", d, f"track_ellipses_0.44{TAG}.csv"))
             for d in MODELOS}
    nom_m = {d: pd.read_csv(os.path.join(ROOT, "models", d, f"track_nominal_0.44{TAG}.csv"))
             for d in MODELOS}
    nom = nom_m[MODELOS[0]]
else:
    ell = pd.read_csv(os.path.join(PIPE, f"track_ellipses_0.44{TAG}.csv"))
    nom = pd.read_csv(os.path.join(PIPE, f"track_nominal_0.44{TAG}.csv"))
    ell_wd = pd.read_csv(os.path.join(PIPE_WD, f"track_ellipses_0.44{TAG}.csv"))
    nom_wd = pd.read_csv(os.path.join(PIPE_WD, f"track_nominal_0.44{TAG}.csv"))
# keep only the ridges present in the outputs (Foundation may be absent)
RIDGE_STYLE = {k: v for k, v in RIDGE_STYLE.items() if k in set(nom.ridge)}

import shapefile
sf = shapefile.Reader(CL_LINE)
coast = []
for sh in sf.shapes():
    pts = np.array(sh.points)
    parts = list(sh.parts) + [len(pts)]
    for a, b in zip(parts[:-1], parts[1:]):
        seg = pts[a:b]
        if ((seg[:, 0] >= EXT[0] - 2) & (seg[:, 0] <= EXT[1] + 2)
                & (seg[:, 1] >= EXT[2] - 2) & (seg[:, 1] <= EXT[3] + 2)).any():
            coast.append(seg)

m = pd.read_csv(MARGIN)
row = m.iloc[(pd.to_numeric(m["age"], errors="coerce") - 0).abs().idxmin()]
nn = sum(1 for c in m.columns if c.startswith("Lon_"))
tr_lon = np.array([row.get(f"Lon_{k}") for k in range(1, nn + 1)], float)
tr_lat = np.array([row.get(f"Lat_{k}") for k in range(1, nn + 1)], float)

por = pd.read_csv(os.path.join(DEP, "Porphyry_Plate_ID_geodyn_v7.csv"))
epi = pd.read_csv(os.path.join(DEP, "Epithermal_Plate_Id_geodyn_v6_thesis.csv"))

EXCLUDED = {m2.strip().lower() for m2 in [
    "K/Ar plagioclase", "K/Ar amphibole", "inference from K/Ar plagioclase",
    "Inference", "inference from K/Ar whole rock", "K/Ar whole rock",
    "Inference from K/Ar (biotite plagioclase whole rock) data",
    "Inference from K/Ar (biotite amphibole) data",
    "Inference from K/Ar biotite data", "K/Ar biotite", "K/Ar feldespar",
    "K/Ar glass", "K/Ar", "N.D.", "Inference from K/Ar",
    "The age of sample S-16 is assigned", "?",
    "The age of sample S-59 is assigned",
    "The age of sample S-H22 (176m) is assigned",
    "The age of sample S-T2 is assigned",
    "The age of sample S-50 is assigned",
    "The age of sample S-31 is assigned",
    "K-Ar feldspar", "K-Ar", "inference with dated units", "K-Ar biotite",
    "K-Ar whole-rock", "K-Ar amphibole", "K-Ar hornblende", "K-Ar plagioclase",
    "radiometric", "K-Ar whole rock",
    "inference from dated samples of same unit with K-Ar"]}
smp = pd.read_excel(GEOCHEM)
_meth = smp["Method"].astype(str).str.strip().str.lower()
_age = pd.to_numeric(smp["Age (Ma)"], errors="coerce")
_bad = _meth.isin(EXCLUDED) | smp["Method"].isna() | (_meth == "") | (_meth == "nan")
smp = smp[(~_bad) | (_age < 1)]

def ellipse_xy(mlo, mla, a_km, b_km, az_deg, n=90):
    th = np.linspace(0, 2 * np.pi, n)
    az = np.radians(az_deg)
    u = np.array([np.sin(az), np.cos(az)])
    v = np.array([-np.cos(az), np.sin(az)])
    E = np.outer(np.cos(th) * a_km, u) + np.outer(np.sin(th) * b_km, v)
    lon = mlo + E[:, 0] / (KM * np.cos(np.radians(mla)))
    lat = mla + E[:, 1] / KM
    return lon, lat

fig, ax = plt.subplots(figsize=(8.4, 7.6))
# Background bathymetry: ETOPO1 Bedrock (Amante and Eakins 2009) at 2', clipped by
# extract_etopo1_F1.py. Muted ramp (greyish-blue ocean, grey continent) with hillshade,
# so that the coloured tracks stay dominant.
from matplotlib.colors import LightSource, LinearSegmentedColormap
_et = np.load(os.path.join(BASE, "data", "etopo1_F1.npz"))
_z = _et["z"].astype(float)
# ramp with its contrast in -5000 to -1500 m, where the aseismic ridges lie
_Z0, _Z1 = -6500.0, 5000.0
_pos = lambda zz: (zz - _Z0) / (_Z1 - _Z0)
_cm = LinearSegmentedColormap.from_list("bat", [
    (0.0, "#4f6b86"), (_pos(-4500), "#7d96ad"), (_pos(-3500), "#a9bccd"),
    (_pos(-2000), "#d3dee8"), (_pos(-1), "#e8eef3"),
    (_pos(0), "#ebe9e4"), (1.0, "#9c9a95")])
_rgb = LightSource(azdeg=315, altdeg=45).shade(
    _z, cmap=_cm, vmin=_Z0, vmax=_Z1, blend_mode="soft", vert_exag=0.05,
    dx=1852 * 2, dy=1852 * 2)
ax.imshow(_rgb, extent=(_et["lon"][0], _et["lon"][-1], _et["lat"][-1], _et["lat"][0]),
          origin="upper", interpolation="bilinear", alpha=0.8, zorder=0)
for seg in coast:
    ax.plot(seg[:, 0], seg[:, 1], color="0.55", lw=0.6, zorder=2)
ax.plot(tr_lon, tr_lat, color="black", lw=1.4, zorder=3)

ax.scatter(smp["Longitude"], smp["Latitude"], s=4, c="0.6", alpha=0.5,
           linewidths=0, zorder=4, rasterized=True)
ax.scatter(por["Longitude"], por["Latitude"], s=42, marker="*", c="gold",
           edgecolor="black", linewidth=0.5, zorder=6)
ax.scatter(epi["Longitude"], epi["Latitude"], s=26, marker="P", c="red",
           edgecolor="black", linewidth=0.4, zorder=6)

def envolvente(polys, n=400):
    """Union of (lon, lat) polygons on a grid: (X, Y, 0/1 mask)."""
    from matplotlib.path import Path
    lo = np.concatenate([p[0] for p in polys]); la = np.concatenate([p[1] for p in polys])
    pad = 0.05
    X, Y = np.meshgrid(np.linspace(lo.min() - pad, lo.max() + pad, n),
                       np.linspace(la.min() - pad, la.max() + pad, n))
    pts = np.column_stack([X.ravel(), Y.ravel()])
    M = np.zeros(len(pts), bool)
    for p in polys:
        M |= Path(np.column_stack(p)).contains_points(pts)
    return X, Y, M.reshape(X.shape).astype(float)


if OCHO:
    for ridge, (col, lab) in RIDGE_STYLE.items():
        for d in MODELOS:
            t = nom_m[d]
            t = t[(t.ridge == ridge) & (t.tau <= 66)].sort_values("tau")
            ax.plot(t.lon, t.lat, color=col, lw=0.8, alpha=0.9, zorder=5)
        e8 = pd.concat([ell_m[d].assign(modelo=d) for d in MODELOS])
        e8 = e8[e8.ridge == ridge]
        for tau, g in e8.groupby("tau"):
            if len(g) < len(MODELOS):
                print(f"  {ridge} tau={tau}: ellipse in {len(g)} of {len(MODELOS)} models")
            polys = [ellipse_xy(r.lon, r.lat, r.a95_km, r.b95_km, r.az_major_deg)
                     for _, r in g.iterrows()]
            X, Y, M = envolvente(polys)
            if not (EXT[0] <= X.mean() <= EXT[1] and EXT[2] <= Y.mean() <= EXT[3]):
                print(f"  out of frame: {ridge} tau={int(tau)}")
                continue
            # 95% faint and 68% strong
            ax.contourf(X, Y, M, levels=[0.5, 1.5], colors=[col], alpha=0.10, zorder=4)
            ax.contour(X, Y, M, levels=[0.5], colors=[col], linewidths=0.5, alpha=0.6,
                       zorder=4.5)
            p68 = [ellipse_xy(r.lon, r.lat, r.a68_km, r.b68_km, r.az_major_deg)
                   for _, r in g.iterrows()]
            X6, Y6, M6 = envolvente(p68)
            ax.contourf(X6, Y6, M6, levels=[0.5, 1.5], colors=[col], alpha=0.38, zorder=4.6)
            ax.contour(X6, Y6, M6, levels=[0.5], colors=[col], linewidths=0.9, zorder=5)
            ax.plot(g.lon, g.lat, "o", color=col, ms=1.6, zorder=6)
            top_ = Y[M > 0.5].max()
            ax.text(X[M > 0.5].mean(), top_ + 0.15, f"{int(tau)}", color=col,
                    fontsize=7.5, fontweight="bold", ha="center", va="bottom",
                    zorder=7, path_effects=_wh, clip_on=True)

for ridge, (col, lab) in ({} if OCHO else RIDGE_STYLE).items():
    t = nom[(nom.ridge == ridge) & (nom.tau <= 66)].sort_values("tau")
    ax.plot(t.lon, t.lat, color=col, lw=1.8, zorder=5)
    e = ell[ell.ridge == ridge]
    for _, r in e.iterrows():
        if not (EXT[0] <= r.lon <= EXT[1] and EXT[2] <= r.lat <= EXT[3]):
            print(f"  out of frame: {ridge} tau={int(r.tau)} "
                  f"({r.lat:.1f}, {r.lon:.1f})")
            continue
        for a, b, alpha in ((r.a95_km, r.b95_km, 0.18),
                            (r.a68_km, r.b68_km, 0.32)):
            lo, la = ellipse_xy(r.lon, r.lat, a, b, r.az_major_deg)
            ax.fill(lo, la, color=col, alpha=alpha, lw=0, zorder=4)
        lo, la = ellipse_xy(r.lon, r.lat, r.a95_km, r.b95_km, r.az_major_deg)
        ax.plot(lo, la, color=col, lw=0.7, zorder=5)
        ax.plot(r.lon, r.lat, "o", color=col, ms=2.4, zorder=6)
        dy = 1.05 * r.a95_km / KM
        ax.text(r.lon, r.lat + dy + 0.25, f"{int(r.tau)}", color=col,
                fontsize=7.5, fontweight="bold", ha="center", va="bottom",
                zorder=7, path_effects=_wh, clip_on=True)

    # ---- revised version (Wilson and DeMets 2025): dashed outline only ----
    t = nom_wd[(nom_wd.ridge == ridge) & (nom_wd.tau <= 66)].sort_values("tau")
    ax.plot(t.lon, t.lat, color=col, lw=1.3, ls=(0, (4, 2)), zorder=5,
            alpha=0.95)
    for _, r in ell_wd[ell_wd.ridge == ridge].iterrows():
        if not (EXT[0] <= r.lon <= EXT[1] and EXT[2] <= r.lat <= EXT[3]):
            continue
        lo, la = ellipse_xy(r.lon, r.lat, r.a95_km, r.b95_km, r.az_major_deg)
        ax.plot(lo, la, color=col, lw=0.9, ls=(0, (3, 1.6)), zorder=6)
        ax.plot(r.lon, r.lat, "s", color="white", markeredgecolor=col,
                markeredgewidth=0.8, ms=3.2, zorder=6)

# with all models the Easter label goes below the young segment (at 44 Ma it falls on the
# Foundation envelope)
_easter = (3, (0, -26)) if OCHO else (44, (0, -14))
for ridge, anchor_tau in (("Easter", _easter[0]), ("TT", 34), ("JFR", 26), ("Foundation", 36)):
    if ridge not in RIDGE_STYLE:
        continue
    col, lab = RIDGE_STYLE[ridge]
    t = nom[(nom.ridge == ridge) & (nom.tau <= 66)].dropna(subset=["lat"]).sort_values("tau")
    j = int(np.argmin(np.abs(t.tau.values - anchor_tau)))
    ax.annotate(lab, (t.lon.values[j], t.lat.values[j]),
                textcoords="offset points",
                xytext=(-10, -4) if ridge == "Foundation" else
                (_easter[1] if ridge == "Easter" else (0, -14)), color=col,
                fontsize=9, fontweight="bold",
                # JFR and San Felix right-aligned: centred, they touch the eastern edge
                ha="right" if ridge in ("Foundation", "JFR", "TT") else "center",
                path_effects=_wh, zorder=7)

ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3])
ax.set_aspect(1.0 / np.cos(np.radians(20)))
ax.set_xlabel("Longitude (°)", fontsize=11)
ax.set_ylabel("Latitude (°)", fontsize=11)
ax.grid(True, color="#e3e3e3", lw=0.5)
ax.set_axisbelow(True)
ax.tick_params(labelsize=9)

handles = [
    Line2D([0], [0], color="black", lw=1.4, label="Peru–Chile trench"),
    Line2D([0], [0], color="0.55", lw=0.8, label="Coastline"),
    Line2D([0], [0], marker="*", color="w", markerfacecolor="gold",
           markeredgecolor="k", markersize=12, label="Porphyry Cu"),
    Line2D([0], [0], marker="P", color="w", markerfacecolor="red",
           markeredgecolor="k", markersize=9, label="Epithermal Au–Ag"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="0.6",
           markeredgecolor="none", markersize=5, label="Igneous sample"),
] + [Line2D([0], [0], color=col, lw=1.8, label=lab)
     for ridge, (col, lab) in RIDGE_STYLE.items()]
if OCHO:
    from matplotlib.patches import Patch
    handles += [
        Line2D([0], [0], color="0.25", lw=0.8, label=f"Nominal track, each of the {len(MODELOS)} models"),
        Patch(facecolor="0.25", alpha=0.38, edgecolor="0.25", lw=0.9,
              label=f"Envelope of the {len(MODELOS)} models' 68% ellipses"),
        Patch(facecolor="0.25", alpha=0.10, edgecolor="0.25", lw=0.5,
              label=f"Envelope of the {len(MODELOS)} models' 95% ellipses"),
    ]
    leg_title = "68 / 95% envelopes at 10-Myr track-age nodes"
else:
    handles += [
        Line2D([0], [0], color="0.25", lw=1.8,
               label="Reference circuit (filled ellipses)"),
        Line2D([0], [0], color="0.25", lw=1.3, ls=(0, (4, 2)),
               label="Revised circuit, W&D 2025 (outlines)"),
    ]
    leg_title = "68 / 95% ellipses at 10-Myr track-age nodes"
if OCHO:   # legend below the map: inside it covers Easter node 0 and Foundation
    fig.set_size_inches(8.4, 8.9)
    leg = ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.09),
                    ncol=3, fontsize=8, framealpha=0.9, edgecolor="0.7",
                    title=leg_title, title_fontsize=8)
else:
    leg = ax.legend(handles=handles, loc="lower left", fontsize=8,
                    framealpha=0.9, edgecolor="0.7",
                    title=leg_title, title_fontsize=8)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(os.environ.get("PAPER_FIG_OUT", BASE), f"F1_track_ellipses{MP.SUFIJO_SALIDA if OCHO else ''}.{ext}"),
                dpi=350, bbox_inches="tight" if OCHO else None)   # no blank space above the map
print("F1_track_ellipses: OK")
