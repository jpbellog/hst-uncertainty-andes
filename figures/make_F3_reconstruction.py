# -*- coding: utf-8 -*-
"""
Figure 3 -- Reconstruction WITH UNCERTAINTY: 8 panels by geological epoch
(tracks reconstructed at the mean age of the interval, trench, restored samples
and deposits), plus the 68/95% CONFIDENCE ELLIPSES of the RECONSTRUCTED position
of each track at the panel age, at the point-age nodes (tau) that are multiples
of 10 Ma and already existed at that age (+ the young end), from the Monte Carlo
of step 10 of the hst_uncertainty pipeline (total variant, f_PAC=0.44, N=400).

The tracks of the end-member models f_PAC=0.50 (dashed) and 0.36 (dash-dot) are
drawn semi-transparent, from step 11 (recon_model_tracks.csv). Every track
segment with tau > 83 Ma (in the three models) is drawn DOTTED AND FAINT: beyond
chron 34y (83 Ma) of Rowan & Rowley the rotation model has no covariances and
uncertainty ellipses CANNOT be computed.

Inputs: margin, igneous samples, deposits and HST_Margin tracks under
PAPER_ROOT/inputs, plus, in the model folder,
  recon_ellipses_0.44.csv
  recon_model_tracks.csv

Output: figures/F3_reconstruction{OUT_SFX}.pdf/.png
"""
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.patheffects as pe

BASE = os.path.dirname(os.path.abspath(__file__))
HST_DIR = (__import__("os").environ["PAPER_ROOT"] + "/inputs/validation/HST_Margin")
MARGIN = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "margin", "Margin_SAM_0.1.csv")  # copy of Margin_SAM (see rotation_models/README.md)
GEOCHEM = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "geochemistry",
          "Full_Data_v11_tt_with_Clusters_v12_geodyn2.xlsx")
DEP = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "deposits")  # copy of the deposit database (see rotation_models/README.md)
POR = os.path.join(DEP, "Porphyry_Plate_ID_geodyn_v7.csv")
EPI = os.path.join(DEP, "Epithermal_Plate_Id_geodyn_v6_thesis.csv")
# Reads the model folder (sigma_omega JFR/San Felix = 6.7 km/Myr, recalculated by
# L.E. Lara). The circuit version is selected with PAPER_PIPE, so that the figure
# is built in the same way for the two kinematic models:
#   py make_F3_reconstruction.py                    -> F3_reconstruction.pdf
#   PAPER_PIPE=WD25_WK08-A py make_F3_...py       -> F3_reconstruction_wd25.pdf
PIPE_NAME = os.environ.get("PAPER_PIPE", "Q22_WK08-A")
OUT_SFX = os.environ.get("PAPER_FIG_SFX", "" if PIPE_NAME == "Q22_WK08-A" else "_wd25")
PIPE = os.path.join(os.path.dirname(BASE), "models", PIPE_NAME)
# PAPER_FIG_TAG = tag of the step 10 and 11 outputs (e.g. _hsb_cc_g13);
# PAPER_FIG_OUT = output folder.
TAG = os.environ.get("PAPER_FIG_TAG", "")
FIG_OUT = os.environ.get("PAPER_FIG_OUT", BASE)
ELL = os.path.join(PIPE, f"recon_ellipses_0.44{TAG}.csv")
MTRK = os.path.join(PIPE, f"recon_model_tracks{TAG}.csv")
TAU_COV_MAX = 83.0   # last node with covariance (R&R chron 34y)

KM = 111.19492664455873
EXT = (-98, -65, -36, -5)
_wh = [pe.withStroke(linewidth=1.6, foreground="white")]

# Colours as in Figs. 1 and 4 (Dark2 + brown), validated for colour blindness with the four
# chains; a red/green pair fails (deutan dE 3.9) and a red Easter track is confused with the
# red of the epithermal deposits.
TRACKS = {14: ("Easter (Nazca R.)", "#d95f02"),
          15: ("San Félix (Taltal R.)", "#7570b3"),
          13: ("Juan Fernández R.", "#1b9e77"),
          16: ("Foundation (Iquique R.)", "#8c510a")}
RIDGE2PID = {"Easter": 14, "TT": 15, "JFR": 13, "Foundation": 16}

INTERVALS = [(72, 56), (56, 41.2), (41.2, 33.9), (33.9, 27.82),
             (27.82, 23.03), (23.03, 13.82), (13.82, 5.33), (5.33, 0)]
TITLES = ["Maastrichtian--Paleocene", "Early Eocene", "Late Eocene",
          "Early Oligocene", "Late Oligocene", "Early Miocene",
          "Late Miocene", "Pliocene--Quaternary"]

# Overlay of the revised circuit version (Wilson and DeMets 2025) on the reference
# one (Wilder 2003 via Quiero et al. 2022). The ELLIPSES are drawn as dashed outlines
# without fill AND the revised f_PAC=0.44 nominal track, also dashed. Step 11 computes
# the track from the rotations for any circuit version.
# Written to F3_reconstruction_overlay.* so that the reference version is NOT overwritten.
PIPE_WD = os.path.join(os.path.dirname(BASE), "models", os.environ.get("PAPER_PIPE_WD", "WD25_WK08-A"))
OVERLAY_WD = os.environ.get("PAPER_F3_OVERLAY", "0") == "1"
if OVERLAY_WD:
    OUT_SFX = "_overlay"
ell_wd = pd.read_csv(os.path.join(PIPE_WD, f"recon_ellipses_0.44{TAG}.csv"))
mtrk_wd = pd.read_csv(os.path.join(PIPE_WD, f"recon_model_tracks{TAG}.csv"),
                      dtype={"pf": str})

ell = pd.read_csv(ELL)
mtrk = pd.read_csv(MTRK, dtype={"pf": str})
PID2RIDGE = {14: "Easter", 15: "TT", 13: "JFR", 16: "Foundation"}
# keep only the ridges present in the outputs (Foundation may be absent)
TRACKS = {k: v for k, v in TRACKS.items() if PID2RIDGE[k] in set(mtrk.ridge)}
EM_STYLE = {"0.5": (2.2, 1.6), "0.36": (4.0, 1.2, 1.0, 1.2)}  # dash patterns
WK_STYLE = {"0.44": (5.0, 2.0), **EM_STYLE}                     # segments without covariance

# With PAPER_FIG_TAG the figure shows all models of the set (paper_models.py) equally, as
# Fig. 1: nominal f_PAC = 0.44 track of each model (thin line) and, at each node, the envelope
# (union) of their 68/95% ellipses. The end-member 0.50/0.36 tracks are not drawn in this
# version (with nine models they would be 27 lines per chain; the sensitivity to f_PAC is in
# Fig. 2).
MULTI = bool(TAG)
LIGERA = os.environ.get("PAPER_F3_LIGERA", "1") == "1"   # PAPER_F3_LIGERA=0: full version (track of each model, 68/95% envelopes)
if MULTI:
    import paper_models as MP
    OUT_SFX = MP.SUFIJO_SALIDA
    ell_m = pd.concat([pd.read_csv(os.path.join(os.path.dirname(BASE), "models", c,
                                                f"recon_ellipses_0.44{TAG}.csv")).assign(modelo=c)
                       for c in MP.CARPETAS])
    mtrk_m = pd.concat([pd.read_csv(os.path.join(os.path.dirname(BASE), "models", c,
                                                 f"recon_model_tracks{TAG}.csv"),
                                    dtype={"pf": str}).assign(modelo=c)
                        for c in MP.CARPETAS])
    # beyond the covariances only WK08-A resolves the track; there the median of the WK08-A
    # models of the set with f_PAC 0.44, 0.36 and 0.50 is drawn, dashed
    mtrk_wk = mtrk_m[mtrk_m.modelo.isin([c for c, _, m in MP.MODELOS if m == "WK08-A"])]
    mtrk_m = mtrk_m[mtrk_m.pf == "0.44"]
    # the single-model drawing (track, end-members and ellipses) is left without data
    mtrk = mtrk.iloc[0:0]
    ell = ell.iloc[0:0]


def envolvente(polys, n=300):
    """Union of (lon, lat) polygons on a grid: (X, Y, 0/1 mask). Same as in Fig. 1."""
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

hst_ages = sorted(float(os.path.splitext(os.path.basename(p))[0])
                  for p in glob.glob(os.path.join(HST_DIR, "*.csv")))
hst_ages = np.array(hst_ages)

def nearest_hst(age):
    return hst_ages[int(np.argmin(np.abs(hst_ages - age)))]

def load_track(age):
    a = nearest_hst(age)
    return pd.read_csv(os.path.join(HST_DIR, f"{a:.1f}.csv"))

mdf = pd.read_csv(MARGIN)
mdf["age"] = pd.to_numeric(mdf["age"], errors="coerce")
lon_cols = [c for c in mdf.columns if c.startswith("Lon_")]
lat_cols = [c for c in mdf.columns if c.startswith("Lat_")]

def margin_at(age):
    row = mdf.loc[(mdf["age"] - age).abs().idxmin()]
    lon = pd.to_numeric(row[lon_cols], errors="coerce").values.astype(float)
    lat = pd.to_numeric(row[lat_cols], errors="coerce").values.astype(float)
    ok = ~np.isnan(lon) & ~np.isnan(lat)
    return lon[ok], lat[ok]

g = pd.read_excel(GEOCHEM)
g["Age"] = pd.to_numeric(g["Age (Ma)"], errors="coerce")
g["Lat"] = pd.to_numeric(g["Latitude_restored"], errors="coerce")
g["Lon"] = pd.to_numeric(g["Longitude_restored"], errors="coerce")
g = g.dropna(subset=["Age", "Lat", "Lon"])

por = pd.read_csv(POR)
epi = pd.read_csv(EPI)
for d in (por, epi):
    d["Age"] = pd.to_numeric(d["Age MY"], errors="coerce")
    d["Lat"] = pd.to_numeric(d["Latitude_ROT"], errors="coerce")
    d["Lon"] = pd.to_numeric(d["Longitude_ROT"], errors="coerce")
por["Ton"] = pd.to_numeric(por["Tonnage Mt"], errors="coerce")
print(f"Porphyries {len(por)}, epithermal {len(epi)}, samples {len(g)}")

def ellipse_xy(mlo, mla, a_km, b_km, az_deg, n=72):
    th = np.linspace(0, 2 * np.pi, n)
    az = np.radians(az_deg)
    u = np.array([np.sin(az), np.cos(az)])
    v = np.array([-np.cos(az), np.sin(az)])
    E = np.outer(np.cos(th) * a_km, u) + np.outer(np.sin(th) * b_km, v)
    lon = mlo + E[:, 0] / (KM * np.cos(np.radians(mla)))
    lat = mla + E[:, 1] / KM
    return lon, lat

fig, axes = plt.subplots(3, 3, figsize=(11.5, 11.0), sharex=True, sharey=True)
axes = axes.flat

def in_iv(s, amax, amin):
    return (s["Age"] <= amax) & (s["Age"] >= amin) if amin == 0 else \
           (s["Age"] <= amax) & (s["Age"] > amin)

for i, (amax, amin) in enumerate(INTERVALS):
    ax = axes[i]
    amid = (amax + amin) / 2.0

    mlon, mlat = margin_at(amid)
    if len(mlon):
        order = np.argsort(mlat)
        ax.plot(mlon[order], mlat[order], "k-", lw=1.3, zorder=4)
        ax.fill_betweenx(mlat[order], -98, mlon[order], color="powderblue",
                         alpha=0.45, zorder=0)

    if MULTI and LIGERA:
        # light version: one track per chain (median of the models at each tau, only
        # tau <= 83 Ma), 68% envelope of the models' ellipses and a label every 20 Myr of
        # track age
        tk = mtrk_m[np.isclose(mtrk_m.t_recon, round(amid, 3), atol=0.02)
                    & (mtrk_m.tau <= TAU_COV_MAX)]
        for ridge, grp in tk.groupby("ridge"):
            col = TRACKS[RIDGE2PID[ridge]][1]
            md = grp.groupby("tau")[["lon", "lat"]].median().sort_index()
            ax.plot(md.lon, md.lat, "-", color=col, lw=1.4, zorder=3)
        tw = mtrk_wk[np.isclose(mtrk_wk.t_recon, round(amid, 3), atol=0.02)
                     & (mtrk_wk.tau >= TAU_COV_MAX)]
        for (ridge, pf), grp in tw.groupby(["ridge", "pf"]):
            col = TRACKS[RIDGE2PID[ridge]][1]
            md = grp.groupby("tau")[["lon", "lat"]].median().sort_index()
            ax.plot(md.lon, md.lat, ls=(0, WK_STYLE[pf]), color=col, lw=1.1, alpha=0.85, zorder=3)
        e = ell_m[np.isclose(ell_m.t_recon, round(amid, 3), atol=0.02)]
        e = e[e.tau % 10 == 0]
        for (ridge, tau), grp in e.groupby(["ridge", "tau"]):
            if not (EXT[0] - 1 <= grp.lon.mean() <= EXT[1] + 1
                    and EXT[2] - 1 <= grp.lat.mean() <= EXT[3] + 1):
                continue
            col = TRACKS[RIDGE2PID[ridge]][1]
            X6, Y6, M6 = envolvente([ellipse_xy(r.lon, r.lat, r.a68_km, r.b68_km,
                                                r.az_major_deg) for _, r in grp.iterrows()])
            ax.contourf(X6, Y6, M6, levels=[0.5, 1.5], colors=[col], alpha=0.22, zorder=2.1)
            ax.contour(X6, Y6, M6, levels=[0.5], colors=[col], linewidths=0.5, zorder=2.5)
            if tau % 20 == 0:
                ax.annotate(f"{int(tau)}", (grp.lon.median(), grp.lat.median()), xytext=(0, 0),
                            textcoords="offset points", color=col, fontsize=7.5,
                            fontweight="bold", ha="center", va="center", zorder=6.5,
                            annotation_clip=True,
                            path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
    elif MULTI:
        tk = mtrk_m[np.isclose(mtrk_m.t_recon, round(amid, 3), atol=0.02)]
        for (ridge, _), grp in tk.groupby(["ridge", "modelo"]):
            col = TRACKS[RIDGE2PID[ridge]][1]
            grp = grp.sort_values("tau")
            g1 = grp[grp.tau <= TAU_COV_MAX]
            ax.plot(g1.lon, g1.lat, "-", color=col, lw=0.7, alpha=0.9, zorder=3)
            g2 = grp[grp.tau >= TAU_COV_MAX]
            ax.plot(g2.lon, g2.lat, ls=(0, (1, 1.6)), color=col, lw=0.8, alpha=0.4, zorder=3)
        e = ell_m[np.isclose(ell_m.t_recon, round(amid, 3), atol=0.02)]
        drop = []
        for ridge, grp in e.groupby("ridge"):
            taus = sorted(set(grp.tau.values))
            for x in taus:
                if x % 10 != 0 and any(abs(x - y) <= 2.5 for y in taus if y != x):
                    drop.append((ridge, x))
        for (ridge, tau), grp in e.groupby(["ridge", "tau"]):
            if (ridge, tau) in drop:
                continue
            if not (EXT[0] - 1 <= grp.lon.mean() <= EXT[1] + 1
                    and EXT[2] - 1 <= grp.lat.mean() <= EXT[3] + 1):
                continue
            col = TRACKS[RIDGE2PID[ridge]][1]
            X, Y, M = envolvente([ellipse_xy(r.lon, r.lat, r.a95_km, r.b95_km, r.az_major_deg)
                                  for _, r in grp.iterrows()])
            ax.contourf(X, Y, M, levels=[0.5, 1.5], colors=[col], alpha=0.13, zorder=2)
            ax.contour(X, Y, M, levels=[0.5], colors=[col], linewidths=0.5, zorder=2.5)
            X6, Y6, M6 = envolvente([ellipse_xy(r.lon, r.lat, r.a68_km, r.b68_km,
                                                r.az_major_deg) for _, r in grp.iterrows()])
            ax.contourf(X6, Y6, M6, levels=[0.5, 1.5], colors=[col], alpha=0.30, zorder=2.1)
            ax.plot(grp.lon, grp.lat, "o", ms=1.6, color=col, zorder=6)
            ax.annotate(f"{int(tau)}", (grp.lon.max(), grp.lat.max()), xytext=(3, 2),
                        textcoords="offset points", color=col, fontsize=7, fontweight="bold",
                        ha="left", va="bottom", zorder=6.5, annotation_clip=True,
                        path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])

    # nominal track of the reference model f_PAC = 0.44 (step 11), computed from the
    # rotations like the end-members, so that the figure can be built with any circuit
    # version. Equivalence with the GPlates grid (HST_DIR) is checked in
    # models/Q22_WK08-A/validate_track_recon.py (max 13 cm over 2090 vertices).
    ref = mtrk[(mtrk.pf == "0.44")
               & np.isclose(mtrk.t_recon, round(amid, 3), atol=0.02)]
    for ridge, grp in ref.groupby("ridge"):
        col = TRACKS[RIDGE2PID[ridge]][1]
        grp = grp.sort_values("tau")
        # segment with covariances (tau <= 83 Ma): solid
        g1 = grp[grp.tau <= TAU_COV_MAX]
        if len(g1):
            ax.plot(g1.lon, g1.lat, "-", color=col, lw=1.6, alpha=0.9, zorder=3)
        # segment WITHOUT covariances (tau > 83 Ma, before R&R chron 34y):
        # dotted and faint = no ellipses can be computed
        g2 = grp[grp.tau >= TAU_COV_MAX]
        if len(g2):
            ax.plot(g2.lon, g2.lat, ls=(0, (1, 1.6)), color=col, lw=1.3,
                    alpha=0.45, zorder=3)

    # nominal tracks of the end-members 0.50 / 0.36 (step 11), faint
    for pf, dashes in EM_STYLE.items():
        mt = mtrk[(mtrk.pf == pf)
                  & np.isclose(mtrk.t_recon, round(amid, 3), atol=0.02)]
        for ridge, grp in mt.groupby("ridge"):
            col = TRACKS[RIDGE2PID[ridge]][1]
            grp = grp.sort_values("tau")
            g1 = grp[grp.tau <= TAU_COV_MAX]
            ax.plot(g1.lon, g1.lat, ls=(0, dashes), color=col, lw=1.1,
                    alpha=0.55, zorder=2.8)
            g2 = grp[grp.tau >= TAU_COV_MAX]
            ax.plot(g2.lon, g2.lat, ls=(0, (1, 1.6)), color=col, lw=1.0,
                    alpha=0.30, zorder=2.8)

    # uncertainty ellipses of the tracks at the panel age (step 10);
    # if the young-end node falls within 2.5 Ma of a multiple of 10 (49 vs
    # 50, 19 vs 20, 38 vs 40), it is skipped to avoid duplicating the ellipse and label
    e = ell[np.isclose(ell.t_recon, round(amid, 3), atol=0.02)].copy()
    drop = []
    for ridge, grp in e.groupby("ridge"):
        taus = sorted(grp.tau.values)
        for x in taus:
            if x % 10 != 0 and any(abs(x - y) <= 2.5 for y in taus if y != x):
                drop.append((ridge, x))
    if len(e):
        e = e[~e.apply(lambda r: (r.ridge, r.tau) in drop, axis=1)]
    for _, r in e.iterrows():
        if not (EXT[0] - 1 <= r.lon <= EXT[1] + 1
                and EXT[2] - 1 <= r.lat <= EXT[3] + 1):
            continue
        col = TRACKS[RIDGE2PID[r.ridge]][1]
        for a, b, alpha in ((r.a95_km, r.b95_km, 0.13),
                            (r.a68_km, r.b68_km, 0.25)):
            lo, la = ellipse_xy(r.lon, r.lat, a, b, r.az_major_deg)
            ax.fill(lo, la, color=col, alpha=alpha, lw=0, zorder=2)
        lo, la = ellipse_xy(r.lon, r.lat, r.a95_km, r.b95_km, r.az_major_deg)
        ax.plot(lo, la, color=col, lw=0.5, zorder=2.5)
        # dot at the node and readable age of the track point (Ma)
        ax.plot(r.lon, r.lat, "o", ms=2.6, color=col, mec="white", mew=0.5, zorder=6)
        ax.annotate(f"{int(r.tau)}", (r.lon, r.lat), xytext=(3, 2), textcoords="offset points",
                    color=col, fontsize=7, fontweight="bold", ha="left", va="bottom",
                    zorder=6.5, path_effects=[pe.withStroke(linewidth=2.2, foreground="white")],
                    annotation_clip=True)

    if OVERLAY_WD:
        # nominal f_PAC = 0.44 track of the revised circuit
        rw = mtrk_wd[(mtrk_wd.pf == "0.44")
                     & np.isclose(mtrk_wd.t_recon, round(amid, 3), atol=0.02)]
        for ridge, grp in rw.groupby("ridge"):
            col = TRACKS[RIDGE2PID[ridge]][1]
            grp = grp.sort_values("tau")
            g1 = grp[grp.tau <= TAU_COV_MAX]
            if len(g1):
                ax.plot(g1.lon, g1.lat, ls=(0, (2.5, 1.5)), color=col, lw=1.4,
                        alpha=0.95, zorder=3.1)
        ew = ell_wd[np.isclose(ell_wd.t_recon, round(amid, 3), atol=0.02)]
        ew = ew[~ew.apply(lambda r: (r.ridge, r.tau) in drop, axis=1)]
        for _, r in ew.iterrows():
            if not (EXT[0] - 1 <= r.lon <= EXT[1] + 1
                    and EXT[2] - 1 <= r.lat <= EXT[3] + 1):
                continue
            col = TRACKS[RIDGE2PID[r.ridge]][1]
            lo, la = ellipse_xy(r.lon, r.lat, r.a95_km, r.b95_km, r.az_major_deg)
            ax.plot(lo, la, color=col, lw=0.7, ls=(0, (2.5, 1.5)), zorder=3)

    # igneous samples from the geochemistry database, faint and below the tracks
    gs = g[in_iv(g, amax, amin)]
    ax.scatter(gs["Lon"], gs["Lat"], c="0.55", s=2.5, alpha=0.35, lw=0, zorder=1.5)

    ps = por[in_iv(por, amax, amin)]
    if len(ps) and MULTI and LIGERA:
        ax.scatter(ps["Lon"], ps["Lat"], marker="*", s=55, c="gold",
                   edgecolor="black", linewidth=0.5, zorder=10)
    elif len(ps):
        sz = 40 + 90 * (np.log10(ps["Ton"].clip(lower=1)) - 2).clip(lower=0)
        ax.scatter(ps["Lon"], ps["Lat"], marker="*", s=sz, c="gold",
                   edgecolor="black", linewidth=0.5, zorder=10)
    es = epi[in_iv(epi, amax, amin)]
    if len(es):
        ax.scatter(es["Lon"], es["Lat"], marker="+", s=60, c="red",
                   linewidth=1.4, zorder=9)

    ax.set_title(f"{chr(ord('A')+i)}) {TITLES[i]} ({amax}--{amin} Ma)\n"
                 f"tracks reconstructed at {amid:.1f} Ma", fontsize=8.5)
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3])
    ax.set_xticks(np.arange(-95, -64, 10))
    ax.grid(True, ls=":", alpha=0.5)
    if i % 3 == 0:
        ax.set_ylabel("Latitude (°)", fontsize=9)
    if i >= 6:
        ax.set_xlabel("Longitude (°)", fontsize=9)

lax = axes[8]; lax.axis("off")
handles = [Line2D([0], [0], color="k", lw=1.3, label="Trench")]
for pid, (name, col) in TRACKS.items():
    handles.append(Line2D([0], [0], color=col, lw=1.8, label=name))
handles += [
    Line2D([0], [0], color="0.35", lw=1.1, ls=(0, EM_STYLE["0.5"]),
           alpha=0.8, label="$f_\\mathrm{PAC}=0.50$ model (dashed)"),
    Line2D([0], [0], color="0.35", lw=1.1, ls=(0, EM_STYLE["0.36"]),
           alpha=0.8, label="$f_\\mathrm{PAC}=0.36$ model (dash-dot)"),
    Line2D([0], [0], color="0.35", lw=1.3, ls=(0, (1, 1.6)), alpha=0.6,
           label="Track older than 83 Ma: no\ncovariance data, no ellipses"),
    Patch(facecolor="0.5", alpha=0.3,
          label="68 / 95% ellipses at 10-Myr\ntrack-age nodes (label = age of\nthe track point, Ma)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="0.55", alpha=0.5,
           markersize=4, label="Restored igneous samples"),
    Line2D([0], [0], marker="*", color="w", markerfacecolor="gold",
           markeredgecolor="k", markersize=14,
           label="Porphyry (size $\\propto$ tonnage)"),
    Line2D([0], [0], marker="+", color="red", lw=0, markeredgecolor="red",
           markersize=11, markeredgewidth=1.5, label="Epithermal"),
]
if OVERLAY_WD:
    handles.insert(4, Line2D([0], [0], color="0.35", lw=1.4,
                             ls=(0, (2.5, 1.5)),
                             label="Revised circuit (Wilson and\n"
                                   "DeMets 2025): track and 95%\n"
                                   "ellipse outlines (dashed)"))
if MULTI:   # no end-members; envelopes of the models
    handles = [h for h in handles if "model (dash" not in h.get_label()]
    handles = [h for h in handles if not h.get_label().startswith("68 / 95% ellipses")]
    handles.insert(1 + len(TRACKS), Line2D([0], [0], color="0.35", lw=0.8,
                   label=f"Nominal track ($f_\\mathrm{{PAC}}=0.44$),\n"
                         f"each of the {len(MP.CARPETAS)} models"))
    handles.insert(3 + len(TRACKS), Patch(facecolor="0.5", alpha=0.3,
                   label=f"Envelope of the {len(MP.CARPETAS)} models'\n"
                         "68 / 95% ellipses at 10-Myr\n"
                         "track-age nodes (label = age\nof the track point, Ma)"))
if MULTI and LIGERA:
    n = len(MP.CARPETAS)
    handles = [Line2D([0], [0], color="k", lw=1.3, label="Trench")]
    for pid, (name, col) in TRACKS.items():
        handles.append(Line2D([0], [0], color=col, lw=1.6, label=name))
    handles += [
        Patch(facecolor="0.5", alpha=0.35, edgecolor="0.4", lw=0.5,
              label=f"Line: median track of the {n} models\n"
                    f"($f_\\mathrm{{PAC}}=0.44$). Patches: envelope\n"
                    f"of their 68% ellipses every 10 Myr\n"
                    "of track age (labelled every 20 Myr)"),
        *[Line2D([0], [0], color="0.35", lw=1.1, ls=(0, WK_STYLE[pf]),
                 label=lab) for pf, lab in (
            ("0.44", f"Track older than {TAU_COV_MAX:g} Ma (no covariances):\n"
                     "median of the WK08-A models, $f_\\mathrm{PAC}=0.44$"),
            ("0.36", "same, $f_\\mathrm{PAC}=0.36$"), ("0.5", "same, $f_\\mathrm{PAC}=0.50$"))],
        Line2D([0], [0], marker="o", color="w", markerfacecolor="0.55", alpha=0.5,
               markersize=4, label="Restored igneous samples"),
        Line2D([0], [0], marker="*", color="w", markerfacecolor="gold",
               markeredgecolor="k", markersize=11, label="Porphyry"),
        Line2D([0], [0], marker="+", color="red", lw=0, markeredgecolor="red",
               markersize=10, markeredgewidth=1.5, label="Epithermal")]
lax.legend(handles=handles, loc="center", fontsize=9, frameon=False)

fig.subplots_adjust(wspace=0.08, hspace=0.28, top=0.97)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG_OUT, f"F3_reconstruction{OUT_SFX}.{ext}"),
                bbox_inches="tight", dpi=300)
plt.close(fig)
print("F5_reconstruction_unc: OK")
