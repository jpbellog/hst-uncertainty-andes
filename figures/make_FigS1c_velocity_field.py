# -*- coding: utf-8 -*-
"""
Figure S1c / S1d (Supplementary) -- convergence VELOCITY field of
Nazca/Farallon-South America, three panels, with the same age and latitude
limits as Fig. S1B:

  (A) absolute magnitude |v|
  (B) component orthogonal to the margin (positive = convergence)
  (C) component parallel to the margin, signed (positive = southward, the
      direction in which the margin nodes run)

It complements Fig. S1B, which shows the uncertainty of the velocity: this
figure shows the velocity itself.

Source of the numbers. They are not read from vc_sigma_*.csv, because those files
cover only the range with published covariances (27.9 or 33.7 Ma) and do not
store the parallel component (vc_ort is stored as an absolute value and the
parallel component is absent). They are computed here directly from the rotations,
which is deterministic and needs no Monte Carlo:

  - stage Euler vector of plate 902003 (f_PAC = 0.44) relative to plate 201
    (South America) over 1 Myr windows [t, t-1], as in step 09;
  - v = w x r evaluated at the nodes of the present-day margin (the same 43 nodes
    between 5 and 36 S used by step 05), so the latitude axis is present-day
    latitude over the whole range and does not change convention within the figure;
  - decomposition with the local azimuth of the margin (node j -> j+10), the same
    convention as steps 05 and 09.

The 1 Myr window is finer than the published stages, which is possible because
there are no covariances to propagate here: it is the nominal velocity.

PAPER_PIPE selects the circuit version:
  py make_FigS1c_velocity_field.py                 -> FigS1c_velocity_field.*
  PAPER_PIPE=WD25_WK08-A py ...                  -> FigS1c_velocity_field_wd25.*

Usage:  py make_FigS1c_velocity_field.py
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm"})

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
TABLES = os.path.join(ROOT, "tables")
PIPE_NAME = os.environ.get("PAPER_PIPE", "Q22_WK08-A")
IS_WD25 = PIPE_NAME != "Q22_WK08-A"
OUT_SFX = "_wd25" if IS_WD25 else ""
PIPE = os.path.join(ROOT, "models", PIPE_NAME)
# With PAPER_FIG_TAG the circuit is taken from PAPER_PIPE among the three of the grid
# (Q22_WK08-A, WD25_WK08-A, SG12_WK08-A) and its _hsb .rot (GTS2020) is used. The
# Nazca-South America velocity does not depend on the absolute frame. The collision lines are
# those of the three frames of the main text with that circuit (paper_models.py), clipped by
# frame_coverage.
CIRC_ROT = {"Q22_WK08-A": ("", "Global_Model_R_BG_2025_hsb.rot"),
            "WD25_WK08-A": ("_wd25", "Global_Model_R_BG_2025_WD25_hsb.rot"),
            "SG12_WK08-A": ("_sg12", "Global_Model_R_BG_2025_SG12_hsb.rot")}
if os.environ.get("PAPER_FIG_TAG", ""):
    OUT_SFX, _rot = CIRC_ROT[PIPE_NAME]
    os.environ["PAPER_ROT_FILE"] = os.path.join(PIPE, "rot", _rot)
elif IS_WD25:
    os.environ.setdefault(
        "PAPER_ROT_FILE",
        os.path.join(PIPE, "rot", "Global_Model_R_BG_2025_WD25_refit.rot"))
sys.path.insert(0, PIPE)
eng = importlib.import_module("02_rot_engine")

MARGIN_CSV = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "margin", "Margin_SAM_0.1.csv")  # Margin_SAM (see rotation_models/README.md)
R_E_CM_YR = 6371e5 / 1e6
D2R = np.pi / 180.0
PLATE, REL, PF = 902003, 201, "0.44"
AGES = np.arange(1.0, 67.0, 1.0)
AGE_MAX = 66.0
LAT_LO, LAT_HI = -36.0, -5.0
_halo = [pe.withStroke(linewidth=2.4, foreground="black")]


def unit(lat, lon):
    la, lo = lat * D2R, lon * D2R
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo),
                     np.sin(la)])


def margin_nodes():
    m = pd.read_csv(MARGIN_CSV)
    row = m.iloc[(pd.to_numeric(m["age"], errors="coerce") - 0).abs().idxmin()]
    nn = sum(1 for c in m.columns if c.startswith("Lon_"))
    lon = np.array([row.get(f"Lon_{k}") for k in range(1, nn + 1)], float)
    lat = np.array([row.get(f"Lat_{k}") for k in range(1, nn + 1)], float)
    k = np.isfinite(lat) & np.isfinite(lon) & (lat <= -5) & (lat >= -36)
    return lat[k], lon[k]


def stage_w(M, t):
    """Euler vector (rad/Myr) of PLATE relative to REL over [t, t-1]."""
    def rel_R(tau):
        R1 = eng.quat_to_matrix(M.total_rotation(REL, tau, 0))
        R2 = eng.quat_to_matrix(M.total_rotation(PLATE, tau, 0))
        return R1.T @ R2
    S = rel_R(t - 1.0) @ rel_R(t).T
    ang = np.arccos(np.clip((np.trace(S) - 1) / 2, -1, 1))
    if ang < 1e-12:
        return np.zeros(3)
    ax = np.array([S[2, 1] - S[1, 2], S[0, 2] - S[2, 0], S[1, 0] - S[0, 1]])
    return ax / np.linalg.norm(ax) * ang


def margin_azimuth(mlat, mlon, j):
    """Local azimuth of the margin at node j, southward (rad).

    Convention of steps 05 and 09: forward difference node j -> j+10.
    Note: those steps use `j2 = min(j+10, n-1)`, which at the last node gives
    j2 == j, i.e. dE = dN = 0 and azimuth 0 (north) instead of the margin azimuth.
    Since both store the orthogonal component as an absolute value, the effect is
    hidden there; here, with sign, that node stands out. It is solved with a
    backward difference when the forward one degenerates.
    """
    n = len(mlat)
    j2 = min(j + 10, n - 1)
    ja, jb = (j, j2) if j2 != j else (max(j - 10, 0), j)
    if ja == jb:
        return 0.0
    dE = (mlon[jb] - mlon[ja]) * np.cos(mlat[ja] * D2R)
    dN = mlat[jb] - mlat[ja]
    return float(np.arctan2(dE, dN))


def decompose(w, mlat, mlon):
    """|v|, orthogonal component and parallel component (signed), cm/yr.

    Signs: v_ort > 0 = convergence (the plate moves toward the continent);
    v_par > 0 = southward, the direction in which the margin nodes run.
    Steps 05 and 09 store abs(v_ort), so they fix no sign; with the normal
    defined as (cos azm, -sin azm) convergence comes out negative (that normal
    points to the ocean), hence the minus sign of n_hat here.
    """
    n = len(mlat)
    vmag = np.empty(n); vort = np.empty(n); vpar = np.empty(n)
    for j in range(n):
        la, lo = mlat[j], mlon[j]
        v3 = np.cross(w, unit(la, lo)) * R_E_CM_YR
        lar, lor = la * D2R, lo * D2R
        e = np.array([-np.sin(lor), np.cos(lor), 0.0])
        nv = np.array([-np.sin(lar) * np.cos(lor), -np.sin(lar) * np.sin(lor),
                       np.cos(lar)])
        vE, vN = v3 @ e, v3 @ nv
        azm = margin_azimuth(mlat, mlon, j)
        t_hat = np.array([np.sin(azm), np.cos(azm)])      # along the margin, southward
        n_hat = np.array([-np.cos(azm), np.sin(azm)])     # normal, toward the continent
        vmag[j] = np.hypot(vE, vN)
        vort[j] = np.array([vE, vN]) @ n_hat
        vpar[j] = np.array([vE, vN]) @ t_hat
    return vmag, vort, vpar


# ------------------------------------------------------------------ field
M = eng.RotModel()
mlat, mlon = margin_nodes()
order = np.argsort(mlat)
MAG = np.empty((len(mlat), len(AGES)))
ORT = np.empty((len(mlat), len(AGES)))
PAR = np.empty((len(mlat), len(AGES)))
for i, t in enumerate(AGES):
    vm, vo, vp = decompose(stage_w(M, t), mlat, mlon)
    MAG[:, i], ORT[:, i], PAR[:, i] = vm[order], vo[order], vp[order]
lat_s = mlat[order]
y_edges = np.concatenate([[lat_s[0] - 0.5 * (lat_s[1] - lat_s[0])],
                          0.5 * (lat_s[:-1] + lat_s[1:]),
                          [lat_s[-1] + 0.5 * (lat_s[-1] - lat_s[-2])]])
x_edges = np.concatenate([AGES - 1.0, [AGES[-1]]])

print(f"[{PIPE_NAME}] |v|  {np.nanmin(MAG):.2f}-{np.nanmax(MAG):.2f} cm/yr "
      f"(median {np.nanmedian(MAG):.2f})")
print(f"[{PIPE_NAME}] v_ort {np.nanmin(ORT):.2f}-{np.nanmax(ORT):.2f} "
      f"(median {np.nanmedian(ORT):.2f})")
print(f"[{PIPE_NAME}] v_par {np.nanmin(PAR):.2f}-{np.nanmax(PAR):.2f} "
      f"(median {np.nanmedian(PAR):.2f}); southward fraction "
      f"{100*np.mean(PAR > 0):.0f}%")
np.savetxt(os.path.join(os.environ.get("PAPER_FIG_OUT", BASE), f"_FigS1c_field{OUT_SFX}.csv"),
           np.column_stack([np.repeat(lat_s, len(AGES)),
                            np.tile(AGES - 0.5, len(lat_s)),
                            MAG.ravel(), ORT.ravel(), PAR.ravel()]),
           delimiter=",", header="lat,age_mid,vc_mag,vc_ort,vc_par",
           comments="", fmt="%.4f")

# ------------------------------------------------------------- lines and deposits
# PAPER_FIG_TAG = tag of the bands (e.g. _hsb_cc_g13), PAPER_TABLES_DIR = folder of the bands,
# PAPER_FIG_OUT = output folder. Without them, the bands of the PAPER_PIPE model are used.
# Foundation (fourth chain) is included only if it is present in the bands.
TAG = os.environ.get("PAPER_FIG_TAG", "")
TABLES = os.environ.get("PAPER_TABLES_DIR", TABLES)
FIG_OUT = os.environ.get("PAPER_FIG_OUT", BASE)
_LABS = [("Easter", "Easter (Nazca R.)", 13),
         ("JFR", "Juan Fernández R.", 21),
         ("TT", "San Félix (Taltal R.)", 16),
         ("Foundation", "Foundation (Iquique R.)", 40)]
if TAG:
    import paper_models as MP
    _circ = {"Q22_WK08-A": "Quiero", "WD25_WK08-A": "WD25", "SG12_WK08-A": "SG12"}[PIPE_NAME]
    _mods = [c for c, ci, _ in MP.MODELOS if ci == _circ]
    _tots = [MP.recortar(pd.read_csv(os.path.join(
                 TABLES, f"mc_collision_bands_total_present{MP.sufijo_bandas(c)}{TAG}.csv"),
                 dtype={"pf": str}), c) for c in _mods]
else:
    _tots = [pd.read_csv(os.path.join(TABLES,
                                      f"mc_collision_bands_total_present_{PIPE_NAME}{TAG}.csv"),
                         dtype={"pf": str})]
# coll[lab] = list of (ages, latitudes, anchor), one per model
coll = {}
for _tot in _tots:
    _tot = _tot[_tot.pf == PF]
    for k, lab, anc in _LABS:
        c = _tot[(_tot.ridge == k) & (_tot.age <= AGE_MAX)].sort_values("age")
        c = c.dropna(subset=["lat_nominal"])
        if len(c):          # some band tables do not include Foundation
            coll.setdefault(lab, []).append((c["age"].values, c["lat_nominal"].values, anc))


def label_on_line(ax, a, la, lab, anchor):
    j = int(np.argmin(np.abs(a - anchor)))
    j0, j1 = max(0, j - 2), min(len(a) - 1, j + 2)
    p0 = ax.transData.transform((a[j0], la[j0]))
    p1 = ax.transData.transform((a[j1], la[j1]))
    ang = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
    if ang > 90: ang -= 180
    if ang < -90: ang += 180
    ax.text(a[j], la[j], lab, color="white", fontsize=6.8, fontweight="bold",
            va="bottom", ha="center", rotation=ang, rotation_mode="anchor",
            zorder=6, path_effects=_halo)


# ------------------------------------------------------------------ figure
_base = plt.cm.jet(np.linspace(0, 1, 256))
_pastel = 0.55 * _base + 0.45 * np.ones_like(_base)
_pastel[:, 3] = 1.0
CMAP = LinearSegmentedColormap.from_list("jet_pastel", _pastel)

vmax = 3 * np.ceil(np.nanmax(MAG) / 3.0)
pmax = np.ceil(max(abs(np.nanmin(PAR)), abs(np.nanmax(PAR))))
PANELS = [("A", MAG, "Absolute rate $|v|$", CMAP,
           dict(vmin=0, vmax=vmax), "$|v|$ (cm/yr)"),
          ("B", ORT, "Trench-orthogonal component $v_\\perp$", CMAP,
           dict(vmin=0, vmax=vmax), "$v_\\perp$ (cm/yr)"),
          ("C", PAR, "Trench-parallel component $v_\\parallel$ "
                     "(positive southward)", "RdBu_r",
           dict(norm=TwoSlopeNorm(vmin=-pmax, vcenter=0.0, vmax=pmax)),
           "$v_\\parallel$ (cm/yr)")]

fig, axes = plt.subplots(3, 1, figsize=(7.4, 10.4), sharex=True, sharey=True)
fig.subplots_adjust(left=0.10, right=0.90, top=0.97, bottom=0.06, hspace=0.16)

for ax, (tag, Z, title, cmap, kw, cbl) in zip(axes, PANELS):
    pcm = ax.pcolormesh(x_edges, y_edges, Z, cmap=cmap, shading="flat",
                        zorder=0, rasterized=True, **kw)
    for lab, lines in coll.items():
        for a, la, anc in lines:
            ax.plot(a, la, "-", color="white", lw=1.6 if len(lines) == 1 else 1.0, zorder=3,
                    path_effects=_halo)
        # label on the line that comes closest to the anchor
        a, la, anc = min(lines, key=lambda x: np.min(np.abs(x[0] - x[2])))
        label_on_line(ax, a, la, lab, anc)
    ax.set_xlim(AGE_MAX, 0); ax.set_ylim(LAT_LO, LAT_HI)
    ax.set_ylabel("Latitude (°)", fontsize=11)
    ax.set_title(f"{tag})  {title}", fontsize=11.5, fontweight="bold", loc="left")
    cb = fig.colorbar(pcm, ax=ax, pad=0.015, fraction=0.05,
                      extend="max" if tag != "C" else "both")
    cb.set_label(cbl, fontsize=10)
axes[-1].set_xlabel("Age (Ma)", fontsize=11)
axes[0].legend(handles=[Line2D([0], [0], color="0.35", lw=2.0,
                               label="Collision latitude" + (", each frame" if TAG else ""))],
               loc="lower left", fontsize=7.5, framealpha=0.9, edgecolor="0.7")

for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG_OUT, f"FigS1c_velocity_field{OUT_SFX}.{ext}"),
                dpi=350)
plt.close(fig)
print(f"FigS1c_velocity_field{OUT_SFX}: OK")
