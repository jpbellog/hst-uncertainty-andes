"""Control figure of the bathymetric test with the 12 models (3 circuits x 4 frames).

Row 1: maps of the NE half of the Nazca Ridge, the Taltal Ridge (San Felix track) and the
Iquique Ridge (Foundation track; its axis joins the highest summits). Background: unsmoothed ETOPO1, to judge by eye whether the axis
follows the chain; contours of the relief smoothed with the FWHM used by each axis. Points: local
maxima of each column (find_peaks, all candidates); black line: accepted axis
(bathymetry_axes.py). Color = frame, line style = circuit. Q22 and SG12 coincide on the
Nazca plate, so the dotted SG12 line falls on the solid Q22 line.

Row 2: median distance from the axis to the nominal track per model (outputs/bathymetry_axes/
resumen.csv; Nazca: segment without ages far from the trench; Taltal and Iquique: the whole axis).

Run from analysis/ after bathymetry_axes.py:  py -I ridge_axes_figure.py
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bathymetry_axes as B  # noqa: E402
from scipy.ndimage import gaussian_filter  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.colors import LightSource  # noqa: E402

COL = {"WK08-A": "#2a78d6", "GMHRF-B": "#eb6834", "GWG20": "#1baf7a", "GMHRF-A'": "#eda100"}
EST = {"Quiero": "-", "WD25": "--", "SG12": ":"}
NOM_C = {"Quiero": "Q22", "WD25": "WD25", "SG12": "SG12"}
PANELES = [  # ridge, summary segment, map window, title
    ("Nazca", "sin_datos_lejos_fosa", (-82.6, -75.0, -24.0, -14.4),
     "Nazca Ridge, NE half (Easter track)"),
    ("Taltal", "todo", (-80.6, -70.6, -28.6, -23.0), "Taltal Ridge (San Félix track)"),
    # window extended to the SW to show the northern segment of the Taltal axis
    ("Iquique", "todo", (-86.5, -70.0, -36.5, -19.2), "Iquique Ridge (Foundation track)"),
]


def main():
    lon, lat, z = B.etopo_box()
    zs = {fw: gaussian_filter(z, sigma=fw / 1.85 / 2.355)
          for fw in sorted({c["fwhm"] for c in B.DORSALES.values()})}
    fla, flo = B.fosa()
    res = pd.read_csv(os.path.join(B.OUT, "resumen.csv"))
    modelos = B.TODOS
    tracks = {m[0]: pd.read_csv(os.path.join(B.ROOT, "models", m[0], f"track_nominal_0.44{B.TAG}.csv"))
              for m in modelos}

    fig = plt.figure(figsize=(16, 10.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[2.3, 1], hspace=0.22, wspace=0.18)
    ls = LightSource(azdeg=315, altdeg=35)
    for k, (nom, tramo, (x0, x1, y0, y1), titulo) in enumerate(PANELES):
        cfg = B.DORSALES[nom]
        # the axis comes from ejes_<ridge>.csv (Taltal already clipped by the Iquique corridor);
        # the candidate maxima are recomputed only to plot them
        e_csv = pd.read_csv(os.path.join(B.OUT, f"ejes_{nom}.csv"))
        e = e_csv
        if cfg["modo"] == "montes":
            crud = pd.read_csv(os.path.join(B.OUT, "cumbres_Iquique.csv"))
        else:
            _, crud = B.eje(lon, lat, zs, cfg, fla, flo)
        ax = fig.add_subplot(gs[0, k])
        ii = (lon >= x0) & (lon <= x1)
        jj = (lat >= y0) & (lat <= y1)
        Zr = z[np.ix_(jj, ii)]
        sh = ls.shade(Zr, cmap=plt.get_cmap("Greys_r"), vmin=-6000, vmax=-1000, vert_exag=0.02,
                      blend_mode="soft")
        ax.imshow(sh, extent=(lon[ii][0], lon[ii][-1], lat[jj][-1], lat[jj][0]), origin="upper",
                  zorder=0, interpolation="bilinear")
        Zs = zs[cfg["fwhm"]][np.ix_(jj, ii)]
        ax.contour(lon[ii], lat[jj], Zs, levels=np.arange(-4500, -1000, 500), colors="0.25",
                   linewidths=0.35, zorder=1)
        ax.plot(flo, fla, color="k", lw=1.0, zorder=2)
        if cfg["modo"] != "montes":
            ax.plot(crud.lon, crud.lat, ".", ms=1.6, color="#c400c4", alpha=0.6, zorder=3)
        seg = (e.col.diff().abs() > B.HUECO_COL * B.D * 1.01).cumsum()
        for _, g in e.groupby(seg):
            ax.plot(g.lon, g.lat, "-", color="white", lw=4.2, zorder=6, solid_capstyle="round")
            ax.plot(g.lon, g.lat, "-", color="k", lw=2.2, zorder=7, solid_capstyle="round")
        if "con_datos" in e_csv and e_csv.con_datos.any():
            d = e_csv[e_csv.con_datos]
            ax.plot(d.lon, d.lat, "o", ms=2.2, color="white", mec="k", mew=0.4, zorder=8)
        if nom == "Iquique":          # other ridges in the window, in grey; corridor
            for o in ("Taltal", "JFR"):
                et = pd.read_csv(os.path.join(B.OUT, f"ejes_{o}.csv"))
                segt = (et.col.diff().abs() > B.HUECO_COL * B.D * 1.01).cumsum()
                for _, g in et.groupby(segt):
                    ax.plot(g.lon, g.lat, "-", color="white", lw=3.6, zorder=6)
                    ax.plot(g.lon, g.lat, "-", color="0.45", lw=1.8, zorder=7)
            cla, clo = B.linea_densa(*B.COPIAPO)
            ax.plot(clo, cla, color="#fdd835", lw=1.8, ls=(0, (5, 2)), zorder=7)
            V = B.vertices_corredor(cfg)
            for A, C in zip(V[:-1], V[1:]):   # corridor edges, segment by segment
                PA, PC = B.xyz(*A), B.xyz(*C)
                n = np.cross(PA, PC); n /= np.linalg.norm(n)
                cla2, clo2 = B.linea_densa(A, C)
                P = B.xyz(cla2, clo2)
                for sg in (-1, 1):
                    ang = sg * cfg["semiancho_km"] / B.R_KM
                    Q = P * np.cos(ang) + n * np.sin(ang)
                    ax.plot(np.degrees(np.arctan2(Q[:, 1], Q[:, 0])), np.degrees(np.arcsin(Q[:, 2])),
                            color="0.15", lw=0.7, ls=(0, (4, 3)), zorder=4)
            fz = pd.read_csv(B.FZ_CHALLENGER, comment="#")
            ax.plot(fz.lon, fz.lat, color="#1e88e5", lw=1.4, ls=(0, (6, 2, 1, 2)), zorder=7)
        if cfg["modo"] == "montes":   # unsmoothed line, summits and excluded summits
            ax.plot(e.lon_bruta, e.lat_bruta, "-", color="0.15", lw=0.6, alpha=0.7, zorder=6.5)
            u = crud[~crud.excluida]; x = crud[crud.excluida]
            ax.plot(x.lon, x.lat, "^", ms=6, mfc="none", mec="#c400c4", mew=1.0, zorder=10)
            ax.plot(u.lon, u.lat, "^", ms=7, color="#c400c4", mec="white", mew=0.8, zorder=10)
        for carpeta, circ, marco in modelos:
            t = tracks[carpeta]
            t = t[(t.ridge == cfg["track"])].dropna(subset=["lat", "lon"]).sort_values("tau")
            w = 1.6 if circ != "SG12" else 1.3
            ax.plot(t.lon, t.lat, EST[circ], color=COL[marco], lw=w, zorder=9 if circ == "SG12" else 8.5)
        cx = cfg["caja"]
        if cfg["modo"] != "montes":
            ax.plot([cx[0], cx[1], cx[1], cx[0], cx[0]], [cx[2], cx[2], cx[3], cx[3], cx[2]],
                    color="0.15", lw=0.7, ls=(0, (4, 3)), zorder=4)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect(1 / np.cos(np.radians((y0 + y1) / 2)))
        ax.set_title(titulo, fontsize=10, loc="left")
        ax.tick_params(labelsize=7)
        ax.grid(color="white", lw=0.3, alpha=0.4)

        # row 2: median per model
        bx = fig.add_subplot(gs[1, k])
        r = res[(res.dorsal == nom) & (res.tramo == tramo)].set_index(["circuito", "marco"])
        ylab, y = [], 0
        for circ in ("Quiero", "WD25", "SG12"):
            for marco in ("WK08-A", "GMHRF-B", "GWG20", "GMHRF-A'"):
                q = r.loc[(circ, marco)]
                bx.plot([0, q.mediana_km], [y, y], color="0.85", lw=1, zorder=1)
                mk = "o" if marco != "GMHRF-A'" else "D"
                bx.plot(q.mediana_km, y, mk, ms=6, color=COL[marco], mec="white", mew=0.8, zorder=2)
                bx.text(q.mediana_km, y, f"  {q.mediana_km:.0f}", va="center", fontsize=6.5,
                        color="0.2")
                ylab.append(f"{NOM_C[circ]} / {marco}")
                y += 1
            y += 0.5
        ys = [i + 0.5 * (i // 4) for i in range(12)]
        bx.set_yticks(ys)
        bx.set_yticklabels(ylab, fontsize=6.5)
        bx.invert_yaxis()
        bx.set_xlim(0, r.mediana_km.max() * 1.18)
        n = int(r.n.iloc[0])
        bx.set_xlabel(f"median distance, track to axis (km); {n} axis points", fontsize=7.5)
        bx.tick_params(axis="x", labelsize=7)
        for s in ("top", "right"):
            bx.spines[s].set_visible(False)

    h = [Line2D([], [], color=COL[m], lw=2, label=m) for m in COL]
    h += [Line2D([], [], color="0.3", ls=EST[c], lw=1.6, label=NOM_C[c] + " circuit") for c in EST]
    h += [Line2D([], [], color="k", lw=2.2, label="ridge axis (this test)"),
          Line2D([], [], ls="", marker=".", color="#c400c4", label="local maxima (all candidates)"),
          Line2D([], [], ls="", marker="o", ms=4, color="white", mec="k", label="axis < 75 km from a dated sample"),
          Line2D([], [], ls="", marker="^", ms=6, color="#c400c4", mec="white", label="Iquique: summits joined by the axis"),
          Line2D([], [], ls="", marker="^", ms=6, mfc="none", mec="#c400c4", label="Iquique: summits excluded (near another ridge)"),
          Line2D([], [], color="0.45", lw=1.8, label="Taltal and Juan Fernández axes (Iquique map)"),
          Line2D([], [], color="#fdd835", lw=1.8, ls=(0, (5, 2)), label="Copiapó Ridge (Iquique map)"),
          Line2D([], [], color="#1e88e5", lw=1.4, ls=(0, (6, 2, 1, 2)), label="Challenger FZ (GSFML)"),
          Line2D([], [], color="k", lw=1, label="trench (margin line, 0 Ma)"),
          Line2D([], [], color="0.15", lw=0.7, ls=(0, (4, 3)), label="search box (Iquique: corridor)")]
    fig.legend(handles=h, loc="upper center", ncol=7, fontsize=7.5, frameon=False,
               bbox_to_anchor=(0.5, 0.995))
    fig.text(0.5, 0.012, "Q22 and SG12 give identical tracks on the Nazca plate (SG12 dotted on top of Q22). "
             "Background: ETOPO1, unsmoothed; contours every 500 m of the smoothed relief used for the axis "
             "(FWHM 40 km). Iquique axis: summits >= 1000 m above the local median depth (100 km), within 75 km of a line from 20.5 S, 71.8 W through 27.5 S, 77.3 W to the Challenger FZ, joined along it and smoothed with a 150 km running median and a 100 km running mean.", ha="center", fontsize=7, color="0.3")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(B.OUT, f"ejes_12_modelos.{ext}"), dpi=170, bbox_inches="tight")
    print("ok", os.path.join(B.OUT, "ejes_12_modelos.png"))


if __name__ == "__main__":
    main()
