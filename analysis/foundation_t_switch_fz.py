# -*- coding: utf-8 -*-
"""Age of the Foundation plate change (t_switch) as the passage of the hotspot from north to south
across the Challenger FZ (Nazca side) / Resolution FZ (Pacific side), criterion of Bello-Gonzalez et
al. (2018, Gondwana Res. 64, p. 9: "passes the Challenger (Resolution) FZ from the north to
south at ~25 Ma").

The position of the hotspot relative to a plate at instant tau is point tau of the track in
the frame of that plate, so the crossing is the intersection, in present-day coordinates of the
plate, of the track (extended to fine tau) with the FZ trace:
- Nazca side: track on 902003 vs GSFML trace FZ_RM_0002-0189 (unnamed in GSFML; it is
  identified as the Challenger FZ because it is the continuous trace that reaches the margin at ~33S,
  where Yanez et al. 2002 place the Challenger; our own identification, to be verified).
- Pacific side: track on 901 vs GSFML trace FZ_KM_0001-0671 "Resolution".

Usage:  py foundation_t_switch_fz.py
Output: outputs/foundation_t_switch_fz.csv and map outputs/iquique_prueba_03b/foundation_fz_cruce.png
"""
import os
import sys
import numpy as np
import shapefile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import track_vs_ridge_axis as r219

GSFML = (__import__("os").environ["PAPER_ROOT"] + "/external/gsfml/SF_shp")
FZ = {"Challenger (Nazca)": ("GSFML_SF_FZ_RM", "FZ_RM_0002-0189", 902003),
      "Resolution (Pacifico)": ("GSFML_SF_FZ_KM", "FZ_KM_0001-0671", 901)}
SEED = (-37.3764, -110.888)
TAU = np.round(np.arange(10.0, 45.0001, 0.05), 3)


def fz_trace(fname, fid):
    r = shapefile.Reader(os.path.join(GSFML, fname))
    for s, rec in zip(r.shapes(), r.records()):
        if rec[0] == fid:
            return np.array(s.points)          # (lon, lat)
    raise KeyError(fid)


def crossings(lo, la, F):
    """Crossings of the polyline (lo, la)[tau] with the polyline F; returns [(tau, lat, lon, sign)]
    with sign = +1 if the hotspot passes from N to S of the FZ as age increases towards the present."""
    out = []
    for i in range(len(lo) - 1):
        p, r = np.array([lo[i], la[i]]), np.array([lo[i + 1] - lo[i], la[i + 1] - la[i]])
        if not np.all(np.isfinite(np.r_[p, r])):
            continue
        for j in range(len(F) - 1):
            q, s = F[j], F[j + 1] - F[j]
            den = r[0] * s[1] - r[1] * s[0]
            if abs(den) < 1e-14:
                continue
            u = ((q - p)[0] * s[1] - (q - p)[1] * s[0]) / den
            v = ((q - p)[0] * r[1] - (q - p)[1] * r[0]) / den
            if 0 <= u < 1 and 0 <= v <= 1:
                x = p + u * r
                out.append((TAU[i] + u * (TAU[i + 1] - TAU[i]), x[1], x[0]))
    return out


def main():
    root = r219.ROOT
    sys.path.insert(0, os.path.join(root, "models", "Q22_WK08-A"))
    eng = r219.load_mod(os.path.join(root, "models", "Q22_WK08-A", "02_rot_engine.py"), "02_rot_engine")
    mc = r219.load_mod(os.path.join(root, "models", "Q22_WK08-A", "04_mc_collisions.py"), "04_mc_collisions")
    mc.TAUS = TAU
    mc.SWITCH.clear()
    traces = {k: fz_trace(f, i) for k, (f, i, _) in FZ.items()}
    rows = []
    fig, ax = plt.subplots(figsize=(11, 6.5))
    for k, F in traces.items():
        ax.plot(F[:, 0], F[:, 1], "k-", lw=2.5)
        ax.annotate(k, (F[0, 0], F[0, 1]), xytext=(0, 6), textcoords="offset points", fontsize=8)
    for d, (rot, nom) in r219.MODELOS.items():
        tot = mc.make_tot(eng.RotModel(os.path.join(root, "models", d, "rot", rot)), None)
        res = {"modelo": nom}
        for k, (f, i, rel) in FZ.items():
            P = mc.track_points(tot, SEED, 16, rel)
            la, lo = mc.project(tot, rel, 0.0, P)       # present-day coordinates of the plate
            # at instant 0 frame 201 coincides with the present; 'project' with t=0 gives the
            # present-day position on plate rel
            c = crossings(lo, la, traces[k])
            res[k] = ";".join(f"{t:.2f}" for t, _, _ in c) if c else ""
            ax.plot(lo, la, lw=1, label=f"{nom} / {rel}" if rel == 902003 else None)
            for t, y, x in c:
                ax.plot(x, y, "o", ms=5, color="red")
        rows.append(res)
        print(f"{nom:14s} " + "  ".join(f"{k}: {res[k] or 'no crossing'}" for k in FZ))
    ax.set_xlim(-140, -70)
    ax.set_ylim(-42, -28)
    ax.set_aspect(1 / np.cos(np.radians(35)))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title("Foundation track (tau 10-45 Ma) in present-day coordinates of Nazca and Pacific vs FZ (GSFML)", fontsize=9)
    fig.savefig(os.path.join(HERE, "outputs", "iquique_prueba_03b", "foundation_fz_cruce.png"), dpi=120,
                bbox_inches="tight")
    import pandas as pd
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "outputs", "foundation_t_switch_fz.csv"), index=False)


if __name__ == "__main__":
    main()
