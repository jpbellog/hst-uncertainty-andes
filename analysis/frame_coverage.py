# -*- coding: utf-8 -*-
"""Scope of each model according to its oldest frame covariance: each version is limited
to the age of its last frame covariance, without filling in with WK08-A.

A collision at age t involves the track segment formed at age
tau_c > t (the one at the trench at t). The row (ridge, pf, t) is valid if
tau_c <= T_cov, with T_cov the age of the oldest node of plate 002 that carries
a covariance. tau_c is computed on the nominal track (_hsb .rot of each model).
It also reports the oldest covariance of the circuit, by plate.

Usage:  py frame_coverage.py            (the twelve models)
Output: outputs/frame_coverage/frame_coverage_{folder}.csv and resumen_frame_coverage.csv
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd

P = (__import__("os").environ["PAPER_ROOT"])
A = os.path.dirname(os.path.abspath(__file__))
MODELOS = {
    "Q22_WK08-A": ("Global_Model_R_BG_2025_hsb.rot", "none", "wk08"),
    "WD25_WK08-A": ("Global_Model_R_BG_2025_WD25_hsb.rot", "none", "wk08"),
    "Q22_GMHRF-A": ("Global_Model_R_BG_2025_GMHRF_gts12_hsb.rot", "none", "rot"),
    "Q22_GMHRF-B": ("Global_Model_R_BG_2025_GMHRF_gts12_hsbgmhrf.rot", "gmhrf", "rot"),
    "WD25_GMHRF-B": ("Global_Model_R_BG_2025_WD25_GMHRF_gts12_hsbgmhrf.rot", "gmhrf", "rot"),
    "Q22_GWG20": ("Global_Model_R_BG_2025_GWG20_hsb.rot", "none", "rot"),
    "WD25_GWG20": ("Global_Model_R_BG_2025_WD25_GWG20_hsb.rot", "none", "rot"),
    "SG12_WK08-A": ("Global_Model_R_BG_2025_SG12_hsb.rot", "none", "wk08"),
    "WD25_GMHRF-A": ("Global_Model_R_BG_2025_WD25_GMHRF_gts12_hsb.rot", "none", "rot"),
    "SG12_GMHRF-A": ("Global_Model_R_BG_2025_SG12_GMHRF_gts12_hsb.rot", "none", "rot"),
    "SG12_GMHRF-B": ("Global_Model_R_BG_2025_SG12_GMHRF_gts12_hsbgmhrf.rot", "gmhrf", "rot"),
    "SG12_GWG20": ("Global_Model_R_BG_2025_SG12_GWG20_hsb.rot", "none", "rot"),
}


def units(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def tau_colision(segs, lat, lon):
    """Formation age of the point (lat, lon) on the track (segs from
    segs_for_age): segment with the smallest excess d1 + d2 - length."""
    p = units(lat, lon)
    a, b = units(segs[:, 1], segs[:, 2]), units(segs[:, 4], segs[:, 5])
    ang = lambda u, v: np.arccos(np.clip(np.sum(u * v, -1), -1, 1))
    d1, d2, L = ang(a, p), ang(b, p), ang(a, b)
    j = int(np.argmin(d1 + d2 - L))
    f = d1[j] / L[j] if L[j] > 0 else 0.0
    return segs[j, 0] + f * (segs[j, 3] - segs[j, 0])


def un_modelo(carpeta):
    rot, prior, cov = MODELOS[carpeta]
    D = os.path.join(P, "models", carpeta)
    os.environ["PAPER_ROT_FILE"] = os.path.join(D, "rot", rot)
    os.environ["PAPER_HS_FIT"] = os.path.join(D, "hotspot_fit.json")
    os.environ["PAPER_HS_PRIOR"] = prior
    os.environ["PAPER_FRAME_COV"] = cov
    os.chdir(D)
    sys.path.insert(0, D)
    for m in ("04_mc_collisions", "03b_hotspot_fit", "02_rot_engine"):
        sys.modules.pop(m, None)
    mc = importlib.import_module("04_mc_collisions")
    sys.path.remove(D)
    M = mc.eng.RotModel()
    if cov == "wk08":
        mc.load_wk08(M)
    mc.clone_9020x_cov(M)
    t_cov = max(r[0] for r in M.table[2] if r[5] is not None)
    circ = {p: max([r[0] for r in M.table[p] if r[5] is not None], default=np.nan)
            for p in mc.CIRCUIT_PLATES if p in M.table}
    mdf, nn = mc.C.load_full_margin(mc.MARGIN_CSV)
    tot = mc.make_tot(M, None)
    rows = []
    for pf, rel in mc.MODELS.items():
        for ridge, (seed, pid) in mc.RIDGES.items():
            Pp = mc.track_points(tot, seed, pid, rel)
            for t in mc.AGES_T:
                lat, lon = mc.project(tot, rel, t, Pp)
                segs = mc.segs_for_age(lat, lon, t)
                r = mc.C.collide_strict(mc.C.margin_segments(mdf, nn, t), segs)
                if r is None:
                    continue
                tc = tau_colision(segs, r[0], r[1])
                rows.append(dict(carpeta=carpeta, ridge=ridge, pf=pf, age=t,
                                 lat_nominal=r[0], tau_c=tc, T_cov_marco=t_cov,
                                 valido=bool(tc <= t_cov)))
    return pd.DataFrame(rows), t_cov, circ


def main():
    out = os.path.join(A, "outputs", "frame_coverage")
    os.makedirs(out, exist_ok=True)
    res = []
    for carpeta in (sys.argv[1:] or MODELOS):
        df, t_cov, circ = un_modelo(carpeta)
        df.to_csv(os.path.join(out, f"frame_coverage_{carpeta}.csv"), index=False,
                  float_format="%.4f")
        for (ridge, pf), g in df.groupby(["ridge", "pf"]):
            inval = g[~g.valido]
            res.append(dict(carpeta=carpeta, ridge=ridge, pf=pf, T_cov_marco=t_cov,
                            ultima_colision=g.age.max(),
                            edad_max_valida=(g.age[g.valido].max() if g.valido.any() else np.nan),
                            primera_invalida=(inval.age.min() if len(inval) else np.nan),
                            tau_c_max=g.tau_c.max()))
        print(carpeta, f"T_cov frame {t_cov}", "circuit (oldest cov):",
              {k: v for k, v in circ.items()}, flush=True)
    R = pd.DataFrame(res)
    R.to_csv(os.path.join(out, "resumen_frame_coverage.csv"), index=False, float_format="%.2f")
    with pd.option_context("display.width", 200):
        print(R.to_string(index=False))


if __name__ == "__main__":
    main()
