# -*- coding: utf-8 -*-
"""Rotated mean of the dated samples compared with 03b. Each dated sample
is rotated back at its age to the hotspot frame (h_i = B(t_i)^T x_i, same B as
03b_hotspot_fit.py) and averaged on the sphere, unweighted and with the fit weights
1/(omega^2 (sigma_t^2 + w_t^2) + w_x^2). The mean is compared with h0 of the
"fijo" (fixed) model of 03b (same set of samples, hotspot_fit.json) and with the dispersion
of the rotated points. Usage: py rotated_mean.py <folder> <rot> (from any directory)."""
import os, sys, json, importlib
import numpy as np, pandas as pd
d, rot = sys.argv[1:3]
os.environ["PAPER_ROT_FILE"] = rot
os.environ.setdefault("PAPER_HS_PRIOR", "none")
os.chdir(d); sys.path.insert(0, d)
fb = importlib.import_module("03b_hotspot_fit")
eng = fb.eng
M = eng.RotModel(rot)
fit = json.load(open("hotspot_fit.json", encoding="utf-8"))
fit = fit.get("tracks", fit)
for ridge, pid, fn, col in fb.TRACKS[:2]:
    r = fit[ridge]
    df = pd.read_csv(os.path.join("data", fn), comment="#")
    excl = {e["sample"] for e in r["excluded"]}
    df = df[(df.use_fit == 1) & ~df["sample"].isin(excl)].reset_index(drop=True)
    tr = fb.Track(M, ridge, pid, df, col)
    h = np.einsum("nji,nj->ni", tr.B0, tr.obs)          # B^T x
    lat0, lon0 = r["models"]["fijo"]["theta"][:2]
    x0 = np.array(eng.latlon_axis(lat0, lon0))
    th = np.array(r["models"]["fijo"]["theta"], float)
    _, _, om = tr.residuals(th, np.ones(len(df), bool))
    wt, wx = r["models"]["fijo"]["w_t"], r["models"]["fijo"]["w_x"]
    w = 1.0 / (om ** 2 * (tr.sig_t ** 2 + wt ** 2) + wx ** 2)
    p0 = np.array(eng.latlon_axis(*fb.SEEDS[pid]))       # present position (0 Ma pseudo-observation)
    hp, wp = np.vstack([h, p0]), np.r_[w, 1.0 / wx ** 2]
    for nom, ww, hh in (("sin pesos", np.ones(len(h)), h), ("con pesos", w, h),
                        ("+ 0 Ma", wp, hp)):
        h_ = hh
        m = (ww[:, None] * h_).sum(0); m /= np.linalg.norm(m)
        la, lo = eng.xyz_to_latlon(m)
        print(f"{ridge:6s} {nom:9s} n={len(h_):2d} media ({la:8.3f}, {lo:8.3f})  "
              f"a {eng.gc_km(la, lo, lat0, lon0):6.1f} km de h0 fijo ({lat0:.3f}, {lon0:.3f})")
    dist = [eng.gc_km(*eng.xyz_to_latlon(v), lat0, lon0) for v in h]
    print(f"{'':6s} dispersion de los puntos rotados respecto a h0: mediana {np.median(dist):.0f} km, max {np.max(dist):.0f} km")
