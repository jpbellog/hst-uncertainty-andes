# -*- coding: utf-8 -*-
"""Nominal collisions (step 04, unperturbed) with two .rot files of the same model:
latitude difference per ridge, f_PAC and age.  Usage: py nominal_sensitivity.py <folder> <rot_a> <rot_b> <csv>"""
import os, sys, importlib
import numpy as np, pandas as pd
d, ra, rb, out = sys.argv[1:5]
os.chdir(d); sys.path.insert(0, d)
mc = importlib.import_module("04_mc_collisions")
mdf, nn = mc.C.load_full_margin(mc.MARGIN_CSV)
ms = {t: mc.C.margin_segments(mdf, nn, t) for t in mc.AGES_T}
res = {}
for tag, fn in (("a", ra), ("b", rb)):
    M = mc.eng.RotModel(fn); tot = mc.make_tot(M, None)
    for pf, rel in mc.MODELS.items():
        for ridge, (seed, pid) in mc.RIDGES.items():
            P = mc.track_points(tot, seed, pid, rel)
            for t in mc.AGES_T:
                lat, lon = mc.project(tot, rel, t, P)
                r = mc.C.collide_strict(ms[t], mc.segs_for_age(lat, lon, t))
                if r is not None:
                    res.setdefault((ridge, pf, t), {})[tag] = r
rows = [dict(ridge=k[0], pf=k[1], age=k[2], lat_a=v.get("a", (np.nan,))[0], lat_b=v.get("b", (np.nan,))[0])
        for k, v in sorted(res.items())]
D = pd.DataFrame(rows); D["dif_km"] = (D.lat_b - D.lat_a) * 111.195
D.to_csv(out, index=False, float_format="%.4f")
g = D[D.pf == "0.44"].assign(tramo=lambda x: pd.cut(x.age, [-1, 27, 45, 66]))
print(g.groupby(["ridge", "tramo"], observed=True).dif_km.agg(lambda s: s.abs().max()).round(1).to_string())
