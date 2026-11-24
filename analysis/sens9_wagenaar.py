# -*- coding: utf-8 -*-
"""Nominal collision with the WVH25 frame (Wagenaar et
al. 2025 with the Quiero circuit; build_rot_wagenaar.py, nominal_collision.py -> coll_WVH25_g13.csv)
in the nine g13 models, compared with the Quiero + WK08-A nominal (pipeline) and with the 68 % envelope of the
nine models (f_PAC 0.44, valid ages of each model).
Usage:  py sens9_wagenaar.py
"""
import os
import sys
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "figures"))
import paper_models as MP  # noqa: E402

W = pd.read_csv(os.path.join(HERE, "outputs", "wagenaar_frame", "coll_WVH25_g13.csv"))
B = []
for d in MP.CARPETAS:
    b = pd.read_csv(os.path.join(ROOT, "models", d, "mc_collision_bands_total_hsb_cc_g13.csv"))
    b = MP.recortar(b[np.isclose(b.pf.astype(float), 0.44)], d)
    b["modelo"] = d
    B.append(b)
B = pd.concat(B)
env = B.groupby(["ridge", "age"]).agg(p16=("lat_p16", "min"), p84=("lat_p84", "max"),
                                       n=("lat_p16", "count")).reset_index()
wk = B[B.modelo == "Q22_WK08-A"][["ridge", "age", "lat_nominal"]].rename(columns={"lat_nominal": "lat_wk08"})
T = W.merge(wk, on=["ridge", "age"]).merge(env, on=["ridge", "age"])
T = T[T.lat.notna() & T.lat_wk08.notna()]
T["dif_deg"] = (T.lat - T.lat_wk08).abs()
T["dentro68"] = (T.lat >= T.p16) & (T.lat <= T.p84)
T.to_csv(os.path.join(HERE, "outputs", "sens9", "wagenaar_9models.csv"), index=False, float_format="%.4f")
print("WVH25 compared with Quiero + WK08-A (|dlat|, degrees) and inside the 68 % envelope of the nine models")
print(T.groupby("ridge").agg(edades=("age", "size"), mediana=("dif_deg", "median"), maximo=("dif_deg", "max"),
                             edad_max=("dif_deg", lambda s: T.loc[s.idxmax(), "age"]),
                             dentro68=("dentro68", "mean")).round(2).to_string())
