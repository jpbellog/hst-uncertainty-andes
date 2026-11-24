# -*- coding: utf-8 -*-
"""Sensitivity: 03b fit without the present reference position as an age-0
observation (PAPER_HS_P0=0), in the twelve models. The seeds have no source in active
volcanism (they are the present positions of the model of Bello-Gonzalez et al. 2018, 2026).

Inputs: rotation_models/fit_hotspots.py and fit_hotspots.py --noP0 (hotspot_fit_rep*.json in each folder and
.rot in rotation_models/outputs/rot and rot_noP0).
Outputs (outputs/noP0/):
  ajustes_noP0.csv   change per model and chain: chosen model, n, drift, w_t, w_x, displacement
                      of the hotspot at 0 Ma (km);
  nominal_noP0.csv   nominal collision with and without the seed (nominal_sensitivity.py), only at ages valid
                      under the frame_coverage of the frame, with the 68% half-width of the g13 band at that age.
Usage:  py sens_noP0.py
"""
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
REP = os.path.join(ROOT, "rotation_models")
OUT = os.path.join(HERE, "outputs", "noP0")
sys.path.insert(0, os.path.join(ROOT, "figures"))
sys.path.insert(0, REP)
import paper_models as MP  # noqa: E402
import fit_hotspots as R  # noqa: E402

KM = 111.195


def gc_km(a, b):
    la1, lo1, la2, lo2 = map(np.radians, (a[0], a[1], b[0], b[1]))
    c = np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2)
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))) * KM)


def ajustes():
    filas = []
    for c in R.MODELOS:
        a = json.load(open(os.path.join(ROOT, "models", c, "hotspot_fit_rep.json"), encoding="utf-8"))["tracks"]
        b = json.load(open(os.path.join(ROOT, "models", c, "hotspot_fit_gauss_conf_age_E_noP0_rep_noP0.json"), encoding="utf-8"))["tracks"]
        for r in a:
            x, y = a[r], b[r]
            filas.append(dict(modelo=c, ridge=r, modelo_con=x["model"], modelo_sin=y["model"],
                              n_con=x["n_used"], n_sin=y["n_used"],
                              deriva_con=np.hypot(*x["theta"][2:4]) if len(x["theta"]) > 2 else 0.0,
                              deriva_sin=np.hypot(*y["theta"][2:4]) if len(y["theta"]) > 2 else 0.0,
                              w_t_con=x["w_t"], w_t_sin=y["w_t"], w_x_con=x["w_x"], w_x_sin=y["w_x"],
                              h0_km=gc_km(x["theta"][:2], y["theta"][:2]),
                              excl_sin=";".join(e["sample"] for e in y["excluded"])))
    return pd.DataFrame(filas)


def nominal():
    out = []
    for c in MP.CARPETAS:
        a = R.salida_rot(c)
        b = a.replace(os.sep + "rot" + os.sep, os.sep + "rot_noP0" + os.sep).replace("/rot/", "/rot_noP0/")
        fn = os.path.join(OUT, f"nominal_{c}.csv")
        # PAPER_HS_FIT: step 04 adds Foundation (and its t_switch) from the json of the fit
        env = dict(os.environ, PAPER_HS_FIT=os.path.join(ROOT, "models", c, "hotspot_fit_rep.json"))
        subprocess.run([sys.executable, os.path.join(HERE, "nominal_sensitivity.py"), os.path.join(ROOT, "models", c), a, b, fn],
                       check=True, stdout=subprocess.DEVNULL, env=env)
        d = pd.read_csv(fn, dtype={"pf": str}).assign(modelo=c)
        d = MP.recortar(d.assign(lat_p50=d.dif_km), c).dropna(subset=["lat_p50"]).drop(columns="lat_p50")
        out.append(d)
    D = pd.concat(out, ignore_index=True)
    B = pd.concat([MP.recortar(pd.read_csv(os.path.join(ROOT, "results",
                   f"mc_collision_bands_total_present{MP.sufijo_bandas(c)}_hsb_cc_g13.csv"), dtype={"pf": str}), c)
                   .assign(modelo=c) for c in MP.CARPETAS])
    B["s68_km"] = (B.lat_p84 - B.lat_p16) / 2 * KM
    B["pf"] = B.pf.astype(float).map(lambda x: f"{x:g}")
    D["pf"] = D.pf.astype(float).map(lambda x: f"{x:g}")
    return D.merge(B[["modelo", "ridge", "pf", "age", "s68_km"]], on=["modelo", "ridge", "pf", "age"], how="left")


def main():
    os.makedirs(OUT, exist_ok=True)
    A = ajustes()
    A.to_csv(os.path.join(OUT, "ajustes_noP0.csv"), index=False, float_format="%.3f")
    pd.set_option("display.width", 200)
    print(A.round(1).to_string(index=False))
    D = nominal()
    D.to_csv(os.path.join(OUT, "nominal_noP0.csv"), index=False, float_format="%.4f")
    D["ad"] = D.dif_km.abs()
    D["exc"] = D.ad > D.s68_km
    g = D[D.pf == "0.44"].groupby("ridge")
    print("\nnine models, f_PAC 0.44, valid ages:")
    print(pd.DataFrame({"mediana_km": g.ad.median(), "max_km": g.ad.max(),
                        "edad_max": g.apply(lambda x: x.age[x.ad.idxmax()], include_groups=False),
                        "n": g.size(), "n_>_semiancho68": g.exc.sum()}).round(1).to_string())
    print(f"\nall f_PAC: {D.exc.sum()} of {len(D)} model-ages with a change larger than the 68% half-width")


if __name__ == "__main__":
    main()
