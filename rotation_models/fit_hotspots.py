# -*- coding: utf-8 -*-
"""Step 03b of the twelve models on the .rot files of build_rotation_models.py, with the environment
of the final run set here, and comparison with the fit of each folder.

Step 03b fits, for each chain (Easter, Juan Fernandez, San Felix, Foundation), the present-day
position of the hotspot in the frame (and its drift, if the BIC prefers it or it is imposed) to the
dated ages, with geological dispersion w_t and w_x and outlier rejection; it writes hotspot_fit*.json
and the rows of plates 13-16 in the .rot.

Common environment: Foundation ages _v3 (default of 03b), w as confidence, fixed or drifting model
chosen by the BIC in each chain (without PAPER_HS_FORCE: without the age-0 observation the BIC
prefers the Easter drift in the twelve models and the fit is identical to the imposed one).
GMHRF-B: PAPER_HS_PRIOR=gmhrf.
The outputs carry the tag _rep (they do not overwrite those of the folders); the .rot with plates
13-16 is written to outputs/rot/<name>_hsb[gmhrf].rot.

The base case does not include the reference present-day position as an age-0 observation
(default of 03b); --conP0 runs the sensitivity with it (PAPER_HS_P0=1, outputs
hotspot_fit_gauss_conf_age_E_conP0_rep_conP0.json and outputs/rot_conP0/).

Usage:  py fit_hotspots.py [folder ...] [--conP0]   (default all twelve; NP=4 processes)
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ROT = os.path.join(HERE, "outputs", "rot")
# --conP0: sensitivity with the present-day position as an age-0 observation (PAPER_HS_P0=1)
CON_P0 = "--conP0" in sys.argv
TAG = "_rep_conP0" if CON_P0 else "_rep"
ROT_OUT = os.path.join(HERE, "outputs", "rot_conP0" if CON_P0 else "rot")

# folder: (.rot of build_rotation_models.py, prior)
MODELOS = {
    "Q22_WK08-A": ("Global_Model_R_BG_2025.rot", "none"),
    "Q22_GMHRF-B": ("Global_Model_R_BG_2025_GMHRF_gts12.rot", "gmhrf"),
    "Q22_GWG20": ("Global_Model_R_BG_2025_GWG20.rot", "none"),
    "WD25_WK08-A": ("Global_Model_R_BG_2025_WD25.rot", "none"),
    "WD25_GMHRF-B": ("Global_Model_R_BG_2025_WD25_GMHRF_gts12.rot", "gmhrf"),
    "WD25_GWG20": ("Global_Model_R_BG_2025_WD25_GWG20.rot", "none"),
    "SG12_WK08-A": ("Global_Model_R_BG_2025_SG12.rot", "none"),
    "SG12_GMHRF-B": ("Global_Model_R_BG_2025_SG12_GMHRF_gts12.rot", "gmhrf"),
    "SG12_GWG20": ("Global_Model_R_BG_2025_SG12_GWG20.rot", "none"),
    "Q22_GMHRF-A": ("Global_Model_R_BG_2025_GMHRF_gts12.rot", "none"),
    "WD25_GMHRF-A": ("Global_Model_R_BG_2025_WD25_GMHRF_gts12.rot", "none"),
    "SG12_GMHRF-A": ("Global_Model_R_BG_2025_SG12_GMHRF_gts12.rot", "none"),
}


def salida_rot(carpeta):
    rot, prior = MODELOS[carpeta]
    return os.path.join(ROT_OUT, rot.replace(".rot", "_hsbgmhrf.rot" if prior == "gmhrf" else "_hsb.rot"))


def correr(carpeta):
    rot, prior = MODELOS[carpeta]
    env = dict(os.environ, PAPER_ROT_FILE=os.path.join(ROT, rot), PAPER_HS_ROT_OUT=salida_rot(carpeta),
               PAPER_HS_PRIOR=prior, PAPER_HS_OUT_TAG=TAG,
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    for k in ("PAPER_HS_W", "PAPER_HS_LAG", "PAPER_HS_EXCLUDE", "PAPER_HS_P0", "PAPER_HS_LIK", "PAPER_HS_FOUNDATION",
              "PAPER_HS_FORCE"):
        env.pop(k, None)
    if CON_P0:
        env["PAPER_HS_P0"] = "1"
    d = os.path.join(ROOT, "models", carpeta)
    with open(os.path.join(d, f"hotspot_fit{TAG}.log"), "w") as log:
        r = subprocess.run([sys.executable, "03b_hotspot_fit.py"], cwd=d, env=env, stdout=log, stderr=subprocess.STDOUT)
    return carpeta, r.returncode


def comparar(carpeta):
    """Maximum difference with hotspot_fit.json of the folder: theta, w_t, w_x, cov and excluded samples."""
    a = json.load(open(os.path.join(ROOT, "models", carpeta, f"hotspot_fit{TAG}.json"), encoding="utf-8"))["tracks"]
    b = json.load(open(os.path.join(ROOT, "models", carpeta, "hotspot_fit.json"), encoding="utf-8"))["tracks"]
    peor, msg = 0.0, []
    for r in b:
        x, y = a[r], b[r]
        if x["model"] != y["model"] or sorted(e["sample"] for e in x["excluded"]) != sorted(e["sample"] for e in y["excluded"]):
            msg.append(f"{r}: model or excluded samples differ")
        for k in ("theta", "w_t", "w_x", "cov"):
            d = np.max(np.abs(np.subtract(x[k], y[k])) / (np.abs(y[k]) + 1e-12))
            peor = max(peor, float(d))
    return peor, msg


def main():
    os.makedirs(ROT_OUT, exist_ok=True)
    pedidas = [a for a in sys.argv[1:] if not a.startswith("--")] or list(MODELOS)
    with ThreadPoolExecutor(int(os.environ.get("NP", "4"))) as ex:
        for c, rc in ex.map(correr, pedidas):
            print(f"{c}: exit {rc}", flush=True)
    for c in pedidas if not CON_P0 else []:
        d, msg = comparar(c)
        print(f"{c}: maximum relative diff. {d:.1e} " + ("; ".join(msg) if msg else "(same model and excluded samples)"))


if __name__ == "__main__":
    main()
