# -*- coding: utf-8 -*-
"""STEP 11 -- NOMINAL traces of the f_PAC=0.5 and 0.36 models per epoch.

For the reconstruction figure with uncertainty (figures/make_F3_reconstruction.py):
reconstructed (unperturbed) tracks of the end-member models 902001 (0.5) and
902002 (0.36) at the mid-age of each of the 8 panels, with the point age
(tau) per vertex so that the tau>83 Ma segment can be marked (no covariances:
Rowan & Rowley ends at chron 34y = 83 Ma).

The reference model 0.44 (plate 902003) is ALSO written. This calculation
reproduces the tracks of an independent GPlates reconstruction (max 13 cm over
2090 vertices).

Salida: recon_model_tracks.csv (pf, t_recon, ridge, tau, lat, lon)
Uso:  py 11_recon_model_tracks.py
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import importlib
eng = importlib.import_module("02_rot_engine")
mc = importlib.import_module("04_mc_collisions")

MODELS = {"0.5": 902001, "0.36": 902002, "0.44": 902003}
INTERVALS = [(72, 56), (56, 41.2), (41.2, 33.9), (33.9, 27.82),
             (27.82, 23.03), (23.03, 13.82), (13.82, 5.33), (5.33, 0)]
T_RECON = [round((a + b) / 2.0, 3) for a, b in INTERVALS]


def main():
    M = eng.RotModel()
    mc.load_circuit_cov(M)   # variant a changes the 802@71.1 pole
    tot0 = mc.make_tot(M, None)
    rows = []
    for pf, rel in MODELS.items():
        for ridge, (seed, pid) in mc.RIDGES.items():
            P = mc.track_points(tot0, seed, pid, rel)
            for t in T_RECON:
                lat, lon = mc.project(tot0, rel, t, P)
                for tau, la, lo in zip(mc.TAUS, lat, lon):
                    if tau < t:      # the point does not exist yet at that age
                        continue
                    if not np.isfinite(la):   # Foundation, tau < t_switch
                        continue
                    rows.append(dict(pf=pf, t_recon=t, ridge=ridge,
                                     tau=tau, lat=la, lon=lo))
    out = pd.DataFrame(rows)
    fn = os.path.join(HERE, f"recon_model_tracks{mc.HS_TAG}.csv")
    out.to_csv(fn, index=False)
    print("wrote:", fn, f"({len(out)} rows)")


if __name__ == "__main__":
    main()
