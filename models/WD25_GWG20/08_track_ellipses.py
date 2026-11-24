# -*- coding: utf-8 -*-
"""STEP 08 -- Uncertainty ellipses of the PRESENT POSITION of the tracks.

For the map figure: uncertainty of the present position (t=0) of each track
point (tau = age of the chain point), sampling the SAME three sources as
step 04 (total variant: circuit + WK08-A + geochronology) with the reference
model f_PAC=0.44.

For each ridge and each tau in {0,10,...,60} Ma it accumulates N positions (lat,lon)
of the track point at t=0 and fits the covariance ellipse in the local plane
(km): 68% (chi2_2 = 2.279) and 95% (chi2_2 = 5.991) semi-axes and azimuth of the
major axis. It also saves the nominal (unperturbed) track to draw the trace.

Outputs: track_ellipses_0.44.csv, track_nominal_0.44.csv
Usage:  py 08_track_ellipses.py [N]   (default N=500)
"""
import os
import sys
import time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import importlib
eng = importlib.import_module("02_rot_engine")
mc = importlib.import_module("04_mc_collisions")

KM = 111.19492664455873
PF, REL = "0.44", 902003
TAUS_SEL = [0, 10, 20, 30, 40, 50, 60]
CHI2_68, CHI2_95 = 2.279, 5.991


def main(N=500):
    M = eng.RotModel()
    # GMHRF and GWG20 frames (PAPER_FRAME_COV=rot): plate 002 already carries in the .rot
    # the covariances of its source; they are not overwritten with those of WK08-A.
    if os.environ.get("PAPER_FRAME_COV", "wk08") == "wk08":
        mc.load_wk08(M)
    mc.clone_9020x_cov(M)
    mc.load_circuit_cov(M)   # PAPER_CIRCUIT_COV
    node_sig = None   # node sigmas of the variant without a hotspot fit; unused with PAPER_HS_FIT
    if os.path.exists(os.path.join(HERE, "node_sigmas.csv")):
        node_sig = pd.read_csv(os.path.join(HERE, "node_sigmas.csv"))
        node_sig["pf"] = node_sig["pf"].astype(str)
    # with PAPER_HS_FIT the geochronology enters through the 03b fit
    hs = mc.load_hotspot_fit(M, mc.HS_FIT) if mc.HS_FIT else None
    rng = np.random.default_rng(20260715)

    idx = [int(np.argmin(np.abs(mc.TAUS - t))) for t in TAUS_SEL]

    # nominal (unperturbed): full trace at t=0
    tot0 = mc.make_tot(M, None)
    nom_rows = []
    for ridge, (seed, pid) in mc.RIDGES.items():
        P = mc.track_points(tot0, seed, pid, REL)
        lat, lon = mc.project(tot0, REL, 0.0, P)
        for tau, la, lo in zip(mc.TAUS, lat, lon):
            nom_rows.append(dict(ridge=ridge, tau=tau, lat=la, lon=lo))
    pd.DataFrame(nom_rows).to_csv(
        os.path.join(HERE, f"track_nominal_{PF}{mc.HS_TAG}.csv"), index=False)

    acc = {(r, t): [] for r in mc.RIDGES for t in TAUS_SEL}
    t0 = time.time()
    for it in range(N):
        pert = mc.sample_perturb(M, rng, "total", node_sig, PF, hs)
        tot = mc.make_tot(M, pert)
        for ridge, (seed, pid) in mc.RIDGES.items():
            P = mc.track_points(tot, seed, pid, REL)
            lat, lon = mc.project(tot, REL, 0.0, P)
            for tsel, j in zip(TAUS_SEL, idx):
                acc[(ridge, tsel)].append((lat[j], lon[j]))
        if (it + 1) % 100 == 0:
            print(f"  {it+1}/{N} ({time.time()-t0:.0f} s)", flush=True)

    rows = []
    for (ridge, tau), pts in sorted(acc.items()):
        A = np.array(pts)
        # Foundation: no track over Nazca/Farallon for tau < t_switch (NaN)
        if not np.isfinite(A).all():
            print(f"{ridge:7s} tau={tau:2d}: no track (tau < t_switch), skipped")
            continue
        mla, mlo = A[:, 0].mean(), A[:, 1].mean()
        dy = (A[:, 0] - mla) * KM                                # north, km
        dx = (A[:, 1] - mlo) * KM * np.cos(np.radians(mla))      # east, km
        C = np.cov(np.vstack([dx, dy]))
        w, V = np.linalg.eigh(C)                                  # w asc
        a_sig, b_sig = np.sqrt(w[1]), np.sqrt(w[0])               # 1 sigma km
        az = (np.degrees(np.arctan2(V[0, 1], V[1, 1]))) % 180     # azimuth of the major axis
        rows.append(dict(ridge=ridge, tau=tau, n=len(A), lat=mla, lon=mlo,
                         a68_km=a_sig * np.sqrt(CHI2_68),
                         b68_km=b_sig * np.sqrt(CHI2_68),
                         a95_km=a_sig * np.sqrt(CHI2_95),
                         b95_km=b_sig * np.sqrt(CHI2_95),
                         az_major_deg=az))
        print(f"{ridge:7s} tau={tau:2d}: 95% ellipse {rows[-1]['a95_km']:.0f} x "
              f"{rows[-1]['b95_km']:.0f} km  az {az:.0f}")
    out = os.path.join(HERE, f"track_ellipses_{PF}{mc.HS_TAG}.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print("wrote:", out)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 500)
