# -*- coding: utf-8 -*-
"""STEP 10 -- Uncertainty ellipses of the RECONSTRUCTED tracks per epoch.

For the reconstruction figure with uncertainty (figures/make_F3_reconstruction.py,
8 panels by geological epoch): 68/95% confidence ellipses of the track position
RECONSTRUCTED at the mid-age of each panel (t_recon), at the track point age
nodes (tau) that are multiples of 10 Ma and ALREADY EXISTED at that age
(tau >= t_recon), plus the young end tau = ceil(t_recon).

Same sampling as step 08 (total variant: circuit + WK08-A + geochronology,
f_PAC = 0.44, with the identity-node fix), but projecting the track to
t_recon in the SAM frame (project from step 04) instead of to t=0.

Output: recon_ellipses_0.44.csv
        (t_recon, ridge, tau, n, lat, lon, a68_km, b68_km, a95_km, b95_km,
         az_major_deg)
Usage:  py 10_recon_ellipses.py [N]   (default 400)
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
CHI2_68, CHI2_95 = 2.279, 5.991

# mid-ages of the 8 panels of the reconstruction figure
INTERVALS = [(72, 56), (56, 41.2), (41.2, 33.9), (33.9, 27.82),
             (27.82, 23.03), (23.03, 13.82), (13.82, 5.33), (5.33, 0)]
T_RECON = [round((a + b) / 2.0, 3) for a, b in INTERVALS]


def taus_for(t):
    """Tau nodes that get an ellipse in a panel of age t: multiples of 10 >= t
    (up to 80) + the young end of the existing track (ceil of t)."""
    base = [x for x in range(10, 81, 10) if x >= t]
    young = int(np.ceil(t))
    if young not in base:
        base = [young] + base
    return base


def main(N=400):
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
    rng = np.random.default_rng(20260717)

    plan = {t: taus_for(t) for t in T_RECON}
    idx = {t: [int(np.argmin(np.abs(mc.TAUS - x))) for x in plan[t]]
           for t in T_RECON}

    acc = {(t, r, x): [] for t in T_RECON for r in mc.RIDGES
           for x in plan[t]}
    t0 = time.time()
    for it in range(N):
        pert = mc.sample_perturb(M, rng, "total", node_sig, PF, hs)
        tot = mc.make_tot(M, pert)
        for ridge, (seed, pid) in mc.RIDGES.items():
            P = mc.track_points(tot, seed, pid, REL)
            for t in T_RECON:
                lat, lon = mc.project(tot, REL, t, P)
                for x, j in zip(plan[t], idx[t]):
                    acc[(t, ridge, x)].append((lat[j], lon[j]))
        if (it + 1) % 100 == 0:
            print(f"  {it+1}/{N} ({time.time()-t0:.0f} s)", flush=True)

    rows = []
    for (t, ridge, x), pts in sorted(acc.items()):
        A = np.array(pts)
        # Foundation: no track over Nazca/Farallon for tau < t_switch (NaN)
        if not np.isfinite(A).all():
            continue
        mla, mlo = A[:, 0].mean(), A[:, 1].mean()
        dy = (A[:, 0] - mla) * KM
        dx = (A[:, 1] - mlo) * KM * np.cos(np.radians(mla))
        C = np.cov(np.vstack([dx, dy]))
        w, V = np.linalg.eigh(C)
        a_sig, b_sig = np.sqrt(w[1]), np.sqrt(w[0])
        az = (np.degrees(np.arctan2(V[0, 1], V[1, 1]))) % 180
        rows.append(dict(t_recon=t, ridge=ridge, tau=x, n=len(A),
                         lat=mla, lon=mlo,
                         a68_km=a_sig * np.sqrt(CHI2_68),
                         b68_km=b_sig * np.sqrt(CHI2_68),
                         a95_km=a_sig * np.sqrt(CHI2_95),
                         b95_km=b_sig * np.sqrt(CHI2_95),
                         az_major_deg=az))
    out = pd.DataFrame(rows)
    fn = os.path.join(HERE, f"recon_ellipses_{PF}{mc.HS_TAG}.csv")
    out.to_csv(fn, index=False)
    print("wrote:", fn)
    for t in T_RECON:
        s = out[out.t_recon == t]
        print(f"  t={t:6.2f} Ma: {len(s)} ellipses, a95 "
              f"{s.a95_km.min():.0f}-{s.a95_km.max():.0f} km")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 400)
