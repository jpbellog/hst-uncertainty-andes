# -*- coding: utf-8 -*-
"""STEP 09 -- VC uncertainty from a Monte Carlo of the CIRCUIT (0-66 Ma).

Extends step 05 (which propagates analytically the covariances of the NZ-SAM
stage vectors of Quiero et al. 2022, available only for 0-28 Ma) to the WHOLE
0-66 Ma range, sampling the nodal covariances of the full circuit in the
.rot ('circuit' variant of step 04):
  - 0-30 Ma: 911 (Wilder via Quiero) + 901/804/802/201 (Croon-REDBACK,
    Granot & Dyment, DeMets & Merkouriev; segments interpolated from Quiero).
  - 28-83 Ma: 902 with the covariances of ROWAN & ROWLEY 2014 (chrons 10y-34y:
    28.278-83 Ma). NOTE: for >28 Ma the OTHER branches of the circuit carry no
    covariance in the .rot, so the sigma for >28 Ma is a LOWER BOUND
    dominated by the FAR-PAC term (R&R), which is by far the largest.

Method: for each realisation the circuit nodes are perturbed within their
covariances (mc.sample_perturb, with the identity-node fix) and the stage
Euler vector of 902003 (f_PAC=0.44) relative to 201 is recomputed over
[t, t-1] Myr windows; the velocity v = w x r is evaluated at the nodes of the
margin RESTORED to age t (Margin_SAM, as for the VC field) and sigma is
accumulated over N realisations (magnitude and orthogonal component, with the
margin azimuth from node j->j+10 as in step 05).

Validation: for 0-28 Ma the comparison with step 05 (analytical, from
Quiero) is printed; both must be of the same order (different time discretisation).

Output: vc_sigma_circuit_mc.csv (age x margin node).
Usage:  py 09_vc_uncertainty_circuit.py [N]   (default 400)
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

MARGIN_CSV = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "margin", "Margin_SAM_0.1.csv")  # copy of Margin_SAM (see rotation_models/README.md)
R_E_CM_YR = 6371e5 / 1e6      # cm/yr per (rad/Myr), as in step 05
D2R = np.pi / 180.0
PLATE, REL, PF = 902003, 201, "0.44"
AGES = np.arange(1.0, 67.0, 1.0)   # window [t, t-1]


def unit(lat, lon):
    la, lo = lat * D2R, lon * D2R
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo),
                     np.sin(la)])


def margin_at(mdf, n_node, t):
    row = mdf.iloc[(mdf["age"] - t).abs().idxmin()]
    lon = np.array([row.get(f"Lon_{k}") for k in range(1, n_node + 1)], float)
    lat = np.array([row.get(f"Lat_{k}") for k in range(1, n_node + 1)], float)
    k = np.isfinite(lat) & np.isfinite(lon) & (lat <= -5) & (lat >= -36)
    return lat[k], lon[k]


def stage_w(tot, t):
    """Euler vector (rad/Myr) of PLATE rel REL over [t, t-1]."""
    def rel_R(tau):
        R1 = eng.quat_to_matrix(tot(REL, tau))
        R2 = eng.quat_to_matrix(tot(PLATE, tau))
        return R1.T @ R2
    S = rel_R(t - 1.0) @ rel_R(t).T
    ang = np.arccos(np.clip((np.trace(S) - 1) / 2, -1, 1))
    if ang < 1e-12:
        return np.zeros(3)
    ax = np.array([S[2, 1] - S[1, 2], S[0, 2] - S[2, 0], S[1, 0] - S[0, 1]])
    ax = ax / np.linalg.norm(ax)
    return ax * ang   # rad/Myr (1 Myr window)


def vc_at_margin(w, mlat, mlon):
    """Magnitude and orthogonal component of v = w x r at the margin nodes."""
    n = len(mlat)
    vmag = np.empty(n); vort = np.empty(n)
    for j in range(n):
        r = unit(mlat[j], mlon[j])
        v3 = np.cross(w, r) * R_E_CM_YR
        la, lo = mlat[j] * D2R, mlon[j] * D2R
        e = np.array([-np.sin(lo), np.cos(lo), 0.0])
        nv = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo),
                       np.cos(la)])
        vE, vN = v3 @ e, v3 @ nv
        j2 = min(j + 10, n - 1)
        dE = (mlon[j2] - mlon[j]) * np.cos(la)
        dN = mlat[j2] - mlat[j]
        azm = np.arctan2(dE, dN)
        n_hat = np.array([np.cos(azm), -np.sin(azm)])   # normal to the margin
        vmag[j] = np.hypot(vE, vN)
        vort[j] = abs(np.array([vE, vN]) @ n_hat)
    return vmag, vort


def main(N=400):
    M = eng.RotModel()
    mc.load_wk08(M)
    mc.clone_9020x_cov(M)
    mc.load_circuit_cov(M)   # PAPER_CIRCUIT_COV, as in steps 08 and 10
    node_sig = None   # node sigmas of the variant without a hotspot fit; unused with PAPER_HS_FIT
    if os.path.exists(os.path.join(HERE, "node_sigmas.csv")):
        node_sig = pd.read_csv(os.path.join(HERE, "node_sigmas.csv"))
        node_sig["pf"] = node_sig["pf"].astype(str)
    rng = np.random.default_rng(20260716)

    mdf = pd.read_csv(MARGIN_CSV)
    mdf["age"] = pd.to_numeric(mdf["age"], errors="coerce")
    n_node = sum(1 for c in mdf.columns if c.startswith("Lon_"))
    margins = {t: margin_at(mdf, n_node, t) for t in AGES}

    # nominal
    tot0 = mc.make_tot(M, None)
    nom = {t: vc_at_margin(stage_w(tot0, t), *margins[t]) for t in AGES}

    # MC
    acc_mag = {t: [] for t in AGES}
    acc_ort = {t: [] for t in AGES}
    t0 = time.time()
    for it in range(N):
        pert = mc.sample_perturb(M, rng, "circuit", node_sig, PF)
        tot = mc.make_tot(M, pert)
        for t in AGES:
            vm, vo = vc_at_margin(stage_w(tot, t), *margins[t])
            acc_mag[t].append(vm)
            acc_ort[t].append(vo)
        if (it + 1) % 50 == 0:
            print(f"  {it+1}/{N} ({time.time()-t0:.0f} s)", flush=True)

    rows = []
    for t in AGES:
        mlat, mlon = margins[t]
        sm = np.std(np.array(acc_mag[t]), axis=0)
        so = np.std(np.array(acc_ort[t]), axis=0)
        vm0, vo0 = nom[t]
        for j in range(len(mlat)):
            rows.append(dict(age=t, lat=mlat[j], lon=mlon[j],
                             vc_mag=vm0[j], sigma_mag=sm[j],
                             vc_ort=vo0[j], sigma_ort=so[j]))
    out = pd.DataFrame(rows)
    fn = os.path.join(HERE, "vc_sigma_circuit_mc.csv")
    out.to_csv(fn, index=False)
    print("wrote:", fn)

    # validation against step 05 (analytical Quiero, 0-28 Ma); step 05 runs only in Q22_WK08-A
    fq = os.path.join(HERE, "vc_sigma_quiero.csv")
    common = []
    if os.path.exists(fq):
        q = pd.read_csv(fq)
        qq = q.groupby(q.age_mid.round(0)).sigma_mag.median()
        mm = out[out.age <= 28].groupby("age").sigma_mag.median()
        common = sorted(set(qq.index) & set(mm.index))
    if common:
        ratio = np.array([mm[a] / qq[a] for a in common if qq[a] > 0])
        print(f"validation vs step 05 (0-28 Ma): median sigma_MC/sigma_Quiero"
              f" = {np.median(ratio):.2f} (p10 {np.percentile(ratio,10):.2f},"
              f" p90 {np.percentile(ratio,90):.2f})")
    for a0, a1 in ((1, 7), (7, 26), (26, 28), (28, 45), (45, 66)):
        s = out[(out.age >= a0) & (out.age < a1)]
        print(f"  {a0}-{a1} Ma: sigma_mag median {s.sigma_mag.median():.2f} "
              f"p90 {s.sigma_mag.quantile(0.9):.2f} cm/yr")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 400)
