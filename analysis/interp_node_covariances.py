# -*- coding: utf-8 -*-
"""The engine (02_rot_engine.finite_rotation) does not interpolate
covariances. It draws h ~ N(0, C) at each NODE, independently between nodes, applies
exp(M(h)) R_node and then interpolates with slerp. Between two nodes the effective
perturbation is ~ (1-f) h0 + f h1, with variance (1-f)^2 C0 + f^2 C1: at the midpoint
the variance is halved if C0 ~ C1 (sd x 0.71).

This script measures the effect on plate 002 (WK08-A frame, nodes every 3-8 Myr), the only
source in the chain with spaced nodes over 0-66 Ma: sd of the perturbation angle
of 002 (relative to its fixed plate) on a 0.25 Myr grid,
compared with the sd at the nodes and with the linear interpolation of the sd between nodes
(which corresponds to full correlation between neighbouring nodes).

Usage:  py interp_node_covariances.py   ->  outputs/interp_node_covariances_002.csv
"""
import importlib
import os
import sys
import numpy as np
import pandas as pd

ROOT = (__import__("os").environ["PAPER_ROOT"] + "/models/Q22_WK08-A")
sys.path.insert(0, ROOT)
os.chdir(ROOT)
mc = importlib.import_module("04_mc_collisions")
eng = mc.eng
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")


def rotvec(A):
    q = eng.matrix_to_quat(A)
    s = np.linalg.norm(q[1:])
    return q[1:] / s * 2 * np.arctan2(s, q[0]) if s > 1e-15 else np.zeros(3)


def main(n=int(os.environ.get("N_INTERP", "3000")), seed=1):
    M = eng.RotModel(os.path.join(ROOT, "rot", "Global_Model_R_BG_2025_hsb.rot"))
    mc.load_wk08(M)
    nodes = [(r[0], r[4], r[5]) for r in M.table[2] if r[5] is not None and r[0] <= 70]
    rng = np.random.default_rng(seed)
    L = {(a, fx): np.linalg.cholesky(C + 1e-18 * np.eye(3)) for a, fx, C in nodes}
    T = np.arange(0.0, 66.01, 0.25)
    nom = [eng.quat_to_matrix(M.finite_rotation(2, t)[0]) for t in T]
    dev = np.zeros((n, T.size))
    for k in range(n):
        pert = {(2, a, fx): L[(a, fx)] @ rng.standard_normal(3) for a, fx in L}
        for j, t in enumerate(T):
            A = eng.quat_to_matrix(M.finite_rotation(2, t, pert)[0])
            dev[k, j] = np.degrees(np.linalg.norm(rotvec(A @ nom[j].T)))
    rms = np.sqrt((dev ** 2).mean(0))                       # deg, rms of the angle
    na = np.array([a for a, _, _ in nodes])
    rn = np.interp(na, T, rms)
    lin = np.interp(T, na, rn)                               # sd with "full correlation"
    D = pd.DataFrame(dict(t=T, rms_mc_gr=rms, rms_lineal_gr=lin, cociente=rms / lin))
    os.makedirs(OUT, exist_ok=True)
    D.to_csv(os.path.join(OUT, f"interp_node_covariances_002{os.environ.get('PAPER_INTERP_COV', '')}.csv"), index=False, float_format="%.5f")
    print("002 nodes with cov:", na.round(2).tolist())
    print(D.iloc[::8].round(4).to_string(index=False))
    print("minimum ratio", D.cociente.min().round(3), "median", D.cociente.median().round(3))


if __name__ == "__main__":
    main()
