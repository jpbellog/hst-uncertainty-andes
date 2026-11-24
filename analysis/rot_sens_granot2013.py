# -*- coding: utf-8 -*-
"""Sensitivity: .rot with the West Antarctica-East Antarctica nodes
(804 rel 802) of Granot et al. (2013, GRL, Table 1)
between 27.859 (last REDBACK node of Granot and Dyment 2018) and 43.8 Ma.

Table 1 (checked against the image of page 282): finite rotations
"East (Fixed) and West Antarctica", cov = (1/k)(a b c; b d e; c e f) 1e-7 rad^2.
The printed angles rotate in the opposite sense to the neighbouring rows of the .rot; the
inverted rotation matches the .rot interpolated to within 13-33 km at Adare
(-71.5, 172 E) and grows monotonically from 27.859 Ma. It is used inverted:
A_rot = A_pub^T, C_rot = A_pub^T C A_pub (as for Tikku and Cande 2000 in circuit_covariance).
The sign was checked with Table S1 of the supplement (2012gl054181ts01.txt,
stored next to the PDF): A.-E.A. 18o (-15.04, -148.29, 23.72) has the convention
of the .rot (801 rel 802 C17o: -14.65, -146.53, 22.88, moving plate first). Closing
the circuit, WANT rel fixed EANT = R(A.-E.A.) R(A.-W.A.)^T gives the Table 1 pole
(within <= 0.4 deg) with a positive angle (12o +1.33, 13o +3.17, 16y +1.39, 18o +4.48);
the inverse product gives the printed negative angle. The printed Table 1 is EANT rel
WANT despite the title "East (Fixed)"; the inversion is correct (g13_sign.py).
16y (k 0.47) is the solution the authors describe as anomalous; it is included.
Option b: node 804 43.8 (Cande et al. 2000 extrapolated, 1.7 deg)
lies 4.33 deg from 18o over 3.7 Ma, and node 52.2 repeats that rotation.
Node 43.8 is replaced by a node at 45.03 (route change to PAC-LHR) equal to 18o, with
its cov and the same draw (SAME_DRAW in 04_mc_collisions.py): this assumes
no motion between 40.13 and 45.03 Ma. After 45.03 the link is not
in the chain.

Base model: with --base the row comments carry the base-model label and the .rot can be
rewritten in place (input rot = output rot). This was applied to the 14 _hsb .rot files of the seven
models (conf and pred), with a backup in <folder>/rot/.

Usage: py rot_sens_granot2013.py <base rot> <output rot> [--base]
"""
import sys
import numpy as np
sys.path.insert(0, (__import__("os").environ["PAPER_ROOT"] + "/models/Q22_WK08-A"))
import importlib
eng = importlib.import_module("02_rot_engine")
sys.path.insert(0, (__import__("os").environ["PAPER_ROOT"] + "/analysis"))
import circuit_covariance as cc

G13 = [  # age, anomaly, lat, lon, ang, kappa, a..f
    (30.94, "12o", -84.78, 190.99, -1.33, 0.81, (19.63, -5.33, 58.44, 2.74, -17.23, 175.88)),
    (33.55, "13o", -81.09, 200.37, -3.17, 0.69, (116.83, 1.27, 301.77, 3.65, -1.43, 786.32)),
    (35.69, "16y", -70.72, 338.78, -1.38, 0.47, (31.15, -7.53, 95.68, 3.20, -23.97, 292.13)),
    (40.13, "18o", -85.87, 220.49, -4.48, 1.02, (103.64, -28.37, 317.60, 10.08, -88.72, 976.84)),
]


def filas(etiqueta="sensibilidad"):
    out = []
    for age, an, la, lo, ang, k, v in G13:
        A = cc.pole_matrix(eng, la, lo, ang)
        C = A.T @ (cc.mat(cc.STD, v) * 1e-7 / k) @ A
        q = eng.matrix_to_quat(A.T)
        if q[0] < 0:
            q = -q
        s = np.linalg.norm(q[1:])
        ax = q[1:] / s
        lat, lon = np.degrees(np.arcsin(ax[2])), np.degrees(np.arctan2(ax[1], ax[0]))
        angd = np.degrees(2 * np.arctan2(s, q[0]))
        c = (C[0, 0], C[0, 1], C[0, 2], C[1, 1], C[1, 2], C[2, 2])
        out.append(f"804 {age} {lat:.4f} {lon:.4f} {angd:.4f} 802 !WANT-EANT {an} @REF Granot_++_2013_GRL "
                   f"Tabla 1 invertida, signo verificado con Tabla S1 ({etiqueta}, rot_sens_granot2013.py) "
                   + " ".join(f"{x:.4E}" for x in c) + "\n")
    # option b: 45.03 = 18o, constant until the route change
    out.append(out[-1].replace("804 40.13 ", "804 45.03 ", 1).replace(
        "!WANT-EANT 18o", "!WANT-EANT 18o constante hasta el cambio de ruta (supuesto nuestro, "
        "opcion b; comparte el sorteo de 40.13)", 1))
    return out


def main(base, salida, etiqueta="sensibilidad"):
    with open(base, encoding="utf-8", errors="surrogateescape") as f:
        L = f.readlines()
    assert not any(l.split()[:2] == ["804", "40.13"] for l in L), "already has the G13 nodes"
    i = next(k for k, l in enumerate(L) if l.split()[:2] == ["804", "43.8"])
    L[i:i + 1] = filas(etiqueta)   # replaces node 43.8 (option b)
    with open(salida, "w", encoding="utf-8", errors="surrogateescape", newline="") as f:
        f.writelines(L)
    M = eng.RotModel(salida)
    for r in M.table[804]:
        if 27 < r[0] < 45:
            print(r[0], r[1:4], "cov" if r[5] is not None else "no cov")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2],
         "modelo base" if "--base" in sys.argv[3:] else "sensibilidad")
