# -*- coding: utf-8 -*-
"""Rotation file with the Nazca (Farallon)-South America motion of Somoza and Ghidella (2012,
EPSL 331-332, Table 1, p. 155) in place of the NZ-PAC-WANT-EANT-AFR-SAM circuit, used as a
separate circuit.

Source (Somoza and Ghidella 2012):
  - Table 1 "Finite reconstructions of Nazca (Farallon) to South America":
    chrons of Cande and Kent (1995), age in Ma, pole (lat N, lon E), negative
    (clockwise) angle and cov = (a b c; b d e; c e f) x 1e-5 rad^2.
  - p. 154: uncertainties "according to the method described by Chang et al.
    (1990, see also Doubrovine et al., 2008b)". Chang et al. (1990, sec. 3.4):
    moving exponential parametrization, R = R_hat Phi(h) (right).
  - Circuit of the authors: NZ-EANT-AFR-SAM up to 16 Ma; FAR-PAC-WANT-EANT-AFR-
    SAM before (p. 154).

What is done:
  - CK95 ages -> GTS2020 with timescale_GTS2020.remap("ck95"). The base .rot
    must already be in GTS2020 (rot_gts2020.py).
  - The .rot carries SAM relative to NZ (201 rel 911, the inverse of Table 1):
    same pole, angle with changed sign. Right-hand source and inverse -> the cov
    enters the engine (left) untransformed: R^T = Phi(-h) R_hat^T.
  - Rows 201 from 0 to 71.2 Ma (GTS2020 71.57) are replaced; after that the
    SAM-AFR of the base .rot follows, with a crossover row at 71.56 Ma equal to the
    interpolated base rotation (no cov). 911 = 902 = 902003 before 21.767 Ma in the .rot,
    so 201 rel 911 covers Nazca and Farallon.
  - Plates 13/14/15 (fit 03b) do not change: the 03b chain uses the circuit only
    at 0 Ma (identity), so the fit of the base model holds as it is.
  - Check: pole of Table 1 against 911 rel 201 of the base .rot at the same
    age (confirms the sense and measures the difference between circuits).

Usage:  py build_rot_sg12.py <base rot> <output rot>
"""
import importlib.util
import os
import sys
import numpy as np

ROOT = (__import__("os").environ["PAPER_ROOT"])
sys.path.insert(0, os.path.join(ROOT, "models", "Q22_WK08-A"))
eng = importlib.import_module("02_rot_engine")
_spec = importlib.util.spec_from_file_location(
    "esc", os.path.join(os.path.dirname(os.path.abspath(__file__)), "timescale_GTS2020.py"))
esc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(esc)

# chron, CK95 age, lat, lon, angle, a..f (x1e-5 rad^2)
SG12 = [
    ("C3n.3n (o)", 4.9, 62.05, 263.45, -3.87, (0.327, 0.287, -0.138, 0.708, -0.226, 0.488)),
    ("C5n.2n (o)", 10.8, 63.93, 262.81, -9.67, (0.383, 0.344, -0.031, 0.953, -0.046, 0.677)),
    ("C5C.1n (y)", 16.0, 63.94, 266.47, -15.04, (0.982, 0.770, -0.356, 1.981, -0.840, 1.263)),
    ("C6n (o)", 20.2, 58.77, 267.57, -20.12, (2.667, -0.519, 0.613, 2.862, -0.128, 0.760)),
    ("C7n (y)", 24.8, 61.45, 266.58, -26.05, (3.074, -0.316, 0.535, 3.506, -0.542, 1.227)),
    ("C10n.1n (y)", 28.3, 66.54, 258.65, -28.98, (3.379, -0.439, 0.927, 3.784, -0.722, 1.652)),
    ("C13n (m)", 33.3, 68.25, 243.93, -32.24, (2.550, -0.492, 1.047, 2.982, 0.112, 0.743)),
    ("C18n.2n (y)", 39.6, 68.82, 226.57, -37.35, (2.689, -0.465, 1.068, 3.377, -0.000, 1.041)),
    ("C21n (m)", 47.1, 68.32, 214.01, -43.65, (2.967, -0.372, 1.085, 4.162, -0.093, 1.989)),
    ("C25n (m)", 56.1, 68.50, 183.87, -49.37, (3.380, 0.238, 0.087, 6.824, -0.240, 6.467)),
    ("C31n (y)", 67.7, 60.64, 171.14, -54.96, (8.856, 3.391, -7.203, 17.447, 7.357, 31.824)),
    ("C32n.1n (m)", 71.2, 58.18, 169.51, -55.75, (17.255, 7.683, -18.627, 24.747, -0.092, 51.398)),
]
REF = '@REF Somoza_&_Ghidella_2012_EPSL @DOI"10.1016/j.epsl.2012.03.003"'


def cov(v):
    a, b, c, d, e, f = v
    return np.array([[a, b, c], [b, d, e], [c, e, f]]) * 1e-5


def rotvec(A):
    q = eng.matrix_to_quat(A)
    s = np.linalg.norm(q[1:])
    return q[1:] / s * 2 * np.arctan2(s, q[0]) if s > 1e-15 else np.zeros(3)


def main(base, salida):
    f = esc.remap("ck95")
    M = eng.RotModel(base)
    with open(base, encoding="utf-8", errors="surrogateescape") as fh:
        L = fh.readlines()
    assert not any(l.split()[:1] == ["201"] and l.split()[5:6] == ["911"] for l in L), "already SG12"
    t_end = float(f(SG12[-1][1]))
    rows = ["201 0.0 90.0 0.0 0.0 911 !SAM-NZ Somoza y Ghidella 2012 (circuito SG12, build_rot_sg12.py)\n"]
    print("chron         CK95   GTS12   misfit vs base .rot (deg)     sd max (km)")
    for cron, a, la, lo, ang, v in SG12:
        t = float(f(a))
        C = cov(v)
        if np.linalg.eigvalsh(C).min() <= 0:
            raise ValueError(f"cov not positive definite: {cron}")
        A_pub = eng.quat_to_matrix(eng.quat_from_pole(la, lo, ang))          # NZ rel SAM
        A_base = eng.quat_to_matrix(M.rel_rotation(911, t, 201))
        mis = np.degrees(np.linalg.norm(rotvec(A_base @ A_pub.T)))
        sd = np.degrees(np.sqrt(np.diag(C)).max()) * 111.195
        print(f"{cron:12s} {a:6.2f} {t:7.3f}   {mis:8.3f}                     {sd:7.1f}")
        c = (C[0, 0], C[0, 1], C[0, 2], C[1, 1], C[1, 2], C[2, 2])
        rows.append(f"201 {t:.3f} {la:.2f} {lo:.2f} {-ang:.2f} 911 !SAM-NZ {REF} @CHRONID\"{cron}\" "
                    f"Tabla 1 invertida (NZ rel SAM publicada), edad CK95 {a} -> GTS2020; cov sin "
                    f"transformar (Chang 1990, derecha; fila inversa) "
                    + " ".join(f"{x:.4E}" for x in c) + "\n")
    # crossover to the SAM-AFR of the base .rot at the end of SG12
    q = M.rel_rotation(201, t_end, 701)
    A = eng.quat_to_matrix(q)
    h = rotvec(A)
    ang = np.degrees(np.linalg.norm(h))
    la, lo = eng.xyz_to_latlon(h / np.linalg.norm(h))
    rows.append(f"201 {t_end:.3f} {la:.4f} {lo:.4f} {ang:.4f} 701 !SAM-AFR cruce SG12: rotacion del "
                f".rot base interpolada a {t_end:.3f} Ma (build_rot_sg12.py)\n")
    out, done = [], False
    for ln in L:
        p = ln.split()
        if p[:1] == ["201"] and len(p) > 5 and p[5] in ("701", "714"):
            if float(p[1]) <= t_end + 1e-9:
                if not done:
                    out += rows
                    done = True
                continue
        out.append(ln)
    with open(salida, "w", encoding="utf-8", errors="surrogateescape", newline="") as fh:
        fh.writelines(out)
    M2 = eng.RotModel(salida)
    for t in (5.0, 20.0, 40.0, 70.0):
        A1 = eng.quat_to_matrix(M2.rel_rotation(911, t, 201))
        print(f"check {t:5.1f} Ma: NZ rel SAM new vs base {np.degrees(np.linalg.norm(rotvec(A1 @ eng.quat_to_matrix(M.rel_rotation(911, t, 201)).T))):.3f} deg")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
