# -*- coding: utf-8 -*-
"""Published covariances for the circuit links that the .rot carries without
covariance (the frame_coverage rule is applied to the circuit; covariances are
added without changing any pole).

Each entry is matched to a .rot node by (plate, age, fixed plate) and only
if the .rot pole coincides with the published one (or with its inverse). Sources and
reading (checked against the original publications):

  201  SAM-AFR   Mueller et al. (1999) Tab. 4 (a..f) and Tab. 2 (kappa). The formula
                 of the paper is (a..f)*g without defining g. g was calibrated with the
                 95% volumes in the text (p. 41: C5 7838, C8 19705, C21 2773
                 km^3): g*kappa = 1.08-1.17e-7, i.e. C = (1/kappa)(a..f)1e-7
                 rad^2 (Royer and Chang convention). Rows 40.1, 47.9 and 67.7 of the
                 .rot are not Mueller poles: the cov of the neighbouring chron is used
                 (18, 21, 30) and they are flagged "vecino".
  802  ANT-SOM   Cande et al. (2010) Tab. 2, 3 and 5, (1/kappa)*1e-7 rad^2. The
                 printed columns a..f are NOT in the order (a b c; b d e; c e f):
                 with that order d<0 in all rows. The only order with all 42
                 matrices positive definite (up to axis permutation) is
                 xx, xy, yy, xz, yz, zz; its structure (corr xy ~ +0.9) matches
                 the independent cov of DeMets et al. (2021) for C13o.
  901  PAC-WANT  Croon et al. (2008) Tab. 1, (1/kappa)*10^-g. Row 33.55 (d=.72)
                 is not positive definite: it is interpolated between 33.06 and 34.66.
                 Row 34.94 printed "738 956":
                 read as .738 .956 (positive definite), flagged "decimales".
  709  SOM-NUB   Lemaux et al. (2002), caption of Fig. 4, 1e-6 rad^2 (node 10.95
                 and node 250, which repeats the same rotation).
  833  LHR-AUS   Gaina et al. (1998) Tab. 2a (kappa) and 2b, (1/kappa)*1e-6 rad^2,
                 Australia-fixed frame. 33o and 34y have no cov (assumed pole).
  801  AUS-EANT  Tikku and Cande (2000) Tab. 1-2, (1/kappa)*w. The table gives the
                 rotation inverse to that of the .rot (antipodal pole, same angle):
                 C_rot = A^T C A with A the published rotation. 32y is not positive
                 definite: excluded.
  804  WANT-EANT DERIVED. Node 43.8 of the .rot has the pole of solution 13o
                 of Cande et al. (2000) (EANT relative to WANT, -0.696 degrees,
                 kappa 0.63, Tab. 2 *10^-5) with the angle extrapolated to C20
                 (1.7 degrees, inverse). It is inverted and scaled by (1.7/0.696)^2.
  901  PAC-LHR   DERIVED. The .rot fixes LHR to the Pacific between 45.03 and 83.5 Ma with
                 the 45 Ma closure rotation; its cov is propagated along the route
                 PAC-WANT-EANT-AUS at 43.79 Ma (LHR-AUS is zero before 52 Ma).

Convention: that of the engine (R_real = exp(M(h)) R_nominal, h in the frame of the
fixed plate). Sources that declare another frame (e.g. DeMets, moving plate)
are pending audit.

Usage:  from circuit_covariance import aplicar;  aplicar(M)   (modifies M.table)
        py circuit_covariance.py                 (check table and csv)
"""
import json
import os
import sys
import numpy as np

R2D = 180.0 / np.pi

# (plate, .rot age, fixed, published pole (lat, lon, ang), kappa, scale,
#  order, (a..f), source, note)
STD = ("xx", "xy", "xz", "yy", "yz", "zz")
ALT = ("xx", "xy", "yy", "xz", "yz", "zz")
MU = "Mueller et al. 1999"
C10 = "Cande et al. 2010"
CR = "Croon et al. 2008"
GA = "Gaina et al. 1998"
TK = "Tikku y Cande 2000"
PUB = [
    # SAM-AFR (Mueller 1999, Tab. 2 poles/kappa, Tab. 4 cov)
    (201, 33.1, 701, (56.17, -33.64, 13.41), .67, 1e-7, STD, (9.525, -2.385, -5.374, 2.117, 1.676, 3.765), MU, "C13"),
    (201, 40.1, 701, (57.10, -33.00, 15.912), 1.44, 1e-7, STD, (31.698, -12.098, -18.390, 5.766, 7.124, 12.671), MU, "vecino C18"),
    (201, 47.9, 701, (56.95, -31.15, 19.107), .89, 1e-7, STD, (11.165, -3.955, -6.408, 2.546, 2.525, 4.334), MU, "vecino C21"),
    (201, 55.9, 701, (61.35, -32.21, 22.273), 1.01, 1e-7, STD, (25.115, -11.568, -13.323, 7.341, 6.681, 8.144), MU, "C25"),
    (201, 67.7, 701, (63.88, -33.61, 24.755), 1.74, 1e-7, STD, (15.017, -6.988, -8.168, 5.207, 4.322, 5.358), MU, "vecino C30"),
    (201, 83.5, 701, (61.88, -34.26, 33.512), .97, 1e-7, STD, (8.711, -4.332, -4.971, 4.014, 3.230, 3.649), MU, "C34"),
    # ANT-SOM (Cande 2010; set 1 = Tab. 2, set 2 = Tab. 3, set 4 = Tab. 5)
    (802, 33.55, 709, (12.69, -44.61, 5.67), .59, 1e-7, ALT, (5.71, 5.12, 5.63, -2.35, -2.53, 3.99), C10, "set1 13o"),
    (802, 40.13, 709, (13.80, -43.75, 7.05), .63, 1e-7, ALT, (6.60, 5.69, 5.56, -0.01, -0.79, 6.80), C10, "set1 18o"),
    (802, 42.54, 709, (11.84, -42.32, 7.53), .93, 1e-7, ALT, (2.36, 1.91, 2.31, -2.28, -2.63, 5.98), C10, "set4 20y"),
    (802, 43.79, 709, (11.84, -42.16, 7.87), .78, 1e-7, ALT, (1.64, 1.23, 1.79, -1.95, -2.42, 4.97), C10, "set4 20o"),
    (802, 46.26, 709, (11.29, -41.54, 8.49), .55, 1e-7, ALT, (1.89, 1.60, 2.24, -2.76, -3.25, 6.86), C10, "set4 21y"),
    (802, 47.91, 709, (9.82, -40.70, 8.83), 1.86, 1e-7, ALT, (2.26, 1.98, 2.67, -3.69, -4.00, 8.38), C10, "set4 21o"),
    (802, 49.71, 709, (9.19, -40.63, 9.21), 1.40, 1e-7, ALT, (2.73, 2.54, 3.70, -4.70, -5.40, 10.84), C10, "set4 22o"),
    (802, 51.74, 709, (9.31, -41.53, 9.61), .74, 1e-7, ALT, (2.78, 2.77, 4.22, -4.41, -5.38, 9.45), C10, "set2 23o"),
    (802, 53.35, 709, (10.16, -43.30, 9.96), .72, 1e-7, ALT, (4.74, 4.34, 5.65, -7.41, -7.62, 13.98), C10, "set2 24o"),
    (802, 55.9, 709, (9.86, -45.24, 10.49), .52, 1e-7, ALT, (11.49, 13.75, 19.32, -20.16, -26.01, 39.54), C10, "set2 25y"),
    (802, 57.55, 709, (10.64, -47.47, 10.78), 1.35, 1e-7, ALT, (12.67, 14.09, 17.47, -18.85, -22.53, 33.30), C10, "set2 26y"),
    (802, 60.92, 709, (7.10, -45.80, 11.08), 1.15, 1e-7, ALT, (16.23, 17.51, 21.93, -24.59, -28.43, 42.88), C10, "set2 27y"),
    (802, 62.5, 709, (4.75, -44.79, 11.39), .66, 1e-7, ALT, (13.59, 15.49, 20.54, -23.54, -28.38, 44.86), C10, "set2 28y"),
    (802, 64.75, 709, (4.79, -45.56, 11.70), .53, 1e-7, ALT, (9.16, 9.94, 12.70, -16.51, -19.13, 32.88), C10, "set2 29o"),
    # PAC-WANT (Croon 2008, Tab. 1; scale 10^-g)
    (901, 28.75, 804, (74.330, -69.543, 22.953), 2.92, 1e-8, STD, (2.33, 1.32, 2.48, 1.36, 1.77, 3.95), CR, "10o"),
    (901, 29.4, 804, (74.378, -68.714, 23.662), 2.20, 1e-8, STD, (1.81, 1.09, 2.14, 1.42, 1.91, 4.16), CR, "11y"),
    (901, 30.1, 804, (74.329, -68.430, 24.232), 1.12, 1e-8, STD, (1.94, .960, 2.34, 1.25, 1.67, 4.75), CR, "11o"),
    (901, 30.48, 804, (74.366, -67.750, 24.708), 2.67, 1e-8, STD, (1.97, 1.00, 2.10, 1.30, 1.74, 3.88), CR, "12y"),
    (901, 30.94, 804, (74.320, -67.593, 25.058), 1.15, 1e-8, STD, (1.92, 1.00, 1.99, 1.30, 1.74, 3.95), CR, "12o"),
    (901, 33.06, 804, (74.444, -64.739, 26.969), 1.83, 1e-8, STD, (3.27, 1.67, 3.13, 1.63, 2.17, 4.92), CR, "13y"),
    (901, 34.66, 804, (74.621, -61.886, 28.395), 1.11, 1e-8, STD, (4.21, 2.70, 4.19, 2.47, 3.27, 5.61), CR, "15y"),
    (901, 34.94, 804, (74.634, -61.485, 28.615), 2.77, 1e-7, STD, (1.22, .805, 1.21, .738, .956, 1.64), CR, "15o decimales"),
    (901, 35.69, 804, (74.702, -60.180, 29.234), 1.46, 1e-8, STD, (4.22, 2.73, 4.36, 2.54, 3.22, 6.03), CR, "16y"),
    (901, 36.34, 804, (74.750, -59.047, 29.808), 2.96, 1e-7, STD, (1.27, .852, 1.30, .801, 1.03, 1.80), CR, "16o"),
    (901, 36.62, 804, (74.763, -58.638, 30.049), 1.15, 1e-8, STD, (4.14, 2.73, 4.31, 2.55, 3.36, 5.90), CR, "17y"),
    (901, 37.47, 804, (74.816, -57.440, 30.678), 2.01, 1e-7, STD, (1.32, .885, 1.37, .825, 1.08, 1.89), CR, "17o"),
    (901, 38.43, 804, (74.855, -56.211, 31.411), 2.13, 1e-8, STD, (4.11, 3.04, 4.83, 3.05, 4.02, 7.09), CR, "18y"),
    (901, 40.13, 804, (74.868, -54.460, 32.620), 2.77, 1e-7, STD, (1.48, 1.04, 1.61, .982, 1.34, 2.24), CR, "18o"),
    (901, 41.39, 804, (74.859, -53.254, 33.532), 1.01, 1e-8, STD, (4.91, 2.75, 5.98, 2.72, 4.00, 9.71), CR, "19"),
    (901, 42.54, 804, (74.856, -52.225, 34.217), 1.85, 1e-7, STD, (3.21, 1.60, 3.94, 1.48, 1.94, 6.15), CR, "20y"),
    (901, 43.79, 804, (74.775, -51.610, 35.288), 1.06, 1e-7, STD, (5.71, .227, 8.84, .258, .295, 14.6), CR, "20o"),
    # SOM-NUB (Lemaux 2002; same rotation at 10.95 and 250 Ma)
    (709, 10.95, 701, (-9.90, -37.0, 0.22), 1.0, 1e-6, STD, (2.081, .779, -3.249, .474, -1.069, 5.863), "Lemaux et al. 2002", "C5"),
    (709, 250.0, 701, (-9.90, -37.0, 0.22), 1.0, 1e-6, STD, (2.081, .779, -3.249, .474, -1.069, 5.863), "Lemaux et al. 2002", "C5 (constante)"),
    # LHR-AUS (Gaina 1998, Australia fixed)
    (833, 55.9, 801, (-15.93, 133.47, -2.112), 2.11, 1e-6, STD, (8.027, -3.805, 6.735, 2.422, -3.710, 6.265), GA, "25y"),
    (833, 67.7, 801, (-9.04, 134.46, -8.83), .99, 1e-6, STD, (2.100, -1.014, 1.981, 5.707, -1.028, 1.954), GA, "31y"),
    (833, 73.6, 801, (-9.53, 137.20, -12.937), .61, 1e-6, STD, (2.838, -1.489, 2.804, 8.307, -1.490, 2.800), GA, "33y"),
    # AUS-EANT (Tikku and Cande 2000; published rotation inverse to that of the .rot)
    (801, 43.79, 802, (15.07, 31.78, 24.55), 1.08, 1e-7, STD, (4.75, -6.85, 5.49, 10.4, -8.84, 17.9), TK, "20o inv"),
    (801, 46.26, 802, (14.00, 33.34, 24.70), 1.06, 1e-7, STD, (5.63, -7.75, 12.1, 11.2, -16.0, 102.7), TK, "21y inv"),
    (801, 53.35, 802, (10.39, 35.59, 25.15), 1.09, 1e-6, STD, (.798, -1.22, 1.43, 1.99, -2.30, 9.50), TK, "24o inv"),
    (801, 60.92, 802, (9.95, 36.52, 25.55), 1.06, 1e-6, STD, (1.87, -3.03, 11.9, 5.15, -20.3, 100.1), TK, "27y inv"),
    (801, 79.08, 802, (5.125, 39.797, 26.572), 1.54, 1e-6, STD, (3.18, -5.14, 10.9, 8.56, -18.2, 62.7), TK, "33o inv"),
    (801, 83.5, 802, (2.054, 40.788, 27.119), 1.00, 1e-5, STD, (1.78, -2.63, 9.55, 3.96, -14.5, 64.7), TK, "34y inv"),
]
# EANT-WANT 13o (Cande et al. 2000, Tab. 1-2), basis of the derived cov of node 43.8
C00_POLE, C00_KAPPA, C00_W = (-18.146, -17.847, -0.696), 0.63, 1e-5
C00_ABC = (2.19, .0039, 5.74, .0041, .0083, 15.1)


def _eng():
    # inside the pipeline the already loaded engine is used (that of the model folder)
    if "02_rot_engine" in sys.modules:
        return sys.modules["02_rot_engine"]
    P = (__import__("os").environ["PAPER_ROOT"] + "/models/Q22_WK08-A")
    if P not in sys.path:
        sys.path.insert(0, P)
    import importlib
    return importlib.import_module("02_rot_engine")


def mat(order, v):
    d = dict(zip(order, v))
    return np.array([[d["xx"], d["xy"], d["xz"]],
                     [d["xy"], d["yy"], d["yz"]],
                     [d["xz"], d["yz"], d["zz"]]])


def pole_matrix(eng, la, lo, ang):
    return eng.quat_to_matrix(eng.quat_from_pole(la, lo, ang))


def rotvec(eng, A):
    q = eng.matrix_to_quat(A)
    s = np.linalg.norm(q[1:])
    return np.zeros(3) if s < 1e-15 else q[1:] / s * 2.0 * np.arctan2(s, q[0])


def _find(M, plate, age, fixed):
    """Node (plate, age, fixed). The ages in this file are on the timescale of
    each source (CK95 or GTS2012); if the age is not found as is, its GTS2020 age is
    searched (rot_gts2020.py): CK95 -> GTS2020, GTS2012 -> GTS2020 or
    top of C34n (83.65). If more than one candidate matches different nodes,
    error."""
    rows = M.table.get(plate, [])
    for i, r in enumerate(rows):
        if abs(r[0] - age) < 1e-6 and r[4] == fixed:
            return i, r
    hit = {}
    for a in _edades_gts20(age):
        for i, r in enumerate(rows):
            if abs(r[0] - a) < 6e-4 and r[4] == fixed:
                hit[i] = r
    if len(hit) > 1:
        raise KeyError(f"ambiguous node in GTS2020: {plate} {age} fixed {fixed}")
    return next(iter(hit.items()), (None, None))


_GTS20 = None


def _edades_gts20(age):
    global _GTS20
    if _GTS20 is None:
        import importlib.util
        here = os.path.dirname(os.path.abspath(__file__))
        s = importlib.util.spec_from_file_location("esc", os.path.join(here, "timescale_GTS2020.py"))
        esc = importlib.util.module_from_spec(s)
        s.loader.exec_module(esc)
        _GTS20 = (esc.remap("ck95"), esc.remap("gts12"), esc.C34_TECHO["gts20"])
    fck, f12, c34 = _GTS20
    out = [float(fck(age)), float(f12(age))]
    if age in (83.0, 83.5):
        out.append(c34)
    return out


# Rows that rot_gts2020.py removes because they are the same chron dated twice
# (C10n.2n(o): Croon 901@28.75 CK95 = Quiero 901@28.278 GTS2012).
BORRADAS_GTS20 = {(901, 28.75, 804)}


def published(eng, M):
    """[(plate, age, fixed, C, source, note, misfit_degrees)]"""
    out = []
    for (p, age, fx, pole, kap, w, order, v, src, note) in PUB:
        C = mat(order, v) * w / kap
        if np.linalg.eigvalsh(C).min() <= 0:
            raise ValueError(f"not positive definite: {p} {age}")
        i, r = _find(M, p, age, fx)
        if r is None and p == 201 and any(x[4] == 911 for x in M.table[201]):
            continue    # .rot SG12 (build_rot_sg12.py): SAM rel NZ replaces SAM-AFR
        if r is None and (p, age, fx) in BORRADAS_GTS20:
            continue
        if r is None:
            raise KeyError(f"node missing from the .rot: {p} {age} fixed {fx}")
        A_rot = pole_matrix(eng, r[1], r[2], r[3])
        A_pub = pole_matrix(eng, *pole)
        if note.endswith("inv"):
            # .rot node = inverse of the published one
            C = A_pub.T @ C @ A_pub
            A_pub = A_pub.T
        mis = np.degrees(np.linalg.norm(rotvec(eng, A_rot @ A_pub.T)))
        out.append((p, age, fx, C, src, note, mis))
    return out


def derived_804(eng, M):
    """Node 804@43.8: 13o of Cande 2000 inverted and scaled to the angle of the .rot."""
    C = mat(STD, C00_ABC) * C00_W / C00_KAPPA
    A = pole_matrix(eng, *C00_POLE)
    C = A.T @ C @ A
    i, r = _find(M, 804, 43.8, 802)
    if r is None:   # Granot 2013 sensitivity .rot (option b): node replaced
        return None
    k = r[3] / -C00_POLE[2]
    return (804, 43.8, 802, k * k * C, "Cande et al. 2000 (derivada)",
            f"13o invertida, escala ({r[3]}/{-C00_POLE[2]})^2", np.nan)


def interp_901_3355(ent):
    """Node 901@33.55 (Croon 13o, printed cov not positive definite): linear
    interpolation between the covs of 13y (33.06) and 15y (34.66), same frame (804 fixed)."""
    C = {a: c for (p, a, fx, c, *_ ) in ent if p == 901 and a in (33.06, 34.66)}
    f = (33.55 - 33.06) / (34.66 - 33.06)
    return (901, 33.55, 804, (1 - f) * C[33.06] + f * C[34.66], CR + " (interpolada)",
            "13o: interpolada entre 13y y 15y", np.nan)


HEL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs", "frame_coverage", "hellinger")


def hellinger_nodes(variante):
    """C32y nodes fitted with our own Hellinger fit (hellinger_circuit.py). The cov
    in the json is already in the frame of the fixed plate and on the sigma^2 scale.
      801@71.07: .rot pole inside the 95% region -> cov assigned to the .rot node.
      802@71.1 : .rot pole (Bernard 2005) outside the 95% region. Variants
                 (both are tested):
                 "a" replaces the pole by the fit and uses its cov;
                 "b" keeps the pole and scales the cov by q/threshold, so that
                     the .rot pole lies on the edge of the 95% region.
    Returns [(plate, age, fixed, C, new_pole or None, source, note)]."""
    out = []
    for caso, age in (("801_71.07", 71.07), ("802_71.1", 71.1)):
        r = json.load(open(os.path.join(HEL, f"fit_{caso}.json")))
        # picks of the African flank in polygon 701; the .rot node uses 709,
        # equal to 701 before 11 Ma under the assumption of a rigid Africa
        r["fija"] = {701: 709}.get(r["fija"], r["fija"])
        C = np.array(r["cov_fija_rad2"])
        src = "Hellinger propio C32y (" + ", ".join(r["refs"]) + ")"
        if r["rot_dentro_95"]:
            out.append((r["movil"], age, r["fija"], C, None, src, "polo del .rot dentro de 95%"))
        elif variante == "a":
            out.append((r["movil"], age, r["fija"], C, tuple(r["polo_ajuste"]), src,
                        "variante a: polo reemplazado por el ajuste"))
        else:
            k = r["chi2_rot_en_region"] / r["umbral_95"]
            out.append((r["movil"], age, r["fija"], k * C, None, src,
                        f"variante b: cov x {k:.2f} (polo del .rot en el borde de 95%)"))
    return out


def rel_cov(eng, M, plate, anchor, t, eps=1e-6):
    """Cov (frame of `anchor`) of the rotation of `plate` relative to `anchor` at
    age t, linearizing over the nodes with covariance (Kirkwood et al. 1999)."""
    q0 = eng.quat_to_matrix(M.rel_rotation(plate, t, anchor))
    C = np.zeros((3, 3))
    for p, rows in M.table.items():
        for (age, la, lo, ang, fx, Cv) in rows:
            if Cv is None:
                continue
            J = np.zeros((3, 3))
            for k in range(3):
                h = np.zeros(3); h[k] = eps
                A = eng.quat_to_matrix(M.rel_rotation(plate, t, anchor, {(p, age, fx): h}))
                J[:, k] = rotvec(eng, A @ q0.T) / eps
            if np.abs(J).max() > 1e-9:
                C += J @ Cv @ J.T
    return C


def aplicar(M, verbose=False, hellinger="b"):
    """Adds to M.table the published covariances, the derived ones and those of
    our own Hellinger fit (variant "a" or "b" for 802@71.1; None to omit them).
    Returns the list of entries (plate, age, fixed, C, source, note, misfit)."""
    eng = _eng()
    ent = published(eng, M)
    ent += [e for e in (derived_804(eng, M), interp_901_3355(ent)) if e is not None]
    for (p, age, fx, C, *_ ) in ent:
        i, r = _find(M, p, age, fx)
        M.table[p][i] = r[:5] + (C,)
    if hellinger:
        for (p, age, fx, C, polo, src, note) in hellinger_nodes(hellinger):
            i, r = _find(M, p, age, fx)
            if polo is not None:
                r = (r[0],) + polo + (r[4],)
            M.table[p][i] = r[:5] + (C,)
            ent.append((p, age, fx, C, src, note, np.nan))
    # PAC-LHR: fixed rotation 45.03-83.5 (LHR welded to the Pacific); cov along the
    # route PAC-WANT-EANT-AUS at 43.79 Ma, with the previous covs already loaded
    Cl = rel_cov(eng, M, 901, 833, _find(M, 901, 43.79, 804)[1][0])   # node C20n(o)
    for age in (45.03, 83.5):
        i, r = _find(M, 901, age, 833)
        M.table[901][i] = r[:5] + (Cl,)
        ent.append((901, age, 833, Cl, "propagada PAC-WANT-EANT-AUS 43.79 Ma",
                    "derivada (LHR soldado a PAC)", np.nan))
    if verbose:
        for (p, age, fx, C, src, note, mis) in ent:
            s = np.sqrt(np.diag(C))
            print(f"{p:4d} {age:7.2f} fixed {fx:4d}  sd_max {np.degrees(s.max())*111.195:7.2f} km"
                  f"  pole misfit {mis:6.3f} deg  {src} ({note})")
    return ent


if __name__ == "__main__":
    import pandas as pd
    eng = _eng()
    M = eng.RotModel((__import__("os").environ["PAPER_ROOT"] + "/models/Q22_WK08-A/rot/Global_Model_R_BG_2025_hsb.rot"))
    ent = aplicar(M, verbose=True)
    rows = [dict(placa=p, edad=a, fija=f, Cxx=C[0, 0], Cxy=C[0, 1], Cxz=C[0, 2],
                 Cyy=C[1, 1], Cyz=C[1, 2], Czz=C[2, 2], fuente=s, nota=n, desajuste_gr=m)
            for (p, a, f, C, s, n, m) in ent]
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs", "frame_coverage")
    os.makedirs(out, exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(out, "circuit_covariance_agregadas.csv"), index=False)
