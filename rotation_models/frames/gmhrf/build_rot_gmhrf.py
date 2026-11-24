# -*- coding: utf-8 -*-
"""Versions of the rotation model with the GMHRF absolute frame of
Doubrovine et al. (2012) in place of WK08-A (Wessel and Kroenke 2008).

Source: supplementary material of Doubrovine, Steinberger and Torsvik (2012),
JGR 117, B09101, doi:10.1029/2011JB009072, in
  external/doubrovine2012/Supporting_Information/
  - jgrb17279-sup-0004-ts03.txt  Table S3, finite rotations in the GMHRF,
                                 block "Pacific" (10-80 Ma, with covariance) and
                                 block "Pacific - ... Wessel and Kroenke
                                 (2008)" (90-120 Ma, without covariance)
  - jgrb17279-sup-0008-ds01.txt  Data Set S1, modelled motion of the
                                 hotspots in the GMHRF (Easter, Juan Fernandez,
                                 San Felix among them)

WHERE THE FRAME ENTERS
In the .rot the hotspots (plates 13 JFR, 14 Easter, 15 San Felix) hang from
plate 002 ("Pacific hotspots"), which in the base model is WK08-A expressed as
hotspots relative to the Pacific (fixed 901). When the tracks are viewed from
South America (201) the global frame 000 cancels: the absolute frame enters ONLY
through plate 002. Changing the frame therefore means rewriting plate 002.

VERSIONS (outputs in rot/)
  A  Global_Model_R_BG_2025_GMHRF.rot
       plate 002 <- inverse of Table S3 "Pacific" (10-80 Ma, with covariance)
       and of the WK08 variant of the same table (90-120 Ma, without covariance).
       The WK08-A nodes older than 120 Ma are kept (no tau of the model scripts
       exceeds 120 Ma). Hotspots FIXED in the GMHRF.
  A' Global_Model_R_BG_2025_GMHRF_refit.rot
       A + refit of anchors 13/14/15 to the geochronology.
  B  Global_Model_R_BG_2025_GMHRF_moving_refit.rot
       A' + hotspot drift from Data Set S1 starting at the age of the
       last node that constrains each anchor (JFR and San Felix 9.26 Ma, Easter
       27.23 Ma). Before that age the geochronology controls (the refit already
       passes the track through the dated edifices); after it, where the anchor
       would otherwise stay frozen, the hotspot moves as the flow model predicts.
       This is the self-consistent version: the GMHRF was fitted assuming that drift.

COVARIANCE CONVENTION (checked in the main text, eq. 1 and the following
paragraph): A_hat = A_0 F(h), perturbation on the RIGHT, h in the frame of the
reconstructed plate. The rotation engine of the models uses R = exp(M(h)) R_nom
(on the left). For plate 002, which is the INVERSE of the published rotation:
    (A F(h))^-1 = F(-h) A^-1   ->   perturbation on the left with h' = -h,
so the covariance passes UNTRANSFORMED: C' = C.
(For WK08-A, which is in the left convention, step 04 applies A^T C A.)

DRIFT (version B): new plates 1013/1014/1015 fixed to 002, identity
up to the start age and then the minimal rotation that takes position
D1(t0) to D1(t), conjugated with the rotation S that takes D1(t0) to the fitted
anchor (S D S^T): in this way the displacement of D1 is applied at the actual
position of the model hotspot, which can be 1-4 degrees away from that of D1.
Plates 13 and 14 hang from 1013/1014; plate 15 hangs from
1015 in 0-2.93 Ma, keeps its segment anchored to 13 in 2.93-9.26 Ma (the rate
borrowed from JFR) and after 9.26 Ma hangs from 1015 with the equivalent
rotation, to follow its own drift.
Repeated D1 rows at the same age (conduit segments) are averaged.
After the last D1 age (plume age, Table S6) the drift is frozen.

Usage:  py -I build_rot_gmhrf.py
"""
import os
import io
import sys
import importlib
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "models", "Q22_WK08-A"))
eng = importlib.import_module("02_rot_engine")

SI = os.path.join(HERE, "..", "..", "..", "external", "doubrovine2012", "Supporting_Information")
S3 = os.path.join(SI, "jgrb17279-sup-0004-ts03.txt")
D1 = os.path.join(SI, "jgrb17279-sup-0008-ds01.txt")
# PAPER_GMHRF_BASE allows building the GMHRF on another version of the
# circuit (e.g. a WD25 .rot whose anchors are already fitted to the
# geochronology with that circuit) and
# PAPER_GMHRF_TAG gives the prefix of the outputs (e.g. "WD25_").
BASE = os.environ.get("PAPER_GMHRF_BASE",
                      __import__("os").environ["PAPER_ROOT"] + "/inputs/rotations/Global_Model_R_BG_2025.rot")
TAG = os.environ.get("PAPER_GMHRF_TAG", "")
WK08_CSV = os.path.join(ROOT, "models", "Q22_WK08-A", "data", "wk08a_covariances.csv")
OUTDIR = os.path.join(HERE, "rot")
OUT_A = os.path.join(OUTDIR, f"Global_Model_R_BG_2025_{TAG}GMHRF.rot")
OUT_AR = os.path.join(OUTDIR, f"Global_Model_R_BG_2025_{TAG}GMHRF_refit.rot")
OUT_B = os.path.join(OUTDIR, f"Global_Model_R_BG_2025_{TAG}GMHRF_moving_refit.rot")

FRAME = 2
REL = 902003                      # f_PAC = 0.44, reference model
SEEDS = {13: (-33.750097, -81.905212),    # Juan Fernandez
         14: (-26.3138, -106.599),        # Easter
         15: (-26.28791, -81.42893)}      # San Felix / Taltal
D1_NAME = {13: "Juan_Fernandez", 14: "Easter", 15: "San_Felix"}
DRIFT_PLATE = {13: 1013, 14: 1014, 15: 1015}
T0_DRIFT = {13: 9.26, 14: 27.23, 15: 9.26}
SENTINEL = 250.0
FRAME_NAME = "GMHRF"   # label in the refit comments (build_rot_gwg20.py changes it)
REF = '@REF Doubrovine_++_2012 @DOI"10.1029/2011JB009072"'


# ----------------------------------------------------------------- reading
def read_s3_block(header):
    """Rows (age, lat, lon, rho, C|None) of the Table S3 block whose
    header is exactly `header`."""
    rows, on = [], False
    for line in io.open(S3, encoding="latin-1"):
        s = line.strip()
        if not s or s.startswith("-"):
            continue
        p = s.split()
        try:
            age = float(p[0])
        except ValueError:
            on = (s == header)
            continue
        if not on:
            continue
        lat, lon, rho = float(p[1]), float(p[2]), float(p[3])
        C = None
        if p[4] != "-":
            c = [float(x) for x in p[4:10]]
            C = np.array([[c[0], c[1], c[2]], [c[1], c[3], c[4]], [c[2], c[4], c[5]]])
        rows.append((age, lat, lon, rho, C))
    if not rows:
        raise ValueError(f"block not found in Table S3: {header!r}")
    return rows


def read_d1(name):
    """Modelled hotspot motion (Data Set S1), with the rows repeated at
    the same age averaged as unit vectors."""
    recs = []
    for line in io.open(D1, encoding="latin-1"):
        p = line.split()
        if len(p) == 4 and p[0] == name:
            recs.append((float(p[1]), float(p[2]), float(p[3])))
    df = pd.DataFrame(recs, columns=["age", "lat", "lon"])
    out = []
    for age, g in df.groupby("age"):
        v = np.mean([eng.latlon_axis(a, b) for a, b in zip(g.lat, g.lon)], axis=0)
        la, lo = eng.xyz_to_latlon(v / np.linalg.norm(v))
        out.append((age, la, lo, len(g)))
    return pd.DataFrame(out, columns=["age", "lat", "lon", "n"])


def d1_at(d1, t):
    """D1 position interpolated at age t (linear in the unit vector)."""
    a = d1.age.values
    if t >= a[-1]:
        return eng.latlon_axis(d1.lat.values[-1], d1.lon.values[-1])
    k = np.searchsorted(a, t)
    if a[k] == t or k == 0:
        return eng.latlon_axis(d1.lat.values[k], d1.lon.values[k])
    f = (t - a[k - 1]) / (a[k] - a[k - 1])
    v = ((1 - f) * eng.latlon_axis(d1.lat.values[k - 1], d1.lon.values[k - 1])
         + f * eng.latlon_axis(d1.lat.values[k], d1.lon.values[k]))
    return v / np.linalg.norm(v)


# ---------------------------------------------------------------- algebra
def minimal_rotation(u, v):
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)
    ax = np.cross(u, v)
    n = np.linalg.norm(ax)
    ang = np.degrees(np.arctan2(n, float(np.dot(u, v))))
    if n < 1e-14:
        return 90.0, 0.0, 0.0
    la, lo = eng.xyz_to_latlon(ax / n)
    return la, lo, ang


def pole_of(R):
    q = eng.matrix_to_quat(R)
    if q[0] < 0:
        q = -q
    s = np.linalg.norm(q[1:])
    if s < 1e-14:
        return 90.0, 0.0, 0.0
    la, lo = eng.xyz_to_latlon(q[1:] / s)
    return la, lo, np.degrees(2 * np.arctan2(s, q[0]))


def inverse_pole(lat, lon, rho):
    """Inverse rotation as the antipodal pole with positive angle."""
    return -lat, (lon + 180.0) if lon <= 0 else (lon - 180.0), rho


def fmt_cov(C):
    return "   ".join(f"{x:.4E}" for x in
                      (C[0, 0], C[0, 1], C[0, 2], C[1, 1], C[1, 2], C[2, 2]))


def gc_deg(u, v):
    return np.degrees(np.arccos(np.clip(np.dot(u, v), -1, 1)))


# ------------------------------------------------------------- version A
def plate2_lines():
    pac = read_s3_block("Pacific")
    pac_wk = read_s3_block("Pacific - calculated using stage rotations of "
                           "Wessel and Kroenke (2008)")
    out = [f"002  0.0   90.0    0.0    0.0  901 !PHS-PAC GMHRF hotspots-Pacific, "
           f"inverse of Pacific in Table S3 {REF}\n"]
    for (age, la, lo, rho, C) in pac + pac_wk:
        ila, ilo, irho = inverse_pole(la, lo, rho)
        tag = ("Table S3 Pacific (Australian circuit 1), cov h'=-h, C'=C"
               if C is not None else
               "Table S3 Pacific via WK08 stage rotations, no covariance")
        cov = f"   {fmt_cov(C)}" if C is not None else ""
        out.append(f"002 {age:6.2f}   {ila:8.4f}  {ilo:9.4f}  {irho:7.3f}  901 "
                   f"!PHS-PAC {REF} @C\"{tag}\"{cov}\n")
    return out, pac, pac_wk


def write_version_a():
    new2, pac, pac_wk = plate2_lines()
    tmax = max(r[0] for r in pac + pac_wk)
    lines, done = [], False
    for line in io.open(BASE, encoding="latin-1"):
        p = line.split()
        if len(p) >= 6 and p[0].isdigit() and int(p[0]) == FRAME:
            if float(p[1]) <= tmax + 1e-9:
                if not done:
                    lines.extend(new2)
                    done = True
                continue
            # WK08-A nodes > 120 Ma: kept, and labelled as such
            line = line.rstrip("\n") + "  [kept from WK08-A beyond GMHRF Table S3]\n"
        lines.append(line)
    io.open(OUT_A, "w", encoding="latin-1", newline="\n").writelines(lines)
    return pac, pac_wk


# ----------------------------------------------------------- refit A'
def present_position(M, plate, tau):
    seed = eng.latlon_axis(*SEEDS[plate])
    p_abs = eng.quat_to_matrix(M.total_rotation(plate, tau, 0)) @ seed
    return eng.quat_to_matrix(M.total_rotation(REL, tau, 0)).T @ p_abs


def fit_node(M, plate, tau, x_obs):
    seed = eng.latlon_axis(*SEEDS[plate])
    R_rel = eng.quat_to_matrix(M.total_rotation(REL, tau, 0))
    R_frame = eng.quat_to_matrix(M.total_rotation(FRAME, tau, 0))
    return minimal_rotation(seed, R_frame.T @ (R_rel @ x_obs))


def write_refit(Mold, src, out):
    Mnew = eng.RotModel(src)
    new_nodes, report = {}, []
    for plate in (13, 14, 15):
        rows = [r for r in Mold.table[plate] if 0.0 < r[0] < SENTINEL and r[4] == FRAME]
        out_rows, errA = [], 0.0
        for (tau, la0, lo0, ang0, fx, C) in rows:
            x_obs = present_position(Mold, plate, tau)
            la_c, lo_c, ang_c = fit_node(Mold, plate, tau, x_obs)
            s = eng.latlon_axis(*SEEDS[plate])
            Rc = eng.quat_to_matrix(eng.quat_from_pole(la_c, lo_c, ang_c))
            R0 = eng.quat_to_matrix(eng.quat_from_pole(la0, lo0, ang0))
            errA = max(errA, float(np.abs(Rc @ s - R0 @ s).max()))
            out_rows.append((tau,) + fit_node(Mnew, plate, tau, x_obs))
        new_nodes[plate] = out_rows
        report.append(f"  plate {plate}: {len(out_rows):2d} nodes, inversion of the "
                      f"old node {errA:.1e}")
    lines, done = [], set()
    for line in io.open(src, encoding="latin-1"):
        p = line.split()
        if len(p) >= 6 and p[0].isdigit() and int(p[0]) in new_nodes:
            plate, tau, fixed = int(p[0]), float(p[1]), int(p[5])
            if fixed == FRAME and 0.0 < tau < SENTINEL:
                if plate not in done:
                    done.add(plate)
                    for (t, la, lo, an) in new_nodes[plate]:
                        lines.append(f"{plate}   {t:.5f}   {la:.6f}   {lo:.6f}   "
                                     f"{an:.6f}   {FRAME}   !HOTSPOT this study, "
                                     f"re-fitted to the {FRAME_NAME} Pacific frame\n")
                continue
            if fixed == FRAME and tau >= SENTINEL:
                t, la, lo, an = new_nodes[plate][-1]
                lines.append(f"{plate}   {tau:.1f}   {la:.6f}   {lo:.6f}   {an:.6f}   "
                             f"{FRAME}   !HOTSPOT this study, sentinel (re-fitted)\n")
                continue
        lines.append(line)
    io.open(out, "w", encoding="latin-1", newline="\n").writelines(lines)
    Mout = eng.RotModel(out)
    worst = 0.0
    for plate, rows in new_nodes.items():
        for (tau, *_r) in rows:
            u = present_position(Mold, plate, tau)
            v = present_position(Mout, plate, tau)
            worst = max(worst, gc_deg(u, v) * np.pi / 180 * 6371.0)
    report.append(f"  the new track passes through the target: max {worst:.1e} km")
    return report


# ------------------------------------------------------------- version B
def write_version_b(src, out):
    M = eng.RotModel(src)
    drift_lines, report = [], []
    for plate in (13, 14, 15):
        d1 = read_d1(D1_NAME[plate])
        t0, pid = T0_DRIFT[plate], DRIFT_PLATE[plate]
        u0 = d1_at(d1, t0)
        tmax = float(d1.age.max())
        # separation D1(t0) vs fitted hotspot position in frame 002
        R13 = eng.quat_to_matrix(M.finite_rotation(plate, t0)[0])
        if plate == 15:
            R13 = eng.quat_to_matrix(M.total_rotation(15, t0, FRAME))
        p_fit = R13 @ eng.latlon_axis(*SEEDS[plate])
        sep = gc_deg(u0, p_fit)
        ages = np.arange(np.ceil(t0), np.floor(tmax) + 1, 1.0)
        drift_lines.append(f"{pid}  0.0   90.0   0.0   0.0   {FRAME}   !{D1_NAME[plate]} "
                           f"drift in GMHRF, identity before {t0} Ma {REF}\n")
        drift_lines.append(f"{pid}  {t0:.2f}   90.0   0.0   0.0   {FRAME}   !drift start\n")
        # S takes D1(t0) to the fitted anchor; the drift is conjugated with S so that the
        # displacement of D1 is applied at the position of the anchor
        S = eng.quat_to_matrix(eng.quat_from_pole(*minimal_rotation(u0, p_fit)))
        last = (90.0, 0.0, 0.0)
        for t in ages:
            D = eng.quat_to_matrix(eng.quat_from_pole(*minimal_rotation(u0, d1_at(d1, t))))
            la, lo, an = pole_of(S @ D @ S.T)
            last = (la, lo, an)
            drift_lines.append(f"{pid}  {t:.2f}   {la:.6f}   {lo:.6f}   {an:.6f}   "
                               f"{FRAME}   !Data Set S1 drift rel. {t0} Ma\n")
        drift_lines.append(f"{pid}  {SENTINEL:.1f}   {last[0]:.6f}   {last[1]:.6f}   "
                           f"{last[2]:.6f}   {FRAME}   !drift frozen after plume age "
                           f"{tmax:.1f} Ma\n")
        report.append(f"  {D1_NAME[plate]:15s} t0={t0:5.2f} Ma, D1 up to {tmax:5.1f} Ma, "
                      f"total drift {last[2]:.2f} deg, D1(t0) vs fitted anchor "
                      f"{sep:.2f} deg (drift moved to the anchor by conjugation)")

    # equivalent rotation of plate 15 at 9.26 Ma relative to 002
    R15eq = eng.quat_to_matrix(M.total_rotation(15, T0_DRIFT[15], FRAME))
    la15, lo15, an15 = pole_of(R15eq)

    lines, inserted = [], False
    for line in io.open(src, encoding="latin-1"):
        p = line.split()
        if len(p) >= 6 and p[0].isdigit() and int(p[0]) in (13, 14, 15):
            plate, tau, fixed = int(p[0]), float(p[1]), int(p[5])
            if plate in (13, 14) and fixed == FRAME:
                line = _refix(line, DRIFT_PLATE[plate])
            elif plate == 15:
                if fixed == FRAME:
                    line = _refix(line, DRIFT_PLATE[15])
                elif fixed == 13 and tau >= SENTINEL:
                    # closes the segment anchored to 13 at 9.26 Ma and switches to 1015
                    t9 = T0_DRIFT[15]
                    C13 = [r for r in M.table[15] if r[4] == 13][0]
                    lines.append(f"15   {t9:.2f}   {C13[1]}   {C13[2]}   {C13[3]}   13   "
                                 f"!SFH-PHS end of JFR-anchored segment\n")
                    for t in (t9, SENTINEL):
                        lines.append(f"15   {t:.2f}   {la15:.6f}   {lo15:.6f}   {an15:.6f}   "
                                     f"{DRIFT_PLATE[15]}   !SFH equivalent at {t9} Ma, "
                                     f"follows its own Data Set S1 drift\n")
                    continue
            if plate == 13 and not inserted:
                lines.extend(drift_lines)
                inserted = True
        lines.append(line)
    io.open(out, "w", encoding="latin-1", newline="\n").writelines(lines)
    return report


def _refix(line, new_fixed):
    p = line.split(None, 6)
    rest = p[6] if len(p) > 6 else "\n"
    return f"{p[0]}   {p[1]}   {p[2]}   {p[3]}   {p[4]}   {new_fixed}   {rest}"


# ------------------------------------------------------------ checks
def checks(pac):
    print("\n[1] Table S3 Pacific vs Pacific in 000 of the base .rot "
          "(GMHRF via Africa + base circuit):")
    Mb = eng.RotModel(BASE)
    for (age, la, lo, rho, C) in pac:
        Rs = eng.quat_to_matrix(eng.quat_from_pole(la, lo, rho))
        Rb = eng.quat_to_matrix(Mb.total_rotation(901, age, 0))
        print(f"    {age:4.0f} Ma: difference {pole_of(Rs.T @ Rb)[2]:6.3f} deg")
    print("\n[2] Maximum standard deviation of the frame (sqrt of max eigenvalue, deg):")
    wk = pd.read_csv(WK08_CSV, comment="#")
    for (age, la, lo, rho, C) in pac:
        r = wk.iloc[(wk.age_ma - age).abs().argmin()]
        Cw = 1e-5 * np.array([[r.a, r.b, r.d], [r.b, r.c, r.e], [r.d, r.e, r.f]])
        print(f"    {age:4.0f} Ma: GMHRF {np.degrees(np.sqrt(np.linalg.eigvalsh(C).max())):5.2f}"
              f" | WK08-A @ {r.age_ma:6.2f} Ma "
              f"{np.degrees(np.sqrt(np.linalg.eigvalsh(Cw).max())):5.2f}")


def track_shift(files):
    print("\n[3] Shift of the present-day track position (f_PAC 0.44) "
          "relative to the base model, km:")
    Mb = eng.RotModel(BASE)
    Ms = {k: eng.RotModel(v) for k, v in files.items()}
    taus = (2, 5, 9, 15, 20, 27, 30, 40, 50, 60)
    print("    " + " " * 22 + "".join(f"{t:>7d}" for t in taus))
    for plate in (13, 14, 15):
        for k, M in Ms.items():
            d = [gc_deg(present_position(Mb, plate, t), present_position(M, plate, t))
                 * np.pi / 180 * 6371.0 for t in taus]
            print(f"    {D1_NAME[plate]:15s}{k:>7s}" + "".join(f"{x:7.0f}" for x in d))


def jump_scan(files):
    """Maximum jump of the present-day track position between consecutive ages
    (step 0.5 Ma). It reveals badly built handover nodes, such as the plate 911 bug
    in WD25_WK08-A."""
    print("\n[4] Jump scan of the nominal track (f_PAC 0.44), step 0.5 Ma, km:")
    taus = np.arange(0.0, 66.01, 0.5)
    for k, fn in files.items():
        M = eng.RotModel(fn)
        for plate in (13, 14, 15):
            P = [present_position(M, plate, t) for t in taus]
            d = np.array([gc_deg(P[i], P[i + 1]) for i in range(len(P) - 1)]) * np.pi / 180 * 6371
            i = int(d.argmax())
            print(f"    {k:>5s} {D1_NAME[plate]:15s} max {d[i]:6.1f} km between "
                  f"{taus[i]:.1f} and {taus[i+1]:.1f} Ma | median {np.median(d):5.1f} km")


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    pac, _ = write_version_a()
    print(f"A  written: {OUT_A}")
    Mold = eng.RotModel(BASE)
    for s in write_refit(Mold, OUT_A, OUT_AR):
        print(s)
    print(f"A' written: {OUT_AR}")
    for s in write_version_b(OUT_AR, OUT_B):
        print(s)
    print(f"B  written: {OUT_B}")
    checks(pac)
    track_shift({"A": OUT_A, "A'": OUT_AR, "B": OUT_B})
    jump_scan({"base": BASE, "A'": OUT_AR, "B": OUT_B})


if __name__ == "__main__":
    main()
