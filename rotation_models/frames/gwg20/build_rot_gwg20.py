# -*- coding: utf-8 -*-
"""Versions of the rotation model with the GWG20 absolute frame (Gaastra,
Gordon and Woodworth 2022) in place of WK08-A (Wessel and Kroenke 2008).

Source: supplementary material of Gaastra KM, Gordon RG, Woodworth DT (2022)
Quantification of Pacific Plate hotspot tracks since 80 Ma. Tectonics 41,
e2021TC006772, doi:10.1029/2021TC006772, in
  external/gaastra2022/Supporting_Information/
  - 2021tc006772-sup-0003-table si-s02.csv  = their Table S6 ("PA_nhotspot_inversion"),
    "Reconstructions of the Pacific plate relative to Pacific hotspots every
    1 Ma between 80 Ma and the present", N-hotspot method of Andrews et al.
    (2006) with Hawaii, Louisville and Rurutu (Rurutu excluded between 10 and 50 Ma),
    FIXED hotspots.
  - 2021tc006772-sup-0001 (pdf), caption of Table S6, note b:
        Cov = (1/kappa_hat) [[a, b, c], [b, d, e], [c, e, f]] x 1e-5 rad^2
    ("dividing by kappa_hat makes chi2 = 1, which corrects for the
    over-dispersion of the problem but is not required if a more conservative
    estimate of uncertainty is desired"). The printed caption also says
    "radians^2 Ma^-2"; for a finite rotation rad^2 is taken. The same label
    is in Table 3 of Andrews et al. (2006); read literally it gives impossible
    sigmas (analysis/gwg20_convention.py).
    The order [[a,b,c],[b,d,e],[c,e,f]] is confirmed because it is the only one of the
    two possible orders that gives positive definite matrices (that of Chang/WK08-A,
    [[a,b,d],[b,c,e],[d,e,f]], gives negative diagonals in this table).

WHERE THE FRAME ENTERS: as in ../gmhrf/build_rot_gmhrf.py, only through
plate 002 (hotspots relative to the Pacific, fixed 901), which is rewritten with the
INVERSE of each rotation of Table S6 (0-80 Ma; the WK08-A nodes older than
80 Ma are kept, and labelled as such).

COVARIANCE CONVENTION: that of Andrews et al. (2006), the source of the N-hotspot
method: u is applied AFTER the PAC->HS rotation (eigenvectors "in fixed hotspot
reference frame", their Fig. 9), A = Phi(u) A_hat. For the inverse
B = A^-1 = Phi(-B_hat u) B_hat, perturbation on the left (that of the rotation engine
of the models) with h' = -B_hat u, so C' = B_hat C B_hat^T.
analysis/gwg20_rotate_covariance.py applies this rotation to the covariances of an
existing .rot.

VERSIONS (outputs in rot/; the hotspots are fixed, there is no version with drift)
  Global_Model_R_BG_2025_{TAG}GWG20.rot        plate 002 <- GWG20, without refit
  Global_Model_R_BG_2025_{TAG}GWG20_refit.rot  + anchors 13/14/15 refitted to the
                                               geochronology (same construction as
                                               build_rot_gmhrf.write_refit)
Empty TAG = circuit of Quiero et al. (2022); TAG "WD25_" = circuit of Wilson
and DeMets (2025) (base models/WD25_WK08-A/rot/Global_Model_R_BG_2025_WD25_refit.rot).

Usage:  py -I build_rot_gwg20.py            (builds the two bases)
"""
import io
import os
import sys
import importlib
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "gmhrf"))
G = importlib.import_module("build_rot_gmhrf")     # shared engine, refit and algebra
eng = G.eng
G.FRAME_NAME = "GWG20"

S6 = os.path.join(HERE, "..", "..", "..", "external", "gaastra2022", "Supporting_Information",
                  "2021tc006772-sup-0003-table si-s02.csv")
OUTDIR = os.path.join(HERE, "rot")
BASES = {"": (__import__("os").environ["PAPER_ROOT"] + "/inputs/rotations/Global_Model_R_BG_2025.rot"),
         "WD25_": os.path.join(ROOT, "models", "WD25_WK08-A", "rot",
                               "Global_Model_R_BG_2025_WD25_refit.rot")}
REF = '@REF Gaastra_++_2022 @DOI"10.1029/2021TC006772"'
FRAME = 2


def read_s6():
    d = pd.read_csv(S6)
    d.columns = ["age", "lat", "lon", "ang", "a", "b", "c", "d", "e", "f", "kappa"]
    rows = []
    for r in d.itertuples():
        C = 1e-5 / r.kappa * np.array([[r.a, r.b, r.c], [r.b, r.d, r.e], [r.c, r.e, r.f]])
        assert np.linalg.eigvalsh(C).min() > 0, f"C not positive definite at {r.age} Ma"
        rows.append((float(r.age), r.lat, r.lon, r.ang, C, r.kappa))
    return rows


def write_version(base, out, rows):
    tmax = max(r[0] for r in rows)
    new2 = [f"002  0.0   90.0    0.0    0.0  901 !PHS-PAC GWG20 hotspots-Pacific, "
            f"inverse of Pacific rel. Pacific hotspots in Table S6 {REF}\n"]
    for (age, la, lo, ang, C, k) in rows:
        if age <= 0:
            continue
        ila, ilo, irho = G.inverse_pole(la, lo, ang)
        B = eng.quat_to_matrix(eng.quat_from_pole(ila, ilo, irho))
        new2.append(f"002 {age:6.2f}   {ila:8.4f}  {ilo:9.4f}  {irho:7.3f}  901 "
                    f"!PHS-PAC {REF} @C\"Table S6, cov/kappa_hat ({k:.3f}), h'=-B h, C'=B C B^T\"   "
                    f"{G.fmt_cov(B @ C @ B.T)}\n")
    lines, done = [], False
    for line in io.open(base, encoding="latin-1"):
        p = line.split()
        if len(p) >= 6 and p[0].isdigit() and int(p[0]) == FRAME:
            if float(p[1]) <= tmax + 1e-9:
                if not done:
                    lines.extend(new2)
                    done = True
                continue
            line = line.rstrip("\n") + "  [kept from WK08-A beyond GWG20 Table S6]\n"
        lines.append(line)
    io.open(out, "w", encoding="latin-1", newline="\n").writelines(lines)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    rows = read_s6()
    print(f"Table S6: {len(rows)} rows, {rows[0][0]:.0f}-{rows[-1][0]:.0f} Ma; "
          f"kappa_hat {min(r[5] for r in rows):.2f}-{max(r[5] for r in rows):.2f}")
    print("Maximum standard deviation of the frame (sqrt of max eigenvalue, deg), GWG20/kappa vs WK08-A:")
    wk = pd.read_csv(os.path.join(ROOT, "models", "Q22_WK08-A", "data", "wk08a_covariances.csv"), comment="#")
    for (age, la, lo, ang, C, k) in rows:
        if age % 10:
            continue
        r = wk.iloc[(wk.age_ma - age).abs().argmin()]
        Cw = 1e-5 * np.array([[r.a, r.b, r.d], [r.b, r.c, r.e], [r.d, r.e, r.f]])
        print(f"   {age:4.0f} Ma: GWG20 {np.degrees(np.sqrt(np.linalg.eigvalsh(C).max())):5.2f}"
              f" (without kappa {np.degrees(np.sqrt(k * np.linalg.eigvalsh(C).max())):5.2f})"
              f" | WK08-A @ {r.age_ma:6.2f} Ma {np.degrees(np.sqrt(np.linalg.eigvalsh(Cw).max())):5.2f}")
    for tag, base in BASES.items():
        out_a = os.path.join(OUTDIR, f"Global_Model_R_BG_2025_{tag}GWG20.rot")
        out_r = os.path.join(OUTDIR, f"Global_Model_R_BG_2025_{tag}GWG20_refit.rot")
        write_version(base, out_a, rows)
        print(f"\n[{tag or 'Quiero'}] written {out_a}")
        for s in G.write_refit(eng.RotModel(base), out_a, out_r):
            print(s)
        print(f"[{tag or 'Quiero'}] written {out_r}")
        Mb, Mr = eng.RotModel(base), eng.RotModel(out_r)
        taus = (2, 5, 9, 15, 20, 27, 30, 40, 50, 60)
        print("    shift of the present-day track relative to its base (km):" +
              "".join(f"{t:>6d}" for t in taus))
        for plate in (13, 14, 15):
            d = [G.gc_deg(G.present_position(Mb, plate, t), G.present_position(Mr, plate, t))
                 * np.pi / 180 * 6371.0 for t in taus]
            print(f"    {G.D1_NAME[plate]:15s}" + " " * 42 + "".join(f"{x:6.0f}" for x in d))
        G.jump_scan({f"{tag}GWG20'": out_r})


if __name__ == "__main__":
    main()
