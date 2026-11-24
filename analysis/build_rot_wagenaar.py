# -*- coding: utf-8 -*-
"""Rotation model with the "minimum continental motion" frame of Wagenaar, Vaes and
van Hinsbergen (2025, JGR Solid Earth 130, e2024JB030430) in place of WK08-A, to see its effect on
the collisions and on the association with deposits.

Source: Wagenaar (2025), Zenodo doi:10.5281/zenodo.17536784, V3.0,
Rotation_files/Continent_frame701.rot = Africa (701) relative to the mantle (000), every 10 Myr
(placed in ../external/wagenaar2025/). mean_poles.csv gives
the Kent ellipse of each pole (not a covariance): the frame enters WITHOUT covariance.

Warning from the Zenodo README: those poles were computed with the relative model of Merdith
et al. (2021) and the authors advise against using them with another one. Here they are combined
with our circuit (Pacific -> Antarctica -> Nubia of the current .rot), so the resulting frame is NOT
exactly that of Wagenaar et al.: it is their motion of Africa relative to the mantle carried to the
Pacific through our circuit. It is a sensitivity test.

Construction: as in build_rot_gwg20.py, only through plate 002 (hotspots relative to the
Pacific, fixed 901): R_901/000(t) = W_701/000(t) * R_901/701(t) (R_901/701 from the base .rot) and
002/901 = inverse of R_901/000. Nodes at 10, 20, ..., 150 Ma (those of the frame); between nodes the
engine interpolates with slerp. Base: models/Q22_WK08-A/rot/Global_Model_R_BG_2025_hsb.rot (circuit of
Quiero et al. 2022, GTS2020, Granot et al. 2013). Plates 13-16 are refitted afterwards with
step 03b in the test folder.

Usage:  py build_rot_wagenaar.py
Output: outputs/wagenaar_frame/Q22_WVH25/ (copy of the scripts and data of models/Q22_WK08-A/ and
        rot/Global_Model_R_BG_2025_WVH25_hsb.rot)
"""
import io
import os
import sys
import glob
import shutil
import importlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "models", "Q22_WK08-A"))
eng = importlib.import_module("02_rot_engine")

W701 = os.path.join(HERE, "..", "external", "wagenaar2025",
                    "Rotation_files", "Continent_frame701.rot")
BASE = os.path.join(ROOT, "models", "Q22_WK08-A", "rot", "Global_Model_R_BG_2025_hsb.rot")
OUT = os.path.join(HERE, "outputs", "wagenaar_frame", "Q22_WVH25")
ROT_OUT = os.path.join(OUT, "rot", "Global_Model_R_BG_2025_WVH25_hsb.rot")
TMAX = 150.0
REF = '@REF Wagenaar_++_2025_JGR @DOI"10.1029/2024JB030430" (Zenodo 10.5281/zenodo.17536784)'


def read_w701():
    rows = []
    for line in io.open(W701, encoding="latin-1"):
        p = line.split()
        if len(p) >= 6 and p[0] == "701" and 0 < float(p[1]) <= TMAX:
            rows.append(tuple(float(x) for x in p[1:5]))
    return rows


def pole(q):
    q = q / np.linalg.norm(q)
    if q[0] < 0:
        q = -q
    ang = 2 * np.degrees(np.arccos(np.clip(q[0], -1, 1)))
    s = np.linalg.norm(q[1:])
    if s < 1e-15:
        return 90.0, 0.0, 0.0
    la, lo = eng.xyz_to_latlon(q[1:] / s)
    return la, lo, ang


def main():
    os.makedirs(os.path.join(OUT, "rot"), exist_ok=True)
    # the pipeline scripts look for rotation_models/ next to their parent folder; in the copy the
    # project root is fixed so that they read the same inputs (deposit database v7, margin, geochemistry)
    rel = "os.path.dirname(os.path.dirname(os.path.abspath(__file__)))"
    for f in glob.glob(os.path.join(ROOT, "models", "Q22_WK08-A", "*.py")):
        s = open(f, encoding="utf-8").read().replace(rel, repr(ROOT))
        open(os.path.join(OUT, os.path.basename(f)), "w", encoding="utf-8").write(s)
    shutil.copytree(os.path.join(ROOT, "models", "Q22_WK08-A", "data"), os.path.join(OUT, "data"),
                    dirs_exist_ok=True)
    M = eng.RotModel(BASE)
    new2 = [f"002  0.0   90.0    0.0    0.0  901 !PHS-PAC WVH25 test: mantle (Wagenaar continent "
            f"frame) relative to Pacific through this circuit {REF}\n"]
    print("  age  pole 901/000 (this circuit + Wagenaar)  |  002/901 new")
    for t, la, lo, ang in read_w701():
        W = eng.quat_from_pole(la, lo, ang)
        q901_701 = M.rel_rotation(901, t, 701)
        q = eng.quat_mult(W, q901_701)                  # 901 -> mantle
        ila, ilo, iang = pole(eng.quat_conj(q))         # 002 (mantle) -> 901
        new2.append(f"002 {t:6.1f}   {ila:9.4f}  {ilo:10.4f}  {iang:8.4f}  901 !PHS-PAC WVH25 "
                    f"test {REF} [no covariance; Kent ellipses only in mean_poles.csv]\n")
        pla, plo, pang = pole(q)
        print(f"{t:5.0f}  {pla:7.2f} {plo:8.2f} {pang:7.2f}  |  {ila:7.2f} {ilo:8.2f} {iang:7.2f}")
    lines, done = [], False
    for line in io.open(BASE, encoding="latin-1"):
        p = line.split()
        if len(p) >= 6 and p[0].isdigit() and int(p[0]) == 2:
            if float(p[1]) <= TMAX + 1e-9:
                if not done:
                    lines.extend(new2)
                    done = True
                continue
            line = line.rstrip("\n") + "  [kept from WK08-A beyond the WVH25 test]\n"
        lines.append(line)
    io.open(ROT_OUT, "w", encoding="latin-1", newline="\n").writelines(lines)
    print(f"written {ROT_OUT}")


if __name__ == "__main__":
    main()
