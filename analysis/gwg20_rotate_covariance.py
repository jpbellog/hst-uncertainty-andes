# -*- coding: utf-8 -*-
"""GWG20 covariance with the perturbation on the left of PAC->HS (Andrews et al. 2006). For
the 002->901 row (B = inverse of Table S6) the engine applies exp(M(h)) B, so C' = B C B^T.

Rewrites only the six covariance numbers of the Gaastra 002 rows in the .rot files of
Q22_GWG20 and WD25_GWG20; ages, poles, angles and the rest of the
file stay unchanged. The original files are backed up in <pipeline>/rot/<BK>/.
Refuses to run if the backup already exists (the script is not idempotent).

Usage:  py -I gwg20_rotate_covariance.py
"""
import io
import os
import re
import shutil
import importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
PIPES = ("Q22_GWG20", "WD25_GWG20")
BK = "backup_pre_gwg20_cov"
spec = importlib.util.spec_from_file_location(
    "eng", os.path.join(ROOT, "models", "Q22_GWG20", "02_rot_engine.py"))
eng = importlib.util.module_from_spec(spec)
spec.loader.exec_module(eng)
TAG_OLD = "h'=-h, C'=C"
TAG_NEW = "h'=-B h, C'=B C B^T"


def fmt_cov(C):
    return "   ".join(f"{x:.4E}" for x in
                      (C[0, 0], C[0, 1], C[0, 2], C[1, 1], C[1, 2], C[2, 2]))


def patch(path):
    out, n = [], 0
    for line in io.open(path, encoding="latin-1"):
        p = line.split()
        if p and p[0] == "002" and "Gaastra" in line and TAG_OLD in line:
            lat, lon, ang = map(float, p[2:5])
            m = re.search(r"((?:\s+-?\d\.\d+E[+-]\d+){6})\s*$", line)
            v = [float(x) for x in m.group(1).split()]
            C = np.array([[v[0], v[1], v[2]], [v[1], v[3], v[4]], [v[2], v[4], v[5]]])
            B = eng.quat_to_matrix(eng.quat_from_pole(lat, lon, ang))
            Cn = B @ C @ B.T
            assert np.linalg.eigvalsh(Cn).min() > 0
            assert abs(np.trace(Cn) - np.trace(C)) < 1e-9 * max(1.0, np.trace(C)) + 1e-12
            line = line[:m.start(1)].replace(TAG_OLD, TAG_NEW) + "   " + fmt_cov(Cn) + "\n"
            n += 1
        out.append(line)
    io.open(path, "w", encoding="latin-1", newline="\n").writelines(out)
    return n


def main():
    for pipe in PIPES:
        rd = os.path.join(ROOT, "models", pipe, "rot")
        bk = os.path.join(rd, BK)
        if os.path.exists(bk):
            raise SystemExit(f"{bk} already exists: not rotating again")
        os.makedirs(bk)
        for f in sorted(os.listdir(rd)):
            if f.endswith(".rot"):
                shutil.copy2(os.path.join(rd, f), bk)
                print(f"{pipe}/{f}: {patch(os.path.join(rd, f))} 002 rows rewritten")


if __name__ == "__main__":
    main()
