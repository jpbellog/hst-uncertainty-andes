# -*- coding: utf-8 -*-
"""Inverse isochron of the Easter samples that the
model uses as isochrons (DRFT-115-2, DRFT-126-1), rebuilt from the ERDA steps
with the York et al. (2004) regression, whose slope and intercept uncertainties
are correct (York 1969 has known errors in the
uncertainty; see Muston et al. 2023 and the 1996 work).

ArArCALC columns (offset o = "Incremental Heating"): x = 39Ar(k)/40Ar(a+r)
in o+29 (2s in o+30), y = 36Ar(a)/40Ar(a+r) in o+31 (2s in o+32), r.i. correlation in
o+35. Steps with 36Ar(a) = 0 (y = 0 without error) are left out: without an error in y the
regression cannot weight them.

40Ar*/39Ar(k) = 1 / x0, x0 = -a/b (intercept with the x axis); age with the J of the sheet and
the Steiger and Jager constants (those used by Ray); then recalibrated as age_A.
"""
import math
import os
import sys
import numpy as np
import xlrd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "external", "erda_ray2012")
LAM_SJ = 5.543e-4                          # /Myr
sys.path.insert(0, HERE)
from ray2012_erda_ages import recal      # noqa: E402  (circular import resolved: recal does not import isocrona)


def york2004(x, sx, y, sy, r, tol=1e-15):
    """York et al. (2004), eq. 13: a, b, sa, sb, cov_ab, mswd (1 sigma)."""
    wx, wy = 1 / sx ** 2, 1 / sy ** 2
    b = np.polyfit(x, y, 1)[0]
    for _ in range(500):
        al = np.sqrt(wx * wy)
        W = wx * wy / (wx + b ** 2 * wy - 2 * b * r * al)
        xb, yb = np.sum(W * x) / np.sum(W), np.sum(W * y) / np.sum(W)
        U, V = x - xb, y - yb
        beta = W * (U / wy + b * V / wx - (b * U + V) * r / al)
        bn = np.sum(W * beta * V) / np.sum(W * beta * U)
        if abs(bn - b) < tol * abs(bn):
            b = bn
            break
        b = bn
    a = yb - b * xb
    xa = xb + beta
    xm = np.sum(W * xa) / np.sum(W)
    u = xa - xm
    sb = math.sqrt(1 / np.sum(W * u ** 2))
    sa = math.sqrt(1 / np.sum(W) + xm ** 2 * sb ** 2)
    mswd = np.sum(W * (y - b * x - a) ** 2) / (len(x) - 2)
    return a, b, sa, sb, -xm * sb ** 2, mswd


def edad(R, J):
    return math.log1p(J * R) / LAM_SJ


def isocrona(erda, solo_plateau=True):
    s = xlrd.open_workbook(os.path.join(SRC, f"{erda}.xls")).sheet_by_index(0)
    o = next(c for c in range(s.ncols) if str(s.cell_value(2, c)).strip().startswith("Incremental"))
    J = float(s.cell_value(5, 51))
    P = []
    for rr in range(5, s.nrows):
        if str(s.cell_value(rr, o + 2)).strip() == "S":
            break
        v = [s.cell_value(rr, o + k) for k in (29, 30, 31, 32, 35)]
        if not all(isinstance(z, float) for z in v[:4]):
            continue
        sel = str(s.cell_value(rr, o + 2)).strip() not in ("", "0.0")
        if v[2] > 0 and v[3] > 0 and (sel or not solo_plateau):
            P.append([v[0], v[1] / 2, v[2], v[3] / 2, v[4] if isinstance(v[4], float) else 0.0])
    P = np.array(P)
    if len(P) < 3:
        return dict(erda=erda, n=len(P))
    x, sx, y, sy, r = P.T
    r = np.clip(r, -0.999, 0.999)
    a, b, sa, sb, cab, mswd = york2004(x, sx, y, sy, r)
    x0 = -a / b
    g = np.array([-1 / b, a / b ** 2])              # d x0 / d(a, b)
    sx0 = math.sqrt(g @ np.array([[sa ** 2, cab], [cab, sb ** 2]]) @ g)
    R, sR = 1 / x0, sx0 / x0 ** 2
    t = edad(R, J)
    st = (edad(R + sR, J) - edad(R - sR, J)) / 2
    ta = recal(t)
    return dict(erda=erda, n=len(P), R=R, t_SJ=t, s2_SJ=2 * st, t_A=ta, s2_A=2 * st * ta / t,
                mswd=mswd, trapped=1 / a)


if __name__ == "__main__":
    for e in (1319, 1323):
        for sp in (True, False):
            print(e, "plateau" if sp else "all", isocrona(e, sp))
