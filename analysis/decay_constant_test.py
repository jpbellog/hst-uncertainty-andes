# -*- coding: utf-8 -*-
"""Effect of the uncertainty of the total 40K decay constant on the 40Ar/39Ar ages of the four
chains, compared with the terms that already enter the Monte Carlo.

Age of a sample u calibrated against a monitor s of age t_s (same irradiation):
    t_u = ln(1 + R (exp(lam t_s) - 1)) / lam,  R fixed by the measurement.
With t_s fixed (FCs 28.201 Ma, Kuiper et al. 2008), lam is perturbed and d t_u is measured.
Values (from ArAr_recalibration.py, read in the original publications):
  lam = 5.463e-10 /a +/- 0.214e-10 (2 sigma; Min et al. 2000 according to Kuiper et al. 2008, p. 3)
  FCs = 28.201 +/- 0.046 Ma (2 sigma; Kuiper et al. 2008)
It is compared with the FCs systematic already used by 03b (meta.systematics.FCs of hotspot_fit.json,
1 sigma relative) and with the geological scatter w_t of each chain (Table S1).
Note: 03b treats the FCs systematic as a relative scale of the age; here the exact sensitivity
is computed with the age equation.
Usage:  py decay_constant_test.py
"""
import json
import os
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
P = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LAM, LAM_2S = 5.463e-10, 0.214e-10
TS = 28.201e6


def edad(R, lam, ts=TS):
    return np.log1p(R * np.expm1(lam * ts)) / lam


def R_de(t, lam=LAM, ts=TS):
    return np.expm1(lam * t) / np.expm1(lam * ts)


meta = json.load(open(os.path.join(P, "models", "Q22_WK08-A", "hotspot_fit.json"), encoding="utf-8"))
s_fcs = meta["meta"]["systematics"]["FCs"]
lam1 = LAM_2S / 2
print(f"lambda 1 sigma relative: {lam1 / LAM * 100:.2f} %   FCs (03b) 1 sigma relative: {s_fcs * 100:.3f} %")
print()
print("age (Ma)  | d t_u per 1 sigma of lambda (kyr) | relative (%) | FCs 1s (kyr) | combined / FCs")
for t in (0.4, 1, 3, 5, 10, 20, 30):
    R = R_de(t * 1e6)
    d = (edad(R, LAM + lam1) - edad(R, LAM - lam1)) / 2
    rel = d / (t * 1e6)
    f = s_fcs * t * 1e6
    print(f"{t:9.1f} | {d / 1e3:31.2f} | {rel * 100:12.4f} | {f / 1e3:12.2f} | {np.hypot(d, f) / f:15.3f}")
print()
# compared with the geological scatter of each chain (model Q22 + WK08-A and range of the nine)
for k, lab in (("Easter", "Easter"), ("JFR", "Juan Fernandez"), ("Foundation", "Foundation")):
    x = meta["tracks"][k]
    tmax = x["t_max"]
    R = R_de(tmax * 1e6)
    d = (edad(R, LAM + lam1) - edad(R, LAM - lam1)) / 2
    print(f"{lab}: oldest sample {tmax:.1f} Ma -> {d / 1e3:.1f} kyr from lambda; "
          f"w_t = {x['w_t']:.2f} Myr ({d / (x['w_t'] * 1e6) * 100:.2f} % of w_t)")
