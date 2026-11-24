# -*- coding: utf-8 -*-
"""Summary of the run_sens9.sh tests in the nine models of the text.

(a) Hotspot-eruption lag and (b) fit without isochrons: nominal collision (step 04 unperturbed,
f_PAC 0.44, same construction as nominal_sensitivity.py) with the .rot of each fit;
differences lag1 - lag0, lag2 - lag0 (all three without the 0 Ma observation) and noYORK2 - final.
They are compared with the 68 % half-width of the total g13 band of the same model, (p84 - p16)/2, only at
the valid ages of the model (paper_models.recortar).
(c) Interpolation: 68 % half-width without the renormalization (_s9_interp0_tot044) compared with the g13 band (with
renormalization), ratio norm / no norm.
Usage:  py sens9_summary.py      Outputs: outputs/sens9/*.csv and an on-screen summary.
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(HERE, "outputs", "sens9")
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.join(ROOT, "figures"))
import paper_models as MP  # noqa: E402

ROT = {"Q22_WK08-A": "Global_Model_R_BG_2025_hsb.rot",
       "Q22_GMHRF-B": "Global_Model_R_BG_2025_GMHRF_gts12_hsbgmhrf.rot",
       "Q22_GWG20": "Global_Model_R_BG_2025_GWG20_hsb.rot",
       "WD25_WK08-A": "Global_Model_R_BG_2025_WD25_hsb.rot",
       "WD25_GMHRF-B": "Global_Model_R_BG_2025_WD25_GMHRF_gts12_hsbgmhrf.rot",
       "WD25_GWG20": "Global_Model_R_BG_2025_WD25_GWG20_hsb.rot",
       "SG12_WK08-A": "Global_Model_R_BG_2025_SG12_hsb.rot",
       "SG12_GMHRF-B": "Global_Model_R_BG_2025_SG12_GMHRF_gts12_hsbgmhrf.rot",
       "SG12_GWG20": "Global_Model_R_BG_2025_SG12_GWG20_hsb.rot"}
VAR = ("final", "_s9_lag0", "_s9_lag1", "_s9_lag2", "_s9_noYORK2")
KM = 111.195
TRAMOS = [-1, 27, 45, 66]


def nominales(d):
    """Nominal collision latitude (f_PAC 0.44) per variant, chain and age."""
    D = os.path.join(ROOT, "models", d)
    sys.path.insert(0, D)
    for k in [k for k in sys.modules if k in ("04_mc_collisions", "02_rot_engine")]:
        del sys.modules[k]
    cwd = os.getcwd()
    os.chdir(D)
    mc = importlib.import_module("04_mc_collisions")
    mdf, nn = mc.C.load_full_margin(mc.MARGIN_CSV)
    ms = {t: mc.C.margin_segments(mdf, nn, t) for t in mc.AGES_T}
    rows = []
    for v in VAR:
        fn = os.path.join(D, "rot", ROT[d]) if v == "final" else os.path.join(D, "rot", "sens9", v + ".rot")
        tot = mc.make_tot(mc.eng.RotModel(fn), None)
        for ridge, (seed, pid) in mc.RIDGES.items():
            P = mc.track_points(tot, seed, pid, 902003)
            for t in mc.AGES_T:
                la, lo = mc.project(tot, 902003, t, P)
                r = mc.C.collide_strict(ms[t], mc.segs_for_age(la, lo, t))
                rows.append(dict(modelo=d, var=v, ridge=ridge, age=t, lat=np.nan if r is None else r[0]))
    os.chdir(cwd)
    sys.path.remove(D)
    return pd.DataFrame(rows)


def bandas(d, nombre):
    B = pd.read_csv(os.path.join(ROOT, "models", d, nombre))
    B = B[np.isclose(B.pf.astype(float), 0.44)]
    B = MP.recortar(B, d)
    B["s68_km"] = (B.lat_p84 - B.lat_p16) / 2 * KM
    return B[["ridge", "age", "s68_km"]]


def main():
    N, I = [], []
    for d in ROT:
        n = nominales(d).pivot_table(index=["modelo", "ridge", "age"], columns="var", values="lat").reset_index()
        n = n.merge(bandas(d, "mc_collision_bands_total_hsb_cc_g13.csv"), on=["ridge", "age"], how="left")
        N.append(n)
        f0 = os.path.join(ROOT, "models", d, "mc_collision_bands_total_s9_interp0_tot044.csv")
        if os.path.exists(f0):
            i = bandas(d, "mc_collision_bands_total_hsb_cc_g13.csv").merge(
                bandas(d, "mc_collision_bands_total_s9_interp0_tot044.csv"), on=["ridge", "age"],
                suffixes=("_norm", "_sin"))
            i["modelo"] = d
            I.append(i)
    N = pd.concat(N)
    N = N[N.s68_km.notna()]
    for a, b, lab in (("_s9_lag1", "_s9_lag0", "lag1"), ("_s9_lag2", "_s9_lag0", "lag2"),
                      ("_s9_noYORK2", "final", "noYORK2")):
        N[f"d_{lab}_km"] = (N[a] - N[b]) * KM
    N["tramo"] = pd.cut(N.age, TRAMOS, labels=["0-27", "27-45", "45-66"])
    N.to_csv(os.path.join(OUT, "nominal_lag_york_9models.csv"), index=False, float_format="%.4f")

    print("== (a) hotspot-eruption lag and (b) without YORK2: max |d| (km) per chain and segment, range "
          "across models; ages with |d| > s68 (sum over the nine models / valid ages)")
    for lab in ("lag1", "lag2", "noYORK2"):
        c = f"d_{lab}_km"
        g = N.assign(ad=N[c].abs(), sobre=N[c].abs() > N.s68_km)
        m = g.groupby(["ridge", "tramo", "modelo"], observed=True).ad.max().groupby(["ridge", "tramo"], observed=True)
        s = g.groupby(["ridge", "tramo"], observed=True).agg(sobre=("sobre", "sum"), n=("sobre", "size"),
                                                           s68=("s68_km", "median"))
        s["max_min"], s["max_max"] = m.min().round(0), m.max().round(0)
        print(f"\n-- {lab}\n{s.to_string()}")
        print(f"   total ages above s68: {int(g.sobre.sum())} of {len(g)}; per chain "
              f"{g.groupby('ridge').sobre.sum().to_dict()}")
    if I:
        I = pd.concat(I).dropna()
        I["cociente"] = I.s68_km_norm / I.s68_km_sin
        I["tramo"] = pd.cut(I.age, TRAMOS, labels=["0-27", "27-45", "45-66"])
        I.to_csv(os.path.join(OUT, "interp_9models.csv"), index=False, float_format="%.4f")
        print("\n== (c) interpolation: ratio s68 norm / no norm (median per segment; range across models)")
        q = I.groupby(["ridge", "tramo", "modelo"], observed=True).cociente.median()
        print(q.groupby(["ridge", "tramo"], observed=True).agg(["min", "median", "max"]).round(2).to_string())
        print(f"   all: median {I.cociente.median():.2f}, p5-p95 {I.cociente.quantile(.05):.2f}-"
              f"{I.cociente.quantile(.95):.2f}, max {I.cociente.max():.2f}; max increase "
              f"{(I.s68_km_norm - I.s68_km_sin).max():.0f} km")


if __name__ == "__main__":
    main()
