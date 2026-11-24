# -*- coding: utf-8 -*-
"""Porphyries vs background test (step 12: one-sided Mann-Whitney, fraction <= 400 km) in the nine
models of the text (g13) and in the following variants:
  - Arriagada: final assignment of step 07 of each folder (restoration of Arriagada et al. 2008);
  - Utrecht: positions and trench of the blocks of Schepers (2017) / Boschman (2021) (run_shortening_step07.sh);
  - "Utrecht, fosa propia" (own trench): Utrecht positions with our trench;
  - WVH25: frame of Wagenaar et al. (2025) with the Quiero circuit (one model, build_rot_wagenaar.py).
Also counts how many porphyries change chain relative to Arriagada in the same model.
Usage:  py sens9_deposits.py     Output: outputs/sens9/test_porphyries_variants.csv
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "models", "Q22_WK08-A"))
sys.path.insert(0, os.path.join(ROOT, "figures"))
s12 = importlib.import_module("12_deposit_track_test")
import paper_models as MP  # noqa: E402

U = os.path.join(HERE, "outputs", "andean_shortening")
W = os.path.join(HERE, "outputs", "wagenaar_frame", "Q22_WVH25")
OUT = os.path.join(HERE, "outputs", "sens9")
os.makedirs(OUT, exist_ok=True)


def main():
    exc = s12.excluded_methods()
    var = {"Arriagada": lambda d: os.path.join(ROOT, "models", d),
           "Utrecht": lambda d: os.path.join(U, d),
           "Utrecht, fosa propia": lambda d: os.path.join(U, d + "_own_trench")}
    rows, cambio = [], []
    for d in MP.CARPETAS:
        base = None
        for v, f in var.items():
            if not os.path.exists(os.path.join(f(d), "deposit_tracks_porphyries.csv")):
                continue
            D = s12.cargar(f(d), exc)
            for pid, cad in list(s12.CAD.items()) + [(0, "Todas")]:
                sel = (lambda t: t) if pid == 0 else (lambda t, p=pid: t[t.hs == p])
                rows.append(dict(variante=v, modelo=d, cadena=cad,
                                 **s12.tests(sel(D["porphyries"]), sel(D["fondo"]))))
            por = pd.read_csv(os.path.join(f(d), "deposit_tracks_porphyries.csv"))
            if v == "Arriagada":
                base = por
            elif base is not None:
                cambio.append(dict(variante=v, modelo=d, cambian=int((por.HS4.values != base.HS4.values).sum()),
                                   n=len(por)))
    if os.path.exists(os.path.join(W, "deposit_tracks_porphyries.csv")):
        D = s12.cargar(W, exc)
        for pid, cad in list(s12.CAD.items()) + [(0, "Todas")]:
            sel = (lambda t: t) if pid == 0 else (lambda t, p=pid: t[t.hs == p])
            rows.append(dict(variante="WVH25", modelo="WVH25", cadena=cad,
                             **s12.tests(sel(D["porphyries"]), sel(D["fondo"]))))
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(OUT, "test_porphyries_variants.csv"), index=False, float_format="%.4g")

    def rango(s):
        s = s.dropna()
        return f"{s.min():.2g}-{s.max():.2g}" if len(s) else ""
    R = T.groupby(["variante", "cadena"], sort=False).agg(
        modelos=("modelo", "size"), p_mw=("mw_p_menor", rango),
        p05=("mw_p_menor", lambda s: int((s < 0.05).sum())),
        f400_por=("f400_dep", rango), f400_bg=("f400_bg", rango),
        n_por=("n_dep", rango), n_bg=("n_bg", rango)).reset_index()
    with pd.option_context("display.width", 200):
        print(R.to_string(index=False))
    if cambio:
        C = pd.DataFrame(cambio)
        print("\nporphyries that change chain relative to Arriagada (range across models):")
        print(C.groupby("variante").cambian.agg(["min", "max"]).to_string())


if __name__ == "__main__":
    main()
