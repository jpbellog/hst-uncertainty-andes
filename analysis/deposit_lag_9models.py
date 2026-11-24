# -*- coding: utf-8 -*-
"""Collision-deposit lag in the nine models of the main text.

Compares the step 07 assignment without lag (deposit_tracks_<type>.csv, run g13) with
PAPER_DEP_LAG = 1, 2 and 5 Myr (deposit_tracks_<type>_lag<L>.csv). A deposit "changes chain"
if its chain within <= 400 km (or none, if it lies farther away) differs between the two runs,
as in the WK08-A figure.
Output: outputs/deposit_lag_9models.csv and a summary on screen.
Usage:  py deposit_lag_9models.py
"""
import os
import sys
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
A = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(A, ".."))
sys.path.insert(0, os.path.join(P, "figures"))
import paper_models as MP

UMBRAL = 400.0


def cadena(df):
    hs = df["HS4"].where(df["Dist4_km"] <= UMBRAL)
    return hs.fillna(-1).to_numpy()


rows = []
for d, circ, marco in MP.MODELOS:
    for tipo in ("porphyries", "epithermal"):
        b = pd.read_csv(os.path.join(P, "models", d, f"deposit_tracks_{tipo}.csv"))
        for L in (1, 2, 5):
            f = os.path.join(P, "models", d, f"deposit_tracks_{tipo}_lag{L}.csv")
            x = pd.read_csv(f)
            assert len(x) == len(b)
            ok = np.isfinite(b["Dist4_km"].to_numpy()) & np.isfinite(x["Dist4_km"].to_numpy())
            cambia = (cadena(b) != cadena(x))[ok]
            dd = np.abs(x["Dist4_km"].to_numpy() - b["Dist4_km"].to_numpy())[ok]
            rows.append(dict(modelo=f"{circ}+{marco}", tipo=tipo, lag=L, n=int(ok.sum()),
                             cambian=int(cambia.sum()),
                             f400_0=float((b["Dist4_km"][ok] <= UMBRAL).mean()),
                             f400_L=float((x["Dist4_km"][ok] <= UMBRAL).mean()),
                             dd_mediana_km=float(np.median(dd))))
R = pd.DataFrame(rows)
R.to_csv(os.path.join(A, "outputs", "deposit_lag_9models.csv"), index=False)
with pd.option_context("display.width", 200):
    print(R.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    print()
    print(R.groupby(["tipo", "lag"]).agg(n=("n", "first"), cambian_min=("cambian", "min"),
                                         cambian_max=("cambian", "max"),
                                         dd_min=("dd_mediana_km", "min"),
                                         dd_max=("dd_mediana_km", "max")).to_string())
