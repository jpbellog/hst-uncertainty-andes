# -*- coding: utf-8 -*-
"""Model sets used by the figures.

The main text shows nine models (circuits of Quiero et al. 2022, Wilson and DeMets 2025 and
Somoza and Ghidella 2012, each with the frames WK08-A, GMHRF-B and GWG20), all weighted equally,
as a common envelope. GMHRF-A' with its three circuits goes to the supplementary material.

PAPER_MODELOS selects the set: "texto" (default, main text), "supl" (GMHRF-A') or "todos" (all twelve).
SUFIJO_SALIDA is appended to the name of each figure ("" for the main text).
"""
import os

# (folder, circuit, frame)
TEXTO = [
    ("Q22_WK08-A", "Quiero", "WK08-A"),
    ("Q22_GMHRF-B", "Quiero", "GMHRF-B"),
    ("Q22_GWG20", "Quiero", "GWG20"),
    ("WD25_WK08-A", "WD25", "WK08-A"),
    ("WD25_GMHRF-B", "WD25", "GMHRF-B"),
    ("WD25_GWG20", "WD25", "GWG20"),
    ("SG12_WK08-A", "SG12", "WK08-A"),
    ("SG12_GMHRF-B", "SG12", "GMHRF-B"),
    ("SG12_GWG20", "SG12", "GWG20"),
]
SUPL = [
    ("Q22_GMHRF-A", "Quiero", "GMHRF-A'"),
    ("WD25_GMHRF-A", "WD25", "GMHRF-A'"),
    ("SG12_GMHRF-A", "SG12", "GMHRF-A'"),
]
CONJUNTOS = {"texto": (TEXTO, ""), "supl": (SUPL, "_supl_gmhrfA"), "todos": (TEXTO + SUPL, "_doce")}

NOMBRE = os.environ.get("PAPER_MODELOS", "texto")
MODELOS, SUFIJO_SALIDA = CONJUNTOS[NOMBRE]
CARPETAS = [m[0] for m in MODELOS]


def sufijo_bandas(carpeta):
    """Suffix of the bands at the present margin (tables/make_present_bands.py, PAPER_BANDS_SUFFIX)."""
    return "_" + carpeta


ALCANCE_DIR = os.path.join(__import__("os").environ["PAPER_ROOT"],
                           "analysis", "outputs", "frame_coverage")
LAT_COLS = ["lat_nominal", "lat_p2_5", "lat_p16", "lat_p50", "lat_p84", "lat_p97_5"]


def recortar(df, carpeta):
    """Each model is valid up to the age of its oldest frame covariance; a collision at age t
    is valid if the colliding segment formed at tau_c <= T_cov (analysis/frame_coverage.py).
    Sets the latitudes to NaN in the (ridge, pf, age) rows that are not valid or not
    evaluated. Returns a copy."""
    import pandas as pd
    a = pd.read_csv(os.path.join(ALCANCE_DIR, f"frame_coverage_{carpeta}.csv"), dtype={"pf": str})
    a["pf"] = a.pf.astype(float).map(lambda x: f"{x:g}")
    ok = set(map(tuple, a.loc[a.valido, ["ridge", "pf", "age"]].values))
    df = df.copy()
    pf = df.pf.astype(float).map(lambda x: f"{x:g}")
    m = [(r, p, float(t)) in ok for r, p, t in zip(df.ridge, pf, df.age)]
    cols = [c for c in LAT_COLS if c in df.columns]
    df.loc[~pd.Series(m, index=df.index), cols] = float("nan")
    return df


def describe():
    """Legend text: how many models and from which set."""
    return f"{len(MODELOS)} models"
