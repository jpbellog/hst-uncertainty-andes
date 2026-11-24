# -*- coding: utf-8 -*-
"""STEP 12 -- Quantitative test of the deposit-track association by hotspot and by model.
Reads the assignment of step 07 (07_deposit_tracks.py) from each of the
model folders (../<folder>/deposit_tracks_*.csv).

For each model and each chain (Easter, Juan Fernandez, San Felix, Foundation) it compares the
deposits assigned to that chain with the magmatic background assigned to the same chain
("background D": age-quality filter v3, 0-63.5 Ma, no latitude
window):
  (1) distance to the track: two-sample KS and one-sided Mann-Whitney (deposits closer);
  (2) side of the track: fraction to the S, one-sided Fisher exact test (deposits more to the S).
Also the pooled set (all chains) and the fraction at <= 400 km (distance threshold).

Figure: one panel per chain with the CDF of the distance to the track for porphyries, epithermal
deposits and background; line = median of the eight models, band = range across models.

Usage:  py 12_deposit_track_test.py   (from models/Q22_WK08-A/, after step 07 in the twelve folders)
Outputs (in models/Q22_WK08-A/): deposit_track_test_by_model.csv, deposit_track_test_resumen.csv,
deposit_track_cdf_por_hotspot.png
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = HERE
MODELOS = {"Q22_WK08-A": "WK08", "WD25_WK08-A": "WD25", "Q22_GMHRF-A": "GMHRF-A'",
           "Q22_GMHRF-B": "GMHRF-B",
           "WD25_GMHRF-B": "WD25+GMHRF-B", "Q22_GWG20": "GWG20",
           "WD25_GWG20": "WD25+GWG20", "SG12_WK08-A": "SG12",
           # full grid of 3 circuits x 4 frames
           "WD25_GMHRF-A": "WD25+GMHRF-A'", "SG12_GMHRF-A": "SG12+GMHRF-A'",
           "SG12_GMHRF-B": "SG12+GMHRF-B", "SG12_GWG20": "SG12+GWG20"}
CAD = {14: "Easter", 13: "Juan Fernandez", 15: "San Felix", 16: "Foundation"}
UMBRAL = 400.0
GRID = np.arange(0, 1201, 5.0)


def excluded_methods():
    return importlib.import_module("07_deposit_tracks").EXCLUDED


def cargar(d, exc):
    por = pd.read_csv(os.path.join(ROOT, "models", d, "deposit_tracks_porphyries.csv"))
    epi = pd.read_csv(os.path.join(ROOT, "models", d, "deposit_tracks_epithermal.csv"))
    g = pd.read_csv(os.path.join(ROOT, "models", d, "deposit_tracks_samples.csv"), low_memory=False)
    age = pd.to_numeric(g["Age (Ma)"], errors="coerce")
    meth = g["Method"].astype(str).str.strip().str.lower()
    bad = meth.isin(exc) | g["Method"].isna() | (meth == "") | (meth == "nan")
    bg = g[((~bad) | (age < 1)) & (age >= 0) & (age <= 63.5)]
    out = {}
    for k, t in (("porphyries", por), ("epithermal", epi), ("fondo", bg)):
        t = t[t.HS4.notna()]
        out[k] = pd.DataFrame({"hs": t.HS4.astype(int).values, "d": t.Dist4_km.values,
                               "s": (t.N_S4 == 0).values})
    return out


def tests(dep, bg):
    r = dict(n_dep=len(dep), n_bg=len(bg))
    if len(dep) < 3 or len(bg) < 3:
        return r
    r["d_med_dep"], r["d_med_bg"] = float(np.median(dep.d)), float(np.median(bg.d))
    r["ks_D"], r["ks_p"] = [float(x) for x in stats.ks_2samp(dep.d, bg.d)]
    r["mw_p_menor"] = float(stats.mannwhitneyu(dep.d, bg.d, alternative="less").pvalue)
    r["f400_dep"], r["f400_bg"] = float(np.mean(dep.d <= UMBRAL)), float(np.mean(bg.d <= UMBRAL))
    r["fS_dep"], r["fS_bg"] = float(dep.s.mean()), float(bg.s.mean())
    tab = [[int(dep.s.sum()), int((~dep.s).sum())], [int(bg.s.sum()), int((~bg.s).sum())]]
    r["fisher_p_masS"] = float(stats.fisher_exact(tab, alternative="greater").pvalue)
    return r


def main():
    exc = excluded_methods()
    rows, cdf = [], {}
    for d, nom in MODELOS.items():
        D = cargar(d, exc)
        for tipo in ("porphyries", "epithermal"):
            for pid, cad in list(CAD.items()) + [(0, "Todas")]:
                sel = (lambda t: t) if pid == 0 else (lambda t, p=pid: t[t.hs == p])
                rows.append(dict(modelo=nom, tipo=tipo, cadena=cad,
                                 **tests(sel(D[tipo]), sel(D["fondo"]))))
        for tipo in ("porphyries", "epithermal", "fondo"):
            for pid, cad in CAD.items():
                x = D[tipo].d[D[tipo].hs == pid].values
                cdf[(nom, tipo, cad)] = (np.searchsorted(np.sort(x), GRID, side="right") / len(x)
                                         if len(x) else np.full(GRID.shape, np.nan))
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(OUT, "deposit_track_test_by_model.csv"), index=False, float_format="%.4g")

    # summary across models
    def rango(s):
        s = s.dropna()
        return f"{s.min():.3g}-{s.max():.3g}" if len(s) else ""
    R = T.groupby(["tipo", "cadena"], sort=False).agg(
        n_dep=("n_dep", rango), n_bg=("n_bg", rango), d_med_dep=("d_med_dep", rango),
        d_med_bg=("d_med_bg", rango), mw_p_menor=("mw_p_menor", rango), ks_p=("ks_p", rango),
        modelos_mw_p05=("mw_p_menor", lambda s: int((s < 0.05).sum())),
        fS_dep=("fS_dep", rango), fS_bg=("fS_bg", rango), fisher_p_masS=("fisher_p_masS", rango),
        modelos_fisher_p05=("fisher_p_masS", lambda s: int((s < 0.05).sum())),
        f400_dep=("f400_dep", rango), f400_bg=("f400_bg", rango)).reset_index()
    R.to_csv(os.path.join(OUT, "deposit_track_test_resumen.csv"), index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(R.to_string(index=False))

    # figure: one panel per chain
    modelos = list(MODELOS.values())
    colores = {"porphyries": "#b2182b", "epithermal": "#2166ac", "fondo": "0.35"}
    etiq = {"porphyries": "Porphyry Cu", "epithermal": "Epithermal", "fondo": "Magmatic background"}
    fig, axs = plt.subplots(1, 4, figsize=(17, 4.6), sharey=True)
    for ax, cad in zip(axs, CAD.values()):
        for tipo in ("fondo", "epithermal", "porphyries"):
            M = np.array([cdf[(m, tipo, cad)] for m in modelos])
            if np.all(np.isnan(M)):
                continue
            n = T[(T.tipo == ("porphyries" if tipo == "fondo" else tipo)) & (T.cadena == cad)]
            nn = n.n_bg if tipo == "fondo" else n.n_dep
            ax.fill_between(GRID, np.nanmin(M, 0), np.nanmax(M, 0), color=colores[tipo], alpha=0.18, lw=0)
            ax.plot(GRID, np.nanmedian(M, 0), color=colores[tipo], lw=2,
                    label=f"{etiq[tipo]} (n = {int(nn.min())}-{int(nn.max())})")
        p = T[(T.tipo == "porphyries") & (T.cadena == cad)].mw_p_menor.dropna()
        ax.axvline(UMBRAL, color="0.5", ls=":", lw=1)
        ax.set_title(f"{cad}\nporphyry vs background, one-sided MW p = "
                     + (f"{p.min():.2g}-{p.max():.2g}" if len(p) else "n/a"), fontsize=10)
        ax.set_xlabel("Distance to nearest track (km)")
        ax.set_xlim(0, 1200)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="lower right")
    axs[0].set_ylabel("Cumulative fraction")
    fig.suptitle("Distance to the nearest of four hotspot tracks; line = median of eight models, "
                 "band = range across models", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "deposit_track_cdf_por_hotspot.png"), dpi=140)
    print(f"\nwritten in {OUT}")


if __name__ == "__main__":
    main()
