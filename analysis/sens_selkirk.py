# -*- coding: utf-8 -*-
"""Foundation and the Selkirk microplate: summary of run_sens_selkirk.sh in the nine
models of the text. The microplate extends from the failed rift at 124 45'W to the eastern boundary at 121 40'W
(O'Connor et al. 1998, EPSL 161, section 2).
(1) Residuals of the final Foundation fit per zone (west of the failed rift, microplate, east).
(2) Fits without the microplate samples (_s9_sinSelkirk: SO100-33DS-1, SO100-38DS-1;
    _s9_sinSelkirkB: also SO100-41DS-1) compared with the final fit: fixed/drift model, drift, w_x and
    hotspot displacement.
(3) Nominal Foundation collision (f_PAC 0.44, same construction as sens9_summary.py) with the .rot
    of each fit, compared with the final fit and with the 68 % half-width of the total g13 band, at the valid ages.
Usage:  py sens_selkirk.py     Outputs: outputs/selkirk/*.csv and an on-screen summary.
"""
import os
import sys
import json
import importlib
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(HERE, "outputs", "selkirk")
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.join(ROOT, "figures"))
sys.path.insert(0, HERE)
import paper_models as MP  # noqa: E402
import sens9_summary as S9  # noqa: E402

KM = 111.195
FR, ESTE = -124.75, -(121 + 40 / 60)
TAGS = ("_s9_sinSelkirk", "_s9_sinSelkirkB")


def gc_km(a, b):
    la1, lo1, la2, lo2 = map(np.radians, (a[0], a[1], b[0], b[1]))
    c = np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2)
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))) * KM)


def zona(lon):
    return np.where(lon < FR, "oeste", np.where(lon <= ESTE, "microplaca", "este"))


def residuos():
    R = []
    for d in S9.ROT:
        r = pd.read_csv(os.path.join(ROOT, "models", d, "hotspot_residuals_Foundation.csv"))
        R.append(r[r.kept == 1].assign(modelo=d, zona=lambda x: zona(x.lon)))
    return pd.concat(R)


def ajustes():
    filas = []
    for d in S9.ROT:
        b = json.load(open(os.path.join(ROOT, "models", d, "hotspot_fit.json"), encoding="utf-8"))["tracks"]["Foundation"]
        for t in TAGS:
            a = json.load(open(os.path.join(ROOT, "models", d, f"hotspot_fit{t}.json"), encoding="utf-8"))["tracks"]["Foundation"]
            filas.append(dict(modelo=d, var=t, modelo_var=a["model"], modelo_final=b["model"],
                              n_var=a["n_used"], n_final=b["n_used"],
                              deriva_var=np.hypot(*a["theta"][2:4]) if len(a["theta"]) > 2 else 0.0,
                              deriva_final=np.hypot(*b["theta"][2:4]) if len(b["theta"]) > 2 else 0.0,
                              w_x_var=a["w_x"], w_x_final=b["w_x"], h0_km=gc_km(a["theta"][:2], b["theta"][:2]),
                              excl_var=";".join(e["sample"] for e in a["excluded"])))
    return pd.DataFrame(filas)


def nominales():
    """With nominal_sensitivity.py (as in sens_conP0.py): Foundation enters step 04 from PAPER_HS_FIT. The
    json of the final fit is used for both rotations; its t_switch differs from that of the variants by <= 0.1 Myr
    and does not cross any integer age."""
    import subprocess
    N = []
    for d in S9.ROT:
        env = dict(os.environ, PAPER_HS_FIT=os.path.join(ROOT, "models", d, "hotspot_fit.json"))
        B = MP.recortar(pd.read_csv(os.path.join(ROOT, "results",
                        f"mc_collision_bands_total_present{MP.sufijo_bandas(d)}_hsb_cc_g13.csv")), d)
        B = B[np.isclose(B.pf.astype(float), 0.44) & (B.ridge == "Foundation")]
        B = B.assign(s68_km=(B.lat_p84 - B.lat_p16) / 2 * KM)[["age", "s68_km"]]
        for t in TAGS:
            fn = os.path.join(OUT, f"nominal_{d}{t}.csv")
            subprocess.run([sys.executable, os.path.join(HERE, "nominal_sensitivity.py"), os.path.join(ROOT, "models", d),
                            os.path.join(ROOT, "models", d, "rot", "sens9", t + ".rot"), os.path.join(ROOT, "models", d, "rot", S9.ROT[d]),
                            fn], check=True, stdout=subprocess.DEVNULL, env=env)
            n = pd.read_csv(fn)
            n = n[(n.ridge == "Foundation") & np.isclose(n.pf.astype(float), 0.44)].assign(modelo=d, var=t)
            N.append(n.merge(B, on="age", how="inner"))
    N = pd.concat(N, ignore_index=True)
    N["ad"] = N.dif_km.abs()
    return N


def main():
    R = residuos()
    R.to_csv(os.path.join(OUT, "residuals_foundation_zone.csv"), index=False, float_format="%.4f")
    print("(1) normalized residuals of the final fit per zone (nine models, samples used):")
    print(R.groupby("zona")[["z_along", "z_across"]].agg(["median", "mean", "count"]).round(2).to_string())
    print(R.groupby(["zona", "sample"])[["z_along", "z_across"]].median().round(2).to_string())
    A = ajustes()
    A.to_csv(os.path.join(OUT, "fits_selkirk.csv"), index=False, float_format="%.4f")
    print("\n(2) Foundation fits without the microplate samples:")
    print(A.round(1).to_string(index=False))
    N = nominales()
    N.to_csv(os.path.join(OUT, "nominal_selkirk.csv"), index=False, float_format="%.4f")
    print("\n(3) nominal Foundation collision (f_PAC 0.44), change relative to the final fit:")
    for t in TAGS:
        s = N[N["var"] == t].reset_index(drop=True)
        k = s.ad.values.argmax()
        print(f"  {t}: median {s.ad.median():.1f} km, max {s.ad.max():.1f} km ({s.modelo[k]}, {s.age[k]:.0f} Ma); "
              f"above the 68% half-width: {(s.ad > s.s68_km).sum()} of {len(s)}")
        print(s.groupby("modelo").ad.agg(["median", "max"]).round(1).to_string())


if __name__ == "__main__":
    main()
