# -*- coding: utf-8 -*-
"""Tables of the Supplementary Material.

Table S1 (tableS1.tex): fit of each hotspot in the twelve models (nine of the main text and
  GMHRF-A'), read from <folder>/hotspot_fit.json (step 03b, confidence fit, g13 run).
  n = ages used / ages that enter the fit, not counting the present-day position (present_obs).
  Drift = norm of the hotspot velocity in the frame (theta[2:], km/Myr); "--" if the hotspot
  is fixed. Shift = distance between the hotspot fitted at 0 Ma and the present-day reference
  position (seed_shift_km).
Table S2 (TableS2_collision_bands.csv): collision latitude bands on the present-day margin, total
  variant, of the twelve models, three f_PAC, 0-66 Ma every 1 Myr, clipped by the frame_coverage of
  the frame covariance (paper_models.recortar); clipped rows have valid = 0.
Usage:  py make_supplementary_tables.py
"""
import json
import os
import sys
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "figures"))
import paper_models as MP

G13 = os.path.join(ROOT, "results")
TAG = "_hsb_cc_g13"
MODELOS = MP.TEXTO + MP.SUPL
CHAINS = [("Easter", "Easter"), ("TT", "San F\\'elix"), ("JFR", "Juan Fern\\'andez"),
          ("Foundation", "Foundation")]
CIRC = {"Quiero": "Q22", "WD25": "WD25", "SG12": "SG12"}


def tabla_s1():
    rows = []
    for carpeta, circ, marco in MODELOS:
        d = json.load(open(os.path.join(ROOT, "models", carpeta, "hotspot_fit.json"), encoding="utf-8"))
        p0 = 1 if d["meta"].get("present_obs") else 0
        for key, lab in CHAINS:
            x = d["tracks"][key]
            th = x["theta"]
            drift = np.hypot(th[2], th[3]) if x["model"] == "deriva" else None
            rows.append(dict(chain=lab, model=f"{CIRC[circ]} + {marco}",
                             n=f"{x['n_used'] - p0}/{x['n_input'] - p0}",
                             hs="drifting" if x["model"] == "deriva" else "fixed",
                             drift="--" if drift is None else f"{drift:.1f}",
                             w_t=f"{x['w_t']:.2f}" + ("$^a$" if x["w_from_jfr"] else ""),
                             w_x=f"{x['w_x']:.0f}" + ("$^a$" if x["w_from_jfr"] else ""),
                             rms=f"{x['rms_along_km']:.0f}/{x['rms_across_km']:.0f}",
                             shift=f"{x['seed_shift_km']:.0f}",
                             excl=", ".join(e["sample"] for e in x["excluded"]) or "--"))
    R = pd.DataFrame(rows)
    out = []
    for lab in [c[1] for c in CHAINS]:
        sub = R[R.chain == lab]
        out.append(r"\midrule" + "\n" + r"\multicolumn{9}{l}{\textit{" + lab + r"}} \\")
        for _, r in sub.iterrows():
            out.append(" & ".join([r.model, r.n, r.hs, r.drift, r.w_t, r.w_x, r.rms, r["shift"],
                                   r.excl.replace("_", r"\_")]) + r" \\")
    open(os.path.join(HERE, "tableS1.tex"), "w", encoding="utf-8").write("\n".join(out) + "\n")
    R.to_csv(os.path.join(HERE, "TableS1_hotspot_fits.csv"), index=False)
    print(R.to_string(index=False))


def tabla_s2():
    dd = []
    for carpeta, circ, marco in MODELOS:
        f = os.path.join(G13, f"mc_collision_bands_total_present{MP.sufijo_bandas(carpeta)}{TAG}.csv")
        b = pd.read_csv(f, dtype={"pf": str})
        r = MP.recortar(b, carpeta)
        b["valid"] = r["lat_p50"].notna().astype(int)
        b.insert(0, "frame", marco)
        b.insert(0, "circuit", circ)
        dd.append(b)
    T = pd.concat(dd)
    T["ridge"] = T.ridge.replace({"TT": "San Felix", "JFR": "Juan Fernandez"})
    T = T[T.age <= 66]
    T.to_csv(os.path.join(HERE, "TableS2_collision_bands.csv"), index=False, float_format="%.3f")
    print("Table S2:", len(T), "rows;", int(T.valid.sum()), "valid")


if __name__ == "__main__":
    tabla_s1()
    tabla_s2()
