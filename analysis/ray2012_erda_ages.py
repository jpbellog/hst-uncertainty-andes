# -*- coding: utf-8 -*-
"""Easter ages rebuilt from the step-heating data of Ray et al. (2012)
(ERDA, see ray2012_erda.py) for step 03b.

For each sample of models/Q22_WK08-A/data/ages_easter_ray2012_recal.csv:
  - plateau: mean of the same steps as the ArArCALC plateau with
    additive overdispersion (step_i ~ N(mu, s_i^2 + om^2), maximum
    likelihood), without multiplying by sqrt(MSWD) and without J (the ArArCALC
    errors of each step are analytical);
  - isochron: inverse isochron rebuilt with York et al. (2004) from the
    steps (ray2012_isochrons.py) when there are at least 3 steps with 36Ar;
    otherwise (DRFT-126-1: ArArCALC "Cannot Calculate"), the published value;
  - total fusion: the published value (not used in the fit except DRFT-100-2);
  - recalibration scenario A (FCT-3 = FCs 28.201 Ma, constants of Min et al.
    2000) with the same equation as ArAr_recalibration.py; the error scales with
    the age.
Adds the irradiation and the J (%1s) of each sample: step 03b treats J as a
systematic error common to each irradiation (one e_k per batch, sigma = mean J of the batch).

Output: columns age_E, err2s_E, irradiation, J_pct1s, src_E added to
ages_easter_ray2012_recal.csv in the seven model folders of CARPETAS (the other
columns do not change) and a copy in outputs/.

Usage:  py -P ray2012_erda_ages.py   (after ray2012_erda.py)
"""
import math
import os
import re
import shutil
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
LAM_SJ, LAM_MIN = 5.543e-10 * 1e6, 5.463e-10 * 1e6
FCS_OLD_RAY, FCS_NEW = 28.03, 28.201
MAT = {"gmass": "groundmass", "plag": "plagioclase", "wr": "whole_rock"}
CARPETAS = ["Q22_WK08-A", "WD25_WK08-A", "Q22_GMHRF-A", "Q22_GMHRF-B",
            "WD25_GMHRF-B", "Q22_GWG20", "WD25_GWG20"]


def recal(t):
    F = lambda a, lam: math.expm1(lam * a)
    R = F(t, LAM_SJ) / F(FCS_OLD_RAY, LAM_SJ)
    return math.log1p(R * F(FCS_NEW, LAM_MIN)) / LAM_MIN


def clave(titulo):
    """'DRFT-121-1 plag' -> ('DRFT-121-1', 'plagioclase'); 'GS7202 DR62-1 wr' ->
    ('GS7202-DR62-1', 'whole_rock')."""
    nom, mat = titulo.rsplit(" ", 1)
    return re.sub(r"\s+", "-", nom.strip()), MAT[mat]


def main():
    import sys
    sys.path.insert(0, HERE)                 # py -P does not add the script folder
    from ray2012_isochrons import isocrona
    E = pd.read_csv(os.path.join(HERE, "outputs", "ray2012_erda_samples.csv"))
    E["key"] = [clave(t) for t in E.erda_title]
    fn = os.path.join(ROOT, "models", "Q22_WK08-A", "data", "ages_easter_ray2012_recal.csv")
    head = [l for l in open(fn, encoding="utf-8") if l.startswith("#")]
    D = pd.read_csv(fn, comment="#")
    D = D.drop(columns=[c for c in ("age_E", "err2s_E", "irradiation", "J_pct1s", "src_E")
                        if c in D.columns])
    base = D["sample"].str.replace(r"p$", "", regex=True)      # DRFT-121-1p = plagioclase
    idx = {k: r for k, r in zip(E.key, E.itertuples())}
    age_E, err_E, irr, jp, src = [], [], [], [], []
    for s, m, r in zip(base, D.material, D.itertuples()):
        e = idx.get((s, m))
        if e is None:
            raise SystemExit(f"{r.sample} ({m}): no ERDA file")
        irr.append(e.irradiation)
        jp.append(e.J_pct1s)
        iso = isocrona(e.erda) if r.age_type == "isochron" else {}
        if r.age_type == "plateau" and e.n_plat > 0:
            a = recal(e.od_age)
            age_E.append(a)
            err_E.append(e.od_err2s * a / e.od_age)
            src.append("plateau ERDA, dispersion aditiva")
        elif iso.get("n", 0) >= 3:
            age_E.append(iso["t_A"])
            err_E.append(iso["s2_A"])
            src.append(f"isocrona inversa ERDA, York 2004 (n {iso['n']}, MSWD {iso['mswd']:.2f})")
        else:
            age_E.append(r.age_A)
            err_E.append(r.err2s_ma * r.age_A / r.age_ma)
            src.append(f"{r.age_type} publicada")
    D["age_E"], D["err2s_E"], D["irradiation"], D["J_pct1s"], D["src_E"] = age_E, err_E, irr, jp, src
    head = [h for h in head if "age_E" not in h] + [
        "# age_E, err2s_E (2 sigma): plateau recomputed from the ERDA steps (ray2012_erda_ages.py),\n",
        "# age_E: additive dispersion without sqrt(MSWD) or J, recalibrated as age_A; inverse isochron (York 2004)\n",
        "# age_E: from the steps if there are >= 3 with 36Ar; otherwise, and for total fusion, the published age. irradiation and J_pct1s (J, %1s) from each ArArCALC file.\n"]
    out = os.path.join(HERE, "outputs", "ages_easter_ray2012_recal.csv")
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.writelines(head)
        D.to_csv(f, index=False, float_format="%.4f", lineterminator="\n")
    for d in CARPETAS:
        shutil.copyfile(out, os.path.join(ROOT, "models", d, "data", "ages_easter_ray2012_recal.csv"))
    u = D[D.use_fit == 1]
    print(u[["sample", "age_type", "age_A", "err2s_ma", "age_E", "err2s_E", "irradiation",
             "J_pct1s"]].to_string(index=False, float_format="%.3f"))
    print("\nbatches:", u.groupby("irradiation").J_pct1s.agg(["count", "mean", "max"]).round(3).to_dict("index"))


if __name__ == "__main__":
    main()
