# -*- coding: utf-8 -*-
"""Step-heating data of Ray et al. (2012) from the EarthRef
Digital Archive (ERDA 1305-1342, one ArArCALC v2.40 .xls per sample,
contributed by R. A. Duncan; Ray et al. 2012, p. 6: "Full measurement data can
be accessed and downloaded from the EarthRef.org Digital Archive"). Copies in
../external/erda_ray2012/ (indice_erda.txt: title of each entry).

Extracts per sample: irradiation, J and its error (%1s), monitor, the steps
(temperature, 39Ar(k), age, 2s error, whether they enter the plateau) and the ArArCALC
results (plateau, normal and inverse isochrons, total fusion) with their MSWD.

Recomputes the plateau with the same steps:
  - weighted mean 1/s^2 without multiplying by sqrt(MSWD) (check: it must give the
    ArArCALC age);
  - mean with additive overdispersion: step_i ~ N(mu, s_i^2 + om^2), mu and om
    by maximum likelihood (in place of sqrt(MSWD)).
The ArArCALC errors of each step are analytical, without J (the plateau gives
"Analytical Error" and "External Error" separately); J is left out and treated
as a systematic error common to each irradiation.

Read with the .xls reader (xlrd); run with py -P (the .xls files are downloaded data).
Usage:  py -P ray2012_erda.py  ->  outputs/ray2012_erda_samples.csv,
                                 outputs/ray2012_erda_escalones.csv
"""
import os
import re
import sys
import numpy as np
import pandas as pd
import xlrd
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "external", "erda_ray2012")
OUT = os.path.join(HERE, "outputs")


def num(v):
    """'± 0.91' or 0.91 -> 0.91; empty -> nan."""
    if isinstance(v, (int, float)):
        return float(v)
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(v))
    return float(m.group()) if m else np.nan


def find(s, col, pat, r0=0):
    for r in range(r0, s.nrows):
        v = s.cell_value(r, col)
        if isinstance(v, str) and re.match(pat, v.strip()):
            return r
    return None


def leer(fn):
    """The ArArCALC "Data Tables" sheet changes its rows and, between versions, its
    columns: everything is located relative to the column o of the header
    "Incremental Heating" (row 2): step o..o+11, results o+4..o+11,
    normal isochron o+22, inverse isochron o+33."""
    s = xlrd.open_workbook(fn).sheet_by_index(0)
    o = next(c for c in range(s.ncols)
             if str(s.cell_value(2, c)).strip().startswith("Incremental"))
    steps = []
    for r in range(5, s.nrows):
        if str(s.cell_value(r, o + 2)).strip() == "S":
            break
        run = str(s.cell_value(r, o)).strip()
        age = num(s.cell_value(r, o + 8))          # in some sheets the age is text
        if not run or not np.isfinite(age) or not np.isfinite(num(s.cell_value(r, o + 9))):
            continue
        steps.append(dict(run=run, temp=s.cell_value(r, o + 1),
                          plateau=str(s.cell_value(r, o + 2)).strip() not in ("", "0.0"),
                          ar39k=num(s.cell_value(r, o + 6)), age=age,
                          err2s=num(s.cell_value(r, o + 9)),
                          pct39=num(s.cell_value(r, o + 11))))
    r5 = 5
    meta = dict(sample=str(s.cell_value(r5, 44)).strip(), material=s.cell_value(r5, 45),
                std_age=num(s.cell_value(r5, 49)), std_pct1s=num(s.cell_value(r5, 50)),
                J=num(s.cell_value(r5, 51)), J_pct1s=num(s.cell_value(r5, 52)),
                irradiation=str(s.cell_value(r5, 63)).strip(),
                standard=str(s.cell_value(r5, 65)).strip())
    c = o + 4
    rp = find(s, c, r"Results") + 3          # plateau row (without a plateau: another label)
    meta.update(plateau_label=str(s.cell_value(rp, c)).strip(), plateau_age=num(s.cell_value(rp, o + 8)),
                plateau_err2s=num(s.cell_value(rp, o + 9)), plateau_mswd=num(s.cell_value(rp, o + 10)),
                plateau_pct39=num(s.cell_value(rp, o + 11)), plateau_n=num(s.cell_value(rp + 1, o + 11)))
    ra = find(s, c, r"Analytical Error", rp)
    re_ = find(s, c, r"External Error", rp)
    meta.update(plateau_err2s_anal=num(s.cell_value(ra, o + 9)),
                plateau_magnif=num(s.cell_value(ra, o + 10)),
                plateau_err2s_ext=num(s.cell_value(re_, o + 9)) if re_ is not None and re_ < ra + 3
                else np.nan)
    for key, ca in (("iso_n", o + 22), ("iso_i", o + 33)):
        if ca + 2 < s.ncols:
            meta.update({f"{key}_age": num(s.cell_value(rp, ca)),
                         f"{key}_err2s": num(s.cell_value(rp, ca + 1)),
                         f"{key}_mswd": num(s.cell_value(rp, ca + 2))})
    rt = find(s, c, r"Total Fusion")
    if rt is not None:
        meta.update(tf_age=num(s.cell_value(rt, o + 8)), tf_err2s=num(s.cell_value(rt, o + 9)))
    return meta, steps


def plateau_ml(a, s1):
    """mu, sd(mu) and om of step_i ~ N(mu, s1_i^2 + om^2) (1 sigma)."""
    def nll(lo):
        v = s1 ** 2 + np.exp(lo) ** 2
        w = 1 / v
        mu = np.sum(w * a) / np.sum(w)
        return 0.5 * np.sum(np.log(v) + (a - mu) ** 2 / v)
    w0 = 1 / s1 ** 2
    mu0 = np.sum(w0 * a) / np.sum(w0)
    mswd = np.sum(w0 * (a - mu0) ** 2) / max(len(a) - 1, 1)
    if mswd <= 1 or len(a) < 3:
        return mu0, np.sqrt(1 / np.sum(w0)), 0.0
    sc = np.std(a) + 1e-6
    r = minimize_scalar(nll, bounds=(np.log(sc * 1e-4), np.log(sc * 10)), method="bounded")
    om = np.exp(r.x)
    w = 1 / (s1 ** 2 + om ** 2)
    return np.sum(w * a) / np.sum(w), np.sqrt(1 / np.sum(w)), om


def main():
    idx = {}
    for ln in open(os.path.join(SRC, "indice_erda.txt"), encoding="utf-8", errors="replace"):
        i, _, t = ln.partition("|")
        if "sample" in t:
            idx[int(i)] = t.split("sample", 1)[1].strip()
    rows, srows = [], []
    for i in sorted(idx):
        meta, steps = leer(os.path.join(SRC, f"{i}.xls"))
        st = pd.DataFrame(steps)
        p = st[st.plateau]
        if len(p) == 0:                       # "Cannot Calculate": no plateau
            meta.update(erda=i, erda_title=idx[i], n_steps=len(st), n_plat=0)
            rows.append(meta)
            srows += [dict(erda=i, sample=meta["sample"], **d) for d in steps]
            continue
        a, s1 = p.age.to_numpy(), p.err2s.to_numpy() / 2
        w = 1 / s1 ** 2
        mu = np.sum(w * a) / np.sum(w)
        mswd = np.sum(w * (a - mu) ** 2) / max(len(a) - 1, 1)
        mu_ml, sd_ml, om = plateau_ml(a, s1)
        meta.update(erda=i, erda_title=idx[i], n_steps=len(st), n_plat=len(p),
                    rec_age=mu, rec_err2s_int=2 * np.sqrt(1 / np.sum(w)), rec_mswd=mswd,
                    rec_err2s_sqrtmswd=2 * np.sqrt(1 / np.sum(w)) * np.sqrt(max(mswd, 1)),
                    od_age=mu_ml, od_err2s=2 * sd_ml, od_omega1s=om)
        rows.append(meta)
        srows += [dict(erda=i, sample=meta["sample"], **d) for d in steps]
    R = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    R.to_csv(os.path.join(OUT, "ray2012_erda_samples.csv"), index=False, float_format="%.5g")
    pd.DataFrame(srows).to_csv(os.path.join(OUT, "ray2012_erda_escalones.csv"), index=False,
                               float_format="%.6g")
    cols = ["erda", "sample", "irradiation", "J_pct1s", "plateau_age", "plateau_err2s", "plateau_mswd",
            "plateau_magnif", "rec_age", "rec_err2s_int", "rec_err2s_sqrtmswd", "od_age", "od_err2s",
            "iso_i_age", "tf_age"]
    with pd.option_context("display.width", 250, "display.max_rows", 100):
        print(R[cols].to_string(index=False, float_format="%.3f"))


if __name__ == "__main__":
    sys.exit(main())
