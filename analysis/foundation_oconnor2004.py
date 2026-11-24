# -*- coding: utf-8 -*-
"""40Ar/39Ar ages of the Foundation chain from the sources published after 1998.
The ages recalculated in 2001 are the base; the script also allows analysing the effect of
adding the 2004 ages and of removing Seamounts 8, 9 and 13a.

Sources:
  - O'Connor et al. (2001, EPSL 192, 633-648, erratum), Dataset 3 (SI3): the 23 ages of
    O'Connor et al. (1998) recalculated with TCR = 28.34 Ma (Renne et al. 1998): plateau +/- 2s,
    MSWD, J +/- and laboratory ID. No step-heating data.
  - O'Connor et al. (2001), Dataset 1 (SI1): step-heating data (age +/- 2s per step, plateau steps
    flagged with 1) and J of the new samples (SO100 74DS-1 to 101GTV-2 and the
    'Hotline' samples FH DR11 to DR15).
  - O'Connor et al. (2004, ch. 11 of Oceanic Hotspots), Table 11.1: compilation of all
    ages with on/off bottom position; it also gives FH DR4-1 (Del Cano), FH DR5-2 and FH DR6-1 (VER,
    ~16 Ma) and FH DR7-16 (failed rift of the Selkirk microplate), from O'Connor et al. (2002), without
    MSWD or steps. Decay constant: total lambda 5.543e-10/a (footnote of Table 11.1; SI3 says
    "5.53 x 10-10", taken as a rounding of that of Steiger and Jager).

What is done:
  - SI1: for each analysis, mean of the plateau steps with additive overdispersion by
    maximum likelihood (plateau_ml of ray2012_erda.py, as for Easter), without sqrt(MSWD) or added
    J. If there are several analyses of the same sample, their weighted mean.
  - SI3: the published error is divided by the MSWD when MSWD > 1 (O'Connor et al. 1998, p. 47:
    "the reported analytical error is multiplied by the MSWD") and multiplied by sqrt(MSWD):
    without steps the additive overdispersion used for Easter cannot be fitted, and for a plateau of
    similar errors sqrt(MSWD) is its approximate equivalent (the pure analytical error, / MSWD,
    would give +/- 0.04 Ma for SO100-50DS-1, MSWD 5.1). Several analyses of one sample: weighted mean.
  - 2002 (Del Cano and VER): value of Table 11.1 as is (the inflation cannot be undone).
    FH DR7-16 (failed rift): use_fit = 0, it does not belong to the hotspot.
  - Recalibration TCR 28.34 Ma (SJ) -> FCs 28.201 Ma (Kuiper 2008) with Min et al. (2000) and
    TCs/FCs = 1.0112 (Renne et al. 1998), the same equation as foundation_oconnor1998.py.
  - On-bottom position of Table 11.1 (2004). SO100-50DS-1: 34 51.492'S (Table 1 of 1998
    gave 35 51.492'S on bottom and 34 51.860'S off bottom; the 2004 table confirms the erratum).
  - J: J_pct1s = 100 x (+/- J) / J, taken as 1 sigma (O'Connor et al. 1998, p. 46: precision
    of J "0.1 to 0.3% (s.e.m.)"). Irradiation batches: OC98a (7 h, Seamounts 1a to 26) and OC98b
    (1.5 h, 27b to 30), O'Connor et al. (1998, p. 46); OC01 for those of SI1 (one batch, ASSUMPTION:
    SI1 does not say how many irradiations there were); OC02 for those of 2002 (J not published: the
    mean J of the others is used, ASSUMPTION).
  - Columns: src (1998, 2001, 2002), set (v1 = the 23 of 1998 recalculated; v2 = all),
    oconnor_excl (Seamounts 8, 9 and 13a).

Outputs: outputs/ages_foundation_oconnor2004_v1.csv (src 1998 only), _v2.csv (all) and _v3.csv
(= v1 with use_fit = 0 for Seamounts 8, 9 and 13a; the base set), copied
to data/ of the twelve folders; outputs/foundation_oconnor2004_control.csv (check against
Table 11.1 and the 1998 positions).
Usage:  py foundation_oconnor2004.py
"""
import importlib.util
import math
import os
import re
import shutil
import numpy as np
import pandas as pd
import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, ".."))
PAP = (__import__("os").environ["PAPER_ROOT"] + "/external/oconnor_pdfs")
SI1 = os.path.join(PAP, "OConnor_etal_2001_En-echelon-volcanic-ridges-Foundation-Chain-"
                        "Pacific-Antarctic-spreading-center_EPSL_SI1.pdf")
CARPETAS = ["Q22_WK08-A", "WD25_WK08-A", "Q22_GMHRF-A", "Q22_GMHRF-B",
            "WD25_GMHRF-B", "Q22_GWG20", "WD25_GWG20",
            "SG12_WK08-A", "WD25_GMHRF-A", "SG12_GMHRF-A",
            "SG12_GMHRF-B", "SG12_GWG20"]

LAM_SJ, LAM_MIN = 5.543e-10 * 1e6, 5.463e-10 * 1e6
FCS_NEW, TCR_2001, R_TCS_FCS = 28.201, 28.34, 1.0112

_s = importlib.util.spec_from_file_location("erda", os.path.join(HERE, "ray2012_erda.py"))
_m = importlib.util.module_from_spec(_s)
_s.loader.exec_module(_m)
plateau_ml = _m.plateau_ml


def F(t, lam):
    return math.expm1(lam * t)


def recal(t_u, t_s_old, lam_old, t_s_new, lam_new):
    return math.log1p(F(t_u, lam_old) / F(t_s_old, lam_old) * F(t_s_new, lam_new)) / lam_new


def dm(d, m):
    return -(d + m / 60.0)


def wmean(a, e):
    w = 1 / np.asarray(e, float) ** 2
    return float(np.sum(w * a) / np.sum(w)), float(np.sqrt(1 / np.sum(w)))


# SI3 (O'Connor 2001, Dataset 3): sample, seamount, material, plateau, 2s, MSWD, J, +/-J, batch,
# on-bottom position (Table 11.1, 2004), excl. O'Connor 1998
SI3 = [
    ("SO100-11DS-1", "1a", "wr", 21.6, 0.2, 0.04, 0.001820, 0.000005, (32, 56.414), (130, 45.459), 0),
    ("FHDR1-3", "Aristotelis", "wr", 18.5, 0.3, 0.1, 0.001880, 0.000009, (32, 30.33), (127, 30.09), 0),
    ("SO100-18DS-1", "5", "wr", 16.6, 0.4, 1.8, 0.001780, 0.000018, (32, 28.871), (126, 0.456), 0),
    ("SO100-17DS-1", "5", "wr", 16.3, 0.3, 0.80, 0.001480, 0.000007, (32, 28.474), (126, 4.069), 0),
    ("SO100-28GTV-2", "8", "plag", 9.0, 0.1, 0.9, 0.001770, 0.000005, (33, 41.735), (124, 54.612), 1),
    ("SO100-25DS-1", "9", "wr", 13.3, 0.1, 1.9, 0.001750, 0.000005, (33, 20.764), (123, 52.525), 1),
    ("SO100-26DS-1", "9", "plag", 13.1, 0.6, 0.30, 0.001740, 0.000005, (33, 31.700), (124, 6.434), 1),
    ("SO100-33DS-1", "10", "plag", 12.7, 0.3, 0.10, 0.001690, 0.000005, (34, 7.384), (122, 22.038), 0),
    ("SO100-38DS-1", "11", "plag", 11.8, 0.2, 0.1, 0.001680, 0.000005, (34, 19.119), (121, 58.729), 0),
    ("SO100-41DS-1", "12b", "plag", 11.8, 0.2, 0.3, 0.001640, 0.000005, (34, 52.352), (121, 33.309), 0),
    ("SO100-45DS-1", "13a", "plag", 13.6, 0.2, 0.5, 0.001620, 0.000005, (35, 3.180), (120, 43.226), 1),
    ("SO100-46DS-2", "16", "wr", 10.5, 0.2, 0.6, 0.001400, 0.000007, (34, 57.336), (120, 24.490), 0),
    ("SO100-50DS-1", "18", "plag", 9.5, 0.2, 0.8, 0.001530, 0.000005, (34, 51.492), (119, 6.724), 0),
    ("SO100-50DS-1", "18", "wr", 9.5, 0.2, 5.1, 0.001460, 0.000007, (34, 51.492), (119, 6.724), 0),
    ("SO100-54DS-1", "19b", "wr", 8.5, 0.1, 1.7, 0.001480, 0.000004, (35, 6.976), (118, 33.132), 0),
    ("SO100-56DS-1", "21a", "plag", 7.8, 0.1, 0.2, 0.001420, 0.000004, (35, 22.943), (118, 5.256), 0),
    ("SO100-63DS-1", "22", "wr", 7.8, 0.1, 1.6, 0.001270, 0.000004, (35, 48.000), (117, 26.305), 0),
    ("SO100-60DS-1", "23", "plag", 7.6, 0.1, 0.2, 0.001310, 0.000004, (35, 27.003), (117, 11.833), 0),
    ("SO100-59DS-1", "24b", "plag", 7.3, 0.2, 0.6, 0.001330, 0.000004, (35, 26.729), (116, 38.813), 0),
    ("SO100-67DS-4", "25", "plag", 6.3, 0.2, 0.6, 0.001560, 0.000005, (36, 1.980), (115, 59.276), 0),
    ("SO100-66DS-1", "26", "plag", 4.9, 0.2, 0.1, 0.001250, 0.000004, (35, 47.484), (115, 39.335), 0),
    ("SO100-66DS-1", "26", "plag", 4.8, 0.2, 0.7, 0.001250, 0.000004, (35, 47.484), (115, 39.335), 0),
    ("SO100-66DS-1", "26", "wr", 4.7, 0.2, 3.2, 0.001590, 0.000005, (35, 47.484), (115, 39.335), 0),
    ("SO100-69DS-1", "27b", "wr", 5.1, 0.1, 0.7, 0.000185, 0.000001, (36, 33.804), (115, 16.627), 0),
    ("SO100-70DS-2", "28", "wr", 3.7, 0.2, 0.8, 0.000190, 0.000001, (36, 20.907), (113, 55.696), 0),
    ("SO100-71DS-1", "30", "wr", 2.1, 0.1, 0.4, 0.000195, 0.000001, (36, 40.812), (113, 28.350), 0),
]
LOTE_B = {"27b", "28", "30"}

# SI1 (2001): name in SI1 -> (sample, edifice, on-bottom position of Table 11.1, plateau
# and 2s of Table 11.1 for checking; mean weight if there are several analyses)
SI1_POS = {
    "SO100 74DS-1": ("SO100-74DS-1", "32 Ohm", (36, 56.987), (112, 12.875), 1.8, 0.3),
    "SO100 75DS-2": ("SO100-75DS-2", "39 Schrodinger", (37, 16.680), (111, 51.521), 1.0, 0.2),
    "SO100 76DS-1": ("SO100-76DS-1", "33 Pascal", (37, 22.190), (112, 6.098), 1.1, 0.1),
    "SO100 77DS-1": ("SO100-77DS-1", "35 Pauling", (37, 3.427), (111, 39.605), 1.1, 0.2),
    "SO100 82DS-1": ("SO100-82DS-1", "North Ridge", (37, 27.530), (111, 12.654), 0.4, 0.1),
    "SO100 87DS-1": ("SO100-87DS-1", "South Ridge", (37, 38.941), (111, 17.397), 0.5, 0.1),
    "SO100 89DS-1": ("SO100-89DS-1", "South Ridge", (37, 31.388), (111, 40.707), 0.7, 0.1),
    "SO100 90DS-1": ("SO100-90DS-1", "38 Rutherford", (38, 21.598), (110, 37.905), 1.1, 0.3),
    "SO100 93DS-1": ("SO100-93DS-1", "36 Planck", (38, 8.063), (111, 32.894), 0.5, 0.1),
    "SO100 94DS-1": ("SO100-94DS-1", "37a Platon", (37, 53.367), (112, 53.847), 0.9, 0.1),
    "SO100 95DS-1": ("SO100-95DS-1", "37b Richter", (37, 53.312), (112, 25.089), 0.9, 0.1),
    "SO100 97DS-2": ("SO100-97DS-2", "33 Pascal", (37, 25.016), (112, 3.419), 1.1, 0.4),
    "SO100 99DS-5": ("SO100-99DS-5", "34 Pasteur", (37, 15.768), (112, 3.044), 1.1, 0.1),
    "SO100 101 GTV-2": ("SO100-101GTV-2", "33 Pascal", (37, 23.149), (112, 5.940), 1.1, 0.3),
    "FH DR11-1": ("FHDR11-1", "31 Newton", (36, 55.477), (113, 4.585), 2.1, 0.2),
    "FH DR12-8": ("FHDR12-8", "Wegener", (37, 55.935), (113, 42.372), 1.9, 0.1),
    "FH DR13-1": ("FHDR13-1", "29 Mendeleiev", (37, 1.600), (114, 2.56), 3.3, 0.1),
    "FH DR14-5": ("FHDR14-5", "Mohorovicic", (37, 27.78), (114, 34.773), 2.8, 0.1),
    "FH DR14-13": ("FHDR14-13", "Mohorovicic", (37, 27.78), (114, 34.773), 0.7, 0.03),
    "FH DR14-14": ("FHDR14-14", "Mohorovicic", (37, 27.78), (114, 34.773), 3.5, 0.1),
    "FH DR14-16": ("FHDR14-16", "Mohorovicic", (37, 27.78), (114, 34.773), 3.6, 0.1),
    "FH DR14-17": ("FHDR14-17", "Mohorovicic", (37, 27.78), (114, 34.773), 3.5, 0.9),
    "FH DR15-4": ("FHDR15-4", "Linne Ridge", (35, 44.983), (114, 21.989), 5.0, 0.2),
}

# 2002 (Table 11.1, 2004): sample, edifice, lat, lon, plateau, 2s, use_fit
T2002 = [
    ("FHDR4-1", "Del Cano", (33, 46.69), (126, 43.83), 16.1, 0.4, 1),
    ("FHDR5-2", "VER Boltzmann", (34, 55.112), (126, 13.798), 16.8, 0.3, 1),
    ("FHDR6-1", "VER Laplace", (34, 34.000), (125, 16.400), 16.6, 0.2, 1),
    ("FHDR7-16", "Failed Rift (Selkirk)", (35, 20.849), (124, 45.723), 20.4, 0.3, 0),
]


def leer_si1():
    """{name: [(J, dJ, ages, err2s), ...]} with the steps flagged 1."""
    with pdfplumber.open(SI1) as p:
        lines = "\n".join(pg.extract_text() or "" for pg in p.pages).splitlines()
    out, cur = {}, None
    hdr = re.compile(r"^((?:SO100|FH) [^()]+?) \(.*\) J = ([\d.]+) ± ([\d.]+)")
    for ln in lines:
        m = hdr.match(ln.strip())
        if m:
            cur = dict(J=float(m.group(2)), dJ=float(m.group(3)), a=[], e=[])
            out.setdefault(m.group(1).strip(), []).append(cur)
            continue
        tk = ln.split()
        if cur is None or len(tk) < 10:
            continue
        # step flag: after "fusion" or after the power (the single fusions of
        # 97DS-2 and 101GTV-2 have one column fewer at the end)
        i = tk.index("fusion") + 1 if "fusion" in tk else 2
        try:
            flag, age, err = int(tk[i]), float(tk[i + 5]), float(tk[i + 6])
        except (ValueError, IndexError):
            continue
        if flag == 1:
            cur["a"].append(age)
            cur["e"].append(err)
    return out


def main():
    tcs = math.log1p(R_TCS_FCS * F(FCS_NEW, LAM_MIN)) / LAM_MIN
    rc = lambda a: recal(a, TCR_2001, LAM_SJ, tcs, LAM_MIN)
    rows, ctrl = [], []

    # 1998 (SI3)
    df3 = pd.DataFrame(SI3, columns=["sample", "smt", "mat", "age", "e2", "mswd", "J", "dJ",
                                     "la", "lo", "ex"])
    df3["e2n"] = np.where(df3.mswd > 1, df3.e2 / np.sqrt(df3.mswd), df3.e2)
    for s, g in df3.groupby("sample", sort=False):
        a, e2 = wmean(g.age, g.e2n) if len(g) > 1 else (float(g.age.iloc[0]), float(g.e2n.iloc[0]))
        r0 = g.iloc[0]
        rows.append(dict(sample=s, edifice=f"Seamount {r0.smt}", lat=dm(*r0.la), lon=dm(*r0.lo),
                         age_ma=a, err2s_ma=e2, err2s_pub=float(g.e2.max()),
                         mswd_max=float(g.mswd.max()), n_analisis=len(g),
                         age_type="plateau", use_fit=1, segment="Foundation",
                         material="/".join(sorted(set(g.mat))), oconnor_excl=int(r0.ex), src=1998,
                         irradiation="OC98b" if r0.smt in LOTE_B else "OC98a",
                         J_pct1s=float(100 * (g.dJ / g.J).mean())))

    # 2001 (SI1, steps)
    si1 = leer_si1()
    for nom, (s, ed, la, lo, a_pub, e_pub) in SI1_POS.items():
        if nom not in si1:
            raise KeyError(f"not in SI1: {nom}")
        mus, sds, oms, jp, ns = [], [], [], [], []
        for an in si1[nom]:
            a, e = np.array(an["a"]), np.array(an["e"])
            if len(a) == 0:
                continue
            mu, sd, om = plateau_ml(a, e / 2.0)
            mus.append(mu); sds.append(sd); oms.append(om)
            jp.append(100 * an["dJ"] / an["J"]); ns.append(len(a))
        a, s1 = wmean(mus, sds) if len(mus) > 1 else (mus[0], sds[0])
        rows.append(dict(sample=s, edifice=ed, lat=dm(*la), lon=dm(*lo), age_ma=a, err2s_ma=2 * s1,
                         err2s_pub=e_pub, mswd_max=np.nan, n_analisis=len(mus),
                         age_type="plateau", use_fit=1, segment="Foundation", material="wr",
                         oconnor_excl=0, src=2001, irradiation="OC01", J_pct1s=float(np.mean(jp))))
        ctrl.append(dict(sample=s, plateau_tabla=a_pub, err2s_tabla=e_pub, plateau_si1=a,
                         err2s_si1=2 * s1, dif=a - a_pub, n_escalones=sum(ns),
                         omega_max=max(oms)))

    jm = float(np.mean([r["J_pct1s"] for r in rows]))
    for s, ed, la, lo, a, e2, uf in T2002:
        rows.append(dict(sample=s, edifice=ed, lat=dm(*la), lon=dm(*lo), age_ma=a, err2s_ma=e2,
                         err2s_pub=e2, mswd_max=np.nan, n_analisis=1, age_type="plateau",
                         use_fit=uf, segment="Foundation", material="wr", oconnor_excl=0,
                         src=2002, irradiation="OC02", J_pct1s=jm))

    df = pd.DataFrame(rows)
    df["age_recal"] = [rc(a) for a in df.age_ma]
    df["set"] = np.where(df.src == 1998, "v1", "v2")

    # check of positions against the 1998 CSV
    old = pd.read_csv(os.path.join(HERE, "outputs", "ages_foundation_oconnor1998_recal.csv"),
                      comment="#")
    m = df.merge(old[["sample", "lat", "lon", "age_recal"]], on="sample", suffixes=("", "_1998"))
    m["dpos_km"] = 111.195 * np.hypot(m.lat - m.lat_1998,
                                      (m.lon - m.lon_1998) * np.cos(np.radians(m.lat)))
    m["dage_recal"] = m.age_recal - m.age_recal_1998

    hdr = ("# Edades 40Ar/39Ar de la cadena Foundation desde O'Connor et al. (2001: SI3 recalculo de 1998\n"
           "# con TCR 28.34 Ma y SI1 por escalon) y O'Connor et al. (2004, Tabla 11.1; Del Cano y VER de\n"
           "# 2002). age_ma: TCR 28.34 Ma; err2s_ma: 2 sigma sin la inflacion por MSWD (SI3: / MSWD si\n"
           "# MSWD > 1; SI1: plateau con sobredispersion aditiva ML). age_recal: FCs 28.201 / Min 2000\n"
           f"# (TCs = {tcs:.4f} Ma). J_pct1s: +/-J/J en %, 1 sigma (supuesto s.e.m.). Generado por\n"
           "# analysis/foundation_oconnor2004.py\n")
    cols = ["sample", "edifice", "lat", "lon", "age_ma", "err2s_ma", "age_type", "use_fit",
            "segment", "material", "oconnor_excl", "src", "set", "irradiation", "J_pct1s",
            "age_recal", "err2s_pub", "mswd_max", "n_analisis"]
    v3 = df[df.src == 1998].copy()
    v3.loc[v3.oconnor_excl == 1, "use_fit"] = 0
    for v, sub in (("v1", df[df.src == 1998]), ("v2", df), ("v3", v3)):
        out = os.path.join(HERE, "outputs", f"ages_foundation_oconnor2004_{v}.csv")
        with open(out, "w", encoding="utf-8", newline="\n") as f:
            f.write(hdr)
            sub[cols].to_csv(f, index=False, float_format="%.6f")
        for d in CARPETAS:
            shutil.copy2(out, os.path.join(P, "models", d, "data"))
    pd.DataFrame(ctrl).to_csv(os.path.join(HERE, "outputs", "foundation_oconnor2004_control.csv"),
                              index=False, float_format="%.4f")
    pd.set_option("display.width", 200)
    print(f"TCs = {tcs:.4f} Ma; factor 2001 -> FCs 28.201: {rc(10.0) / 10.0:.6f}")
    print("\nCheck of SI1 against Table 11.1 (2004):")
    print(pd.DataFrame(ctrl).round(3).to_string(index=False))
    print("\nChange relative to the 1998 CSV (position and recalibrated age):")
    print(m[["sample", "dpos_km", "dage_recal", "err2s_ma", "mswd_max"]].round(3).to_string(index=False))
    print(f"\nrows: v1 {int((df.src == 1998).sum())}, v2 {len(df)} (use_fit {int(df.use_fit.sum())})")
    print(df.groupby("irradiation").J_pct1s.mean().round(3).to_string())


if __name__ == "__main__":
    main()
