# -*- coding: utf-8 -*-
"""Numbers of the main text and of the supplementary material.

Reads the saved outputs of the nine models (g13 collision bands, ellipses, hotspot fits,
sensitivity tests) and writes each number that appears in main.tex and supl.tex next to the
number written in the text, with the .tex line. Nothing is computed by hand: a number of the
text that does not come out of this script is not supported.

Usage:  py manuscript_numbers.py            (all sections)
        py manuscript_numbers.py bandas     (one section)
Output: outputs/manuscript_numbers.csv (ref, descripcion, calculado, texto, coincide) and a summary on screen.

Conventions: f_PAC = 0.44 unless stated otherwise; each model is clipped at the age of its
frame covariance (paper_models.recortar); 68% half-width = (p84 - p16)/2; 95% width =
p97.5 - p2.5; latitudes in degrees, mapped to the present margin.
"""
import csv
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "figures"))
os.environ.setdefault("PAPER_MODELOS", "texto")
import paper_models as MP  # noqa: E402

G13 = os.path.join(ROOT, "results")
MAIN = os.path.join(HERE, "table2_rows.tex")   # rows of Table 2 of main.tex (the manuscript text is not published)
SAL = os.path.join(HERE, "outputs")
TAG = "_hsb_cc_g13"
PF = "0.44"
CADENAS = {"Easter": "Easter", "TT": "San Felix", "JFR": "Juan Fernandez", "Foundation": "Foundation"}

FILAS = []


def chk(ref, desc, calc, texto):
    """Records one number: calc is the string formatted as in the text."""
    FILAS.append(dict(ref=ref, descripcion=desc, calculado=calc, texto=texto, coincide=calc == texto))


def red(x, nd=1):
    """Round half up (convention of the text; format() rounds half to even)."""
    from decimal import Decimal, ROUND_HALF_UP
    q = Decimal(1).scaleb(-nd)
    return str(Decimal(repr(float(x))).quantize(q, rounding=ROUND_HALF_UP))


def rango(v, fmt="{:.1f}"):
    nd = int(re.search(r"\.(\d)f", fmt).group(1))
    a, b = red(np.nanmin(v), nd), red(np.nanmax(v), nd)
    return a if a == b else f"{a}--{b}"


def bandas(variante="total", modelos=None, pf=PF):
    """Bands at the present margin of the models, clipped by the frame coverage."""
    out = []
    for c in (modelos or MP.CARPETAS):
        fn = os.path.join(G13, f"mc_collision_bands_{variante}_present{MP.sufijo_bandas(c)}{TAG}.csv")
        d = MP.recortar(pd.read_csv(fn, dtype={"pf": str}), c)
        if pf is not None:
            d = d[d.pf.astype(float) == float(pf)]
        out.append(d.assign(modelo=c))
    d = pd.concat(out, ignore_index=True)
    d["hw68"] = (d.lat_p84 - d.lat_p16) / 2
    d["w95"] = d.lat_p97_5 - d.lat_p2_5
    return d


# ------------------------------------------------------------------ Results: collision bands
def sec_bandas():
    d = bandas().dropna(subset=["lat_p50"])
    med = d.groupby("ridge").hw68.median()
    for r, t in (("TT", "0.8"), ("Foundation", "0.9"), ("Easter", "1.2"), ("JFR", "1.2")):
        chk("main:163", f"median 68% half-width, single model, {CADENAS[r]}", f"{med[r]:.1f}", t)
    mx = d.groupby("ridge").hw68.max()
    chk("main:163", "maximum 68% half-width per chain (range across chains)", rango(mx), "2.1--3.6")
    chk("main:45/281", "median 68% half-width, range across chains", rango(med), "0.8--1.2")
    chk("main:45", "maximum 68% half-width", f"{d.hw68.max():.1f}", "3.6")
    # 95% envelope of all models and the three f_PAC
    a = bandas(pf=None).dropna(subset=["lat_p50"])
    env = a.groupby(["ridge", "age"]).agg(lo=("lat_p2_5", "min"), hi=("lat_p97_5", "max"))
    env["w"] = env.hi - env.lo
    em = env.groupby("ridge").w.median()
    for r, t in (("Foundation", "6.4"), ("TT", "8.0"), ("Easter", "9.4"), ("JFR", "10.0")):
        chk("main:163", f"median width of the 95% envelope, {CADENAS[r]}", f"{em[r]:.1f}", t)
    ex = env.groupby("ridge").w.max()
    chk("main:163", "maximum 95% envelope, Easter/Foundation",
        rango(ex[["Easter", "Foundation"]], "{:.0f}"), "15")
    am = env.reset_index().loc[lambda q: q.groupby("ridge").w.idxmax()].set_index("ridge").age
    chk("main:163", "age of the maximum 95% envelope, Easter/Foundation (Ma)",
        rango(am[["Easter", "Foundation"]], "{:.0f}"), "16")
    chk("main:163", "maximum 95% envelope JFR and its age", f"{ex['JFR']:.0f} a {am['JFR']:.0f}", "13 a 61")
    chk("main:281", "median 95% envelope, range across chains", rango(em, "{:.0f}"), "6--10")
    for r in CADENAS:
        e = env.loc[r]
        print(f"   95% envelope {r}: maximum {e.w.max():.2f} at {e.w.idxmax():.0f} Ma")


def fmt_lat(x):
    return "--" if not np.isfinite(x) else f"{x:.1f}"


def sec_tabla1():
    """Table 1: median of the medians (f_PAC 0.44) with the number of valid models, and the
    95% envelope of all models and f_PAC. The values of the text are read from main.tex."""
    tex = open(MAIN, encoding="utf8").read().split("\n")
    d = bandas().dropna(subset=["lat_p50"])
    a = bandas(pf=None).dropna(subset=["lat_p50"])
    orden = ["Easter", "TT", "JFR", "Foundation"]
    for i, l in enumerate(tex):
        m = re.match(r"^(\d+) & (.*)\\\\$", l.strip())
        if not m or "scriptsize" not in tex[i + 1]:
            continue
        age = float(m.group(1))
        celdas = [c.strip() for c in m.group(2).split("&")]
        # below each age, row of the 68% union (f_PAC 0.44) and row of the
        # 95% envelope (all f_PAC)
        u68_txt = re.findall(r"\[([^\]]*)\]", tex[i + 1])
        env_txt = re.findall(r"\[([^\]]*)\]", tex[i + 2])
        for k, r in enumerate(orden):
            s = d[(d.ridge == r) & (d.age == age)]
            calc = "--" if s.empty else f"{s.lat_p50.median():.1f} ({len(s)})"
            txt = re.sub(r"[\$\s]", "", celdas[k]).replace("(", " (")
            chk(f"main:{i + 1}", f"Table 1, {CADENAS[r]} {age:.0f} Ma, median (n)", calc, txt)
            if s.empty:
                continue
            calc = f"{s.lat_p16.min():.1f}, {s.lat_p84.max():.1f}"
            txt = re.sub(r"[\$\s]", "", u68_txt[k]).replace(",", ", ") if k < len(u68_txt) else ""
            chk(f"main:{i + 2}", f"Table 1, {CADENAS[r]} {age:.0f} Ma, 68% union", calc, txt)
            e = a[(a.ridge == r) & (a.age == age)]
            calc = f"{e.lat_p2_5.min():.1f}, {e.lat_p97_5.max():.1f}"
            txt = re.sub(r"[\$\s]", "", env_txt[k]).replace(",", ", ") if k < len(env_txt) else ""
            chk(f"main:{i + 3}", f"Table 1, {CADENAS[r]} {age:.0f} Ma, 95% envelope", calc, txt)


def sec_frame_coverage():
    """Oldest valid collision age per model and chain (clipping by frame coverage)."""
    MP2 = __import__("importlib").import_module("paper_models")
    todos = [m for m in MP2.CONJUNTOS["todos"][0]]
    filas = []
    for c, ci, fr in todos:
        d = bandas(modelos=[c], pf=None).dropna(subset=["lat_p50"])
        for (r, pf), g in d.groupby(["ridge", "pf"]):
            filas.append(dict(c=c, marco=fr, ridge=r, pf=float(pf), tmax=g.age.max()))
    A = pd.DataFrame(filas)
    nw = A[A.marco != "WK08-A"]
    chk("main:147", "oldest valid age GMHRF-B and GWG20, all chains and f_PAC",
        rango(nw[nw.marco.isin(["GMHRF-B", "GWG20"])].tmax, "{:.0f}"), "15--45")
    a = A[A.pf == 0.44]
    wk = a[a.marco == "WK08-A"]
    chk("main:147/supl:74", "WK08-A, oldest age (without Foundation)",
        rango(wk[wk.ridge != "Foundation"].tmax, "{:.0f}"), "66")
    chk("main:147/supl:74", "WK08-A, oldest age Foundation",
        rango(wk[wk.ridge == "Foundation"].tmax, "{:.0f}"), "44--46")
    txt = {"GMHRF-B": ["20--22", "23--25", "33--35", "17--18"],
           "GMHRF-A'": ["25--26", "27--31", "38--41", "17--18"],
           "GWG20": ["29--33", "31--35", "38--41", "20--21"]}
    for fr, tt in txt.items():
        for r, t in zip(["Easter", "JFR", "TT", "Foundation"], tt):
            chk("supl:74", f"{fr}, {CADENAS[r]}: oldest valid age",
                rango(a[(a.marco == fr) & (a.ridge == r)].tmax, "{:.0f}"), t)


def sec_comparacion():
    """Fig. 3 (model comparison): reads F5_model_comparison_numbers.csv
    (figures/make_F5_model_comparison.py). Full 68% widths (p84 - p16)."""
    f = pd.read_csv(os.path.join(G13, "F5_model_comparison_numbers.csv")).set_index("cadena")
    for r, t in (("Foundation", "1.1"), ("Easter", "1.7"), ("TT", "2.6"), ("JFR", "2.6")):
        chk("main:206", f"range of medians (median over ages), {CADENAS[r]}",
            f"{f.rango_med_deg[r]:.1f}", t)
    chk("main:206", "median 68% width, range across chains", rango(f.w68_med_deg), "1.8--2.7")
    chk("main:206", "% ages with range <= 68% width, Easter/JFR/Foundation",
        rango(100 * f.frac_rango_le_w68[["Easter", "JFR", "Foundation"]], "{:.0f}"), "61--78")
    chk("main:206", "% ages with range <= 68% width, San Felix",
        f"{100 * f.frac_rango_le_w68['TT']:.0f}", "18")
    for r, t in (("Foundation", "0.6"), ("Easter", "1.3"), ("TT", "1.3"), ("JFR", "1.9")):
        chk("main:206", f"range across circuits with fixed frame, {CADENAS[r]}",
            f"{f.circ_med_deg[r]:.1f}", t)
    chk("main:206", "range across frames with fixed circuit", rango(f.marco_med_deg), "0.7--1.4")


NOMBRE12 = {"Q22_WK08-A": "WK08", "WD25_WK08-A": "WD25", "SG12_WK08-A": "SG12",
            "Q22_GMHRF-B": "GMHRF-B", "WD25_GMHRF-B": "WD25+GMHRF-B",
            "SG12_GMHRF-B": "SG12+GMHRF-B", "Q22_GWG20": "GWG20",
            "WD25_GWG20": "WD25+GWG20", "SG12_GWG20": "SG12+GWG20"}


def sec_depositos():
    """Section 5.2 and Fig. 6: test of step 12 (models/Q22_WK08-A/deposit_track_test_by_model.csv), numbers
    of Fig. 6 (F4_deposit_distance_numbers.csv) and collision-deposit lag
    (outputs/deposit_lag_9models.csv, deposit_lag_9models.py)."""
    t = pd.read_csv(os.path.join(ROOT, "models", "Q22_WK08-A", "deposit_track_test_by_model.csv"))
    t = t[t.modelo.isin(NOMBRE12.values())]
    assert t.modelo.nunique() == 9
    p = t[(t.tipo == "porphyries") & (t.cadena == "Todas")]
    e10 = lambda x: f"{x / 10 ** np.floor(np.log10(x)):.0f}e{int(np.floor(np.log10(x)))}"
    chk("main:261/supl:164", "Mann-Whitney porphyries vs background, all chains (min--max)",
        f"{e10(p.mw_p_menor.min())}--{e10(p.mw_p_menor.max())}", "3e-9--2e-7")
    chk("main:261", "median distance porphyries (km)", rango(p.d_med_dep, "{:.0f}"), "119--216")
    chk("main:261", "median distance background (km)", rango(p.d_med_bg, "{:.0f}"), "275--348")
    chk("main:261", "% porphyries at <= 400 km", rango(100 * p.f400_dep, "{:.0f}"), "79--90")
    chk("main:261", "% background at <= 400 km", rango(100 * p.f400_bg, "{:.0f}"), "54--60")
    q = t[t.tipo == "porphyries"]
    for cad, nt, dt in (("Juan Fernandez", "9", "41--45"), ("San Felix", "3", "13--18"),
                        ("Foundation", "7", "4--13"), ("Easter", "0", "4--8")):
        s = q[q.cadena == cad]
        chk("main:261", f"{cad}: models with p<0.05", f"{int((s.mw_p_menor < 0.05).sum())}", nt)
        chk("main:261", f"{cad}: porphyries assigned", rango(s.n_dep, "{:.0f}"), dt)
    chk("main:261", "Easter: background assigned", rango(q[q.cadena == "Easter"].n_bg, "{:.0f}"), "2--10")
    ep = t[(t.tipo == "epithermal") & (t.cadena == "Todas")]
    print(f"   epithermal, models with p<0.05 (all chains): {(ep.mw_p_menor < 0.05).sum()}/9")
    chk("main:261", "Fisher porphyries farther south, minimum p (all chains and models)",
        f"{q.fisher_p_masS.min():.2f}", "0.29")
    f = pd.read_csv(os.path.join(G13, "F4_deposit_distance_numbers.csv"))
    fp = f[f.clase == "por"]
    chk("main:263", "median of the median distances of porphyries (km)",
        f"{fp.dist_mediana_km.median():.0f}", "165")
    chk("main:263", "median distances of porphyries, range of the nine models (km)",
        rango(fp.dist_mediana_km, "{:.0f}"), "119--216")
    chk("main:263", "% porphyries within the 68% half-width (median of models)",
        f"{100 * fp.frac_dentro_68.median():.0f}", "31")
    chk("main:268", "median 68% half-width at the age of the porphyries (km)",
        f"{fp.semiancho68_mediano_km.median():.0f}", "88")
    L = pd.read_csv(os.path.join(SAL, "deposit_lag_9models.csv"))
    for tipo, ref, tc in (("porphyries", "main:263/supl:168", ["1--5", "3--11", "14--30"]),
                          ("epithermal", "supl:168", ["0--5", "1--7", "9--14"])):
        for lag, tt in zip((1, 2, 5), tc):
            s = L[(L.tipo == tipo) & (L.lag == lag)]
            chk(ref, f"{tipo} that change chain, L={lag} Myr", rango(s.cambian, "{:.0f}"), tt)
    for lag, tt in zip((1, 2, 5), ("18--46", "38--92", "90--171")):
        s = L[(L.tipo == "porphyries") & (L.lag == lag)]
        chk("main:263", f"change of median distance, L={lag} Myr (km)", rango(s.dd_mediana_km, "{:.0f}"), tt)


def sec_descomposicion():
    """Section 4.2 and Fig. S2: 95% width with a single sampled source (variants circuit, wk08 =
    frame of each model, hotspot = hotspot fit); dominant source = the one with the largest width."""
    w = []
    for v in ("circuit", "wk08", "hotspot"):
        d = bandas(v).dropna(subset=["lat_p50"])
        w.append(d[["modelo", "ridge", "age", "w95"]].assign(fuente=v))
    W = pd.concat(w).pivot_table(index=["modelo", "ridge", "age"], columns="fuente", values="w95").dropna()
    W["lider"] = W[["circuit", "wk08", "hotspot"]].idxmax(axis=1)
    fr = W.lider.value_counts(normalize=True) * 100
    chk("main:45/218/281", "% of ages with the frame as dominant source", f"{fr.get('wk08', 0):.0f}", "97")
    chk("main:45/218", "% of ages with the circuit dominant", f"{fr.get('circuit', 0):.0f}", "3")
    chk("main:218", "% of ages with the hotspot fit dominant", f"{fr.get('hotspot', 0):.1f}", "0.4")
    pr = W.groupby("ridge").lider.apply(lambda s: 100 * (s == "wk08").mean())
    for r, t in (("Easter", "96"), ("TT", "96"), ("JFR", "98"), ("Foundation", "98")):
        chk("main:218", f"% frame dominant, {CADENAS[r]}", f"{pr[r]:.0f}", t)
    med = W.groupby("ridge")[["wk08", "circuit", "hotspot"]].median()
    chk("main:218", "median 95% width of the frame, range across chains", rango(med.wk08), "3.2--4.4")
    chk("main:218", "median 95% width of the circuit, range across chains", rango(med.circuit), "1.0--1.5")
    chk("main:218", "median 95% width of the fit, range across chains", rango(med.hotspot), "0.4--1.4")
    circ = dict((c, ci) for c, ci, _ in MP.MODELOS)
    W2 = W.reset_index()
    W2["circ"] = W2.modelo.map(circ)
    print("   % frame dominant per circuit:",
          W2.groupby("circ").lider.apply(lambda s: round(100 * (s == "wk08").mean(), 1)).to_dict())


def dist_km(la1, lo1, la2, lo2):
    p1, p2, l1, l2 = map(np.radians, (la1, la2, lo1, lo2))
    c = np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(l2 - l1)
    return 6371.0 * np.arccos(np.clip(c, -1, 1))


def sec_elipses():
    """Section 4.3 and Fig. 1: <carpeta>/track_ellipses_0.44_hsb_cc_g13.csv (step 04, ellipses every
    10 Myr of track age) and track_vs_ridge_axis (Nazca Ridge)."""
    E = pd.concat([pd.read_csv(os.path.join(ROOT, "models", c, f"track_ellipses_{PF}{TAG}.csv")).assign(modelo=c)
                   for c in MP.CARPETAS])
    chk("main:225", "a95 at the dated ends (tau 0), km", rango(E[E.tau == 0].a95_km, "{:.0f}"), "16--172")
    chk("main:225", "a95 at tau 10 Ma, km", rango(E[E.tau == 10].a95_km, "{:.0f}"), "230--414")
    chk("main:225", "a95 at tau 20-60 Ma, km", rango(E[E.tau.between(20, 60)].a95_km, "{:.0f}"), "117--500")
    chk("main:225", "maximum axial ratio a95/b95", f"{(E.a95_km / E.b95_km).max():.0f}:1", "8:1")
    chk("main:225", "median azimuth of the major axis", f"{(E.az_major_deg % 180).median():.0f}", "84")
    az = E.az_major_deg % 180
    chk("main:225", "interquartile range of the azimuth of the major axis", f"{az.quantile(.25):.0f}--{az.quantile(.75):.0f}", "51--108")

    def sep(g):
        v = g[["lat", "lon"]].values
        return max(dist_km(a[0], a[1], b[0], b[1]) for a in v for b in v)
    S = E.groupby(["ridge", "tau"]).apply(sep, include_groups=False).rename("sep").reset_index()
    chk("main:225", "maximum separation of centres between models, tau 0 (range across chains)",
        rango(S[S.tau == 0].sep.round(-1), "{:.0f}"), "70--230")
    s3 = S[S.tau.between(30, 60)].groupby("ridge").sep.max()
    chk("main:225", "maximum separation of centres, tau 30-60 (range across chains)",
        rango(s3.round(-1), "{:.0f}"), "310--650")
    print("   separation per chain and tau:\n", S.pivot(index="tau", columns="ridge", values="sep").round(0))
    R = pd.read_csv(os.path.join(SAL, "track_vs_ridge_axis_g13.csv"))   # track_vs_ridge_axis.py
    chk("main:225", "Nazca Ridge: track ages compared (Myr)", rango(R.tau, "{:.0f}"), "24--37")
    m = R.groupby("modelo").track_menos_centro_km.median()
    nw = m[m.index.str.contains("WK08-A|GMHRF-B")]
    chk("main:225", "Easter NW of the axis, WK08-A and GMHRF-B (median, km)", rango(nw, "{:.0f}"), "64--77")
    chk("main:225", "Easter relative to the axis, GWG20 (max |median|, km)",
        f"{m[m.index.str.contains('GWG20')].abs().max():.0f}", "18")
    chk("main:225", "median width of the ridge (km, rounded to 10)", f"{round(R.ancho_km.median(), -1):.0f}", "140")


def sec_flatslab():
    """Section 5.1 and Fig. 5: flat_slab_tracks_test.csv (flat_slab_tracks.py: overlap intervals
    of the 68% union, the 95% envelope and all the medians with the latitude strip of each
    segment of wilson2025 Table 3) and the bands for the arrival of JFR."""
    F = pd.read_csv(os.path.join(G13, "flat_slab_tracks_test.csv")).fillna("")
    F["cadena"] = F.cadena.str.replace("é", "e").str.replace("á", "a")
    g = lambda c, s, col: F[(F.cadena == c) & (F.segmento == s)][col].iloc[0]
    chk("main:243", "JFR in Pampean, 68% union (Ma)", g("Juan Fernandez", "Pampean", "solape_68"), "0-19; 45-50")
    chk("main:245", "Easter in Peruvian, 68% union (Ma)", g("Easter", "Peruvian", "solape_68"), "0-22; 39-50")
    chk("main:245", "Easter in Peruvian, all medians (Ma)", g("Easter", "Peruvian", "todas_medianas_dentro"), "4-13; 50")
    chk("main:245", "Foundation in Peruvian, all medians (Ma)", g("Foundation", "Peruvian", "todas_medianas_dentro"), "9-16; 32-42")
    chk("main:247", "San Felix below Altiplano, all medians (Ma)", g("San Felix", "Altiplano", "todas_medianas_dentro"), "19-29; 34-43")
    chk("main:247", "San Felix in Puna, all medians (Ma)", g("San Felix", "Puna", "todas_medianas_dentro"), "0-14")
    chk("main:247", "JFR in Puna, 68% union (Ma)", g("Juan Fernandez", "Puna", "solape_68"), "15-50")
    pay = F[F.segmento == "Payenia"]
    chk("main:247", "Payenia, 68% union of any chain", ";".join(pay.solape_68) or "ninguna", "ninguna")
    print("   Payenia, 95% envelope:", dict(zip(pay.cadena, pay.solape_95)))
    d = bandas().dropna(subset=["lat_p50"])
    j = d[d.ridge == "JFR"]
    lleg = j[j.lat_p50 <= -31].groupby("modelo").age.max()
    chk("main:243", "JFR: age at which each median reaches <= 31 S (range across models)",
        rango(lleg, "{:.0f}"), "12--15")
    s = pd.concat([j[(j.modelo == m) & (j.age <= a)] for m, a in lleg.items()])
    chk("main:243", "JFR: latitude of the medians after arrival", rango(-s.lat_p50, "{:.0f}"), "31--33")


def ajustes(carpetas=None):
    """Hotspot fits (step 03b): <carpeta>/hotspot_fit.json."""
    import json
    filas = []
    for c in (carpetas or MP.CARPETAS):
        d = json.load(open(os.path.join(ROOT, "models", c, "hotspot_fit.json"), encoding="utf-8"))
        p0 = 1 if d["meta"].get("present_obs") else 0
        for r, x in d["tracks"].items():
            m = x["models"]
            bic = (m["deriva"]["bic"] < m["fijo"]["bic"]) if "deriva" in m else False
            filas.append(dict(modelo=c, ridge=r, n_used=x["n_used"] - p0, n_input=x["n_input"] - p0,
                              deriva=x["model"] == "deriva", bic_deriva=bic, w_t=x["w_t"], w_x=x["w_x"],
                              shift=x["seed_shift_km"], from_jfr=x["w_from_jfr"],
                              rms_al=x["rms_along_km"], rms_ac=x["rms_across_km"]))
    return pd.DataFrame(filas)


def sec_tabla2():
    """Table 2 (hotspot fits in the nine models), from hotspot_fit.json."""
    A = ajustes()
    for r, nombre, tn, twt, twx, tdr in (("Easter", "Easter", "27--29/30", "1.21--2.59", "63--85", "9/9"),
                                         ("JFR", "Juan Fernandez", "15--16/16", "0.24--0.44", "3--23", "9/9"),
                                         ("TT", "San Felix", "2/2", None, None, "0/9"),
                                         ("Foundation", "Foundation", "19/19", "0.48--0.61", "39--59", "6/9")):
        s = A[A.ridge == r]
        chk("main:107-110", f"Table 2, {nombre}: n", f"{rango(s.n_used, '{:.0f}')}/{rango(s.n_input, '{:.0f}')}", tn)
        if twt:
            chk("main:107-110", f"Table 2, {nombre}: w_t", rango(s.w_t, "{:.2f}"), twt)
            chk("main:107-110", f"Table 2, {nombre}: w_x", rango(s.w_x, "{:.0f}"), twx)
        chk("main:107-110", f"Table 2, {nombre}: models with drifting hotspot", f"{s.deriva.sum()}/9", tdr)
        if r == "Easter":
            print(f"   Easter: drift preferred by the BIC in {s.bic_deriva.sum()}/9")
    chk("main:218", "Easter: w_x (km)", rango(A[A.ridge == "Easter"].w_x, "{:.0f}"), "63--85")
    ea = A[A.ridge == "Easter"]
    for fr, patron, t in (("WK08-A", "WK08-A$", "229--232"), ("GMHRF-B", "GMHRF-B", "126--148"),
                          ("GWG20", "GWG20", "280--296")):
        chk("supl:104", f"{fr}: distance of the fitted Easter hotspot to the reference position (km)",
            rango(ea[ea.modelo.str.contains(patron)]["shift"], "{:.0f}"), t)


def sec_velocidad():
    """Section 4.2 and Methods: 1 sigma of the MAGNITUDE of the convergence velocity from the circuit
    Monte Carlo (step 09, <carpeta>/vc_sigma_circuit_mc.csv), median over the margin nodes
    and the ages in [a0, a1); ratio MC / analytic of Quiero et al. 2022 (step 05,
    models/Q22_WK08-A/vc_sigma_quiero.csv) as in the validation of step 09."""
    med = {}
    for c in ("Q22_WK08-A", "WD25_WK08-A", "SG12_WK08-A"):
        o = pd.read_csv(os.path.join(ROOT, "models", c, "vc_sigma_circuit_mc.csv"))
        med[c] = [o[(o.age >= a) & (o.age < b)].sigma_mag.median() for a, b in ((1, 7), (7, 26))]
        if c == "Q22_WK08-A":
            q = pd.read_csv(os.path.join(ROOT, "models", "Q22_WK08-A", "vc_sigma_quiero.csv"))
            qq = q.groupby(q.age_mid.round(0)).sigma_mag.median()
            mm = o[o.age <= 28].groupby("age").sigma_mag.median()
            com = sorted(set(qq.index) & set(mm.index))
            chk("main:143", "median ratio sigma MC / analytic (Q22, 0-28 Ma)",
                f"{np.median([mm[a] / qq[a] for a in com]):.1f}", "1.3")
    m = pd.DataFrame(med, index=["1-7", "7-26"]).T
    print("   median sigma_mag (cm/yr):\n", m.round(2))
    chk("main:220", "1 sigma of the velocity, 1-7 Ma, range across circuits (cm/yr)", rango(m["1-7"]), "0.3--0.5")
    chk("main:220", "1 sigma of the velocity, 7-26 Ma, range across circuits (cm/yr)", rango(m["7-26"]), "0.7--2.2")
    chk("main:220", "circuit with the smallest value, 7-26 Ma", m["7-26"].idxmin(), "SG12_WK08-A")


S9 = os.path.join(SAL, "sens9")


def deriva_kmmyr(carpeta, fn):
    """Drift rate of each hotspot (|theta[2:4]|, km/Myr) in a 03b fit."""
    import json
    d = json.load(open(os.path.join(ROOT, "models", carpeta, fn), encoding="utf-8"))
    return {r: (np.hypot(x["theta"][2], x["theta"][3]) if x["model"] == "deriva" else 0.0)
            for r, x in d["tracks"].items()}


def sec_sensibilidad():
    """Section S4 of the supplement and Limitations: sensitivity tests of the nine models
    (run_sens9.sh, sens9_summary.py, sens9_deposits.py, sens9_wagenaar.py; outputs/sens9/)."""
    lag = {L: pd.DataFrame([deriva_kmmyr(c, f"hotspot_fit_s9_lag{L}.json")
                            for c in MP.CARPETAS]) for L in (0, 1, 2)}
    for L, t in ((0, "17--30"), (1, "10--19"), (2, "7--22")):
        chk("supl:152", f"JFR, drift with L={L} Myr (km/Myr)", rango(lag[L].JFR, "{:.0f}"), t)
    N = pd.read_csv(os.path.join(S9, "nominal_lag_york_9models.csv"))
    for L, t in ((1, "481"), (2, "602")):
        c = f"d_lag{L}_km"
        k = N[c].abs().idxmax()
        chk("supl:152", f"maximum change of the nominal collision, L={L} (km)", f"{N[c].abs().max():.0f}", t)
        chk("supl:152", f"chain and age of the maximum, L={L}", f"{N.ridge[k]} {N.age[k]:.0f}", "JFR 14")
    for r, t in (("Easter", "211"), ("TT", "227")):
        chk("supl:152", f"{CADENAS[r]}: maximum change L=1 (km)", f"{N[N.ridge == r].d_lag1_km.abs().max():.0f}", t)
    fijo1 = [c for c in MP.CARPETAS if json.load(open(os.path.join(ROOT, "models", c, "hotspot_fit_s9_lag1.json"),
             encoding="utf-8"))["tracks"]["Easter"]["model"] == "fijo"]
    ea = N[N.ridge == "Easter"].groupby("modelo").d_lag1_km.apply(lambda x: x.abs().max())
    chk("supl:152", "Easter L=1: maximum change in the models with fixed hotspot (km)", rango(ea[ea.index.isin(fijo1)], "{:.0f}"),
        "141--211")
    chk("supl:152", "Easter L=1: maximum change in the other models (km, <=)", f"{ea[~ea.index.isin(fijo1)].max():.0f}", "50")
    v = N[N.ridge.isin(["Easter", "JFR", "TT"])].dropna(subset=["s68_km"])
    chk("supl:152", "valid model-ages Easter/JFR/San Felix", f"{len(v)}", "1164")
    for L, t, tj in ((1, "72", "43"), (2, "162", "134")):
        ex = v[v[f"d_lag{L}_km"].abs() > v.s68_km]
        chk("supl:152", f"model-ages with change > 68% half-width, L={L}", f"{len(ex)}", t)
        chk("supl:152", f"of which in JFR, L={L}", f"{(ex.ridge == 'JFR').sum()}", tj)
    y = pd.DataFrame([deriva_kmmyr(c, "hotspot_fit_s9_noYORK2.json") for c in MP.CARPETAS])
    fin = pd.DataFrame([deriva_kmmyr(c, "hotspot_fit.json") for c in MP.CARPETAS])
    chk("supl:154", "Easter: change of drift without YORK2 relative to the final fit (km/Myr, max)",
        f"{(y.Easter - fin.Easter).abs().max():.1f}", "1.6")
    e = N[N.ridge == "Easter"].dropna(subset=["d_noYORK2_km"])
    for fr, patron in (("WK08-A", "WK08-A$"), ("GWG20", "GWG20")):
        nb = e[e.modelo.str.contains(patron)]
        print(f"   without YORK2, {fr}: max {nb.d_noYORK2_km.abs().max():.1f} km; before 27 Ma "
              f"{nb[nb.age < 27].d_noYORK2_km.abs().max():.1f} km")
    nb = e[e.modelo.str.contains("WK08-A$")]
    chk("supl:154", "Easter without YORK2, WK08-A, before 27 Ma (km, max)",
        f"{nb[nb.age < 27].d_noYORK2_km.abs().max():.0f}", "6")
    k = nb.d_noYORK2_km.abs().idxmax()
    chk("supl:154", "Easter without YORK2, WK08-A (km, max) and age", f"{nb.d_noYORK2_km.abs().max():.0f} a {nb.age[k]:.0f}",
        "10 a 55")
    for fr, patron, t in (("GWG20", "GWG20", "152 a 29"), ("GMHRF-B", "GMHRF", "69 a 14")):
        gb = e[e.modelo.str.contains(patron)]
        mx = gb.d_noYORK2_km.abs().max()
        ed = gb[gb.d_noYORK2_km.abs().round() == round(mx)].age   # ages at which the rounded maximum is reached
        chk("supl:154", f"Easter without YORK2, {fr} (km, max) and ages", f"{mx:.0f} a {rango(ed, '{:.0f}')}", t)
    chk("supl:154", "model-ages without YORK2 with change > 68% half-width",
        f"{(e.d_noYORK2_km.abs() > e.s68_km).sum()}", "1")
    I = pd.read_csv(os.path.join(S9, "interp_9models.csv"))
    q = 100 * (I.cociente - 1)
    chk("supl:160", "interpolation: median increase of the 68% half-width (%)", f"{q.median():.0f}", "18")
    chk("supl:160", "interpolation: p5--p95 (%)", f"{q.quantile(.05):.0f}--{q.quantile(.95):.0f}", "4--33")
    chk("supl:160", "interpolation: maximum (%)", f"{q.max():.0f}", "44")
    chk("supl:160", "interpolation: absolute maximum (km)", f"{(I.s68_km_norm - I.s68_km_sin).max():.0f}", "98")
    W = pd.read_csv(os.path.join(S9, "wagenaar_9models.csv"))
    for r, t, td in (("Easter", "2.1", "55"), ("JFR", "0.6", "94"), ("TT", "1.0", "82"), ("Foundation", "0.8", "84")):
        s = W[W.ridge == r]
        chk("main:275/supl:164", f"Wagenaar, {CADENAS[r]}: median difference (degrees)", f"{s.dif_deg.abs().median():.1f}", t)
        chk("supl:164", f"Wagenaar, {CADENAS[r]}: % ages within the 68% envelope", f"{100 * s.dentro68.mean():.0f}", td)
    s = W[W.ridge == "Easter"]
    k = s.dif_deg.abs().idxmax()
    chk("supl:164", "Wagenaar, Easter: maximum (degrees) and age", f"{s.dif_deg.abs().max():.1f} a {s.age[k]:.0f}", "10.0 a 64")
    T = pd.read_csv(os.path.join(S9, "test_porphyries_variants.csv"))
    print("   deposit variants:", T.variante.unique().tolist())


def sec_utrecht():
    """Andean shortening: andean_shortening_utrecht.py (restored positions,
    outputs/andean_shortening/position_summary.csv), step 07 with the blocks (run_shortening_step07.sh)
    and test of step 12 per variant (sens9_deposits.py, outputs/sens9/test_porphyries_variants.csv)."""
    U = os.path.join(SAL, "andean_shortening")
    R = pd.read_csv(os.path.join(U, "position_summary.csv")).set_index("conjunto")
    chk("main:275/supl:166", "porphyries: median displacement Utrecht - Arriagada (km)",
        f"{R.dif_posicion_km_med['porphyries']:.0f}", "104")
    chk("supl:166", "porphyries: maximum displacement (km)", f"{R.dif_posicion_km_max['porphyries']:.0f}", "178")
    cam = {"Utrecht": [], "Utrecht, fosa propia": []}
    for c in MP.CARPETAS:
        base = pd.read_csv(os.path.join(ROOT, "models", c, "deposit_tracks_porphyries.csv")).HS4.values
        for v, sfx in (("Utrecht", ""), ("Utrecht, fosa propia", "_own_trench")):
            u = pd.read_csv(os.path.join(U, c + sfx, "deposit_tracks_porphyries.csv")).HS4.values
            cam[v].append(int((u != base).sum()))
    chk("supl:166", "porphyries that change chain, trench of the blocks", rango(cam["Utrecht"], "{:.0f}"), "1--6")
    chk("supl:166", "porphyries that change chain, own trench", rango(cam["Utrecht, fosa propia"], "{:.0f}"), "1--3")
    chk("main:275", "porphyries that change chain, both trenches",
        rango(cam["Utrecht"] + cam["Utrecht, fosa propia"], "{:.0f}"), "1--6")
    T = pd.read_csv(os.path.join(S9, "test_porphyries_variants.csv"))
    T = T[T.modelo.isin(MP.CARPETAS)]
    for v in ("Utrecht", "Utrecht, fosa propia"):
        s = T[T.variante == v]
        n = s[s.cadena != "Todas"].groupby("cadena").mw_p_menor.apply(lambda x: int((x < 0.05).sum()))
        print(f"   {v}: models with p<0.05 per chain {n.to_dict()}; "
              f"all chains p max {s[s.cadena == 'Todas'].mw_p_menor.max():.1e}")
    u = T[T.variante.str.startswith("Utrecht")]
    chk("supl:166", "Utrecht: Mann-Whitney all chains and JFR, maximum p < 1e-5",
        str(bool(u[u.cadena.isin(["Todas", "Juan Fernandez"])].mw_p_menor.max() < 1e-5)), "True")
    for v, t in (("Utrecht", "1"), ("Utrecht, fosa propia", "1")):
        s = T[(T.variante == v) & (T.cadena == "San Felix")]
        chk("supl:166", f"San Felix, models with p<0.05, {v}", f"{(s.mw_p_menor < 0.05).sum()}", t)
    for v in ("Utrecht", "Utrecht, fosa propia"):
        s = T[(T.variante == v) & (T.cadena == "Foundation")]
        chk("supl:166", f"Foundation, models with p<0.05, {v}", f"{(s.mw_p_menor < 0.05).sum()}", "9")
    s = u[u.cadena == "Todas"]
    chk("supl:166", "Utrecht: % porphyries at <= 400 km", rango(100 * s.f400_dep, "{:.0f}"), "77--92")
    chk("supl:166", "Utrecht: % background at <= 400 km", rango(100 * s.f400_bg, "{:.0f}"), "51--62")
    w = pd.read_csv(os.path.join(S9, "test_porphyries_variants.csv"))
    w = w[(w.variante == "WVH25") & (w.cadena == "Todas")]
    chk("supl:164", "Wagenaar: Mann-Whitney all chains", f"{w.mw_p_menor.iloc[0]:.0e}", "5e-10")


def sec_geocrono():
    """Section S3: decay constant (decay_constant_test.py, run here) and
    mean of the rotated samples (rotated_mean.py -> outputs/rotated_mean_9models.txt)."""
    import subprocess
    o = subprocess.run([sys.executable, os.path.join(HERE, "decay_constant_test.py")],
                       capture_output=True, text=True, encoding="utf-8").stdout
    tab = {float(m.group(1)): (float(m.group(2)), float(m.group(3)), float(m.group(5)))
           for m in re.finditer(r"^\s*([\d.]+) \|\s+(-?[\d.]+) \|\s+(-?[\d.]+) \|\s+([\d.]+) \|\s+([\d.]+)", o, re.M)}
    chk("supl:96", "maximum age change for 1 sigma of lambda (kyr)",
        f"{max(abs(v[0]) for v in tab.values()):.0f}", "1")
    for a, t in ((0.4, "0.015"), (10.0, "0.010")):
        chk("supl:96", f"relative change at {a} Ma (%)", f"{tab[a][1]:.3f}", t)
    chk("supl:96", "maximum widening from the calibration term (%)",
        f"{100 * (max(v[2] for v in tab.values()) - 1):.1f}", "1.7")
    pw = [abs(float(x)) for x in re.findall(r"w_t = [\d.]+ Myr \((-?[\d.]+) % of w_t\)", o)]
    chk("supl:96", "maximum effect as % of w_t (bound)", str(max(pw) < 0.25), "True")
    txt = open(os.path.join(SAL, "rotated_mean_9models.txt"), encoding="utf-8").read()
    bloques = re.split(r"^=== ", txt, flags=re.M)[1:]
    filas = []
    for b in bloques:
        c = b.split("\n")[0].strip()
        for r in ("Easter", "JFR"):
            m = re.search(rf"^{r}\s+con pesos .* a\s+([\d.]+) km", b, re.M)
            s = re.search(rf"^{r}\s+con pesos.*?\n.*?\n\s+dispersion .*?mediana ([\d.]+) km", b, re.M | re.S)
            if m:
                filas.append(dict(modelo=c, ridge=r, dist=float(m.group(1)), disp=float(s.group(1))))
    P = pd.DataFrame(filas)
    P = P[P.modelo.isin(MP.CARPETAS)]
    g = P.modelo.str.contains("GMHRF")
    P = P.reset_index(drop=True)
    g = P.modelo.str.contains("GMHRF")
    modelos = {r: {json.load(open(os.path.join(ROOT, "models", c, "hotspot_fit.json"), encoding="utf-8"))["tracks"][r]["model"]
                   for c in MP.CARPETAS} for r in ("Easter", "JFR")}
    chk("supl:106", "model chosen by the BIC in the nine (Easter; JFR)",
        f"{'/'.join(sorted(modelos['Easter']))}; {'/'.join(sorted(modelos['JFR']))}", "deriva; deriva")
    s = P[P.ridge == "JFR"]
    q = s.modelo.isin(["Q22_WK08-A", "SG12_WK08-A"])
    chk("supl:106", "JFR: rotated mean to the fixed hotspot, except Q22 and SG12 + WK08-A (km)", rango(s[~q].dist, "{:.0f}"), "4--45")
    chk("supl:106", "JFR: same, Q22 and SG12 + WK08-A (km)", rango(s[q].dist, "{:.0f}"), "159")
    s = P[P.ridge == "Easter"]
    g = s.modelo.str.contains("GMHRF")
    chk("supl:106", "Easter: rotated mean to the fixed hotspot, WK08-A and GWG20 (km)", rango(s[~g].dist, "{:.0f}"), "66--73")
    chk("supl:106", "Easter: same, GMHRF-B (km)", rango(s[g].dist, "{:.0f}"), "181--201")
    chk("supl:106", "JFR: median dispersion of the rotated samples (km)", rango(P[P.ridge == "JFR"].disp, "{:.0f}"), "42--167")
    chk("supl:106", "Easter: same (km)", rango(P[P.ridge == "Easter"].disp, "{:.0f}"), "141--253")
    chk("supl:106", "in each model and chain, distance < median dispersion", str(bool((P.dist < P.disp).all())), "True")


# ------------------------------------------------------------------ Suppl. S2: circuit and frames
def sec_circuito():
    """Own Hellinger fits (hellinger_circuit.py) and unit of the GWG20 covariance (Table S6)."""
    H = os.path.join(SAL, "frame_coverage", "hellinger")
    a1 = json.load(open(os.path.join(H, "fit_801_71.07.json"), encoding="utf-8"))
    a2 = json.load(open(os.path.join(H, "fit_802_71.1.json"), encoding="utf-8"))
    chk("supl:58", "AUS-EANT 71.1 Ma: picks and segments", f"{a1['n']} picks, {a1['s']} segmentos", "48 picks, 6 segmentos")
    chk("supl:58", "AUS-EANT: kappa", f"{a1['kappa_hat']:.2f}", "0.16")
    chk("supl:58", "AUS-EANT: .rot pole inside the 95% region", str(a1["rot_dentro_95"]), "True")
    chk("supl:58", "EANT-AFR 71.1 Ma: picks and segments", f"{a2['n']} picks, {a2['s']} segmentos", "49 picks, 5 segmentos")
    chk("supl:58", "EANT-AFR: .rot pole outside the 95% region", str(not a2["rot_dentro_95"]), "True")
    chk("supl:58", "EANT-AFR: distance from the .rot pole to the fit (degrees)", f"{a2['dif_rot_grados']:.2f}", "0.39")
    # inflation: the 95% region reaches the .rot pole when chi2 / f^2 = threshold
    chk("supl:58", "EANT-AFR: factor on the standard deviations",
        f"{np.sqrt(a2['chi2_rot_en_region'] / a2['umbral_95']):.2f}", "2.39")
    V = pd.read_csv(os.path.join(SAL, "frame_coverage", "variants_802.csv"))   # variants_802.py
    V = V.assign(lat_p50=V.dlat_km)   # recortar() sets NaN in the lat_* columns outside the frame coverage
    V = pd.concat([MP.recortar(V[V.modelo == c], c) for c in MP.CARPETAS]).dropna(subset=["lat_p50"])
    chk("supl:58", "fitted pole instead of the .rot one: maximum change of the nominal collision, 9 models and 3 f_PAC (km)",
        f"{V.dlat_km.abs().max():.0f}", "12")
    S6 = os.path.join(ROOT, "external", "gaastra2022", "Supporting_Information",
                      "2021tc006772-sup-0003-table si-s02.csv")
    g = pd.read_csv(S6)
    g.columns = ["age", "lat", "lon", "ang", "a", "b", "c", "d", "e", "f", "k"]
    # Cov = (1/kappa) [[a,b,c],[b,d,e],[c,e,f]] x 1e-5 (caption of Table S6)
    s = np.array([np.degrees(np.sqrt(np.linalg.eigvalsh(
        1e-5 / r.k * np.array([[r.a, r.b, r.c], [r.b, r.d, r.e], [r.c, r.e, r.f]])).max())) for r in g.itertuples()])
    chk("supl:72", "GWG20: maximum sigma read as rad^2 (degrees)", rango(s), "0.3--1.7")
    chk("supl:72", "GWG20: maximum sigma read as rad^2 Ma^-2 (x t, degrees, rounded to 10)",
        f"{round((s * g.age.values).max(), -1):.0f}", "130")


# ------------------------------------------------------------------ Suppl. S4: GMHRF-A'
def sec_gmhrfA():
    """GMHRF-A' versus GMHRF-B in the three circuits (same clipping by frame coverage)."""
    PARES = (("Q22_GMHRF-A", "Q22_GMHRF-B"),
             ("WD25_GMHRF-A", "WD25_GMHRF-B"),
             ("SG12_GMHRF-A", "SG12_GMHRF-B"))
    A = ajustes([c for p in PARES for c in p]).set_index(["modelo", "ridge"])
    e = {c: A.loc[(c, "Easter")] for p in PARES for c in p}
    q, s = PARES[0], PARES[2]
    chk("supl:162", "Easter rms along, B -> A' (Q22 and SG12, km)",
        f"{e[q[1]].rms_al:.0f} -> {e[q[0]].rms_al:.0f}; {e[s[1]].rms_al:.0f} -> {e[s[0]].rms_al:.0f}", "211 -> 225; 211 -> 225")
    chk("supl:162", "Easter: distance of the fitted hotspot to the reference position, B -> A' (Q22, km)",
        f"{e[q[1]]['shift']:.0f} -> {e[q[0]]['shift']:.0f}", "148 -> 235")
    dif, raz = [], []
    for a, b in PARES:
        m = bandas(modelos=[a]).merge(bandas(modelos=[b]), on=["ridge", "age"], suffixes=("_a", "_b"))
        m = m.dropna(subset=["lat_p50_a", "lat_p50_b"])
        m["d_km"] = (m.lat_p50_b - m.lat_p50_a).abs() * 111.195
        m["r68"] = (m.lat_p84_b - m.lat_p16_b) / (m.lat_p84_a - m.lat_p16_a)
        dif.append(m.assign(par=a))
    # Foundation has no motion in Doubrovine et al. (2012): A' and B coincide
    D = pd.concat(dif, ignore_index=True)
    D = D[D.ridge != "Foundation"]
    g = D.groupby(["par", "ridge"])
    med = g.d_km.median()
    chk("supl:162", "|difference of medians| A' - B, median per circuit and chain (km)", rango(med, "{:.0f}"), "15--198")
    mx = g.d_km.max()
    chk("supl:162", "maximum |difference of medians|, range across circuits and chains (km)", rango(mx, "{:.0f}"), "347--702")
    ed = D.loc[D.groupby(["par", "ridge"]).d_km.idxmax()]
    chk("supl:162", "age of the maximum (Ma)", rango(ed.age, "{:.0f}"), "14--16")
    print("   maxima:", ed[["par", "ridge", "age", "d_km"]].round(0).to_string(index=False))
    r68 = g.r68.median()
    chk("supl:162", "ratio of 68% widths B/A', median per circuit and chain", rango(r68, "{:.2f}"), "0.98--1.02")
    E = [pd.read_csv(os.path.join(ROOT, "models", c, f"track_ellipses_0.44{TAG}.csv")) for p in PARES for c in p]
    q = [E[i].merge(E[i + 1], on=["ridge", "tau"], suffixes=("_a", "_b")) for i in (0, 2, 4)]
    q = pd.concat(q, ignore_index=True)
    chk("supl:162", "track ellipses: ratio of 95% semi-axes B/A' (median per chain, range)",
        rango(pd.concat([(q.a95_km_b / q.a95_km_a), (q.b95_km_b / q.b95_km_a)]).groupby(
            pd.concat([q.ridge, q.ridge]).values).median(), "{:.2f}"), "0.99--1.00")
    print("   valid ages per pair and chain:", g.size().to_dict())


# ------------------------------------------------------------------ Suppl. S3 and Methods: input data
def sec_entradas():
    """Easter ages recomputed from ERDA (ray2012_erda.py, ray2012_erda_ages.py), San Felix and
    counts of deposits and background samples (step 12)."""
    R = pd.read_csv(os.path.join(SAL, "ray2012_erda_samples.csv"))
    d = (R.rec_age - R.plateau_age).abs().dropna()
    chk("supl:84", "weighted mean without dispersion = ArArCALC within 0.005 Ma (samples)",
        f"{(d <= 0.005).sum()}/{len(d)}", "34/35")
    E = pd.read_csv(os.path.join(ROOT, "models", "Q22_WK08-A", "data", "ages_easter_ray2012_recal.csv"), comment="#")
    chk("supl:84", "irradiation batches", f"{E.irradiation.nunique()}", "8")
    chk("supl:84", "J per sample (%1s)", rango(E.J_pct1s, "{:.2f}"), "0.13--0.38")
    chk("supl:84", "DRFT-126-1: published error (Myr)", f"{E.set_index('sample').loc['DRFT-126-1', 'err2s_ma']:.2f}", "8.06")
    chk("supl:84", "recalibration FCT-3 -> FCs (%)", f"+{red(100 * (E.age_A / E.age_ma - 1).median(), 2)}", "+0.61")
    T = pd.read_csv(os.path.join(ROOT, "models", "Q22_WK08-A", "data", "ages_taltal_recal.csv"), comment="#").set_index("sample")
    chk("supl:86", "San Ambrosio: age and error used (Ma, 2s)",
        f"{T.loc['SAMBROSIO293', 'age_ma']:.2f} +/- {T.loc['SAMBROSIO293', 'err2s_ma']:.2f}", "2.93 +/- 0.30")
    chk("supl:86", "position error of the islands (km)", rango(T.pos_err_km), "3.9")
    J = pd.read_csv(os.path.join(ROOT, "models", "Q22_WK08-A", "data", "ages_jfr_lara2018.csv"), comment="#")
    J = J[J.get("use_fit", 1) == 1] if "use_fit" in J else J
    chk("supl:152", "JFR: age range (Ma)", rango(J.age_ma), "0.6--9.3")
    D = pd.read_csv(os.path.join(ROOT, "models", "Q22_WK08-A", "deposit_track_test_by_model.csv"))
    t = D[D.cadena == "Todas"]
    for tipo, txt in (("porphyries", "71"), ("epithermal", "24")):
        chk("main:154", f"{tipo} with assigned track (all models)", rango(t[t.tipo == tipo].n_dep, "{:.0f}"), txt)
    chk("main:154", "background samples", rango(t.n_bg, "{:.0f}"), "1150")


def sec_bic():
    """Suppl. S3: Easter BIC in the nine models. The 03b fit does not force the drift
    (no PAPER_HS_FORCE) and has no 0 Ma observation; the unforced fit (tag _s9_noFORCE,
    run_bic_noforce.sh) is identical to the one in hotspot_fit.json."""
    import json
    mod, dif = [], []
    for c in MP.CARPETAS:
        li = json.load(open(os.path.join(ROOT, "models", c, "hotspot_fit_s9_noFORCE.json"), encoding="utf-8"))["tracks"]["Easter"]
        mod.append(li["model"]); dif.append(li["models"]["fijo"]["bic"] - li["models"]["deriva"]["bic"])
    chk("supl:104", "Easter unforced: models with drift preferred", str(mod.count("deriva")), "9")
    chk("supl:104", "Easter: BIC fixed - drift (range, nine models)", rango(pd.Series(dif), "{:.1f}"), "3.0--28.1")


def sec_conP0():
    """Seeded variant (fit_hotspots.py --conP0) versus the base case, from sens_conP0.py.
    Not in SECCIONES: it is not reported in the text. Kept to redo the comparison."""
    A = pd.read_csv(os.path.join(SAL, "conP0", "ajustes_conP0.csv"))
    A = A[A.modelo.isin(MP.CARPETAS)]
    for r, t in (("Easter", "43--118"), ("JFR", "5--49"), ("TT", "31--87"), ("Foundation", "6--42")):
        chk("supl:154", f"seeded: hotspot displacement, {CADENAS[r]} (km)",
            rango(A[A.ridge == r].h0_km, "{:.0f}"), t)
    cam = A[A.modelo_con != A.modelo_sin]
    chk("supl:154", "seeded: fixed/drift changes (model, chain, seeded)",
        "; ".join(f"{m} {r} {c}" for m, r, c in zip(cam.modelo, cam.ridge, cam.modelo_con)),
        "WD25_WK08-A Easter fijo")
    ex = A[A.excl_con.fillna("") != A.excl_sin.fillna("")]
    chk("supl:154", "seeded: chain-models with different exclusions",
        "; ".join(f"{r} {m}" for m, r in zip(ex.modelo, ex.ridge)),
        "JFR Q22_GWG20; Easter WD25_WK08-A; JFR SG12_GWG20")
    D = pd.read_csv(os.path.join(SAL, "conP0", "nominal_conP0.csv"))
    D = D[(D.pf == 0.44) & D.modelo.isin(MP.CARPETAS)].assign(ad=lambda x: x.dif_km.abs())
    g = D.groupby("ridge")
    chk("supl:154", "seeded: median change Easter/JFR/San Felix/Foundation (km)",
        ", ".join(f"{g.ad.median()[r]:.0f}" for r in ("Easter", "JFR", "TT", "Foundation")), "14, 7, 42, 6")
    k = g.ad.idxmax()
    chk("supl:154", "seeded: maximum change and age Easter/JFR/San Felix/Foundation",
        "; ".join(f"{D.ad[k[r]]:.0f} a {D.age[k[r]]:.0f}" for r in ("Easter", "JFR", "TT", "Foundation")),
        "171 a 29; 100 a 32; 228 a 12; 68 a 21")
    e = D[D.ad > D.s68_km]
    chk("supl:154", "seeded: model-ages above the 68% half-width (total; San Felix; Easter; n)",
        f"{len(e)}; {(e.ridge == 'TT').sum()}; {(e.ridge == 'Easter').sum()}; {len(D)}", "33; 28; 5; 1418")


def sec_selkirk():
    """Suppl. S4, Selkirk microplate: run_sens_selkirk.sh and sens_selkirk.py (outputs/selkirk/)."""
    SK = os.path.join(SAL, "selkirk")
    R = pd.read_csv(os.path.join(SK, "residuals_foundation_zone.csv"))
    m = R[R.zona == "microplaca"]
    chk("supl:S4 Selkirk", "samples of the fit on the microplate", ", ".join(sorted(m["sample"].unique())),
        "SO100-33DS-1, SO100-38DS-1")
    chk("supl:S4 Selkirk", "minimum age of the samples west of the failed rift (Ma)",
        f"{R[R.zona == 'oeste'].age_ma.min():.0f}", "16")
    med = m.groupby("sample")[["z_along", "z_across"]].median().abs().max().max()
    chk("supl:S4 Selkirk", "normalized residuals: maximum median < 0.6 and maximum <= 2",
        str(bool(med < 0.6 and m[["z_along", "z_across"]].abs().max().max() <= 2)), "True")
    A = pd.read_csv(os.path.join(SK, "fits_selkirk.csv"))
    chk("supl:S4 Selkirk", "models with fixed/drift change", str((A.modelo_var != A.modelo_final).sum()), "0")
    chk("supl:S4 Selkirk", "maximum hotspot displacement (km)", f"{A.h0_km.max():.0f}", "13")
    chk("supl:S4 Selkirk", "maximum change of the drift (km/Myr)",
        f"{(A.deriva_var - A.deriva_final).abs().max():.1f}", "0.6")
    N = pd.read_csv(os.path.join(SK, "nominal_selkirk.csv"))
    g = N.groupby("var").ad.median()
    chk("supl:S4 Selkirk", "median change of the nominal collision (km, range of variants)", rango(g, "{:.0f}"), "4--8")
    e = N[N.ad > N.s68_km]
    chk("supl:S4 Selkirk", "model-ages above the 68% half-width (per variant) and total",
        f"{e.groupby('var').size().max()} de {N.groupby('var').size().max()}", "2 de 414")
    chk("supl:S4 Selkirk", "models and ages of those cases",
        "; ".join(sorted({f'{a} {b:.0f}' for a, b in zip(e.modelo, e.age)})),
        "Q22_GWG20 21; WD25_GWG20 20")


def sec_tswitch():
    """Suppl. S3 (Foundation): age of the Challenger FZ crossing (t_switch of 03b, hotspot_fit.json;
    trace in models/Q22_WK08-A/data/fz_challenger_gsfml.csv)."""
    T = [json.load(open(os.path.join(ROOT, "models", c, "hotspot_fit.json"), encoding="utf-8"))["tracks"]["Foundation"]
         for c in MP.CARPETAS]
    chk("supl:90", "t_switch of Foundation in the nine models (Ma)", rango([t["t_switch"] for t in T], "{:.1f}"),
        "23.2--26.0")
    chk("supl:90", "t_switch rule (FZ crossing in all)", "/".join(sorted({t["t_switch_rule"] for t in T})), "fz")


def sec_jdoble():
    """Suppl. S3 (systematic terms): common 0.5 % J in Lara (2018), whose errors already include it."""
    d = pd.read_csv(os.path.join(ROOT, "models", "Q22_WK08-A", "data", "ages_jfr_lara2018.csv"), comment="#")
    J = 0.005 * d.age_ma.max()
    W = min(json.load(open(os.path.join(ROOT, "models", c, "hotspot_fit.json"), encoding="utf-8"))["tracks"]["JFR"]["w_t"]
            for c in MP.CARPETAS)
    chk("supl:94", "maximum J 0.5 % (1s) in JFR (Myr)", f"{J:.2f}", "0.05")
    chk("supl:94", "minimum w_t of JFR in the nine models (Myr)", f"{W:.2f}", "0.24")
    chk("supl:94", "J max / w_t min < 1/5", str(J / W < 0.2), "True")


SECCIONES = {"bandas": sec_bandas, "tabla1": sec_tabla1, "frame_coverage": sec_frame_coverage, "elipses": sec_elipses,
             "selkirk": sec_selkirk, "tswitch": sec_tswitch, "jdoble": sec_jdoble,
             "geocrono": sec_geocrono,
             "utrecht": sec_utrecht,
             "sensibilidad": sec_sensibilidad,
             "velocidad": sec_velocidad,
             "tabla2": sec_tabla2,
             "flatslab": sec_flatslab,
             "comparacion": sec_comparacion, "depositos": sec_depositos,
             "descomposicion": sec_descomposicion,
             "circuito": sec_circuito, "gmhrfA": sec_gmhrfA, "entradas": sec_entradas, "bic": sec_bic}


def main():
    pedidas = sys.argv[1:] or list(SECCIONES)
    for s in pedidas:
        print(f"== {s}")
        SECCIONES[s]()
    out = os.path.join(HERE, "outputs", "manuscript_numbers.csv")
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ref", "descripcion", "calculado", "texto", "coincide"])
        w.writeheader()
        w.writerows(FILAS)
    for r in FILAS:
        print(f"{'OK ' if r['coincide'] else 'DIF'} {r['ref']:<12} {r['descripcion'][:62]:<62} "
              f"calc {r['calculado']:<10} text {r['texto']}")
    n = sum(r["coincide"] for r in FILAS)
    print(f"{n}/{len(FILAS)} match -> {out}")


if __name__ == "__main__":
    main()
