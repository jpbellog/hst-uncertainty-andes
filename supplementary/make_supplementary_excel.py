# -*- coding: utf-8 -*-
"""Excel data workbook of the Supplementary Material.

Collects the long data tables of the Supplementary Material as sheets of a single file,
TablesS_data.xlsx. Requires make_supplementary_tables.py to have been run (Tables S1 and S2 as CSV).
Sheets:
  README              description of each sheet and column
  S1_hotspot_fits     fit of each hotspot in the twelve models (TableS1_hotspot_fits.csv)
  S2_collision_bands  collision latitude bands, twelve models, three f_PAC, 0-66 Ma (TableS2)
  S3_geochronology    age compilation of the four chains: published age and age used
                      (recalibrated) with its 2 sigma, type, use in the fit and in how many of the
                      twelve models the outlier rule removes it (hotspot_fit.json of each folder)
  S4_hotspot_seeds    present-day reference position of each hotspot (SEEDS of 03b_hotspot_fit.py)
  S5_convergence      Nazca-South America convergence rate and its 1 sigma by stage at the margin
                      nodes closest to 20 and 30 S, from the published stage covariances
                      (Quiero et al. 2022: models/Q22_WK08-A/vc_sigma_quiero.csv; Wilson and DeMets 2025:
                      models/WD25_WK08-A/vc_sigma_wd25.csv)
Usage:  py make_supplementary_excel.py
"""
import json
import os
import sys
import importlib.util
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "figures"))
import paper_models as MP  # noqa: E402

DATA = os.path.join(ROOT, "models", "Q22_WK08-A", "data")
CIRC = {"Quiero": "Q22", "WD25": "WD25", "SG12": "SG12"}
CADENAS = [  # (key in hotspot_fit.json, name, file, age column used, 2s column used, source)
    ("Easter", "Easter", "ages_easter_ray2012_recal.csv", "age_E", "err2s_E", "Ray et al. (2012); EarthRef step data"),
    ("JFR", "Juan Fernandez", "ages_jfr_lara2018.csv", "age_ma", None, "Lara et al. (2018)"),
    ("TT", "San Felix", "ages_taltal_recal.csv", "age_recal", None, "Haase et al. (2000); Gonzalez-Ferran (1987)"),
    ("Foundation", "Foundation", "ages_foundation_oconnor2004_v3.csv", "age_recal", None,
     "O'Connor et al. (1998, 2001, 2004)"),
]


def geocrono():
    fits = {}
    for carpeta, circ, marco in MP.TEXTO + MP.SUPL:
        d = json.load(open(os.path.join(ROOT, "models", carpeta, "hotspot_fit.json"), encoding="utf-8"))
        fits[f"{CIRC[circ]} + {marco}"] = d["tracks"]
    out = []
    for key, nom, fn, col, ecol, src in CADENAS:
        df = pd.read_csv(os.path.join(DATA, fn), comment="#")
        a_pub = df.age_ma.astype(float)
        a_use = df[col].astype(float)
        e_use = df[ecol].astype(float) if ecol else df.err2s_ma.astype(float) * (a_use / a_pub).where(a_pub > 0, 1.0)
        rem = {s: [m for m, t in fits.items() if s in {e["sample"] for e in t[key]["excluded"]}]
               for s in df["sample"]}
        t = pd.DataFrame({
            "chain": nom, "sample": df["sample"],
            "edifice_or_material": df["edifice"] if "edifice" in df else df.get("material"),
            "lat": df.lat, "lon": df.lon,
            "age_published_Ma": a_pub, "err2s_published_Myr": df.err2s_ma,
            "age_used_Ma": a_use.round(4), "err2s_used_Myr": e_use.round(4),
            "age_type": df.age_type, "in_fit": df.use_fit.astype(int),
            "removed_by_outlier_rule_in_n_of_12_models": [len(rem[s]) for s in df["sample"]],
            "removed_in_models": ["; ".join(rem[s]) for s in df["sample"]],
            "source": src})
        if "pos_err_km" in df:
            t["position_err_km"] = df.pos_err_km
        t["note"] = ""
        if "assumed_err" in df:
            # San Ambrosio (Gonzalez-Ferran 1987): published as 2.93 +/- 0.15 Ma without confidence level;
            # taken as 1 sigma (0.30 at 2 sigma). In the CSV err2s_ma is already 0.30 and age_type "unverified".
            a = (df.assumed_err == 1).values
            t.loc[a, "err2s_published_Myr"] = 0.15
            t.loc[a, "age_type"] = "K-Ar"
            t.loc[a, "note"] = ("published as 2.93 +/- 0.15 Ma without confidence level; "
                                "taken as 1 sigma (err2s_used = 0.30)")
        out.append(t)
    return pd.concat(out, ignore_index=True)


def semillas():
    spec = importlib.util.spec_from_file_location("hb", os.path.join(ROOT, "models", "Q22_WK08-A", "03b_hotspot_fit.py"))
    txt = open(spec.origin, encoding="utf-8").read()
    i = txt.index("SEEDS = {")
    seeds = eval(txt[i + len("SEEDS = "):txt.index("}", i) + 1])
    nom = {13: "Juan Fernandez", 14: "Easter", 15: "San Felix", 16: "Foundation"}
    return pd.DataFrame([dict(chain=nom[k], plate_id=k, lat=v[0], lon=v[1]) for k, v in seeds.items()])


def convergencia():
    out = []
    for lab, fn in (("Quiero et al. (2022) stage covariances", os.path.join(ROOT, "models", "Q22_WK08-A", "vc_sigma_quiero.csv")),
                    ("Wilson and DeMets (2025) Nazca-South America stage covariances",
                     os.path.join(ROOT, "models", "WD25_WK08-A", "vc_sigma_wd25.csv"))):
        d = pd.read_csv(fn)
        d["stage_young_Ma"] = d[["t_ini", "t_fin"]].min(axis=1)
        d["stage_old_Ma"] = d[["t_ini", "t_fin"]].max(axis=1)
        lats = sorted(d.lat.unique())
        for q in (-20.0, -30.0):
            la = min(lats, key=lambda x: abs(x - q))
            s = d[d.lat == la].sort_values("stage_young_Ma")
            out.append(pd.DataFrame({"source": lab, "node_lat": s.lat.round(3), "node_lon": s.lon.round(3),
                                     "stage_young_Ma": s.stage_young_Ma, "stage_old_Ma": s.stage_old_Ma,
                                     "Vc_cm_yr": s.vc_mag.round(2), "sigma_Vc_cm_yr": s.sigma_mag.round(2),
                                     "Vort_cm_yr": s.vc_ort.round(2), "sigma_Vort_cm_yr": s.sigma_ort.round(2)}))
    return pd.concat(out, ignore_index=True)


README = [
    ("General", "Data tables of the Supplementary Material. Latitudes and longitudes in degrees, negative south and west. "
                "Models: circuits Q22 (Quiero et al. 2022), WD25 (Wilson and DeMets 2025), SG12 (Somoza and Ghidella 2012); "
                "hotspot frames WK08-A, GMHRF-B, GWG20 (main text) and GMHRF-A' (Supplementary Section S4)."),
    ("S1_hotspot_fits", "Table S1 in tabular form: fit of each hotspot in the twelve models (see the caption of Table S1)."),
    ("S2_collision_bands", "Monte Carlo collision latitude on the present-day margin, all sources sampled, N = 1000, "
                           "at 1 Myr steps for f_PAC = 0.50, 0.44 and 0.36. lat_nominal: unperturbed model; "
                           "lat_p2_5 to lat_p97_5: percentiles (68% interval p16-p84, 95% interval p2.5-p97.5). "
                           "valid = 0: the collision requires a frame rotation older than the frame covariances "
                           "(Supplementary Section S2); these rows are not used in the paper."),
    ("S3_geochronology", "Ages of the four chains. age_published / err2s_published: as published. age_used / err2s_used: "
                         "recalibrated to Fish Canyon sanidine 28.201 Ma (Kuiper et al. 2008) with the decay constants of "
                         "Min et al. (2000) (Supplementary Section S3); for Easter, recomputed from the step data. "
                         "err2s_published for San Ambrosio: published without confidence level (see note column). in_fit = 0: set aside before the fit by the selection rules. removed_by_outlier_rule: number of "
                         "the twelve models in which the outlier rule removes the sample, and which."),
    ("S4_hotspot_seeds", "Present-day reference position from which each track is computed (plate ID of the hotspot)."),
    ("S5_convergence", "Nazca (Farallon)-South America convergence rate Vc and its trench-orthogonal component Vort "
                       "with 1-sigma uncertainties, by stage, at the margin grid nodes closest to 20 and 30 S, "
                       "propagated from the published stage covariances (Eq. 6 of the main text). Stage limits are those "
                       "of each source."),
]


def main():
    s1 = pd.read_csv(os.path.join(HERE, "TableS1_hotspot_fits.csv"))
    s2 = pd.read_csv(os.path.join(HERE, "TableS2_collision_bands.csv"))
    fn = os.path.join(HERE, "TablesS_data.xlsx")
    with pd.ExcelWriter(fn, engine="openpyxl") as xw:
        pd.DataFrame(README, columns=["sheet", "description"]).to_excel(xw, sheet_name="README", index=False)
        s1.to_excel(xw, sheet_name="S1_hotspot_fits", index=False)
        s2.to_excel(xw, sheet_name="S2_collision_bands", index=False)
        geocrono().to_excel(xw, sheet_name="S3_geochronology", index=False)
        semillas().to_excel(xw, sheet_name="S4_hotspot_seeds", index=False)
        convergencia().to_excel(xw, sheet_name="S5_convergence", index=False)
        for ws in xw.book.worksheets:
            ws.freeze_panes = "A2"
            for col in ws.columns:
                w = max(len(str(c.value)) if c.value is not None else 0 for c in col[:200])
                ws.column_dimensions[col[0].column_letter].width = min(max(10, w + 2), 60 if ws.title != "README" else 120)
    print(f"written {fn}")


if __name__ == "__main__":
    main()
