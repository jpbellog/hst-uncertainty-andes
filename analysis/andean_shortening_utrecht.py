# -*- coding: utf-8 -*-
"""Restoration of Andean shortening with the block model of Schepers et al. (2017), updated by
Boschman (2021), as implemented in the Utrecht model (van Hinsbergen et al. 2026, Supplementary
Files 1: Utrecht_Paleogeography_15Dec25.rot and Paleogeography/Continental polygons.gpml), in
place of our blocks (Arriagada et al. 2008).

What changes (upper plate only; tracks and circuit unchanged, model WK08):
  - each porphyry, epithermal deposit and sample is assigned to the Utrecht polygon that contains
    it (if several, the one with the smallest area), and its position at age t is R_{pid/201}(t) x0
    with the Utrecht .rot (relative to stable South America, 201). Validated with our .rot: the
    same computation with PLATEID1 reproduces Latitude_ROT (max 0.01 km).
  - the 69 trench points (Margin_SAM.shp, present day) are assigned to the nearest polygon and
    restored in the same way, at the same ages as Margin_SAM_0.1.csv (every 0.1 Myr).
The assignment of step 07 (WK08) and the test of step 12 are then repeated with those positions
and that trench, together with the nominal collision (f_PAC 0.44).

Usage:  py andean_shortening_utrecht.py
Outputs: outputs/andean_shortening/
"""
import os
import sys
import json
import importlib
import numpy as np
import pandas as pd
import shapefile
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
PIPE = os.path.join(ROOT, "models", "Q22_WK08-A")
OUT = os.path.join(HERE, "outputs", "andean_shortening")
UROT = os.path.join(HERE, "..", "external", "vanhinsbergen2026", "Utrecht_Paleogeography_15Dec25.rot")
MSAM = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "margin")  # copy of Margin_SAM (see rotation_models/README.md)
os.environ.setdefault("PAPER_ROT_FILE", os.path.join(PIPE, "rot", "Global_Model_R_BG_2025_hsb.rot"))
os.environ.setdefault("PAPER_HS_FIT", os.path.join(PIPE, "hotspot_fit.json"))
sys.path.insert(0, PIPE)
eng = importlib.import_module("02_rot_engine")
s07 = importlib.import_module("07_deposit_tracks")
mc = s07.mc


def polys():
    P = json.load(open(os.path.join(OUT, "poligonos_andes.json"), encoding="utf-8"))
    out = []
    for p in P:
        # South American plates only (2xx, 2xxx, 2xxxx); excludes polygons that cross the date
        # line (e.g. Fiji, 8005), which the planar inclusion test handles incorrectly
        if p["pid"] is None or str(p["pid"])[0] != "2" or (max(p["lon"]) - min(p["lon"])) > 90 \
                or p["t_end"] not in ("0", "0.0", None) and float(p["t_end"]) > 0:
            continue
        lo, la = np.array(p["lon"]), np.array(p["lat"])
        area = 0.5 * abs(np.dot(lo, np.roll(la, 1)) - np.dot(la, np.roll(lo, 1)))
        out.append((p["pid"], p["name"], Path(np.c_[lo, la]), area, lo, la))
    return out


def pid_of(lo, la, PP, nearest=False):
    best = None
    for pid, name, path, area, plo, pla in PP:
        if path.contains_point((lo, la)) and (best is None or area < best[3]):
            best = (pid, name, 0.0, area)
    if best is None and nearest:
        for pid, name, path, area, plo, pla in PP:
            d = np.min(np.hypot(plo - lo, (pla - la) * 1.0) * np.cos(np.radians(la)) ** 0) * 111.2
            if best is None or d < best[2]:
                best = (pid, name, d, area)
    return best


def restore(M, pid, t, la, lo):
    q = M.rel_rotation(int(pid), float(t), 201)
    return eng.xyz_to_latlon(eng.quat_to_matrix(q) @ eng.latlon_axis(la, lo))


def main():
    PP = polys()
    MU = eng.RotModel(UROT)
    # deposits and samples
    por = pd.read_csv(os.path.join(s07.DEP, "Porphyry_Plate_ID_geodyn_v7.csv"))
    epi = pd.read_csv(os.path.join(s07.DEP, "Epithermal_Plate_Id_geodyn_v6_thesis.csv"))
    g = pd.read_excel(s07.GEOCHEM)
    resumen = []
    for nom, d, age, lat, lon, lr, lor in (
            ("porphyries", por, "Age MY", "Latitude", "Longitude", "Latitude_ROT", "Longitude_ROT"),
            ("epithermal", epi, "Age MY", "Latitude", "Longitude", "Latitude_ROT", "Longitude_ROT"),
            ("samples", g, "Age (Ma)", "Latitude", "Longitude", "Latitude_restored", "Longitude_restored")):
        a = pd.to_numeric(d[age], errors="coerce").values
        la0 = pd.to_numeric(d[lat], errors="coerce").values
        lo0 = pd.to_numeric(d[lon], errors="coerce").values
        newla, newlo, upid = np.full(len(d), np.nan), np.full(len(d), np.nan), np.full(len(d), np.nan)
        for i in range(len(d)):
            if not (np.isfinite(a[i]) and np.isfinite(la0[i]) and np.isfinite(lo0[i])):
                continue
            b = pid_of(lo0[i], la0[i], PP, nearest=True)
            upid[i] = b[0]
            newla[i], newlo[i] = restore(MU, b[0], a[i], la0[i], lo0[i])
        old_la = pd.to_numeric(d[lr], errors="coerce").values
        old_lo = pd.to_numeric(d[lor], errors="coerce").values
        dk = np.hypot(newla - old_la, (newlo - old_lo) * np.cos(np.radians(old_la))) * 111.2
        sh_u = np.hypot(newla - la0, (newlo - lo0) * np.cos(np.radians(la0))) * 111.2
        sh_a = np.hypot(old_la - la0, (old_lo - lo0) * np.cos(np.radians(la0))) * 111.2
        d = d.copy()
        d["PID_UTR"] = upid
        d[lr], d[lor] = newla, newlo
        d.to_csv(os.path.join(OUT, f"{nom}_utrecht.csv"), index=False)
        ok = np.isfinite(dk)
        resumen.append(dict(conjunto=nom, n=int(ok.sum()),
                            restauracion_arriagada_km_med=float(np.nanmedian(sh_a)),
                            restauracion_utrecht_km_med=float(np.nanmedian(sh_u)),
                            dif_posicion_km_med=float(np.nanmedian(dk)),
                            dif_posicion_km_p90=float(np.nanpercentile(dk, 90)),
                            dif_posicion_km_max=float(np.nanmax(dk))))
        if nom == "samples":
            gx = os.path.join(OUT, "samples_utrecht.xlsx")
            d.to_excel(gx, index=False)
    pd.DataFrame(resumen).to_csv(os.path.join(OUT, "position_summary.csv"), index=False)
    print(pd.DataFrame(resumen).round(1).to_string(index=False))

    # trench
    r = shapefile.Reader(os.path.join(MSAM, "Margin_SAM.shp"))
    pts = [s.points[0] for s in r.shapes()]                         # (lon, lat), present day
    tp = [pid_of(lo, la, PP, nearest=True) for lo, la in pts]
    old = pd.read_csv(os.path.join(MSAM, "Margin_SAM_0.1.csv"))
    new = old.astype(float)
    for k, ((lo, la), b) in enumerate(zip(pts, tp), start=1):
        for j, t in enumerate(old.age.values):
            new.loc[j, f"Lat_{k}"], new.loc[j, f"Lon_{k}"] = restore(MU, b[0], t, la, lo)
    new.to_csv(os.path.join(OUT, "Margin_SAM_0.1_utrecht.csv"), index=False)
    pd.DataFrame([dict(punto=k, lon=lo, lat=la, pid=b[0], poligono=b[1], dist_km=round(b[2], 1))
                  for k, ((lo, la), b) in enumerate(zip(pts, tp), start=1)]).to_csv(
        os.path.join(OUT, "fosa_bloques.csv"), index=False)
    dd = []
    for k in range(1, len(pts) + 1):
        for t in (10, 20, 30, 40, 50, 60):
            j = int(np.argmin(np.abs(old.age.values - t)))
            a1, o1 = old.loc[j, f"Lat_{k}"], old.loc[j, f"Lon_{k}"]
            a2, o2 = new.loc[j, f"Lat_{k}"], new.loc[j, f"Lon_{k}"]
            if pts[k - 1][1] < -5 and pts[k - 1][1] > -45:
                dd.append(dict(punto=k, lat0=pts[k - 1][1], t=t,
                               dif_km=np.hypot(a2 - a1, (o2 - o1) * np.cos(np.radians(a1))) * 111.2))
    dd = pd.DataFrame(dd)
    print("\ntrench 5-45S, difference Utrecht - Arriagada (km), median and max by age:")
    print(dd.groupby("t").dif_km.agg(["median", "max"]).round(0).to_string())


if __name__ == "__main__":
    main()
