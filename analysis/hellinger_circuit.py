# -*- coding: utf-8 -*-
"""Our own Hellinger fits for the circuit nodes without a usable published
covariance (Hellinger fits with GSFML picks).

Nodes that can be fitted with picks (chron C32n.1n y, 71.07 Ma in CK95):
  801@71.07 AUS-EANT  (Tikku and Cande 2000 32y: published cov not positive definite)
  802@71.1  EANT-AFR  (Bernard et al. 2005 32y: FZ pole, no cov)
Nodes 833@79.1 and 833@86.0 (LHR-AUS) and 802@75.5 cannot be fitted: GSFML
has no Tasman Sea picks older than C33y (Gaina et al. 1998) and
75.5 Ma does not correspond to any chron (Bernard et al. 2005, FZ pole).

Data: GSFML.global.picks.gmt (Seton et al. 2014). The plate of each pick comes from a
point-in-polygon test with the static polygons of Mueller et al. (2019).
Duplicates between references (< 8 km apart, e.g. Segoufin et al. 2004 and Cande and
Patriat 2015, which digitize the same anomalies) are counted once. sd = 3.5 km per pick (scheme U, Cande et al.
2010), as in our own two-plate fitting script.

Sections: GSFML has no fracture zones that reach these isochrons, so
sections are defined by windows on the position reconstructed with the .rot pole
(in the fixed-plate frame), cut at the visible offsets of the isochron.
A section that includes an unrecognized offset lowers kappa and enlarges the cov.

Convention: Hellinger gives A = A0 exp(M(h)) (h in the moving-plate frame).
The engine uses exp(M(h')) R (h' in the fixed-plate frame): h' = A h,
C' = A (cov_h / kappa) A^T.

Usage: py hellinger_circuit.py   (writes outputs/frame_coverage/hellinger/)
"""
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from matplotlib.path import Path

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))), "tools", "hellinger_chang"))  # copy of Chang's port
from chang import hellinger, region  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(AQUI, "outputs", "frame_coverage", "hellinger")
GSFML = (__import__("os").environ["PAPER_ROOT"] + "/external/gsfml/GSFML.global.picks.gmt")
POLY = (__import__("os").environ["PAPER_ROOT"] + "/external/muller2019/Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons_2019_v1.shp")
OSGEO = __import__("os").environ.get("OSGEO4W_ROOT", r"C:\OSGeo4W64")   # GDAL (ogr2ogr)
SD = 3.5
R_KM = 6371.0

# node: (moving plate, fixed plate, .rot pole, chron, end, pick box, references,
#        sections [(lon_min, lon_max, lat_min, lat_max)] in the fixed-plate frame)
CASOS = {
    "801_71.07": dict(
        movil=801, fija=802, polo=(-9.48, -142.98, 26.13), cron="C32n.1n", extremo="y",
        caja=(70, 160, -68, -30), refs=["Tikku+Cande_1999_JGR"],
        secciones=[(103.5, 108.9, -63, -61.8), (108.9, 113.8, -63, -61.8),
                   (113.8, 119.0, -63, -61.8), (119.0, 124.5, -63, -61.8),
                   (124.5, 129.0, -63, -61.8), (129.0, 133.5, -63, -61.8)]),
    "802_71.1": dict(
        movil=802, fija=701, polo=(-1.2, -42.4, 12.38), cron="C32n.1n", extremo="y",
        caja=(0, 70, -72, -25), refs=None,
        secciones=[(17.0, 23.5, -48.2, -46.2), (24.5, 26.0, -48.2, -46.8),
                   (26.5, 28.5, -49.0, -46.8), (34.5, 38.6, -40.0, -38.3),
                   (38.7, 41.9, -37.6, -36.6)]),
}


def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def latlon(v):
    v = np.atleast_2d(v)
    return np.degrees(np.arcsin(np.clip(v[:, 2], -1, 1))), np.degrees(np.arctan2(v[:, 1], v[:, 0]))


def pole_matrix(lat, lon, ang):
    k = xyz(lat, lon)
    a = np.radians(ang)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K


def rotvec(A):
    a = np.arccos(np.clip((np.trace(A) - 1) / 2, -1, 1))
    if a < 1e-12:
        return np.zeros(3)
    return a / (2 * np.sin(a)) * np.array([A[2, 1] - A[1, 2], A[0, 2] - A[2, 0], A[1, 0] - A[0, 1]])


def leer_picks():
    rows, meta = [], None
    for l in open(GSFML, encoding="utf-8", errors="replace"):
        if l.startswith("# @D"):
            meta = l[4:].strip().split("|")
        elif l[0] not in "#>" and meta:
            x = l.split()
            if len(x) >= 2:
                rows.append(dict(lon=float(x[0]), lat=float(x[1]), chron=meta[0], end=meta[1],
                                 ref=meta[4]))
    return pd.DataFrame(rows)


def poligonos():
    gj = os.path.join(OUT, "static_poly_sur.geojson")
    if not os.path.exists(gj):
        env = dict(os.environ, PATH=os.path.join(OSGEO, "bin") + ";" + os.environ["PATH"],
                   GDAL_DATA=os.path.join(OSGEO, "share", "gdal"),
                   PROJ_LIB=os.path.join(OSGEO, "share", "proj"))
        subprocess.run([os.path.join(OSGEO, "bin", "ogr2ogr.exe"), "-f", "GeoJSON", "-overwrite",
                        "-spat", "0", "-75", "180", "-20", gj, POLY], env=env, check=True,
                       capture_output=True)
    out = []
    for ft in json.load(open(gj))["features"]:
        g, p = ft["geometry"], ft["properties"]
        for poly in ([g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"]):
            out.append((Path(np.array(poly[0])[:, :2]), [Path(np.array(h)[:, :2]) for h in poly[1:]],
                        p.get("PLATEID1")))
    return out


def placa(polys, lon, lat):
    for path, holes, pid in polys:
        if path.contains_point((lon, lat)) and not any(h.contains_point((lon, lat)) for h in holes):
            return pid
    return None


PREFERENCIA = ["Cande+Patriat_2015_GJI", "Tikku+Cande_1999_JGR"]   # first in the deduplication


def sin_duplicados(d, km=8.0):
    """A pick is discarded if another, already accepted reference has one within km (the
    same anomaly digitized twice); picks within one reference are kept."""
    d = d.assign(_o=[PREFERENCIA.index(r) if r in PREFERENCIA else 9 for r in d.ref]).sort_values("_o")
    keep, X, refs = [], xyz(d.lat.values, d.lon.values), d.ref.values
    for i in range(len(d)):
        if all(refs[i] == refs[j] or np.degrees(np.arccos(np.clip(X[i] @ X[j], -1, 1))) * 111.195 > km
               for j in keep):
            keep.append(i)
    return d.iloc[keep].drop(columns="_o")


def preparar(caso, picks, polys):
    c = CASOS[caso]
    a, b, la0, la1 = c["caja"]
    d = picks[(picks.chron == c["cron"]) & (picks.end == c["extremo"])
              & picks.lon.between(a, b) & picks.lat.between(la0, la1)].copy()
    if c["refs"]:
        d = d[d.ref.isin(c["refs"])]
    d["placa"] = [placa(polys, x, y) for x, y in zip(d.lon, d.lat)]
    d = d[d.placa.isin([c["movil"], c["fija"]])]
    d = pd.concat([sin_duplicados(g) for _, g in d.groupby("placa")])
    A = pole_matrix(*c["polo"])
    rl, ro = latlon(xyz(d.lat.values, d.lon.values) @ A.T)
    mov = (d.placa == c["movil"]).values
    d["lat_f"] = np.where(mov, rl, d.lat)
    d["lon_f"] = np.where(mov, ro, d.lon)
    d["lado"] = np.where(mov, 1, 2)
    d["seccion"] = 0
    for k, (x0, x1, y0, y1) in enumerate(c["secciones"], 1):
        m = d.lon_f.between(x0, x1) & d.lat_f.between(y0, y1)
        d.loc[m, "seccion"] = k
    # only sections with picks on both sides
    ok = d[d.seccion > 0].groupby("seccion").lado.nunique()
    d = d[d.seccion.isin(ok[ok == 2].index)]
    return d


def ajustar(caso, d):
    c = CASOS[caso]
    cr = [(int(r.lado), int(r.seccion), float(r.lat), float(r.lon), SD) for r in d.itertuples()]
    res = hellinger.fit_hellinger(cr, c["polo"])
    kap = res["kappa_hat"]
    A = pole_matrix(res["lat"], res["lon"], res["rho"])
    C_mov = res["cov_h"] / kap                 # moving-plate frame (A0 exp(M(h)))
    C_fix = A @ C_mov @ A.T                    # fixed-plate frame (engine convention)
    A_rot = pole_matrix(*c["polo"])
    dif = rotvec(A_rot @ A.T)                  # .rot vs fit difference, fixed-plate frame
    q = float(dif @ np.linalg.solve(C_fix, dif))
    # region 95%: h' C'^-1 h' <= 3 F(3, df) (= xchi * kappa, with C' already divided by kappa)
    xchi = kap * region.critical_value(0.95, df=res["df"], kappa_hat=kap, known_kappa=False)
    sd_km = np.sqrt(np.linalg.eigvalsh(C_fix)) * R_KM
    out = dict(caso=caso, movil=c["movil"], fija=c["fija"], cron=c["cron"] + c["extremo"],
               n=res["ndat"], s=res["nsect"], df=res["df"], kappa_hat=kap, r_min=res["r_min"],
               polo_ajuste=(res["lat"], res["lon"], res["rho"]), polo_rot=c["polo"],
               dif_rot_grados=float(np.degrees(np.linalg.norm(dif))),
               chi2_rot_en_region=q, umbral_95=float(xchi),
               rot_dentro_95=bool(q <= xchi), sd_km_ejes=sd_km.tolist(),
               cov_fija_rad2=C_fix.tolist(), cov_movil_rad2=C_mov.tolist(),
               refs=sorted(d.ref.unique().tolist()),
               picks_por_seccion=d.groupby(["seccion", "lado"]).size().unstack(fill_value=0)
               .to_dict(orient="index"))
    return out


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    picks = leer_picks()
    polys = poligonos()
    for caso in CASOS:
        d = preparar(caso, picks, polys)
        d.to_csv(os.path.join(OUT, f"picks_{caso}.csv"), index=False)
        r = ajustar(caso, d)
        json.dump(r, open(os.path.join(OUT, f"fit_{caso}.json"), "w"), indent=1)
        print(f"\n{caso}: n {r['n']}, s {r['s']}, df {r['df']}, kappa {r['kappa_hat']:.3f}")
        print("  picks per section (side 1 moving / 2 fixed):", r["picks_por_seccion"])
        print("  fitted pole", np.round(r["polo_ajuste"], 3), " .rot pole", r["polo_rot"])
        print(f"  difference {r['dif_rot_grados']:.3f} deg; chi2 of the .rot {r['chi2_rot_en_region']:.2f}"
              f" vs 95% threshold {r['umbral_95']:.2f} -> inside: {r['rot_dentro_95']}")
        print("  sd per axis (km):", np.round(r["sd_km_ejes"], 1))
