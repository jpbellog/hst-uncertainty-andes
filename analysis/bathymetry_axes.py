# -*- coding: utf-8 -*-
"""Bathymetric test of the tracks: axis of each ridge from ETOPO1 and distance from the nominal
tracks of the twelve models to that axis.

Method: in a lat-lon box that contains the ridge, the local maxima of the bathymetry are found in
each longitude column; they are joined into a line (the axis), and the maxima are checked to belong
to the chain under analysis and not to other relief.

  1. ETOPO1 (bedrock, 1') smoothed with a Gaussian of 20 km full width at half maximum (as in
     track_vs_ridge_axis.py).
  2. In each column of the box, local maxima with prominence >= PROM_M.
  3. Tracking from a seed placed on the chain: in each column the maximum closest in latitude to
     the last accepted one is taken, if it lies within SALTO_GRADOS; if there is none, the column
     stays empty (gap). Only maxima with relief >= RELIEVE_M above the column median within the
     box are accepted.
  4. The axis is the running median (MEDIANA_COL columns) of the accepted latitudes.
  5. Points closer than FOSA_KM to the trench (present-day margin line from
     inputs/margin/Margin_SAM_0.1.csv, age 0) are discarded, to avoid picking the inner wall.
  6. Iquique runs nearly N-S: the main axis comes from the maxima per latitude row and the
     version per longitude is kept as a control.
  7. Parameters per ridge. JFR: 20 km smoothing and a 0.12 degree corridor. Taltal (scattered
     seamounts): 40 km and 0.4 degrees. Iquique (seamounts on a swell): 40 km and 0.3 degrees.
     Nazca (ridge ~140 km wide): 40 km, 0.7 degrees and a running median of 21 columns, because
     in a N-S section the crest is wide and flat and the maximum jumps from one edge to the other.
  8. Points closer than CERCA_FOSA_KM to the margin line are flagged (cerca_fosa): there the
     outer bulge can be confused with the ridge, and the summary separates them.

For each axis point and each model: distance to the nominal track (f_PAC 0.44, output
track_nominal_0.44_hsb_cc_g13.csv of step 04, densified), signed (+ = axis to the left of the
direction of increasing tau), tau of the nearest point, and whether the axis falls inside the 68%
and 95% ellipses (track_ellipses_0.44_hsb_cc_g13.csv, nodes every 10 Myr, linearly interpolated
in tau; outside the node range the end node is used and flagged). Also the distance from the axis
point to the nearest dated sample of the chain (con_datos if < DATOS_KM).

Outputs in outputs/bathymetry_axes/: ejes_<ridge>.csv, distancias.csv, resumen.csv, map .png/.pdf.
Usage:  py -I bathymetry_axes.py
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(HERE, "outputs", "bathymetry_axes")
# ETOPO1 bedrock grid (GeoTIFF)
ETOPO = (__import__("os").environ["PAPER_ROOT"] + "/external/etopo1/ETOPO1_Bed_g_geotiff.tif")
NX, NY, OFF0 = 21601, 10801, 87470          # as in track_vs_ridge_axis.py
LON0, LAT0, D = -180.0, 90.0, 1 / 60.0
BOX = (-87.0, -70.0, -37.0, -14.0)   # extended to the SW for the Iquique corridor
R_KM = 6371.0
PROM_M = 150.0
RELIEVE_M = 300.0
HUECO_COL = 15          # columns without axis that break the line on the map
SALTO_GRADOS = 0.12
MEDIANA_COL = 11
FOSA_KM = 60.0
DATOS_KM = 75.0
CERCA_FOSA_KM = 120.0   # from the margin line; there the outer bulge can pass for the ridge
TAG = "_hsb_cc_g13"

sys.path.insert(0, os.path.join(ROOT, "figures"))
import paper_models as MP  # noqa: E402

TODOS = MP.TEXTO + MP.SUPL
# ridge: model track, box (lon_min, lon_max, lat_min, lat_max), seed (lon, lat),
# sweep direction, age file
DORSALES = {
    "Nazca": dict(track="Easter", caja=(-82.5, -76.0, -25.5, -14.6), semilla=(-80.0, -19.5),
                  modo="lon", fwhm=40.0, salto=0.7, mediana=21, edades="ages_easter_ray2012_recal.csv"),
    # Iquique (Bello-Gonzalez et al. 2018, sec. 2.2): seamounts on a swell, from north of the
    # Taltal Ridge to the trench at the Arica bend. It runs nearly N-S.
    # The Iquique axis joins the summits of the highest seamounts (mode "montes", eje_montes)
    # within a +/-75 km corridor around a polyline from the northern summit (20.52 S, 71.82 W)
    # to the Challenger FZ (see vertices_corredor); summits closer than
    # 50 km to the axes of Taltal (complete, unclipped), JFR and Copiapo are excluded.
    "Iquique": dict(track="Foundation", modo="montes", fwhm=40.0, edades=None,
                    # polyline: northern summit, summit at 27.47 S where the belt turns S, and the
                    # intersection of that direction (towards 32.5 S, 80.0 W) with the Challenger FZ
                    corredor=((-20.5167, -71.8167), (-27.4667, -77.2667)), hacia=(-32.5, -80.0),
                    semiancho_km=75.0,
                    excluir=("Taltal", "JFR", "Copiapo"), caja=(-87.0, -70.0, -37.0, -18.5)),
    "Taltal": dict(track="TT", caja=(-80.4, -71.6, -27.2, -24.6), semilla=(-80.1, -26.3),
                   modo="lon", fwhm=40.0, salto=0.4, edades="ages_taltal_recal.csv"),
    "JFR": dict(track="JFR", caja=(-81.2, -72.2, -34.4, -32.6), semilla=(-78.85, -33.62),
                modo="lon", fwhm=20.0, edades="ages_jfr_lara2018.csv"),
}


def etopo_box():
    with open(ETOPO, "rb") as f:
        assert f.read(2) == b"II"
    z = np.memmap(ETOPO, dtype="<i2", mode="r", offset=OFF0, shape=(NY, NX))
    i0, i1 = int(round((BOX[0] - LON0) / D)), int(round((BOX[1] - LON0) / D))
    j0, j1 = int(round((LAT0 - BOX[3]) / D)), int(round((LAT0 - BOX[2]) / D))
    sub = np.array(z[j0:j1 + 1, i0:i1 + 1], dtype=float)
    return LON0 + D * np.arange(i0, i1 + 1), LAT0 - D * np.arange(j0, j1 + 1), sub


def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def gc_km(la1, lo1, la2, lo2):
    a, b = xyz(la1, lo1), xyz(la2, lo2)
    return R_KM * np.arccos(np.clip(np.sum(a * b, axis=-1), -1, 1))


def azimut(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2)
    dl = np.radians(lo2 - lo1)
    return np.degrees(np.arctan2(np.sin(dl) * np.cos(p2),
                                 np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)))


def fosa():
    m = pd.read_csv(os.path.join(ROOT, "inputs", "margin", "Margin_SAM_0.1.csv"))
    r = m[m.age == 0].iloc[0]
    n = (len(m.columns) - 1) // 2
    lo = np.array([r[f"Lon_{i}"] for i in range(1, n + 1)])
    la = np.array([r[f"Lat_{i}"] for i in range(1, n + 1)])
    # densify every 0.05 degrees
    L, A = [], []
    for k in range(n - 1):
        s = max(2, int(np.hypot(lo[k + 1] - lo[k], la[k + 1] - la[k]) / 0.05))
        L.append(np.linspace(lo[k], lo[k + 1], s, endpoint=False))
        A.append(np.linspace(la[k], la[k + 1], s, endpoint=False))
    return np.concatenate(A), np.concatenate(L)


def eje(lon, lat, zs, cfg, fla, flo):
    zs = zs[cfg["fwhm"]]
    x0, x1, y0, y1 = cfg["caja"]
    ii = np.where((lon >= x0) & (lon <= x1))[0]
    jj = np.where((lat >= y0) & (lat <= y1))[0]
    Z = zs[np.ix_(jj, ii)]                      # rows = lat (decreasing), columns = lon
    L, A = lon[ii], lat[jj]
    if cfg["modo"] == "lat":                    # sweep by rows: transpose
        Z, L, A = Z.T, A, L                     # now "columns" = lat, coordinate = lon
    crudos, acept = [], []
    sem_col, sem_coord = (cfg["semilla"] if cfg["modo"] == "lon" else cfg["semilla"][::-1])
    c0 = int(np.argmin(np.abs(L - sem_col)))
    picos = {}
    for c in range(len(L)):
        prof = Z[:, c]
        pk, pr = find_peaks(prof, prominence=PROM_M)
        med = np.median(prof)
        picos[c] = [(A[p], prof[p], prof[p] - med) for p in pk]
        for a, zz, rel in picos[c]:
            crudos.append((L[c], a, zz, rel))
    # the seed snaps to the highest maximum within 0.8 degrees in its column
    cs = [q for q in picos[c0] if abs(q[0] - sem_coord) <= 0.8 and q[2] >= RELIEVE_M]
    if cs:
        sem_coord = max(cs, key=lambda q: q[1])[0]
    for paso in (1, -1):
        ref, ultima = sem_coord, (c0 if paso == 1 else c0 - 1) - paso
        c = c0 if paso == 1 else c0 - 1
        while 0 <= c < len(L):
            # corridor: SALTO_GRADOS, plus 0.3 times the width of the columns without an accepted
            # maximum since the last one (gaps between seamounts); within the corridor the highest is taken
            tol = cfg.get("salto", SALTO_GRADOS) + 0.3 * abs(c - ultima - paso) * D
            cand = [(abs(a - ref), a, zz, rel) for a, zz, rel in picos[c] if rel >= RELIEVE_M]
            cand = [q for q in cand if q[0] <= tol]
            if cand:
                _, a, zz, rel = max(cand, key=lambda q: q[2])
                acept.append((L[c], a, zz, rel))
                ref, ultima = a, c
            c += paso
    acept.sort()
    crudos = pd.DataFrame(crudos, columns=["col", "coord", "z", "relieve"])
    e = pd.DataFrame(acept, columns=["col", "coord", "z", "relieve"]).drop_duplicates("col")
    e["coord_eje"] = e.coord.rolling(cfg.get("mediana", MEDIANA_COL), center=True, min_periods=3).median()
    if cfg["modo"] == "lon":
        e["lon"], e["lat"], e["lat_max"] = e.col, e.coord_eje, e.coord
        crudos["lon"], crudos["lat"] = crudos.col, crudos.coord
    else:
        e["lat"], e["lon"], e["lon_max"] = e.col, e.coord_eje, e.coord
        crudos["lat"], crudos["lon"] = crudos.col, crudos.coord
    # distance to the trench
    dmin = np.array([gc_km(la, lo, fla, flo).min() for la, lo in zip(e.lat, e.lon)])
    e["dist_fosa_km"] = dmin
    e = e[e.dist_fosa_km >= FOSA_KM].reset_index(drop=True)
    return e, crudos


MONTE_RELIEVE_M = 1000.0   # summit above the local median (unsmoothed ETOPO1)
MONTE_RADIO_KM = 100.0     # radius of the local median
MONTE_VENTANA = 21         # 1' cells (~35 km): one summit per window
MONTE_MEDIANA_KM = 150.0   # window of the running median of the Iquique axis
MONTE_PROMEDIO_KM = 100.0  # running mean applied after the median
EXCLUIR_KM = 50.0          # summits closer than this to another ridge are not used for Iquique
# Copiapo Ridge: fixed segment, used only to exclude its seamounts
COPIAPO = ((-28.09, -80.24), (-27.5, -72.21))


FZ_CHALLENGER = os.path.join(ROOT, "models", "Q22_WK08-A", "data", "fz_challenger_gsfml.csv")


def vertices_corredor(cfg):
    """Corridor vertices: those of cfg and, at the end, the intersection of the extension of the
    last segment (towards cfg["hacia"]) with the trace of the Challenger FZ (GSFML)."""
    if "_vert" in cfg:
        return cfg["_vert"]
    f = pd.read_csv(FZ_CHALLENGER, comment="#")
    fl, fo = [], []
    for k in range(len(f) - 1):
        a, o = linea_densa((f.lat[k], f.lon[k]), (f.lat[k + 1], f.lon[k + 1]))
        fl.append(a); fo.append(o)
    fl, fo = np.concatenate(fl), np.concatenate(fo)
    P0, P1 = xyz(*cfg["corredor"][-1]), xyz(*cfg["hacia"])
    om = np.arccos(np.clip(P0 @ P1, -1, 1))
    w = np.linspace(0, 3.0, 6000)
    Q = (np.sin((1 - w) * om)[:, None] * P0 + np.sin(w * om)[:, None] * P1) / np.sin(om)
    la, lo = np.degrees(np.arcsin(Q[:, 2])), np.degrees(np.arctan2(Q[:, 1], Q[:, 0]))
    d = np.array([gc_km(a, o, fl, fo).min() for a, o in zip(la, lo)])
    k = int(np.argmin(d))
    assert d[k] < 5, "the extension does not cross the Challenger FZ"
    cfg["_vert"] = list(cfg["corredor"]) + [(float(la[k]), float(lo[k]))]
    return cfg["_vert"]


def corredor_coords(cfg, la, lo):
    """(s, d, L): distance along the corridor polyline from its NE end (of the nearest node,
    every 1'), minimum distance to it and total length, in km."""
    V = vertices_corredor(cfg)
    pla, plo = [], []
    for k in range(len(V) - 1):
        a, o = linea_densa(V[k], V[k + 1])
        pla.append(a if k == 0 else a[1:]); plo.append(o if k == 0 else o[1:])
    pla, plo = np.concatenate(pla), np.concatenate(plo)
    paso = np.r_[0.0, gc_km(pla[:-1], plo[:-1], pla[1:], plo[1:])]
    sacum = np.cumsum(paso)
    la, lo = np.atleast_1d(np.asarray(la, float)), np.atleast_1d(np.asarray(lo, float))
    s, d = np.empty(len(la)), np.empty(len(la))
    for i, (a, o) in enumerate(zip(la, lo)):
        dd = gc_km(a, o, pla, plo)
        k = int(np.argmin(dd))
        s[i], d[i] = sacum[k], dd[k]
    return s, d, sacum[-1]


def linea_densa(p0, p1, paso=D):
    P0, P1 = xyz(*p0), xyz(*p1)
    om = np.arccos(np.clip(P0 @ P1, -1, 1))
    w = np.linspace(0, 1, max(2, int(np.ceil(np.degrees(om) / paso))))
    Q = (np.sin((1 - w) * om)[:, None] * P0 + np.sin(w * om)[:, None] * P1) / np.sin(om)
    return np.degrees(np.arcsin(Q[:, 2])), np.degrees(np.arctan2(Q[:, 1], Q[:, 0]))


def eje_montes(lon, lat, z, cfg, fla, flo, otros):
    """Iquique axis: summits of unsmoothed ETOPO1 that are a maximum in a window of
    MONTE_VENTANA' and stand >= MONTE_RELIEVE_M above the median within MONTE_RADIO_KM, inside the
    cfg corridor, >= FOSA_KM from the trench and >= EXCLUIR_KM from the ridges in `otros`
    ({name: (lat, lon)}); they are joined in order along the corridor by great circles,
    densified every 1'."""
    from scipy.ndimage import maximum_filter
    # the threshold is measured above the local median; since that median does not fall below
    # ~-6000 m outside the trench (which is excluded), a maximum below -4500 m cannot reach it,
    # so it is discarded first
    pk = np.argwhere((z == maximum_filter(z, size=MONTE_VENTANA)) & (z >= -4500.0))
    la_p, lo_p = lat[pk[:, 0]], lon[pk[:, 1]]
    s, d, L = corredor_coords(cfg, la_p, lo_p)
    ok = (d <= cfg["semiancho_km"]) & (s >= -10) & (s <= L + 10)
    filas = []
    for (j, i), la, lo, ss, dd in zip(pk[ok], la_p[ok], lo_p[ok], s[ok], d[ok]):
        dj = int(np.ceil(MONTE_RADIO_KM / 111.19 / D))
        di = int(np.ceil(dj / np.cos(np.radians(la))))
        jj = slice(max(j - dj, 0), j + dj + 1)
        ii = slice(max(i - di, 0), i + di + 1)
        LO, LA = np.meshgrid(lon[ii], lat[jj])
        disco = gc_km(LA, LO, la, lo) <= MONTE_RADIO_KM
        med = np.median(z[jj, ii][disco])
        filas.append(dict(lat=la, lon=lo, z=z[j, i], relieve=z[j, i] - med, s_km=ss, d_km=dd))
    m = pd.DataFrame(filas)
    m = m[m.relieve >= MONTE_RELIEVE_M].copy()
    m["dist_fosa_km"] = [gc_km(a, o, fla, flo).min() for a, o in zip(m.lat, m.lon)]
    for nom, (ola, olo) in otros.items():
        m[f"dist_{nom}_km"] = [gc_km(a, o, ola, olo).min() for a, o in zip(m.lat, m.lon)]
    excl = np.zeros(len(m), bool)
    for nom in otros:
        excl |= m[f"dist_{nom}_km"].values < EXCLUIR_KM
    m["excluida"] = excl | (m.dist_fosa_km < FOSA_KM)
    m = m.drop_duplicates(["lat", "lon"]).sort_values("s_km").reset_index(drop=True)
    m.to_csv(os.path.join(OUT, "cumbres_Iquique.csv"), index=False)
    m = m[~m.excluida].reset_index(drop=True)
    la_e, lo_e = [], []
    for k in range(len(m) - 1):
        a, o = linea_densa((m.lat[k], m.lon[k]), (m.lat[k + 1], m.lon[k + 1]))
        la_e.append(a[:-1]); lo_e.append(o[:-1])
    la_e.append(m.lat.values[-1:]); lo_e.append(m.lon.values[-1:])
    e = pd.DataFrame({"lat": np.concatenate(la_e), "lon": np.concatenate(lo_e)})
    # running median of lat and lon along the axis: removes the zigzag between neighbouring
    # summits; window of MONTE_MEDIANA_KM, with the mean step of the densified line
    paso = np.mean(gc_km(e.lat.values[:-1], e.lon.values[:-1], e.lat.values[1:], e.lon.values[1:]))
    nv = max(3, int(round(MONTE_MEDIANA_KM / paso)) | 1)
    e["lat_bruta"], e["lon_bruta"] = e.lat, e.lon
    # median (removes jumps to neighbouring summits) and then running mean of MONTE_PROMEDIO_KM
    # (removes the steps left by the separate median in lat and lon)
    npm = max(3, int(round(MONTE_PROMEDIO_KM / paso)) | 1)
    for c in ("lat", "lon"):
        m_ = e[c + "_bruta"].rolling(nv, center=True, min_periods=1).median()
        e[c] = m_.rolling(npm, center=True, min_periods=1).mean()
    e["col"] = D * np.arange(len(e))           # no gaps: the line is continuous
    e["dist_fosa_km"] = [gc_km(a, o, fla, flo).min() for a, o in zip(e.lat, e.lon)]
    e = e[e.dist_fosa_km >= FOSA_KM].reset_index(drop=True)
    crudos = m.assign(col=m.lon, coord=m.lat)
    return e, crudos


def track_denso(t):
    t = t.dropna(subset=["lat", "lon"]).sort_values("tau")
    P = xyz(t.lat.values, t.lon.values)
    tau = t.tau.values
    out_t, out_p = [], []
    for k in range(len(tau) - 1):
        if tau[k + 1] - tau[k] > 1.01:          # segment with a gap (plate change)
            continue
        w = np.linspace(0, 1, 21)[:-1]
        om = np.arccos(np.clip(P[k] @ P[k + 1], -1, 1))
        if om < 1e-12:
            continue
        seg = (np.sin((1 - w) * om)[:, None] * P[k] + np.sin(w * om)[:, None] * P[k + 1]) / np.sin(om)
        out_p.append(seg)
        out_t.append(tau[k] + w * (tau[k + 1] - tau[k]))
    out_p.append(P[-1:])
    out_t.append(tau[-1:])
    return np.concatenate(out_t), np.concatenate(out_p)


def elipse_en(ell, tau):
    e = ell.sort_values("tau")
    ts = e.tau.values.astype(float)
    fuera = bool(tau < ts.min() or tau > ts.max())
    tt = np.clip(tau, ts.min(), ts.max())
    az = np.unwrap(np.radians(2 * e.az_major_deg.values)) / 2   # undirected axis: period 180
    return (np.interp(tt, ts, e.a68_km), np.interp(tt, ts, e.b68_km),
            np.interp(tt, ts, e.a95_km), np.interp(tt, ts, e.b95_km),
            np.degrees(np.interp(tt, ts, az)), fuera)


def main():
    os.makedirs(OUT, exist_ok=True)
    lon, lat, z = etopo_box()
    zs = {fw: gaussian_filter(z, sigma=fw / 1.85 / 2.355)
          for fw in sorted({c["fwhm"] for c in DORSALES.values()})}
    fla, flo = fosa()
    ejes, crudos = {}, {}
    orden = [n for n, c in DORSALES.items() if c["modo"] != "montes"] +             [n for n, c in DORSALES.items() if c["modo"] == "montes"]
    iq = DORSALES["Iquique"]
    ejes_completos = {}
    for nom in orden:
        cfg = DORSALES[nom]
        if cfg["modo"] == "montes":
            # Taltal: excluded with its complete axis, before clipping
            otros = {o: (ejes_completos.get(o, ejes[o]).lat.values, ejes_completos.get(o, ejes[o]).lon.values)
                     for o in cfg["excluir"] if o in ejes}
            otros["Copiapo"] = linea_densa(*COPIAPO)
            e, c = eje_montes(lon, lat, z, cfg, fla, flo, otros)
        else:
            e, c = eje(lon, lat, zs, cfg, fla, flo)
            if nom == "Taltal":   # drop points inside the Iquique corridor
                ejes_completos["Taltal"] = e
                _, dd, _ = corredor_coords(iq, e.lat.values, e.lon.values)
                e = e[dd > iq["semiancho_km"]].reset_index(drop=True)
        if cfg["edades"]:
            ed = pd.read_csv(os.path.join(ROOT, "models", "Q22_WK08-A", "data", cfg["edades"]), comment="#")
            e["dist_dato_km"] = [gc_km(la, lo, ed.lat.values, ed.lon.values).min()
                                 for la, lo in zip(e.lat, e.lon)]
        else:
            e["dist_dato_km"] = np.inf
        e["con_datos"] = e.dist_dato_km < DATOS_KM
        e["cerca_fosa"] = e.dist_fosa_km < CERCA_FOSA_KM
        e.to_csv(os.path.join(OUT, f"ejes_{nom}.csv"), index=False)
        ejes[nom], crudos[nom] = e, c
        print(f"{nom:12s} axis: {len(e)} columns, lon {e.lon.min():.2f} to {e.lon.max():.2f}, "
              f"lat {e.lat.min():.2f} to {e.lat.max():.2f}, with data {int(e.con_datos.sum())}")
    filas = []
    for carpeta, circ, marco in TODOS:
        nomv = pd.read_csv(os.path.join(ROOT, "models", carpeta, f"track_nominal_0.44{TAG}.csv"))
        ellv = pd.read_csv(os.path.join(ROOT, "models", carpeta, f"track_ellipses_0.44{TAG}.csv"))
        texto = (carpeta, circ, marco) in MP.TEXTO
        for nom, cfg in DORSALES.items():
            t = nomv[nomv.ridge == cfg["track"]]
            tau, P = track_denso(t)
            ell = ellv[ellv.ridge == cfg["track"]]
            e = ejes[nom]
            Q = xyz(e.lat.values, e.lon.values)
            for q, (_, r) in zip(Q, e.iterrows()):
                cosd = P @ q
                k = int(np.argmax(cosd))
                d = R_KM * np.arccos(np.clip(cosd[k], -1, 1))
                kk = min(max(k, 1), len(P) - 2)
                tdir = P[kk + 1] - P[kk - 1]
                signo = np.sign(np.cross(P[k], tdir) @ q) or 1.0
                pla = np.degrees(np.arcsin(P[k][2]))
                plo = np.degrees(np.arctan2(P[k][1], P[k][0]))
                a68, b68, a95, b95, az, fuera = elipse_en(ell, tau[k])
                th = np.radians(azimut(pla, plo, r.lat, r.lon) - az)
                x, y = d * np.cos(th), d * np.sin(th)
                filas.append(dict(dorsal=nom, carpeta=carpeta, circuito=circ, marco=marco,
                                  texto=texto, lon_eje=r.lon, lat_eje=r.lat, con_datos=r.con_datos,
                                  cerca_fosa=r.cerca_fosa,
                                  tau=tau[k], extremo=k in (0, len(P) - 1),
                                  dist_km=d, dist_signo_km=signo * d,
                                  dentro68=(x / a68) ** 2 + (y / b68) ** 2 <= 1,
                                  dentro95=(x / a95) ** 2 + (y / b95) ** 2 <= 1,
                                  elipse_extrapolada=fuera))
    df = pd.DataFrame(filas)
    df.to_csv(os.path.join(OUT, "distancias.csv"), index=False)
    v = df[~df.extremo]
    res = []
    for (nom, carpeta, circ, marco, texto), g in v.groupby(
            ["dorsal", "carpeta", "circuito", "marco", "texto"], sort=False):
        for tramo, h in (("todo", g), ("sin_datos", g[~g.con_datos]), ("con_datos", g[g.con_datos]),
                         ("sin_datos_lejos_fosa", g[~g.con_datos & ~g.cerca_fosa])):
            if len(h) == 0:
                continue
            res.append(dict(dorsal=nom, tramo=tramo, carpeta=carpeta, circuito=circ, marco=marco,
                            texto=texto, n=len(h), n_total=int((df.dorsal == nom).sum() / len(TODOS)),
                            tau_min=h.tau.min(), tau_max=h.tau.max(),
                            mediana_km=h.dist_km.median(), p90_km=h.dist_km.quantile(0.9),
                            mediana_signo_km=h.dist_signo_km.median(),
                            frac68=h.dentro68.mean(), frac95=h.dentro95.mean(),
                            frac_elipse_extrap=h.elipse_extrapolada.mean()))
    R = pd.DataFrame(res)
    R.to_csv(os.path.join(OUT, "resumen.csv"), index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 500):
        print(R[R.tramo != "con_datos"].round(2).to_string(index=False))
    mapa(lon, lat, zs, ejes, crudos)


def mapa(lon, lat, zss, ejes, crudos):
    zs = zss[20.0]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col_m = {"WK08-A": "#1f4e79", "GMHRF-B": "#b5651d", "GWG20": "#2e7d32", "GMHRF-A'": "#7b3294"}
    est_c = {"Quiero": "-", "WD25": "--", "SG12": ":"}
    paneles = [("Nazca", (-83, -75.5, -26, -14.3)), ("Iquique", (-86.5, -70, -36.5, -19)),
               ("Taltal", (-81, -71, -27.6, -24.2)), ("JFR", (-82, -71.5, -34.8, -32.2))]
    fig, axs = plt.subplots(2, 2, figsize=(14, 11))
    for ax, (nom, (x0, x1, y0, y1)) in zip(axs.ravel(), paneles):
        ii = (lon >= x0) & (lon <= x1)
        jj = (lat >= y0) & (lat <= y1)
        ax.pcolormesh(lon[ii], lat[jj], zs[np.ix_(jj, ii)], cmap="Greys_r", vmin=-5500, vmax=-1500,
                      shading="auto")
        ax.contour(lon[ii], lat[jj], zs[np.ix_(jj, ii)], levels=[-4000, -3500, -3000, -2500],
                   colors="0.45", linewidths=0.3)
        c = crudos[nom]
        ax.plot(c.lon, c.lat, ".", ms=1.2, color="#4fc3f7", zorder=2)
        e = ejes[nom]
        tramo = (e.col.diff().abs() > HUECO_COL * D * 1.01).cumsum()
        for k, (_, g) in enumerate(e.groupby(tramo)):
            ax.plot(g.lon, g.lat, "-", color="red", lw=2.2, zorder=4,
                    label="ridge axis (this test)" if k == 0 else None)
        ax.plot(e.lon[e.con_datos], e.lat[e.con_datos], "o", ms=2.5, color="yellow", zorder=5,
                label=f"axis < {DATOS_KM:.0f} km from a dated sample")
        cfg = DORSALES[nom]
        for carpeta, circ, marco in MP.TEXTO:
            t = pd.read_csv(os.path.join(ROOT, "models", carpeta, f"track_nominal_0.44{TAG}.csv"))
            t = t[t.ridge == cfg["track"]].dropna()
            ax.plot(t.lon, t.lat, est_c[circ], color=col_m[marco], lw=1.1, zorder=3,
                    label=f"{circ}+{marco}")
        x_0, x_1, y_0, y_1 = cfg["caja"]
        ax.plot([x_0, x_1, x_1, x_0, x_0], [y_0, y_0, y_1, y_1, y_0], "c-", lw=0.8)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect(1 / np.cos(np.radians((y0 + y1) / 2)))
        ax.set_title(nom if nom != "Taltal" else "Taltal (San Félix track)")
        ax.tick_params(labelsize=7)
    axs[0, 0].legend(fontsize=6, loc="lower right")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"mapa_ejes.{ext}"), dpi=150)


if __name__ == "__main__":
    main()
