# -*- coding: utf-8 -*-
"""Position of the nominal Easter track (f_PAC 0.44, present position) relative to the
Nazca Ridge in the ETOPO1 bathymetry, in the models listed in MODELOS (current _hsb .rot files, w conf).

For each track point on the ridge (west of the trench) a profile of +-300 km perpendicular
to the track is taken in smoothed ETOPO1 (Gaussian of ~20 km), and the script measures:
  - the signed distance from the track to the crest (profile maximum) and to the ridge center
    (centroid of the relief above half the height between the crest and the regional floor); positive =
    track to the NW (to the left of the direction of increasing tau, which points to the trench in the NE);
  - the ridge edges, where the profile crosses half the relief between the crest and the
    regional floor (median of the outermost 60 km of the profile), and the track position
    relative to them.

Outputs: outputs/track_vs_ridge_axis_g13.csv (per model and tau) and an on-screen summary.
Usage:  py -I track_vs_ridge_axis.py
"""
import os
import sys
import importlib.util
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, map_coordinates

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
# ETOPO1 bedrock grid (GeoTIFF)
ETOPO = (__import__("os").environ["PAPER_ROOT"] + "/external/etopo1/ETOPO1_Bed_g_geotiff.tif")
# the nine models of the text (g13, conf .rot files)
MODELOS = {  # folder: (conf rot, name in the text)
    "Q22_WK08-A": ("Global_Model_R_BG_2025_hsb.rot", "Quiero+WK08-A"),
    "Q22_GMHRF-B": ("Global_Model_R_BG_2025_GMHRF_gts12_hsbgmhrf.rot", "Quiero+GMHRF-B"),
    "Q22_GWG20": ("Global_Model_R_BG_2025_GWG20_hsb.rot", "Quiero+GWG20"),
    "WD25_WK08-A": ("Global_Model_R_BG_2025_WD25_hsb.rot", "WD25+WK08-A"),
    "WD25_GMHRF-B": ("Global_Model_R_BG_2025_WD25_GMHRF_gts12_hsbgmhrf.rot",
                                         "WD25+GMHRF-B"),
    "WD25_GWG20": ("Global_Model_R_BG_2025_WD25_GWG20_hsb.rot", "WD25+GWG20"),
    "SG12_WK08-A": ("Global_Model_R_BG_2025_SG12_hsb.rot", "SG12+WK08-A"),
    "SG12_GMHRF-B": ("Global_Model_R_BG_2025_SG12_GMHRF_gts12_hsbgmhrf.rot",
                                         "SG12+GMHRF-B"),
    "SG12_GWG20": ("Global_Model_R_BG_2025_SG12_GWG20_hsb.rot", "SG12+GWG20"),
}
R_KM = 6371.0
# ETOPO1 "g" (grid registered): 21601 x 10801 int16, uncompressed, one row per strip
NX, NY, OFF0 = 21601, 10801, 87470
LON0, LAT0, D = -180.0, 90.0, 1 / 60.0
BOX = (-90.0, -70.0, -32.0, -10.0)   # lon_min, lon_max, lat_min, lat_max


def load_mod(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def etopo_box():
    with open(ETOPO, "rb") as f:
        assert f.read(2) == b"II"
    z = np.memmap(ETOPO, dtype="<i2", mode="r", offset=OFF0, shape=(NY, NX))
    i0, i1 = int(round((BOX[0] - LON0) / D)), int(round((BOX[1] - LON0) / D))
    j0, j1 = int(round((LAT0 - BOX[3]) / D)), int(round((LAT0 - BOX[2]) / D))
    sub = np.array(z[j0:j1 + 1, i0:i1 + 1], dtype=float)
    lon = LON0 + D * np.arange(i0, i1 + 1)
    lat = LAT0 - D * np.arange(j0, j1 + 1)
    return lon, lat, sub


def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def ll(v):
    v = v / np.linalg.norm(v)
    return np.degrees(np.arcsin(v[2])), np.degrees(np.arctan2(v[1], v[0]))


def main():
    lon, lat, z = etopo_box()
    # Gaussian smoothing, 20 km full width at half maximum (1' ~ 1.85 km in latitude)
    zs = gaussian_filter(z, sigma=20 / 1.85 / 2.355)

    def sample(la, lo):
        jj = (lat[0] - np.asarray(la)) / D
        ii = (np.asarray(lo) - lon[0]) / D
        return map_coordinates(zs, [jj, ii], order=1, mode="nearest")

    sys.path.insert(0, os.path.join(ROOT, "models", "Q22_WK08-A"))
    eng = load_mod(os.path.join(ROOT, "models", "Q22_WK08-A", "02_rot_engine.py"), "02_rot_engine")
    mc = load_mod(os.path.join(ROOT, "models", "Q22_WK08-A", "04_mc_collisions.py"), "04_mc_collisions")
    seed, pid = mc.RIDGES["Easter"]
    s_km = np.linspace(-300, 300, 241)
    rows = []
    for d, (rot, nombre) in MODELOS.items():
        M = eng.RotModel(os.path.join(ROOT, "models", d, "rot", rot))
        tot = mc.make_tot(M, None)
        P = mc.track_points(tot, seed, pid, 902003)
        la, lo = mc.project(tot, 902003, 0.0, P)
        for k in range(1, len(mc.TAUS) - 1):
            tau = mc.TAUS[k]
            if not (-26.5 <= la[k] <= -14.5 and -84.0 <= lo[k] <= -76.0):
                continue
            p, pa, pb = xyz(la[k], lo[k]), xyz(la[k - 1], lo[k - 1]), xyz(la[k + 1], lo[k + 1])
            t = pb - pa
            t -= p * (p @ t)
            t /= np.linalg.norm(t)          # increasing tau
            n = np.cross(p, t)              # left of the direction of increasing tau
            pts = [ll(np.cos(s / R_KM) * p + np.sin(s / R_KM) * n) for s in s_km]
            prof = sample([q[0] for q in pts], [q[1] for q in pts])
            z0 = sample([la[k]], [lo[k]])[0]
            if z0 > -500:                    # already on the continent or the inner trench
                continue
            ic = int(np.argmax(prof))
            fondo = np.median(np.r_[prof[:24], prof[-24:]])
            half = fondo + 0.5 * (prof[ic] - fondo)
            iL = ic
            while iL > 0 and prof[iL] > half:
                iL -= 1
            iR = ic
            while iR < len(prof) - 1 and prof[iR] > half:
                iR += 1
            w = np.clip(prof[iL:iR + 1] - half, 0, None)
            centro = float(np.sum(w * s_km[iL:iR + 1]) / np.sum(w))
            rows.append(dict(modelo=nombre, carpeta=d, tau=tau, lat=la[k], lon=lo[k],
                             z_track=z0, z_cresta=prof[ic], fondo=fondo,
                             track_menos_cresta_km=-s_km[ic], track_menos_centro_km=-centro,
                             borde_SE_km=s_km[iL], borde_NO_km=s_km[iR],
                             ancho_km=s_km[iR] - s_km[iL]))
    df = pd.DataFrame(rows)
    out = os.path.join(HERE, "outputs", "track_vs_ridge_axis_g13.csv")
    df.to_csv(out, index=False)
    # distance from the track (s=0) to the NW edge: positive if the track lies outside, NW of the edge
    df["fuera_NO_km"] = -df["borde_NO_km"]
    print("track - ridge center (km, + = track to the NW); center = centroid of the relief above half height")
    for nombre, g in df.groupby("modelo", sort=False):
        print(f"  {nombre:13s} n={len(g):2d} tau {g.tau.min():4.0f}-{g.tau.max():4.0f} Ma | "
              f"track-center median {g.track_menos_centro_km.median():5.0f} "
              f"[{g.track_menos_centro_km.min():5.0f}, {g.track_menos_centro_km.max():5.0f}] | "
              f"points outside the NW edge {int((g.borde_NO_km < 0).sum())}/{len(g)} | "
              f"median width {g.ancho_km.median():4.0f} km")
    print(f"written {out}")


if __name__ == "__main__":
    main()
