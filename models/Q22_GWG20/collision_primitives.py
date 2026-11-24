# -*- coding: utf-8 -*-
"""Geodesic primitives and track-trench intersection algorithm.

Python port of the MATLAB script Inter_Coor_Margin_TT.m, which computes the
collision latitude of each aseismic ridge (Easter/Nazca, San Felix/Taltal,
Juan Fernandez) with the trench, by age, by intersecting the great circle of the
margin segment with that of the hotspot track segment.

Two routes:
  - run_algorithm(): EXACT replica of the .m, reads an already assembled
    *_margin.xlsx table (preselected margin-track pairs). Used to VALIDATE the
    port against the existing inter_fosa_*_{0.5,0.36}.xlsx.
  - collide_strict() + margin_segments(): GENERAL method that does not depend on
    the preselection in _margin.xlsx: it intersects the reconstructed track with
    the FULL MARGIN (Margin_SAM, independent of the f_PAC model) with strict
    selection (the intersection must fall INSIDE the margin segment and the
    track segment). Used to build the line of the Pf=0.44 model, whose
    inter_fosa tables were not precomputed. It is NOT assumed that 0.44 falls
    between 0.5 and 0.36.

Column convention of the *_margin.xlsx tables (same as in the .m):
  edad, edad_track1, edad_track2, lat_m1, lon_m1, lat_m2, lon_m2,
  lat_t1, lon_t1, lat_t2, lon_t2
(edad = age)
"""
import os
import numpy as np
import pandas as pd

D2R = np.pi / 180.0

# ----------------------------------------------------------------- primitives
def azimuth(lat1, lon1, lat2, lon2):
    """Initial great-circle azimuth P1->P2 (degrees), as in the Mapping Toolbox."""
    la1, la2 = lat1 * D2R, lat2 * D2R
    dlon = (lon2 - lon1) * D2R
    y = np.sin(dlon) * np.cos(la2)
    x = np.cos(la1) * np.sin(la2) - np.sin(la1) * np.cos(la2) * np.cos(dlon)
    return np.degrees(np.arctan2(y, x)) % 360.0

def _unit(lat, lon):
    la, lo = lat * D2R, lon * D2R
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])

def _gc_normal(lat, lon, az):
    """Normal (pole) of the great circle through (lat,lon) with azimuth az."""
    la, lo, a = lat * D2R, lon * D2R, az * D2R
    p = _unit(lat, lon)
    east = np.array([-np.sin(lo), np.cos(lo), 0.0])
    north = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    t = np.sin(a) * east + np.cos(a) * north
    n = np.cross(p, t)
    return n / np.linalg.norm(n)

def gcxgc(lat1, lon1, az1, lat2, lon2, az2):
    """Intersection of two great circles (point+azimuth each).
    Returns (lat_a, lat_b, lon_a, lon_b) as MATLAB gcxgc (antipodal pair)."""
    n1 = _gc_normal(lat1, lon1, az1)
    n2 = _gc_normal(lat2, lon2, az2)
    d = np.cross(n1, n2)
    nd = np.linalg.norm(d)
    if nd < 1e-12:
        return (np.nan,) * 4
    d = d / nd
    lat_a = np.degrees(np.arcsin(d[2]))
    lon_a = np.degrees(np.arctan2(d[1], d[0]))
    return (lat_a, -lat_a, lon_a, (lon_a + 180) % 360 - 180)

def distance(lat1, lon1, lat2, lon2):
    """Great-circle distance (degrees)."""
    la1, la2 = lat1 * D2R, lat2 * D2R
    dlon = (lon2 - lon1) * D2R
    d = np.arccos(np.clip(np.sin(la1) * np.sin(la2) +
                          np.cos(la1) * np.cos(la2) * np.cos(dlon), -1, 1))
    return np.degrees(d)

# ------------------------------------------------ EXACT algorithm (validation)
COLS = ["edad", "et1", "et2", "lat_m1", "lon_m1", "lat_m2", "lon_m2",
        "lat_t1", "lon_t1", "lat_t2", "lon_t2"]

def run_algorithm(coord):
    """Exact replica of the loop of Inter_Coor_Margin_TT.m over an already
    assembled table (columns COLS). Returns {age: (lat, lon)}."""
    c = coord[COLS].to_numpy(float)
    out, mind = {}, {}
    for i in range(len(c)):
        ed = c[i, 0]
        azm = azimuth(c[i, 5], c[i, 6], c[i, 3], c[i, 4])   # margin (m2)->(m1)
        azt = azimuth(c[i, 7], c[i, 8], c[i, 9], c[i, 10])  # track  (t1)->(t2)
        la, lb, lna, lnb = gcxgc(c[i, 5], c[i, 6], azm, c[i, 7], c[i, 8], azt)
        if not np.isfinite(la):
            continue
        lat = la if abs(la - c[i, 3]) < abs(lb - c[i, 3]) else lb
        lon = min(lna, lnb)
        between = (c[i, 8] < lon < c[i, 10]) or (c[i, 10] < lon < c[i, 8])
        if between:
            d1 = distance(c[i, 7], c[i, 8], lat, lon)
            d2 = distance(c[i, 9], c[i, 10], lat, lon)
            if (d1 + d2) < mind.get(ed, 1e9):
                mind[ed] = d1 + d2
                out[ed] = (lat, lon)
    return out

# --------------------------------------- GENERAL method (full margin, 0.44)
def load_full_margin(margin_csv):
    """Loads Margin_SAM_0.1.csv and returns (df, n_nodes)."""
    m = pd.read_csv(margin_csv)
    m["age"] = pd.to_numeric(m["age"], errors="coerce")
    nn = sum(1 for c in m.columns if c.startswith("Lon_"))
    return m, nn

def margin_segments(mdf, nn, age):
    """Consecutive segments of the FULL margin at the given age."""
    row = mdf.iloc[(mdf["age"] - age).abs().idxmin()]
    lon = np.array([row.get(f"Lon_{k}") for k in range(1, nn + 1)], float)
    lat = np.array([row.get(f"Lat_{k}") for k in range(1, nn + 1)], float)
    segs = []
    for j in range(len(lat) - 1):
        if np.isfinite(lat[j]) and np.isfinite(lat[j + 1]):
            segs.append((lat[j], lon[j], lat[j + 1], lon[j + 1]))
    return np.array(segs)

def _units(lat, lon):
    """Unit vectors (N,3) for lat/lon arrays (degrees)."""
    la = np.asarray(lat, float) * D2R
    lo = np.asarray(lon, float) * D2R
    return np.stack([np.cos(la) * np.cos(lo),
                     np.cos(la) * np.sin(lo),
                     np.sin(la)], axis=-1)

def collide_strict(msegs, tsegs, tol_lat=0.3):
    """Best pair (margin segment x track segment) per age with STRICT
    selection: the intersection must fall within the lat range of the margin
    segment (+-tol_lat) and within the lon range of the track segment.
    Vectorised version (numpy): the plane of each great circle is the cross
    product of the unit vectors of its two endpoints (identical to
    point+azimuth). Returns (lat, lon) or None.

    msegs: (Nm,4) [lat_m1,lon_m1,lat_m2,lon_m2]
    tsegs: (Nt,6) [et1,lat_t1,lon_t1,et2,lat_t2,lon_t2]"""
    msegs = np.asarray(msegs, float); tsegs = np.asarray(tsegs, float)
    if msegs.size == 0 or tsegs.size == 0:
        return None
    m1 = _units(msegs[:, 0], msegs[:, 1]); m2 = _units(msegs[:, 2], msegs[:, 3])
    t1 = _units(tsegs[:, 1], tsegs[:, 2]); t2 = _units(tsegs[:, 4], tsegs[:, 5])
    nm = np.cross(m1, m2); nm /= np.linalg.norm(nm, axis=1, keepdims=True)
    nt = np.cross(t1, t2); nt /= np.linalg.norm(nt, axis=1, keepdims=True)
    D = np.cross(nm[:, None, :], nt[None, :, :])          # (Nm,Nt,3)
    nrm = np.linalg.norm(D, axis=2, keepdims=True)
    nrm[nrm < 1e-12] = np.nan
    D = D / nrm
    m_lo = np.minimum(msegs[:, 0], msegs[:, 2]); m_hi = np.maximum(msegs[:, 0], msegs[:, 2])
    t_lo = np.minimum(tsegs[:, 2], tsegs[:, 5]); t_hi = np.maximum(tsegs[:, 2], tsegs[:, 5])
    best, bestd = None, np.inf
    for sign in (1.0, -1.0):
        P = sign * D
        lat = np.degrees(np.arcsin(np.clip(P[..., 2], -1, 1)))    # (Nm,Nt)
        lon = np.degrees(np.arctan2(P[..., 1], P[..., 0]))
        in_m = (lat >= (m_lo[:, None] - tol_lat)) & (lat <= (m_hi[:, None] + tol_lat))
        in_t = (lon > t_lo[None, :]) & (lon < t_hi[None, :])
        ok = in_m & in_t & np.isfinite(lat)
        if not ok.any():
            continue
        d1 = np.degrees(np.arccos(np.clip(np.einsum('mnk,nk->mn', P, t1), -1, 1)))
        d2 = np.degrees(np.arccos(np.clip(np.einsum('mnk,nk->mn', P, t2), -1, 1)))
        dsum = np.where(ok, d1 + d2, np.inf)
        idx = np.unravel_index(np.argmin(dsum), dsum.shape)
        if dsum[idx] < bestd:
            bestd = dsum[idx]; best = (float(lat[idx]), float(lon[idx]))
    return best

# -------------------------------------------------------------- track readers
def track_segs_from_marginxlsx(fn):
    """Track segments per age, extracted from a *_margin.xlsx table.
    Returns {age: ndarray[(et1,lat1,lon1,et2,lat2,lon2)]}."""
    m = pd.read_excel(fn, header=None, skiprows=1).iloc[:, :11]
    m.columns = COLS
    m = m.apply(pd.to_numeric, errors="coerce").dropna(subset=["edad"])
    out = {}
    for ed, sub in m.groupby("edad"):
        out[ed] = sub[["et1", "lat_t1", "lon_t1", "et2", "lat_t2", "lon_t2"]] \
            .round(4).drop_duplicates().to_numpy()
    return out

def track_segs_from_xy(folder, age):
    """Segments of the reconstructed track read from motion_path_output_{age}Ma.xy
    (GPlates format: data lines 'Lon Lat SeamountAge', headers with '>')."""
    fn = os.path.join(folder, f"motion_path_output_{age:.2f}Ma.xy")
    if not os.path.exists(fn):
        return None
    lon, lat, sa = [], [], []
    for line in open(fn):
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        p = line.split()
        if len(p) >= 3:
            lon.append(float(p[0])); lat.append(float(p[1])); sa.append(float(p[2]))
    lon, lat, sa = np.array(lon), np.array(lat), np.array(sa)
    segs = []
    for k in range(len(lat) - 1):
        segs.append((sa[k], lat[k], lon[k], sa[k + 1], lat[k + 1], lon[k + 1]))
    return np.array(segs)
