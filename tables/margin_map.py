# -*- coding: utf-8 -*-
"""Mapping of positions on the margin restored to age t onto the PRESENT-DAY
margin (age 0), through the node correspondence of Margin_SAM_0.1.csv.

Step 04 of the pipeline computes the track-trench collision in the frame of South
America fixed in its present-day position, but against the margin RESTORED to the
collision age (same node indexing for every age in the CSV). This module
re-restores to the present: it locates the collision latitude between two margin
nodes at age t and returns the position interpolated between the SAME nodes at
age 0. Since the mapping is monotonic along the margin, applied to latitude
percentiles it gives the percentiles of the mapped distribution.
"""
import os
import numpy as np
import pandas as pd

MARGIN_CSV = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "margin", "Margin_SAM_0.1.csv")  # copy of Margin_SAM (see rotation_models/README.md)

_cache = {}

def _load():
    if "m" not in _cache:
        m = pd.read_csv(MARGIN_CSV)
        m["age"] = pd.to_numeric(m["age"], errors="coerce")
        nn = sum(1 for c in m.columns if c.startswith("Lon_"))
        _cache["m"], _cache["nn"] = m, nn
    return _cache["m"], _cache["nn"]

def _nodes(age):
    m, nn = _load()
    row = m.iloc[(m["age"] - age).abs().idxmin()]
    lat = np.array([row.get(f"Lat_{k}") for k in range(1, nn + 1)], float)
    lon = np.array([row.get(f"Lon_{k}") for k in range(1, nn + 1)], float)
    ok = np.isfinite(lat) & np.isfinite(lon)
    return lat[ok], lon[ok]

def check_monotonic(age):
    """True if the latitude of the margin nodes changes monotonically."""
    lat, _ = _nodes(age)
    return bool(np.all(np.diff(lat) < 0) or np.all(np.diff(lat) > 0))

def map_lat_to_present(age, lat_t):
    """Latitude (and longitude) on the margin at age `age` -> corresponding
    position on the present-day margin. Returns (lat_0, lon_0).
    NaN if lat_t falls outside the range of the margin at that age."""
    if not np.isfinite(lat_t):
        return np.nan, np.nan
    la_t, lo_t = _nodes(age)
    la_0, lo_0 = _nodes(0.0)
    n = min(len(la_t), len(la_0))
    la_t, lo_t, la_0, lo_0 = la_t[:n], lo_t[:n], la_0[:n], lo_0[:n]
    for k in range(n - 1):
        a, b = la_t[k], la_t[k + 1]
        if (a >= lat_t >= b) or (b >= lat_t >= a):
            f = 0.0 if a == b else (lat_t - a) / (b - a)
            return (la_0[k] + f * (la_0[k + 1] - la_0[k]),
                    lo_0[k] + f * (lo_0[k + 1] - lo_0[k]))
    return np.nan, np.nan
