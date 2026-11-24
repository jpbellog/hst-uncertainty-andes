# -*- coding: utf-8 -*-
"""Writes the collision bands RE-RESTORED TO THE PRESENT.

Reads <model folder>/mc_collision_bands_{variant}.csv
(step 04: collision computed against the margin restored to age t, with
South America fixed in its present-day position) and maps the LATITUDE
percentiles of each row to the present-day margin through the node
correspondence of Margin_SAM_0.1.csv (see margin_map.py). Since the mapping is
monotonic along the margin in the range of collision latitudes (the only
non-monotonic segment of the polygon, ~40.9 S, is out of range), the mapped
percentiles are the percentiles of the mapped distribution.

Output: tables/mc_collision_bands_{variant}_present.csv with columns
ridge, pf, age, n, lat_nominal, lon_nominal, lat_p2_5, lat_p16, lat_p50,
lat_p84, lat_p97_5, lon_p50 (lon_* = longitude of the present-day margin at the
corresponding mapped latitude).

Usage: py make_present_bands.py
"""
import os
import numpy as np
import pandas as pd
import margin_map as MM

BASE = os.path.dirname(os.path.abspath(__file__))
# Model folder (sigma_omega of JFR/San Felix = 6.7 km/Myr, recalculated by L.E. Lara).
# The paper compares TWO kinematic models. The environment variables PAPER_PIPE
# (source model folder) and PAPER_BANDS_SUFFIX (suffix of the output CSVs) also
# produce the Wilson and DeMets (2025) version without touching the base model:
#   PAPER_PIPE=WD25_WK08-A PAPER_BANDS_SUFFIX=_wd25 py make_present_bands.py
PIPE = os.path.join(os.path.dirname(BASE), "models", os.environ.get("PAPER_PIPE", "Q22_WK08-A"))
SUFFIX = os.environ.get("PAPER_BANDS_SUFFIX", "")
# PAPER_FIG_TAG is the tag of the step 04 outputs (e.g. _hsb_cc_g13); it is appended
# to the input and output names. PAPER_TABLES_DIR changes the output folder.
TAG = os.environ.get("PAPER_FIG_TAG", "")
OUT_DIR = os.environ.get("PAPER_TABLES_DIR", BASE)
VARIANTS = ("total", "circuit", "wk08", "geochron", "hotspot")
LAT_COLS = ["lat_nominal", "lat_p2_5", "lat_p16", "lat_p50", "lat_p84", "lat_p97_5"]

for variant in VARIANTS:
    src = os.path.join(PIPE, f"mc_collision_bands_{variant}{TAG}.csv")
    if not os.path.exists(src):          # e.g. geochron does not exist with PAPER_HS_FIT
        print(f"{src} not found, skipped")
        continue
    d = pd.read_csv(src, dtype={"pf": str})
    out = []
    for _, r in d.iterrows():
        row = dict(ridge=r.ridge, pf=r.pf, age=r.age, n=r.n)
        for c in LAT_COLS:
            la, lo = MM.map_lat_to_present(r.age, r[c])
            row[c] = la
            if c == "lat_nominal":
                row["lon_nominal"] = lo
            if c == "lat_p50":
                row["lon_p50"] = lo
        out.append(row)
    out = pd.DataFrame(out)[["ridge", "pf", "age", "n", "lat_nominal", "lon_nominal",
                             "lat_p2_5", "lat_p16", "lat_p50", "lat_p84",
                             "lat_p97_5", "lon_p50"]]
    fn = os.path.join(OUT_DIR, f"mc_collision_bands_{variant}_present{SUFFIX}{TAG}.csv")
    out.to_csv(fn, index=False)
    nna = out[LAT_COLS].isna().sum().sum() - int(d[LAT_COLS].isna().sum().sum())
    print(f"written {fn} ({len(out)} rows; new NaN from mapping: {nna})")
