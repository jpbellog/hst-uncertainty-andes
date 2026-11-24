# -*- coding: utf-8 -*-
"""Nominal track-trench collision latitude (f_PAC 0.44, unperturbed) per chain and age
(0-66 Ma), with the same construction as step 04 (track_points, project, collide_strict).
Run it with the environment of a model (PAPER_ROT_FILE, PAPER_HS_FIT) and the folder of its scripts.

Usage:  py nominal_collision.py <model_folder> <output.csv>
"""
import os
import sys
import importlib
import pandas as pd

D, OUTCSV = sys.argv[1], sys.argv[2]
sys.path.insert(0, D)
eng = importlib.import_module("02_rot_engine")
mc = importlib.import_module("04_mc_collisions")
C = mc.C

tot = mc.make_tot(eng.RotModel(), None)
mdf, nn = C.load_full_margin(mc.MARGIN_CSV)
rows = []
for ridge, (seed, pid) in mc.RIDGES.items():
    P = mc.track_points(tot, seed, pid, 902003)
    for t in mc.AGES_T:
        la, lo = mc.project(tot, 902003, t, P)
        r = C.collide_strict(C.margin_segments(mdf, nn, t), mc.segs_for_age(la, lo, t))
        rows.append(dict(ridge=ridge, age=t, lat=None if r is None else r[0]))
pd.DataFrame(rows).to_csv(OUTCSV, index=False)
print(f"written {OUTCSV}")
