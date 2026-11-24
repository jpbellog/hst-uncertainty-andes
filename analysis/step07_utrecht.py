# -*- coding: utf-8 -*-
"""Step 07 of a model with the restored positions and the trench of the Utrecht block
model (andean_shortening_utrecht.py). Uses the 07_deposit_tracks.py of the model folder
without modifying it: it only changes DEP, GEOCHEM, MARGIN_CSV and the output folder.

Usage:  PAPER_ROT_FILE=... PAPER_HS_FIT=... py step07_utrecht.py <model_folder> [own_trench]
  own_trench: keeps our trench (Margin_SAM_0.1.csv) and changes only the positions; the trench
  restored with the Utrecht blocks has steps of up to ~130 km between blocks (e.g. at
  ~24.3S, between 2527 and 2516), which shift the nearby collisions.
Output: outputs/andean_shortening/<folder>[_own_trench]/deposit_tracks_*.csv
"""
import os
import sys
import shutil
import importlib

HERE = os.path.dirname(os.path.abspath(__file__))
U = os.path.join(HERE, "outputs", "andean_shortening")
D = sys.argv[1]
sys.path.insert(0, D)
s07 = importlib.import_module("07_deposit_tracks")

dep = os.path.join(U, "dep")
os.makedirs(dep, exist_ok=True)
for a, b in (("porphyries_utrecht.csv", "Porphyry_Plate_ID_geodyn_v7.csv"),
             ("epithermal_utrecht.csv", "Epithermal_Plate_Id_geodyn_v6_thesis.csv")):
    shutil.copy2(os.path.join(U, a), os.path.join(dep, b))
s07.DEP = dep
s07.GEOCHEM = os.path.join(U, "samples_utrecht.xlsx")
FOSA_PROPIA = len(sys.argv) > 2 and sys.argv[2] == "own_trench"
if not FOSA_PROPIA:
    s07.mc.MARGIN_CSV = os.path.join(U, "Margin_SAM_0.1_utrecht.csv")
s07.HERE = os.path.join(U, os.path.basename(os.path.normpath(D)) + ("_own_trench" if FOSA_PROPIA else ""))
os.makedirs(s07.HERE, exist_ok=True)
s07.main()
