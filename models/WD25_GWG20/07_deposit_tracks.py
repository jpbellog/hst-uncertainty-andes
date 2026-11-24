# -*- coding: utf-8 -*-
"""STEP 07 -- Assignment of deposits and samples to the nearest hotspot track.

Criterion: distance to the nearest reconstructed hotspot track, with the four chains of the
03b fit (Easter, Juan Fernandez, San Felix and Foundation/Iquique): at the age t of each
deposit or sample, at its restored position (Latitude_ROT / Latitude_restored), the tracks
are reconstructed with South America fixed (f_PAC 0.44, plate 902003), clipped at their
intersection with the margin of that age (Margin_SAM_0.1.csv), and the subducted segment
is kept. The distance is computed by brute force over the densified track vertices
(step <= 1 km). HS4 = plate of the nearest track, Dist4_km, N_S4 = 1 if the point lies N of
its nearest track point, lat_inter4 = nominal collision of the assigned track. Foundation
only over Nazca/Farallon (tau >= t_switch of 03b).
Validation: with the tracks used to build deposit database v7, the method reproduces the
database assignment (71/71 same HS, max 0.01 km).

Collision-to-deposit lag: with PAPER_DEP_LAG = L (Myr) a deposit of age t is compared with
the tracks and the margin at t + L; the outputs carry the suffix _lag<L>.

Environment: PAPER_ROT_FILE (_hsb .rot of the model) and PAPER_HS_FIT (hotspot_fit.json), as in
steps 04-11. The test across the eight models is step 12 (12_deposit_track_test.py).
Usage:  py 07_deposit_tracks.py
Outputs: deposit_tracks_porphyries.csv, deposit_tracks_epithermal.csv, deposit_tracks_samples.csv
"""
import os
import sys
import json
import importlib
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
eng = importlib.import_module("02_rot_engine")
mc = importlib.import_module("04_mc_collisions")

R_KM = 6371.0
STEP_KM = 1.0
REL = 902003
LAG = float(os.environ.get("PAPER_DEP_LAG", "0"))
SUF = f"_lag{LAG:g}" if LAG else ""
DEP = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "deposits")  # copy of the deposit database (see rotation_models/README.md)
GEOCHEM = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "geochemistry",
          "Full_Data_v11_tt_with_Clusters_v12_geodyn2.xlsx")
# age-quality filter (version 3): age methods that are excluded
EXCLUDED = {m.strip().lower() for m in [
    "K/Ar plagioclase", "K/Ar amphibole", "inference from K/Ar plagioclase",
    "Inference", "inference from K/Ar whole rock", "K/Ar whole rock",
    "Inference from K/Ar (biotite plagioclase whole rock) data",
    "Inference from K/Ar (biotite amphibole) data",
    "Inference from K/Ar biotite data", "K/Ar biotite", "K/Ar feldespar",
    "K/Ar glass", "K/Ar", "N.D.", "Inference from K/Ar",
    "The age of sample S-16 is assigned", "?",
    "The age of sample S-59 is assigned",
    "The age of sample S-H22 (176m) is assigned",
    "The age of sample S-T2 is assigned",
    "The age of sample S-50 is assigned",
    "The age of sample S-31 is assigned",
    "K-Ar feldspar", "K-Ar", "inference with dated units", "K-Ar biotite",
    "K-Ar whole-rock", "K-Ar amphibole", "K-Ar hornblende", "K-Ar plagioclase",
    "radiometric", "K-Ar whole rock",
    "inference from dated samples of same unit with K-Ar"]}


def units(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def ll(v):
    return np.degrees(np.arcsin(np.clip(v[..., 2], -1, 1))), np.degrees(np.arctan2(v[..., 1], v[..., 0]))


def densify(la, lo):
    """Vertices (n, 3) every <= STEP_KM along the polyline."""
    V = units(la, lo)
    out = [V[:1]]
    for a, b in zip(V[:-1], V[1:]):
        ang = np.arccos(np.clip(a @ b, -1, 1))
        n = max(1, int(np.ceil(ang * R_KM / STEP_KM)))
        if ang < 1e-12:
            continue
        f = np.linspace(0, 1, n + 1)[1:]
        s = np.sin(ang)
        out.append((np.sin((1 - f) * ang)[:, None] * a + np.sin(f * ang)[:, None] * b) / s)
    return np.vstack(out)


def clip_at(la, lo, c):
    """Subducted segment: from the collision point c = (lat, lon) to the oldest end
    (the non-subducted segment does not affect the continental lithosphere)."""
    if c is None:
        return la, lo
    V = units(la, lo)
    pc = units(np.array([c[0]]), np.array([c[1]]))[0]
    best, bj = np.inf, None
    for j in range(len(la) - 1):
        L = np.arccos(np.clip(V[j] @ V[j + 1], -1, 1))
        d = np.arccos(np.clip(V[j] @ pc, -1, 1)) + np.arccos(np.clip(pc @ V[j + 1], -1, 1)) - L
        if d < best:
            best, bj = d, j
    return np.r_[c[0], la[bj + 1:]], np.r_[c[1], lo[bj + 1:]]


def main():
    assert mc.HS_FIT, "run with PAPER_HS_FIT (and PAPER_ROT_FILE) from the 03b fit"
    C = mc.C
    tot = mc.make_tot(eng.RotModel(), None)
    with open(mc.HS_FIT, encoding="utf-8") as f:
        fit = json.load(f)["tracks"]
    assert "Foundation" in fit, f"{mc.HS_FIT} does not contain Foundation"
    P = {r: (pid, mc.track_points(tot, seed, pid, REL)) for r, (seed, pid) in mc.RIDGES.items()}
    mdf, nn = C.load_full_margin(mc.MARGIN_CSV)
    cache = {}

    def tracks_at(t):
        t = round(float(t), 2)
        if t in cache:
            return cache[t]
        ms = C.margin_segments(mdf, nn, t)
        res = []
        for r, (pid, Pr) in P.items():
            la, lo = mc.project(tot, REL, t, Pr)
            k = (mc.TAUS >= t - 1e-9) & np.isfinite(la)
            la, lo, ta = la[k], lo[k], mc.TAUS[k]
            if len(la) < 2:
                continue
            segs = np.array([(ta[j], la[j], lo[j], ta[j + 1], la[j + 1], lo[j + 1])
                             for j in range(len(ta) - 1)])
            c = C.collide_strict(ms, segs)
            cl, co = clip_at(la, lo, c)
            res.append((r, pid, densify(cl, co), np.nan if c is None else c[0]))
        cache[t] = res
        return res

    def assign(df, lat_col, lon_col, age_col):
        df = df.copy()
        a = pd.to_numeric(df[age_col], errors="coerce")
        la = pd.to_numeric(df[lat_col], errors="coerce")
        lo = pd.to_numeric(df[lon_col], errors="coerce")
        cols = {k: np.full(len(df), np.nan) for k in
                ("HS4", "Dist4_km", "N_S4", "lat_inter4", "Lat_track4", "Lon_track4")}
        ok = np.flatnonzero(np.isfinite(a) & np.isfinite(la) & np.isfinite(lo) & (a >= 0) & (a <= 66))
        for t in np.unique(a.values[ok].round(2)):
            idx = ok[np.isclose(a.values[ok].round(2), t)]
            X = units(la.values[idx], lo.values[idx])
            best = np.full(len(idx), np.inf)
            for r, pid, V, lat_c in tracks_at(t + LAG):
                dots = X @ V.T
                j = np.argmax(dots, axis=1)
                d = np.arccos(np.clip(dots[np.arange(len(idx)), j], -1, 1)) * R_KM
                better = d < best
                if not better.any():
                    continue
                best[better] = d[better]
                tl, tn = ll(V[j[better]])
                ii = idx[better]
                cols["HS4"][ii] = pid
                cols["Dist4_km"][ii] = d[better]
                cols["Lat_track4"][ii], cols["Lon_track4"][ii] = tl, tn
                cols["N_S4"][ii] = (la.values[ii] > tl).astype(float)
                cols["lat_inter4"][ii] = lat_c
        for k, v in cols.items():
            df[k] = v
        return df

    por = pd.read_csv(os.path.join(DEP, "Porphyry_Plate_ID_geodyn_v7.csv"))
    epi = pd.read_csv(os.path.join(DEP, "Epithermal_Plate_Id_geodyn_v6_thesis.csv"))
    g = pd.read_excel(GEOCHEM)
    jobs = [("porphyries", por, "Latitude_ROT", "Longitude_ROT", "Age MY"),
            ("epithermal", epi, "Latitude_ROT", "Longitude_ROT", "Age MY"),
            ("samples", g, "Latitude_restored", "Longitude_restored", "Age (Ma)")]
    names = {13: "JFR", 14: "Easter", 15: "TT", 16: "Foundation"}
    for nom, df, lac, loc, agc in jobs:
        out = assign(df, lac, loc, agc)
        out.to_csv(os.path.join(HERE, f"deposit_tracks_{nom}{SUF}.csv"), index=False)
        hs_db = pd.to_numeric(out["HS"], errors="coerce")
        k = np.isfinite(out.HS4)
        print(f"\n== {nom}: {int(k.sum())} with age 0-66 Ma and restored position (of {len(out)})")
        print("   HS4: " + ", ".join(f"{names[int(p)]} {int((out.HS4 == p).sum())}" for p in (14, 13, 15, 16)))
        same = k & (hs_db == out.HS4)
        print(f"   HS equal to the database HS: {int(same.sum())}; moved to Foundation: "
              f"{int((k & (out.HS4 == 16)).sum())}; changed among the three database chains: "
              f"{int((k & (out.HS4 != 16) & (hs_db != out.HS4)).sum())}")
        dd = (out.Dist4_km - pd.to_numeric(out["Dist_Sample_Track"], errors="coerce"))[same]
        print(f"   with the same HS, Dist4 - Dist_db: median {dd.median():.1f} km, "
              f"|.| p90 {dd.abs().quantile(0.9):.1f} km")
        if nom != "samples":
            f = out[k & (out.HS4 == 16)]
            name_col = "Name Depos" if "Name Depos" in f.columns else f.columns[1]
            for _, r in f.iterrows():
                print(f"     Foundation: {str(r[name_col])[:28]:28s} {float(r[agc]):6.2f} Ma  "
                      f"lat_ROT {float(r[lac]):7.2f}  d {r.Dist4_km:6.0f} km  (HS db {r.HS})")
    print(f"\nwrote deposit_tracks_*{SUF}.csv in {HERE}")


if __name__ == "__main__":
    main()
