# -*- coding: utf-8 -*-
"""Assignment of deposits and samples to the nearest track with the four chains (Easter,
Juan Fernandez, San Felix and Foundation / Iquique Ridge), with the criterion of the thesis
(Ch. 4, "distance to the nearest reconstructed hotspot track"): at the age of each deposit
or sample the tracks are reconstructed with South America fixed, clipped at their intersection
with the margin at that age (the subducted segment is kept, as in the thesis), and the minimum
distance over all tracks is taken. HS is the plate of the nearest track and N_S = 1 if the point
lies north of its nearest track point.

Differences from the thesis:
- Brute-force distance over the densified track vertices (step <= 1 km), in place of the
  projection onto each segment.
- Nominal tracks of the chosen model (current _hsb .rot, f_PAC 0.44, plate 902003),
  with the anchors of 03b; Foundation only on Nazca/Farallon (tau >= t_switch, Iquique
  Ridge). The margin is that of the thesis (Margin_SAM_0.1.csv, row with the nearest age).
- lat_inter = latitude of the nominal collision of the assigned track at that age (NaN if that
  track does not cross the margin).

Inputs: the databases read by step 07 (restored positions Latitude_ROT/Longitude_ROT and
Latitude_restored/Longitude_restored). They are not modified.
Outputs (outputs/asignacion_4tracks/): porphyries_4tracks.csv, epithermal_4tracks.csv,
samples_4tracks.csv, with the database columns plus HS4, Dist4_km, N_S4, lat_inter4,
Lat_track4, Lon_track4; and a summary on screen (HS changes relative to the thesis database).

Usage:  py assign_nearest_track.py [model_folder]   (default pipeline = WK08)
Outputs go to outputs/asignacion_4tracks/<model_folder>/.
"""
import os
import sys
import json
import importlib.util
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
PIPE = os.path.join(ROOT, "models", "Q22_WK08-A")      # modules (identical in the eight folders)
# Development version; the pipeline step is models/Q22_WK08-A/07_deposit_tracks.py.
sys.path.insert(0, HERE)
import track_vs_ridge_axis as r219       # noqa: E402  (folder -> trusted .rot)
MOD = sys.argv[1] if len(sys.argv) > 1 else "Q22_WK08-A"
ROT = os.path.join(ROOT, "models", MOD, "rot", r219.MODELOS[MOD][0])
FIT = os.path.join(ROOT, "models", MOD, "hotspot_fit.json")
# Collision -> deposit lag. With PAPER_DEP_LAG = L (Myr) the deposit of age t (at its
# position restored to t) is compared with the tracks and the margin at t + L.
LAG = float(os.environ.get("PAPER_DEP_LAG", "0"))
OUT = os.path.join(HERE, "outputs", "asignacion_4tracks", MOD + (f"_lag{LAG:g}" if LAG else ""))
R_KM = 6371.0
STEP_KM = 1.0
REL = 902003


def load_mod(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


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
    (as in the thesis: the unsubducted segment does not affect the continental lithosphere)."""
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
    sys.path.insert(0, PIPE)
    os.environ["PAPER_ROT_FILE"] = ROT
    os.environ["PAPER_HS_FIT"] = FIT
    load_mod(os.path.join(PIPE, "02_rot_engine.py"), "02_rot_engine")
    mc = load_mod(os.path.join(PIPE, "04_mc_collisions.py"), "04_mc_collisions")
    C = mc.C
    M = mc.eng.RotModel(ROT)
    tot = mc.make_tot(M, None)
    with open(FIT, encoding="utf-8") as f:
        fit = json.load(f)["tracks"]
    assert "Foundation" in fit, f"the hotspot_fit.json of {MOD} does not include Foundation"
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

    os.makedirs(OUT, exist_ok=True)
    dep = load_mod(os.path.join(PIPE, "07_impact_geochem_ESI.py"), "07_impact_geochem_ESI")
    por = pd.read_csv(os.path.join(dep.DEP, "Porphyry_Plate_ID_geodyn_v5.csv"))
    epi = pd.read_csv(os.path.join(dep.DEP, "Epithermal_Plate_Id_geodyn_v4_thesis_filtered.csv"))
    g = pd.read_excel(dep.GEOCHEM)
    jobs = [("porphyries", por, "Latitude_ROT", "Longitude_ROT", "Age MY"),
            ("epithermal", epi, "Latitude_ROT", "Longitude_ROT", "Age MY"),
            ("samples", g, "Latitude_restored", "Longitude_restored", "Age (Ma)")]
    names = {13: "JFR", 14: "Easter", 15: "TT", 16: "Foundation"}
    for nom, df, lac, loc, agc in jobs:
        out = assign(df, lac, loc, agc)
        out.to_csv(os.path.join(OUT, f"{nom}_4tracks.csv"), index=False)
        hs_db = pd.to_numeric(out["HS"], errors="coerce")
        k = np.isfinite(out.HS4)
        print(f"\n== {nom}: {int(k.sum())} with age 0-66 Ma and restored position (of {len(out)})")
        print("   HS4: " + ", ".join(f"{names[int(p)]} {int((out.HS4 == p).sum())}" for p in (14, 13, 15, 16)))
        same = k & (hs_db == out.HS4)
        print(f"   HS equal to the thesis database: {int(same.sum())}; move to Foundation: "
              f"{int((k & (out.HS4 == 16)).sum())}; change among the three thesis chains: "
              f"{int((k & (out.HS4 != 16) & (hs_db != out.HS4)).sum())}")
        dd = (out.Dist4_km - pd.to_numeric(out["Dist_Sample_Track"], errors="coerce"))[same]
        print(f"   with the same HS, Dist4 - Dist_thesis: median {dd.median():.1f} km, "
              f"|.| p90 {dd.abs().quantile(0.9):.1f} km")
        if nom != "samples":
            f = out[k & (out.HS4 == 16)]
            name_col = "Name Depos" if "Name Depos" in f.columns else f.columns[1]
            for _, r in f.iterrows():
                print(f"     Foundation: {str(r[name_col])[:28]:28s} {float(r[agc]):6.2f} Ma  "
                      f"lat_ROT {float(r[lac]):7.2f}  d {r.Dist4_km:6.0f} km  (HS thesis {r.HS})")
    print(f"\nwritten to {OUT}")


if __name__ == "__main__":
    main()
