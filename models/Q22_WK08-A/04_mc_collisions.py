# -*- coding: utf-8 -*-
"""STEP 04 -- Monte Carlo of the track-trench collision points (error bands).

Monte Carlo simulation (N=1000): each realisation samples
  (a) the nodal rotations of the plate circuit with their covariances from
      Global_Model_R_BG_2025.rot (911 NZ-PAC and 902 FAR-PAC from Wilder 2003 /
      Rowan & Rowley 2014 via Quiero et al. 2022; 901 PAC-WANT from Croon 2008
      filtered with REDBACK; 804/802 Granot & Dyment 2018; 201 DeMets &
      Merkouriev 2019),
  (b) the rotations of the WK08-A hotspot frame (plate 002) with the
      covariances of Table 3 of Wessel & Kroenke (2008),
  (c) the position of the hotspots (plates 13/14/15, fitted to the published
      geochronology) with the along/across sigmas of step 03
      (eq. 2 of O'Neill et al. 2005; Ray 2012, Lara 2018, Haase 2000, C&L 2015).
With the perturbed model the track is reconstructed and the track-trench
intersection is recomputed (collide_strict, collision_lines pipeline) per age
and f_PAC model.

Outputs: mc_collision_bands_{variant}.csv with percentiles of collision lat/lon
per (ridge, pf, age); variant: total | circuit | wk08 | geochron.

Usage:  py 04_mc_collisions.py [N] [comma-separated variants]
"""
import os
import sys
import time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# collision_primitives.py is in this same folder
import importlib
eng = importlib.import_module("02_rot_engine")
import collision_primitives as C

MARGIN_CSV = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "margin", "Margin_SAM_0.1.csv")  # copy of Margin_SAM (see rotation_models/README.md)
MODELS = {"0.5": 902001, "0.36": 902002, "0.44": 902003}
RIDGES = {"Easter": ((-26.3138, -106.599), 14),
          "TT":     ((-26.28791, -81.42893), 15),
          "JFR":    ((-33.750097, -81.905212), 13)}
CIRCUIT_PLATES = [911, 902, 901, 804, 802, 201, 90201, 90202, 902001]

# With PAPER_CIRCUIT_COV=b (or a) the covariances of the links that the .rot
# carries without them are loaded (circuit_covariance.aplicar; variant of node
# 802@71.1), Somalia, Lord Howe and Australia are also sampled, and the
# perturbations are indexed by (plate, age, fixed) so that the two rows of a
# change of fixed plate do not share the draw. Without the variable none of
# this is applied.
CIRCUIT_COV = os.environ.get("PAPER_CIRCUIT_COV", "").strip().lower() or None
if CIRCUIT_COV not in (None, "a", "b"):
    raise ValueError(f"PAPER_CIRCUIT_COV={CIRCUIT_COV!r}: use a, b or empty")
COV_CIRCUITO_DIR = (__import__("os").environ["PAPER_ROOT"] + "/analysis")
if CIRCUIT_COV:
    CIRCUIT_PLATES = CIRCUIT_PLATES + [709, 833, 801]
# Rows that repeat the same rotation (same pole and cov) of another plate: they
# are the same measurement and share the draw. 902 21.767 901 is the copy of
# NZ-PAC 6Bn.1 that opens FAR-PAC at the break-up; with independent draws the
# propagated cov of the NZ-FAR identity rows (911 21.767 902, 902 21.767 911)
# would no longer be zero (see the comment in the .rot). Only with
# PAPER_CIRCUIT_COV.
# 804 45.03 802 exists only in the Granot 2013 sensitivity .rot (constant 18o).
# The .rot files of the eight models carry instead the node of
# Cande et al. (2000) at 43.462 (GTS2020); the pair 44.846 is unused.
# Same pairs in GTS2020 (rot_gts2020.py): 21.767 -> 21.806,
# 45.03 -> 44.846 and 40.13 -> 40.073 Ma.
SAME_DRAW = {(902, 21.767, 901): (911, 21.767, 901),
             (804, 45.03, 802): (804, 40.13, 802),
             (902, 21.806, 901): (911, 21.806, 901),
             (804, 44.846, 802): (804, 40.073, 802)}
AGES_T = np.arange(0, 67, 1.0)
TAUS = np.concatenate([np.arange(0, 84, 1.0), [93.0, 103.0, 113.0, 120.0]])

WK08_CONV = os.environ.get("PAPER_WK08_CONV", "right")

# Hotspot fit (03b_hotspot_fit.py): with PAPER_HS_FIT=<hotspot_fit*.json> the
# hotspot position is sampled from its parametric fit (variant "hotspot",
# which replaces "geochron" also within "total"). Requires running with
# PAPER_ROT_FILE pointing to the _hsb .rot written by the same fit.
HS_FIT = os.environ.get("PAPER_HS_FIT")
hsfit = importlib.import_module("03b_hotspot_fit") if HS_FIT else None
PLATE_RIDGE = {14: "Easter", 13: "JFR", 15: "TT"}


def hs_tag():
    """Output suffix with the hotspot fit (steps 04, 07, 08, 10, 11):
    "_hsb" with w as confidence and "_hsb_pred" with w as prediction; empty
    without PAPER_HS_FIT. With PAPER_CIRCUIT_COV "_cc" (variant b) or "_cca"
    (variant a) is added, and at the end PAPER_TAG_EXTRA (e.g. "_g13", nodes of
    Granot et al. 2013 in the .rot)."""
    cc = {None: "", "b": "_cc", "a": "_cca"}[CIRCUIT_COV] + os.environ.get("PAPER_TAG_EXTRA", "")
    if not HS_FIT:
        return cc
    import json
    with open(HS_FIT, encoding="utf-8") as f:
        w_mode = json.load(f)["meta"].get("w_mode", "conf")
    return ("_hsb" if w_mode == "conf" else f"_hsb_{w_mode}") + cc


HS_TAG = hs_tag()

# Fourth chain: the chains of the fit that are not in RIDGES (Foundation,
# plate 16) are added from the json. With a plate change (t_switch in the json)
# the track exists on the relative plate of the model (Nazca/Farallon, Iquique
# ridge) only for tau >= t_switch; the younger part lies on the Pacific plate
# (Foundation chain), far from the margin, and is not computed.
SWITCH = {}


def add_fit_ridges(path):
    import json
    with open(path, encoding="utf-8") as f:
        tracks = json.load(f)["tracks"]
    for ridge, r in tracks.items():
        pid = r["plate"]
        if ridge not in RIDGES:
            RIDGES[ridge] = (tuple(hsfit.SEEDS[pid]), pid)
            PLATE_RIDGE[pid] = ridge
        if "t_switch" in r:
            SWITCH[pid] = float(r["t_switch"])


if HS_FIT:
    add_fit_ridges(HS_FIT)


def load_circuit_cov(M):
    """Added circuit covariances (see CIRCUIT_COV). Call after
    load_wk08 / clone_9020x_cov; does nothing without PAPER_CIRCUIT_COV."""
    if not CIRCUIT_COV:
        return
    if COV_CIRCUITO_DIR not in sys.path:
        sys.path.append(COV_CIRCUITO_DIR)
    import circuit_covariance
    circuit_covariance.aplicar(M, hellinger=CIRCUIT_COV)


def load_hotspot_fit(M, path):
    """Reads the fit and checks that nodes 13/14/15 of the loaded .rot are
    those of the fit (otherwise the nominal and the realisations would not match)."""
    import json
    with open(path, encoding="utf-8") as f:
        js = json.load(f)
    fit = js["tracks"]
    prior = js["meta"].get("prior", "none")
    if prior != hsfit.PRIOR:
        raise RuntimeError(f"{path} was fitted with PAPER_HS_PRIOR={prior}; "
                           f"run step 04 with the same value")
    for ridge, r in fit.items():
        th = np.array(r["theta"])
        for (age, la, lo, ang, fx, Cv) in M.table[r["plate"]]:
            A = eng.quat_to_matrix(eng.quat_from_pole(la, lo, ang))
            B = eng.quat_to_matrix(eng.quat_from_pole(
                *hsfit.node_pole(r["plate"], th, age, r["t_max"])))
            if fx != 2 or np.abs(A - B).max() > 1e-6:
                raise RuntimeError(f"plate {r['plate']} of the .rot does not match "
                                   f"{path}: run with PAPER_ROT_FILE=..._hsb.rot")
        r["L"] = np.linalg.cholesky(np.array(r["cov"]))
    return fit


def rotvec(A):
    """Rotation vector h with exp(M(h)) = A."""
    q = eng.matrix_to_quat(A)
    v = q[1:]
    s = np.linalg.norm(v)
    if s < 1e-15:
        return np.zeros(3)
    return v / s * 2.0 * np.arctan2(s, q[0])


def sample_hotspot(M, rng, fit, pert):
    """One realisation of plates 13/14/15: theta = theta_hat + L z +
    sum_k g_k e_k, with e_k common to the three tracks (calibration)."""
    eps = {k: s * rng.standard_normal()
           for k, s in hsfit.SYSTEMATICS.items()}
    for ridge, r in fit.items():
        th = np.array(r["theta"]) + r["L"] @ rng.standard_normal(len(r["theta"]))
        for k, g in r["g"].items():
            th = th + np.array(g) * eps[k]
        pid = r["plate"]
        for (age, la, lo, ang, fx, Cv) in M.table[pid]:
            A0 = eng.quat_to_matrix(eng.quat_from_pole(la, lo, ang))
            A1 = eng.quat_to_matrix(eng.quat_from_pole(
                *hsfit.node_pole(pid, th, age, r["t_max"])))
            pert[(pid, age)] = rotvec(A1 @ A0.T)

def load_wk08(M):
    """Attaches the WK08-A covariances (Table 3, [[a,b,d],[b,c,e],[d,e,f]]*1e-5 rad^2,
    PAC rel hotspots) to plate 002 of the .rot (which is the INVERSE: hotspots rel
    PAC). The engine perturbs on the left (exp(M(h)) R_node). Wessel and Kroenke
    (2008) do not state their convention:
      right (Chang 1988; Kirkwood et al. 1999): W_real = W exp(M(h)), and the
        inverse is exp(M(-h)) W^T  ->  C' = C;
      left: W_real = exp(M(h)) W, and the inverse is exp(M(-A h)) A  ->  C' = A C A^T
        (A = matrix of 002).
    PAPER_WK08_CONV=right|left, default right (as in W&D 2025, GMHRF and
    GWG20)."""
    if WK08_CONV not in ("right", "left"):
        raise ValueError(f"PAPER_WK08_CONV={WK08_CONV!r}: use right or left")
    wk = pd.read_csv(os.path.join(HERE, "data", "wk08a_covariances.csv"), comment="#")
    new = []
    for (age, la, lo, ang, fx, Cv) in M.table[2]:
        m = wk[(wk.age_ma - age).abs() < 0.02]
        if len(m):
            r = m.iloc[0]
            a, b, c, d, e, f = [float(r[k]) for k in "abcdef"]
            Cw = 1e-5 * np.array([[a, b, d], [b, c, e], [d, e, f]])
            A = eng.quat_to_matrix(eng.quat_from_pole(la, lo, ang))
            Cv = Cw if WK08_CONV == "right" else A @ Cw @ A.T
        new.append((age, la, lo, ang, fx, Cv))
    M.table[2] = new

def clone_9020x_cov(M):
    """90202 (Pf=0.36, 53-83 Ma) has no covariances; it is given those of the
    equivalent nodes of 90201 (same chrons of Rowan & Rowley 2014)."""
    cov01 = {age: Cv for (age, la, lo, ang, fx, Cv) in M.table[90201]
             if Cv is not None}
    new = []
    for (age, la, lo, ang, fx, Cv) in M.table[90202]:
        if Cv is None and age in cov01 and fx == 901:
            Cv = cov01[age]
        new.append((age, la, lo, ang, fx, Cv))
    M.table[90202] = new

def same_rotation(M, a, b):
    """True if rows (plate, age, fixed) a and b have the same pole and cov."""
    ra, rb = ([r for r in M.table.get(k[0], []) if r[0] == k[1] and r[4] == k[2]]
              for k in (a, b))
    if len(ra) != 1 or len(rb) != 1 or ra[0][5] is None or rb[0][5] is None:
        return False
    return (np.allclose(ra[0][1:4], rb[0][1:4], atol=1e-9)
            and np.allclose(ra[0][5], rb[0][5], rtol=1e-9, atol=0))

def sample_perturb(M, rng, variant, node_sig, pf, hs=None):
    """Builds the dict of nodal perturbations of one realisation."""
    pert = {}
    if hs is not None and variant in ("total", "hotspot"):
        sample_hotspot(M, rng, hs, pert)
    if variant in ("total", "circuit"):
        for p in CIRCUIT_PLATES:
            if p not in M.table:
                continue
            for (age, la, lo, ang, fx, Cv) in M.table[p]:
                if Cv is None:
                    continue
                # IDENTITY node with covariance: in the .rot only
                # 902001@83 (f_PAC=0.5) exists, which REPEATS the covariance of
                # the 83 Ma node of 90201 (same eigenvalues, verified).
                # Sampling it would double (~sqrt(2)) that term only for the
                # 0.5 model; it is skipped.
                if ang == 0.0 and la == 90.0:
                    continue
                try:
                    L = np.linalg.cholesky(Cv + 1e-18 * np.eye(3))
                except np.linalg.LinAlgError:
                    w, V = np.linalg.eigh(Cv)
                    L = V @ np.diag(np.sqrt(np.clip(w, 0, None)))
                key = (p, age, fx) if CIRCUIT_COV else (p, age)
                pert[key] = L @ rng.standard_normal(3)
        if CIRCUIT_COV:
            for dst, src in SAME_DRAW.items():
                if dst in pert and src in pert and same_rotation(M, dst, src):
                    pert[dst] = pert[src]
    if variant in ("total", "wk08"):
        for (age, la, lo, ang, fx, Cv) in M.table[2]:
            if Cv is None:
                continue
            L = np.linalg.cholesky(Cv + 1e-18 * np.eye(3))
            pert[(2, age)] = L @ rng.standard_normal(3)
    if variant in ("total", "geochron") and hs is None:
        for ridge, (seed, pid) in RIDGES.items():
            sig = node_sig[(node_sig.ridge == ridge) & (node_sig.pf == pf)]
            for (age, la, lo, ang, fx, Cv) in M.table[pid]:
                if age > 200:
                    continue
                s = sig[(sig.age_node - age).abs() < 0.01]
                if not len(s):
                    continue
                s_al = np.deg2rad(float(s.iloc[0].sigma_along_deg))
                s_ac = np.deg2rad(float(s.iloc[0].sigma_across_deg))
                # nominal absolute position of the hotspot at the node time
                q = M.total_rotation(pid, age, 0)
                p = eng.quat_to_matrix(q) @ eng.latlon_axis(*seed)
                # along-track direction: derivative of the hotspot position
                q2 = M.total_rotation(pid, min(age + 0.5, 120.0), 0)
                p2 = eng.quat_to_matrix(q2) @ eng.latlon_axis(*seed)
                t_hat = p2 - p - np.dot(p2 - p, p) * p
                nt = np.linalg.norm(t_hat)
                if nt < 1e-9:
                    # node without local drift: use local east as along
                    la0, lo0 = eng.xyz_to_latlon(p)
                    t_hat = np.array([-np.sin(np.radians(lo0)),
                                      np.cos(np.radians(lo0)), 0.0])
                else:
                    t_hat = t_hat / nt
                n_hat = np.cross(p, t_hat)
                dp = (s_al * rng.standard_normal() * t_hat
                      + s_ac * rng.standard_normal() * n_hat)
                pert[(pid, age)] = np.cross(p, dp)
    return pert

def make_tot(M, pert):
    """total_rotation memoised within one realisation (the plate tree is
    traversed only once per (plate, age))."""
    cache = {}
    def tot(plate, t):
        k = (plate, t)
        if k not in cache:
            cache[k] = M.total_rotation(plate, t, 0, pert)
        return cache[k]
    return tot

def track_points(tot, seed, pid, rel):
    """Track positions in the frame of the relative plate (p_plate per tau)."""
    P = np.empty((len(TAUS), 3))
    s = eng.latlon_axis(*seed)
    for i, tau in enumerate(TAUS):
        p_abs = eng.quat_to_matrix(tot(pid, tau)) @ s
        P[i] = eng.quat_to_matrix(tot(rel, tau)).T @ p_abs
    if pid in SWITCH:
        P[TAUS < SWITCH[pid]] = np.nan
    return P

def project(tot, rel, t, P):
    R = eng.quat_to_matrix(tot(201, t)).T @ eng.quat_to_matrix(tot(rel, t))
    V = P @ R.T
    lat = np.degrees(np.arcsin(np.clip(V[:, 2], -1, 1)))
    lon = np.degrees(np.arctan2(V[:, 1], V[:, 0]))
    return lat, lon

def segs_for_age(lat, lon, t):
    k = (TAUS >= t - 1e-9) & np.isfinite(lat)
    la, lo, ta = lat[k], lon[k], TAUS[k]
    segs = []
    for j in range(len(ta) - 1):
        segs.append((ta[j], la[j], lo[j], ta[j + 1], la[j + 1], lo[j + 1]))
    return np.array(segs)

def main(N=1000, variants=("total", "circuit", "wk08", "geochron"), pfs=None):
    """pfs: subset of f_PAC models to recompute (e.g. ("0.5",)). If given,
    the rows of those pf are MERGED into the existing CSV files (the others
    are kept)."""
    models = {k: v for k, v in MODELS.items() if pfs is None or k in pfs}
    M = eng.RotModel()
    # GMHRF and GWG20 frames (PAPER_FRAME_COV=rot): plate 002 already carries in
    # the .rot the covariances of its source; they are not overwritten by WK08-A.
    if os.environ.get("PAPER_FRAME_COV", "wk08") == "wk08":
        load_wk08(M)
    clone_9020x_cov(M)
    load_circuit_cov(M)
    hs = load_hotspot_fit(M, HS_FIT) if HS_FIT else None
    if hs is not None and "geochron" in variants:
        raise ValueError("with PAPER_HS_FIT the variant is 'hotspot', not 'geochron'")
    node_sig = None   # node sigmas of the variant without a hotspot fit; unused with PAPER_HS_FIT
    if os.path.exists(os.path.join(HERE, "node_sigmas.csv")):
        node_sig = pd.read_csv(os.path.join(HERE, "node_sigmas.csv"))
        node_sig["pf"] = node_sig["pf"].astype(str)
    mdf, nn = C.load_full_margin(MARGIN_CSV)
    msegs = {t: C.margin_segments(mdf, nn, t) for t in AGES_T}
    rng = np.random.default_rng(20260714)

    # nominal (no perturbation) + margin windows per ridge
    nominal = {}
    tot0 = make_tot(M, None)
    for pf, rel in models.items():
        for ridge, (seed, pid) in RIDGES.items():
            P = track_points(tot0, seed, pid, rel)
            for t in AGES_T:
                lat, lon = project(tot0, rel, t, P)
                r = C.collide_strict(msegs[t], segs_for_age(lat, lon, t))
                if r is not None:
                    nominal[(pf, ridge, t)] = r
    # margin window per (ridge, t): nominal lat +- 8 degrees (5-6x faster)
    win = {}
    for ridge in RIDGES:
        lats = [v[0] for (pf2, rg, t2), v in nominal.items() if rg == ridge]
        if not lats:   # no nominal collision (e.g. a track that does not reach the margin)
            lats = [-90.0, 10.0]
        lo_all, hi_all = min(lats) - 8.0, max(lats) + 8.0
        for t in AGES_T:
            ls = [v[0] for (pf2, rg, t2), v in nominal.items()
                  if rg == ridge and t2 == t]
            lo, hi = (min(ls) - 8.0, max(ls) + 8.0) if ls else (lo_all, hi_all)
            ms = msegs[t]
            k = (np.maximum(ms[:, 0], ms[:, 2]) >= lo) &                 (np.minimum(ms[:, 0], ms[:, 2]) <= hi)
            win[(ridge, t)] = ms[k]

    for variant in variants:
        t0 = time.time()
        acc = {}   # (pf, ridge, t) -> [lats], [lons]
        for pf, rel in models.items():
            for it in range(N):
                pert = sample_perturb(M, rng, variant, node_sig, pf, hs)
                tot = make_tot(M, pert)
                for ridge, (seed, pid) in RIDGES.items():
                    P = track_points(tot, seed, pid, rel)
                    for t in AGES_T:
                        lat, lon = project(tot, rel, t, P)
                        r = C.collide_strict(win[(ridge, t)],
                                             segs_for_age(lat, lon, t))
                        if r is None:
                            continue
                        acc.setdefault((pf, ridge, t), []).append(r)
            print(f"[{variant}] Pf={pf}: {N} realisations "
                  f"({time.time()-t0:.0f} s)", flush=True)
        rows = []
        for (pf, ridge, t), pts in sorted(acc.items()):
            A = np.array(pts)
            la, lo = A[:, 0], A[:, 1]
            nom = nominal.get((pf, ridge, t), (np.nan, np.nan))
            rows.append(dict(ridge=ridge, pf=pf, age=t, n=len(la),
                             lat_nominal=nom[0], lon_nominal=nom[1],
                             lat_p2_5=np.percentile(la, 2.5),
                             lat_p16=np.percentile(la, 16),
                             lat_p50=np.percentile(la, 50),
                             lat_p84=np.percentile(la, 84),
                             lat_p97_5=np.percentile(la, 97.5),
                             lat_sd=np.std(la),
                             lon_p2_5=np.percentile(lo, 2.5),
                             lon_p50=np.percentile(lo, 50),
                             lon_p97_5=np.percentile(lo, 97.5),
                             lon_sd=np.std(lo)))
        out = pd.DataFrame(rows)
        # PAPER_SUFFIX runs one pf per process, merged afterwards
        # (analysis/merge_mc_parts_tag.py), without write races.
        # with the hotspot fit the outputs carry "_hsb" and do not overwrite
        # the bands without it
        sfx = os.environ.get("PAPER_SUFFIX", HS_TAG)
        fn = os.path.join(HERE, f"mc_collision_bands_{variant}{sfx}.csv")
        if sfx:
            out.to_csv(fn, index=False)
            print(f"[{variant}] written {fn} ({len(out)} rows, "
                  f"{time.time()-t0:.0f} s)", flush=True)
            continue
        if pfs is not None and os.path.exists(fn):
            old = pd.read_csv(fn, dtype={"pf": str})
            old = old[~old.pf.isin(pfs)]
            out = pd.concat([old, out], ignore_index=True)
            out = out.sort_values(["ridge", "pf", "age"])
        out.to_csv(fn, index=False)
        print(f"[{variant}] written {fn} ({len(out)} rows, "
              f"{time.time()-t0:.0f} s)", flush=True)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    if len(sys.argv) > 2:
        vs = tuple(sys.argv[2].split(","))
    else:
        vs = ("total", "circuit", "wk08", "hotspot" if HS_FIT else "geochron")
    pfs = tuple(sys.argv[3].split(",")) if len(sys.argv) > 3 else None
    main(N, vs, pfs)
