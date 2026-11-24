# -*- coding: utf-8 -*-
"""STEP 02 -- Rotation engine for the plate model + covariances from the .rot file.

Reads Global_Model_R_BG_2025.rot (the base rotation model), including the
covariance matrices (Cxx,Cxy,Cxz,Cyy,Cyz,Czz in rad^2) stored at the end of the
comment of each line (taken from Quiero et al. 2022 / Wilder 2003 /
Rowan & Rowley 2014 / Croon 2008-REDBACK / DeMets & Merkouriev 2019 /
Granot & Dyment 2018).

Provides:
  - RotModel: .rot parsing, interpolated finite rotation per plate
    (pole/angle interpolation as in GPlates: slerp of the stage quaternion),
    total plate->anchor rotation by recursion over the plate tree, and Monte
    Carlo sampling of the nodal rotations with their covariances.
  - motion_path(seed, plate_seed, plate_rel, anchor, tau_list, t): reproduces the
    GPlates motion path (validated against the exported .xy files).

VALIDATION (run this file): reproduces motion_path_output_{t}Ma.xy for the
three tracks / three f_PAC models and prints the maximum error to the console.

Covariance convention (Chang/Kirkwood/Hellinger, as in Quiero et al. 2022):
R_real = exp(M(h)) * R_nominal, with h ~ N(0, C) in radians (frame fixed to the
reference plate of that line). When rotations are composed the covariances are
transported: C_total = C_1 + A_1 C_2 A_1^T (Kirkwood et al. 1999, eq. 19).
"""
import os
import re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# The environment variable PAPER_ROT_FILE runs the whole pipeline against an
# alternative version of the model (e.g. the .rot with the Nazca-Pacific branch
# of Wilson and DeMets 2025, see build_rot_wd25.py) without other changes.
ROT_FILE = os.environ.get(
    "PAPER_ROT_FILE",
    os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "rotations",
                 "Global_Model_R_BG_2025.rot"))  # base rotation model
D2R = np.pi / 180.0
# PAPER_INTERP_COV=norm preserves the variance between nodes. By default
# (empty) each node is sampled, perturbed and interpolated with slerp: the
# effective perturbation between nodes is ~ (1-f) h0 + f h1, with variance
# (1-f)^2 C0 + f^2 C1 (half at the midpoint if C0 ~ C1). With "norm", if both
# nodes of the interval are perturbed, the nominal rotation is interpolated and
# h = ((1-f) h0 + f h1) / sqrt((1-f)^2 + f^2) is applied: equal to C at the nodes,
# continuous in t, and with variance ~ C between nodes (see analysis/
# interp_node_covariances.py). h0 and h1 are in the frame of the same fixed plate.
INTERP_COV = os.environ.get("PAPER_INTERP_COV", "").strip().lower()
if INTERP_COV not in ("", "norm"):
    raise ValueError(f"PAPER_INTERP_COV={INTERP_COV!r}: use norm or empty")

# ------------------------------------------------------------------ algebra
def latlon_axis(lat, lon):
    la, lo = lat * D2R, lon * D2R
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])

def quat_from_pole(lat, lon, angle_deg):
    ax = latlon_axis(lat, lon)
    a = 0.5 * angle_deg * D2R
    return np.concatenate([[np.cos(a)], np.sin(a) * ax])

def quat_mult(q, r):
    w1, x1, y1, z1 = q; w2, x2, y2, z2 = r
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2])

def quat_conj(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])

def quat_slerp(q0, q1, f):
    d = float(np.dot(q0, q1))
    if d < 0:
        q1, d = -q1, -d
    d = min(1.0, max(-1.0, d))
    th = np.arccos(d)
    if th < 1e-12:
        out = (1 - f) * q0 + f * q1
        return out / np.linalg.norm(out)
    return (np.sin((1 - f) * th) * q0 + np.sin(f * th) * q1) / np.sin(th)

def quat_to_matrix(q):
    w, x, y, z = q
    return np.array([
        [1 - 2*(y*y + z*z), 2*(x*y - w*z),     2*(x*z + w*y)],
        [2*(x*y + w*z),     1 - 2*(x*x + z*z), 2*(y*z - w*x)],
        [2*(x*z - w*y),     2*(y*z + w*x),     1 - 2*(x*x + y*y)]])

def exp_M(h):
    """Small rotation exp(M(h)) for h in radians (Rodrigues)."""
    th = np.linalg.norm(h)
    if th < 1e-15:
        return np.eye(3)
    k = h / th
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)

def xyz_to_latlon(v):
    return (np.degrees(np.arcsin(np.clip(v[2], -1, 1))),
            np.degrees(np.arctan2(v[1], v[0])))

# ---------------------------------------------------------------- .rot parsing
_NUM = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"

def parse_rot(fn=ROT_FILE):
    """Returns {moving: [(age, lat, lon, angle, fixed, C(3x3)|None), ...]}."""
    table = {}
    with open(fn, errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("!"):
                continue
            m = re.match(rf"^\s*(\d+)\s+({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})\s+(\d+)(.*)$",
                         line)
            if not m:
                continue
            mov = int(m.group(1)); age = float(m.group(2))
            lat = float(m.group(3)); lon = float(m.group(4))
            ang = float(m.group(5)); fix = int(m.group(6))
            rest = m.group(7)
            C = None
            nums = re.findall(_NUM, rest)
            if len(nums) >= 6:
                c = [float(x) for x in nums[-6:]]
                # discards false covariances (comments with 6+ numbers that are
                # not matrices: requires the small magnitudes typical of rad^2)
                if all(abs(v) < 1e-2 for v in c):
                    C = np.array([[c[0], c[1], c[2]],
                                  [c[1], c[3], c[4]],
                                  [c[2], c[4], c[5]]])
            table.setdefault(mov, []).append((age, lat, lon, ang, fix, C))
    for mov in table:
        table[mov].sort(key=lambda r: r[0])
    return table

class RotModel:
    def __init__(self, fn=ROT_FILE):
        self.table = parse_rot(fn)

    def _segments(self, plate):
        """Splits the sequence of a plate into intervals with a constant fixed plate."""
        rows = self.table[plate]
        segs = []
        for r in rows:
            if segs and segs[-1][0] == r[4]:
                segs[-1][1].append(r)
            else:
                segs.append((r[4], [r]))
        return segs

    def finite_rotation(self, plate, t, perturb=None):
        """Quaternion of the finite rotation of `plate` relative to its fixed
        plate at age t, and the fixed plate used. Interpolates with slerp (as
        GPlates). `perturb`: dict {(plate, age_nodo): h(3,)} of MC perturbations
        applied to the NODES (exp(M(h)) * R_node) before interpolation. The key
        (plate, age_nodo, fixed) takes precedence: it distinguishes the two rows
        of a change of fixed plate at the same age."""
        rows = self.table[plate]
        ages = [r[0] for r in rows]
        if t <= ages[0]:
            i0 = i1 = 0
        elif t >= ages[-1]:
            i0 = i1 = len(rows) - 1
        else:
            i1 = next(k for k, a in enumerate(ages) if a >= t)
            i0 = i1 - 1 if ages[i1] > t else i1
        r0, r1 = rows[i0], rows[i1]
        # if the interval crosses a change of fixed plate, use the nearest node
        if r0[4] != r1[4]:
            r0 = r1 if abs(r1[0] - t) <= abs(t - r0[0]) else r0
            r1 = r0
        def node_h(r):
            if not perturb:
                return None
            h = perturb.get((plate, r[0], r[4]))
            return perturb.get((plate, r[0])) if h is None else h

        def node_quat(r, h):
            q = quat_from_pole(r[1], r[2], r[3])
            if h is not None:
                A = exp_M(h) @ quat_to_matrix(q)
                q = matrix_to_quat(A)
            return q
        h0, h1 = node_h(r0), node_h(r1)
        if r0[0] == r1[0]:
            return node_quat(r0, h0), r0[4]
        f = (t - r0[0]) / (r1[0] - r0[0])
        if INTERP_COV == "norm" and h0 is not None and h1 is not None:
            q = quat_slerp(node_quat(r0, None), node_quat(r1, None), f)
            h = ((1 - f) * h0 + f * h1) / np.sqrt((1 - f) ** 2 + f ** 2)
            return matrix_to_quat(exp_M(h) @ quat_to_matrix(q)), r0[4]
        return quat_slerp(node_quat(r0, h0), node_quat(r1, h1), f), r0[4]

    def total_rotation(self, plate, t, anchor=0, perturb=None, _depth=0):
        """Quaternion of the total rotation plate->absolute frame (000), or up to
        `anchor` if it is found in the chain. R_total = R_fixed_total * R_rel."""
        if plate == anchor or _depth > 30:
            return np.array([1.0, 0, 0, 0])
        if plate not in self.table:
            return np.array([1.0, 0, 0, 0])
        q, fixed = self.finite_rotation(plate, t, perturb)
        if fixed == anchor:
            return q
        qf = self.total_rotation(fixed, t, anchor, perturb, _depth + 1)
        return quat_mult(qf, q)

    def rel_rotation(self, plate, t, anchor, perturb=None):
        """R plate->anchor = R_anchor^-1 * R_plate (both to the absolute frame)."""
        qp = self.total_rotation(plate, t, 0, perturb)
        qa = self.total_rotation(anchor, t, 0, perturb)
        return quat_mult(quat_conj(qa), qp)

    def motion_path(self, seed_latlon, plate_seed, plate_rel, anchor, taus, t,
                    perturb=None):
        """Reproduces the GPlates motion path: position of the seed (fixed to
        plate_seed) 'deposited' on plate_rel at each age tau, seen in the
        `anchor` frame at age t. Returns an array (len(taus), 2) lat, lon."""
        s = None
        out = []
        for tau in taus:
            # position of the seed in the absolute frame at tau
            qs = self.total_rotation(plate_seed, tau, 0, perturb)
            p_abs = quat_to_matrix(qs) @ latlon_axis(*seed_latlon)
            # attached to the relative plate (plate coordinates at tau)
            qr = self.total_rotation(plate_rel, tau, 0, perturb)
            p_plate = quat_to_matrix(qr).T @ p_abs
            # seen in the anchor frame at age t
            qrt = self.total_rotation(plate_rel, t, 0, perturb)
            qat = self.total_rotation(anchor, t, 0, perturb)
            p = quat_to_matrix(qat).T @ (quat_to_matrix(qrt) @ p_plate)
            out.append(xyz_to_latlon(p))
        return np.array(out)

def matrix_to_quat(A):
    w = np.sqrt(max(0.0, 1 + A[0, 0] + A[1, 1] + A[2, 2])) / 2
    if w > 1e-8:
        x = (A[2, 1] - A[1, 2]) / (4 * w)
        y = (A[0, 2] - A[2, 0]) / (4 * w)
        z = (A[1, 0] - A[0, 1]) / (4 * w)
    else:
        d = np.diag(A)
        k = int(np.argmax(d))
        x = y = z = 0.0  # degenerate case (does not occur with small rotations)
    return np.array([w, x, y, z])

# ------------------------------------------------------------- validation
TRACK_ROOT = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "validation",
                          "motion_paths")  # GPlates motion paths used in the validation (S5)
SEEDS = {"Easter": (-26.3138, -106.599, 14),
         "San_Felix": (-26.28791, -81.42893, 15),
         "Juan Fernandez": (-33.750097, -81.905212, 13)}
MODELS = {"0.5": 902001, "0.36": 902002, "0.44": 902003}
FOLDER = {("Easter", "0.36"): "Easter_Pf=0.36",
          ("Easter", "0.5"):  "Easter_Pf=0.5",
          ("Easter", "0.44"): "Easter_Pf=0.44/Easter_Pf=0.44",
          ("San_Felix", "0.36"): "TT_Pf=0.36",
          ("San_Felix", "0.5"):  "TT_Pf=0.5",
          ("San_Felix", "0.44"): "TT_Pf=0.44/TT_Pf=0.44",
          ("Juan Fernandez", "0.36"): "JFR_Pf=0.36",
          ("Juan Fernandez", "0.5"):  "JFR_Pf=0.5",
          ("Juan Fernandez", "0.44"): "JFR_Pf=0.44/JFR_Pf=0.44"}

def read_xy(fn):
    lon, lat, age = [], [], []
    for line in open(fn):
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        p = line.split()
        if len(p) >= 3:
            lon.append(float(p[0])); lat.append(float(p[1])); age.append(float(p[2]))
    return np.array(lon), np.array(lat), np.array(age)

def gc_km(lat1, lon1, lat2, lon2):
    la1, la2 = np.radians(lat1), np.radians(lat2)
    dlon = np.radians(lon2 - lon1)
    d = np.arccos(np.clip(np.sin(la1) * np.sin(la2)
                          + np.cos(la1) * np.cos(la2) * np.cos(dlon), -1, 1))
    return d * 6371.0

if __name__ == "__main__":
    M = RotModel()
    print("plates parsed:", len(M.table))
    ncov = sum(1 for rows in M.table.values() for r in rows if r[5] is not None)
    print("nodes with covariance:", ncov)
    for t in (10.0, 36.0, 56.0):
        for name, (slat, slon, pid) in SEEDS.items():
            for pf, rel in MODELS.items():
                fn = os.path.join(TRACK_ROOT, FOLDER[(name, pf)],
                                  f"motion_path_output_{t:.2f}Ma.xy")
                if not os.path.exists(fn):
                    continue
                lon, lat, age = read_xy(fn)
                sim = M.motion_path((slat, slon), pid, rel, 201, age, t)
                err = gc_km(sim[:, 0], sim[:, 1], lat, lon)
                print(f"t={t:5.1f} {name:14s} Pf={pf:4s}: n={len(age):3d} "
                      f"err median {np.median(err):8.3f} km  max {err.max():8.3f} km")
