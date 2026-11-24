"""
Validation of triple.py DIRECTLY against hellinger3.exe (original Fortran, compiled).

Uses ../reference/dataT.dat (gen_reference_triple.py) and the output ../reference/fitT_run.log
of hellinger3.exe, and compares: the three rotations (1->2, 1->3, 2->3), kappa_hat, df and the
H11.2 matrices of the rotations 1->2 and 1->3.

hellinger3.exe was compiled with MinGW gfortran 6.3.0 (`-std=legacy -w`, hellinger3.f dtrans.f
numlib.f bingham.f). Run command (see ../reference/README.md):
  printf 'dataT.dat\n33 -23 20\n-28 63 20\nn\n0.95\ny\nn\nn\nn\nn\n' | ./bin/hellinger3.exe > fitT_run.log
"""
import sys, os, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chang import triple, geometry as geo

REF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reference")
NUM = r"[-+]?\d+\.\d*(?:[EeDd][-+]?\d+)?"


def read_dat(path):
    """Reads a Hellinger-format .dat -> list of (iside, isect, lat, lon, sd_km)."""
    out = []
    with open(path) as f:
        f.readline()  # nsect
        for ln in f:
            p = ln.split()
            if len(p) >= 5:
                out.append((int(p[0]), int(p[1]), float(p[2]), float(p[3]), float(p[4])))
    return out


def parse_log(path):
    lines = open(path).read().splitlines()
    out = {"rot": {}, "H112": {}}
    for i, l in enumerate(lines):
        m = re.search(r"fitted rotation side\s+(\d)\s+to side\s+(\d)", l)
        if m:
            out["rot"][(int(m.group(1)), int(m.group(2)))] = [float(x) for x in re.findall(NUM, lines[i + 1])]
        if "kappahat, degrees of freedom" in l:
            out["kappa"], out["df"] = [float(x) for x in re.findall(NUM, l)][:2]
        m = re.search(r"simultaneous confidence region for rotationside\s+(\d)\s+to\s+(\d)", l)
        if m and "matrix h11.2" in lines[i + 1]:
            v = [float(x) for x in re.findall(NUM, lines[i + 2])]
            H = np.zeros((3, 3)); k = 0
            for a in range(3):          # lower-triangular packing by rows
                for b in range(a + 1):
                    H[a, b] = H[b, a] = v[k]; k += 1
            out["H112"].setdefault((int(m.group(1)), int(m.group(2))), H)
    return out


def latlonrho(A):
    axis, ang = geo.matrix_to_axis_angle(A)
    lat, lon = geo.xyz_to_latlon(axis)
    return lat, lon, np.rad2deg(ang)


def axis_angle_diff(a, b):
    """Angular difference (deg) between two rotations given as (lat, lon, rho)."""
    A = geo.rotation_matrix(*a); B = geo.rotation_matrix(*b)
    return np.rad2deg(geo.matrix_to_axis_angle(A.T @ B)[1])


data = read_dat(os.path.join(REF, "dataT.dat"))
f = parse_log(os.path.join(REF, "fitT_run.log"))
r = triple.fit_triple(data, (33, -23, 20), (-28, 63, 20))

print("=== triple.py vs hellinger3.exe (dataT.dat) ===")
py_rot = {(1, 2): (r["lat1"], r["lon1"], r["rho1"]),
          (1, 3): (r["lat2"], r["lon2"], r["rho2"]),
          (2, 3): latlonrho(r["A3"])}
for k in ((1, 2), (1, 3), (2, 3)):
    fo, pe = f["rot"][k], py_rot[k]
    print(f"rot {k[0]}->{k[1]}  Fortran {fo[0]:9.4f} {fo[1]:10.4f} {fo[2]:8.4f} | "
          f"Python {pe[0]:9.4f} {pe[1]:10.4f} {pe[2]:8.4f} | diff {axis_angle_diff(fo, pe):.2e} deg")
print(f"kappa_hat Fortran {f['kappa']:.6f}  Python {r['kappa_hat']:.6f}  "
      f"rel diff {abs(f['kappa'] - r['kappa_hat']) / f['kappa']:.2e}")
print(f"df        Fortran {f['df']:.0f}  Python {r['df']}")
for k, key in (((1, 2), "H112_1"), ((1, 3), "H112_2")):
    Hf, Hp = f["H112"][k], r[key]
    rel = np.linalg.norm(Hf - Hp) / np.linalg.norm(Hf)
    print(f"H11.2 {k[0]}->{k[1]}: relative difference (Frobenius norm) {rel:.2e}")
