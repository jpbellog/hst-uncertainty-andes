"""
Generates synthetic data in Hellinger format for two INDEPENDENT rotations with
UNEQUAL kappas, to produce a Fortran reference (hellinger1.exe + add2.exe) of the
combination pipeline with unequal kappas (Welch method of Kirkwood and Chang 1998).

Hellinger data file format (input of hellinger1):
    line 1: number of sections
    lines : iside(1/2)  isect  lat  lon  sd_km
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from chang import geometry as geo

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20260714)


def make_dataset(A_true, nsect, mpts, assign_sd_km, actual_sd_km, fname):
    """Generates a two-plate dataset (side1 -> side2 = A_true) and writes it."""
    sd_rad = actual_sd_km / geo.EARTH_KM_PER_RAD
    lines = [f"{nsect}"]
    ndat = 0
    for i in range(1, nsect + 1):
        eta = rng.normal(size=3); eta /= np.linalg.norm(eta)
        for side, A in ((1, np.eye(3)), (2, A_true)):
            for _ in range(mpts):
                v = rng.normal(size=3); v -= (v @ eta) * eta; v /= np.linalg.norm(v)
                p = A @ v + rng.normal(scale=sd_rad, size=3); p /= np.linalg.norm(p)
                lat, lon = geo.xyz_to_latlon(p)
                lines.append(f"{side} {i} {lat:.5f} {lon:.5f} {assign_sd_km:.2f}")
                ndat += 1
    with open(os.path.join(HERE, fname), "w") as f:
        f.write("\n".join(lines) + "\n")
    return ndat


# Rotation A: kappa ~ 1 (assigned sd = actual noise)
A_true = geo.rotation_matrix(20.0, -30.0, 15.0)
nA = make_dataset(A_true, 6, 6, assign_sd_km=8.0, actual_sd_km=8.0, fname="dataA.dat")

# Rotation B: kappa ~ 4 (assigned sd overestimated relative to the actual noise)
B_true = geo.rotation_matrix(50.0, 10.0, 20.0)
nB = make_dataset(B_true, 6, 6, assign_sd_km=8.0, actual_sd_km=4.0, fname="dataB.dat")

print("Datasets generated:")
print(f"  dataA.dat: rotation A=(20,-30,15), {nA} points, assigned sd 8 km, "
      "actual noise 8 km (kappa~1)")
print(f"  dataB.dat: rotation B=(50,10,20), {nB} points, assigned sd 8 km, "
      "actual noise 4 km (kappa~4)")
print("\nTrue combined rotation C = A * B:")
C_true = A_true @ B_true
lat, lon = geo.xyz_to_latlon(geo.matrix_to_axis_angle(C_true)[0])
ang = np.rad2deg(geo.matrix_to_axis_angle(C_true)[1])
print(f"  C=(lat {lat:.3f}, lon {lon:.3f}, rho {ang:.3f})")
