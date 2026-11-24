"""
Generates dataT.dat: synthetic triple-junction data in Hellinger format
(line 1 = n sections; then `iside isect lat lon sd_km`), to compare
triple.py against hellinger3.exe.

Same as validate_triple_and_tests.py: true A1 (plate1->2) and A2 (plate1->3);
one normal eta_i per section; side 1 on the great circle perpendicular to eta_i,
side 2 = A1 @ points, side 3 = A2 @ points, with Gaussian noise of sd_km.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from chang import geometry as geo

rng = np.random.default_rng(2024)
A1t = geo.rotation_matrix(30.0, -20.0, 18.0)
A2t = geo.rotation_matrix(-25.0, 60.0, 22.0)
sd_km = 5.0
sd_rad = sd_km / geo.EARTH_KM_PER_RAD
nsect, mpts = 6, 6


def perp(eta):
    v = rng.normal(size=3); v -= (v @ eta) * eta
    return v / np.linalg.norm(v)


def noisy(p):
    p = p + rng.normal(scale=sd_rad, size=3)
    return p / np.linalg.norm(p)


here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "dataT.dat"), "w") as f:
    f.write(f"{nsect}\n")
    for i in range(1, nsect + 1):
        eta = rng.normal(size=3); eta /= np.linalg.norm(eta)
        for side, A in ((1, np.eye(3)), (2, A1t), (3, A2t)):
            for _ in range(mpts):
                lat, lon = geo.xyz_to_latlon(noisy(A @ perp(eta)))
                f.write(f"{side} {i} {lat:.6f} {lon:.6f} {sd_km:.3f}\n")
print("dataT.dat written")
