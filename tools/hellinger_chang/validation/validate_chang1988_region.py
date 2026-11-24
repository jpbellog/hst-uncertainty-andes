"""
Validation of region.py (bingham_Q, endmembers) against Chang (1988) and the worked
example of Hanna and Chang (1990) for anomaly 13.

Data from the paper:
- H11.2 = [[1.420,-0.106,-0.281],[-0.106,0.265,-0.097],[-0.281,-0.097,0.113]] x 1e6
- qhat = 0.9706 + 0.2407*(0.1410 i - 0.2250 j + 0.9641 k)
- Q = 4 M(qhat)^T H11.2 M(qhat), printed by Hanna (x 1e4):
    [[  2.76,   1.09,  17.1,  -7.72],
     [  1.09, 571.0,   73.3, -71.0],
     [ 17.1,   73.3,  116.0, -55.0],
     [ -7.72, -71.0,  -55.0,  29.8]]
- Chang 1988: eigenvalues of H11.2 = 2644, 3.113e5, 1.483e6; eigenvector of the smallest
  (2644) at 62.6 S, 116.7 W.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from chang import geometry as geo
from chang import region as reg

H112 = np.array([[1.420, -0.106, -0.281],
                 [-0.106, 0.265, -0.097],
                 [-0.281, -0.097, 0.113]]) * 1e6
qhat = np.array([0.9706, 0.2407 * 0.1410, 0.2407 * -0.2250, 0.2407 * 0.9641])

Q_paper = np.array([[2.76, 1.09, 17.1, -7.72],
                    [1.09, 571.0, 73.3, -71.0],
                    [17.1, 73.3, 116.0, -55.0],
                    [-7.72, -71.0, -55.0, 29.8]]) * 1e4

print("=== Validation of region.bingham_Q vs Hanna and Chang (1988/1990) anomaly 13 ===")
print("|qhat| =", np.linalg.norm(qhat), "(must be ~1)")
Q = reg.bingham_Q(qhat, H112)
print("\ncomputed Q (x1e4):\n", np.round(Q / 1e4, 2))
print("\nQ of the paper (x1e4):\n", np.round(Q_paper / 1e4, 2))
rel = np.max(np.abs(Q - Q_paper)) / np.max(np.abs(Q_paper))
print(f"\nmax|Q - Q_paper| / max|Q_paper| = {rel:.4f}  "
      "(expected small; the Q of the paper is rounded to 3 figures)")

print("\n=== Eigenvalues of H11.2 (Chang 1988) ===")
evals, evecs = np.linalg.eigh(H112)
print("eigenvalues :", [f"{e:.4g}" for e in evals])
print("expected    : 2644, 3.113e5, 1.483e6")
z = evecs[:, 0]  # smallest eigenvalue
lat, lon = geo.xyz_to_latlon(z)
# normalize the hemisphere to compare with 62.6 S, 116.7 W
if lat > 0:
    lat, lon = -lat, lon + 180.0 if lon < 0 else lon - 180.0
print(f"eigenvector of the smallest eigenvalue: lat {lat:.1f}, lon {lon:.1f} "
      "(expected 62.6 S, 116.7 W)")

print("\n=== Maximum rotation angles (95% region, known kappa=1) ===")
xchi = reg.critical_value(0.95, known_kappa=True)  # chi2_0.95(3)
for r in reg.confidence_region_axes(H112, xchi):
    print(f"  axis lat {r['lat_u']:7.2f} lon {r['lon_u']:8.2f}  "
          f"lambda={r['eigenvalue']:.4g}  max_angle={r['max_angle_deg']:.3f} deg")
