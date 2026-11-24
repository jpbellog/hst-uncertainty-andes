"""
Validation of geometry.py against the worked example of O'Neill et al. (2005) A7.

O'Neill A7: u = [1, -2, 1]  =>  F(u) is a right-handed rotation of 2.4495 rad
(140.3 degrees) about the axis [0.4082, -0.8165, 0.4082] = (24.09 N, 63.43 W), with
matrix (printed in the paper):
    [[-0.4749, -0.8505, -0.2261],
     [-0.3294,  0.4100, -0.8504],
     [ 0.8160, -0.3294, -0.4749]]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from chang import geometry as geo

print("=== Validation of geometry.py vs O'Neill A7 ===")
u = np.array([1.0, -2.0, 1.0])
print("u =", u, " |u| =", np.linalg.norm(u), "(expected 2.4495 rad = 140.3 deg)")
axis = u / np.linalg.norm(u)
print("axis =", np.round(axis, 4), "(expected [0.4082, -0.8165, 0.4082])")
lat, lon = geo.xyz_to_latlon(axis)
print(f"axis lat/lon = {lat:.2f}, {lon:.2f} (expected 24.09 N, 63.43 W = -63.43)")

F_paper = np.array([[-0.4749, -0.8505, -0.2261],
                    [-0.3294, 0.4100, -0.8504],
                    [0.8160, -0.3294, -0.4749]])
F_rod = geo.rodrigues(u)
F_quat = geo.quat_to_matrix(geo.exp_quat(u))
print("\nF via rodrigues:\n", np.round(F_rod, 4))
print("max|rodrigues - quat_to_matrix(exp)| =", np.max(np.abs(F_rod - F_quat)),
      "(must agree: two routes to the same map)")
print("max|rodrigues - F_paper| =", np.max(np.abs(F_rod - F_paper)),
      "(expected ~1e-4, rounding of the paper)")

print("\n=== Round-trips quaternion <-> matrix <-> axis/angle ===")
q = geo.rotation_to_quat(30.0, 45.0, 27.0)
lat, lon, rho = geo.quat_to_rotation(q)
print(f"rotation_to_quat(30,45,27) -> quat_to_rotation = {lat:.4f}, {lon:.4f}, {rho:.4f}")
A = geo.quat_to_matrix(q)
h, r = geo.matrix_to_axis_angle(A)
print("orthogonal matrix? max|A^T A - I| =", np.max(np.abs(A.T @ A - np.eye(3))))
print("det(A) =", np.linalg.det(A), "(expected 1)")
h_log = geo.log_quat(q)
print("does log_quat(q) recover rho*axis? |exp(log q) - q| =",
      np.max(np.abs(geo.exp_quat(h_log) - q)))
