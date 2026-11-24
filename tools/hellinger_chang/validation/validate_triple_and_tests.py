"""
Validation of triple.py (triple junction) and tests.py with synthetic data.

Triple junction: two true rotations A1 (plate1->2) and A2 (plate1->3) are chosen.
For each section i a normal eta_i is chosen; side 1 is points on the great circle
perpendicular to eta_i, side 2 = A1 @ (those points), side 3 = A2 @ (those points),
with noise. The fit must recover A1 and A2, and df = ndat - 2s - 6.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from chang import geometry as geo
from chang import triple, tests

rng = np.random.default_rng(2024)

lat1, lon1, rho1 = 30.0, -20.0, 18.0     # true A1
lat2, lon2, rho2 = -25.0, 60.0, 22.0     # true A2
A1t = geo.rotation_matrix(lat1, lon1, rho1)
A2t = geo.rotation_matrix(lat2, lon2, rho2)

sd_km = 5.0
sd_rad = sd_km / geo.EARTH_KM_PER_RAD
nsect, mpts = 6, 6


def perp(eta):
    v = rng.normal(size=3); v -= (v @ eta) * eta
    return v / np.linalg.norm(v)


def noisy(p):
    p = p + rng.normal(scale=sd_rad, size=3)
    return p / np.linalg.norm(p)


crossings = []
for i in range(1, nsect + 1):
    eta = rng.normal(size=3); eta /= np.linalg.norm(eta)
    for side, A in ((1, np.eye(3)), (2, A1t), (3, A2t)):
        for _ in range(mpts):
            p = A @ perp(eta)
            crossings.append((side, i, *geo.xyz_to_latlon(noisy(p)), sd_km))

print("=== Validation of triple.py (triple junction) ===")
print(f"true A1: {lat1}, {lon1}, {rho1}   true A2: {lat2}, {lon2}, {rho2}")
res = triple.fit_triple(crossings, (lat1 + 3, lon1 - 3, rho1 + 2),
                        (lat2 - 3, lon2 + 3, rho2 - 2))
print(f"fitted A1 : {res['lat1']:.3f}, {res['lon1']:.3f}, {res['rho1']:.3f}")
print(f"fitted A2 : {res['lat2']:.3f}, {res['lon2']:.3f}, {res['rho2']:.3f}")
for name, A_t, A_f in (("A1", A1t, res['A1']), ("A2", A2t, res['A2'])):
    _, ang = geo.matrix_to_axis_angle(A_t.T @ A_f)
    print(f"  angular error {name}: {np.rad2deg(ang):.4f} deg")
print(f"df = {res['df']} (expected {res['ndat']}-2*{res['nsect']}-6 = "
      f"{res['ndat'] - 2 * res['nsect'] - 6}),  kappa_hat = {res['kappa_hat']:.3f}")
print("H11.2(A1) positive definite:",
      np.all(np.linalg.eigvalsh(res['H112_1']) > 0),
      " H11.2(A3=A2 A1^T) positive definite:",
      np.all(np.linalg.eigvalsh(res['H112_3']) > 0))

print("\n=== Validation of tests.py ===")
# identity test on A1 (rotation of 18 deg): must REJECT the identity (p ~ 0)
idt = tests.identity_test(res['A1'], res['H112_1'], kappa_hat=res['kappa_hat'],
                          df=res['df'], known_kappa=False)
print(f"identity_test(A1): stat={idt['statistic']:.1f}, p={idt['p_value']:.2e}, "
      f"rejects identity={idt['reject_at_95']} (expected True: A1 != I)")
# identity test on a rotation close to the identity
A_small = geo.rotation_matrix(10, 20, 0.02)
idt2 = tests.identity_test(A_small, res['H112_1'], kappa_hat=res['kappa_hat'],
                           df=res['df'], known_kappa=False)
print(f"identity_test(rot 0.02 deg): p={idt2['p_value']:.3f}, "
      f"rejects={idt2['reject_at_95']} (expected False: ~identity)")
# pre-test of kappa equality
ket = tests.kappa_equality_test(1.0, 50, 1.1, 48)
print(f"kappa_equality_test(1.0/1.1): F={ket['F']:.3f}, p={ket['p_value']:.3f}, "
      f"equal={ket['equal']}")
