"""
Validation of hellinger.py (fit by sections) with synthetic data.

The papers do not publish tables of raw crossings in Hellinger format, so the test
is the recovery of a KNOWN rotation from noisy synthetic data:
- A true rotation A_true is chosen.
- For each section i a normal eta_i is chosen; side 1 (moving) is points on the
  great circle perpendicular to eta_i; side 2 (fixed) is those points rotated by
  A_true (if q ⊥ eta_i then A q ⊥ A eta_i). Gaussian noise is added.
- The fit must recover A_true; H11.2 (positive definite), kappa_hat and the
  confidence region are checked.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from chang import geometry as geo
from chang import hellinger as hel
from chang import region as reg

rng = np.random.default_rng(12345)

# true rotation
lat_t, lon_t, rho_t = 32.0, -40.0, 24.0
A_true = geo.rotation_matrix(lat_t, lon_t, rho_t)

sd_km = 6.0
sd_rad = sd_km / geo.EARTH_KM_PER_RAD
nsect = 6
mpts = 6  # points per side per section


def random_perp(eta, rng):
    """random unit vector perpendicular to eta."""
    v = rng.normal(size=3)
    v = v - (v @ eta) * eta
    return v / np.linalg.norm(v)


def add_noise(p, sd_rad, rng):
    p2 = p + rng.normal(scale=sd_rad, size=3)
    return p2 / np.linalg.norm(p2)


crossings = []
for i in range(1, nsect + 1):
    eta = rng.normal(size=3)
    eta = eta / np.linalg.norm(eta)
    for _ in range(mpts):
        p = random_perp(eta, rng)          # point on the great circle (side 1)
        crossings.append((1, i, *geo.xyz_to_latlon(add_noise(p, sd_rad, rng)), sd_km))
    for _ in range(mpts):
        p = random_perp(eta, rng)
        q = A_true @ p                      # rotated to side 2
        crossings.append((2, i, *geo.xyz_to_latlon(add_noise(q, sd_rad, rng)), sd_km))

print("=== Validation of hellinger.py with synthetic data ===")
print(f"true rotation      : lat {lat_t}, lon {lon_t}, rho {rho_t}")
print(f"data: {len(crossings)} crossings, {nsect} sections, sd={sd_km} km")

# initial guess close to (but not equal to) the true one
res = hel.fit_hellinger(crossings, (lat_t + 3, lon_t - 4, rho_t + 2))
print(f"\nfitted rotation    : lat {res['lat']:.3f}, lon {res['lon']:.3f}, "
      f"rho {res['rho']:.3f}")
print(f"r_min (Phi_0)      : {res['r_min']:.4f}")
print(f"df, kappa_hat      : {res['df']}, {res['kappa_hat']:.3f}")

# angular error between A_true and A_fit
A_err = A_true.T @ res['A']
_, ang = geo.matrix_to_axis_angle(A_err)
print(f"angular error |A_true^-1 A_fit| = {np.rad2deg(ang):.4f} deg "
      f"(expected << rho; ~noise/sqrt(N))")

evals = np.linalg.eigvalsh(res['H112'])
print(f"\nH11.2 eigenvalues  : {np.round(evals, 1)} "
      f"(positive definite: {np.all(evals > 0)})")

xchi = reg.critical_value(0.95, df=res['df'], kappa_hat=res['kappa_hat'],
                          known_kappa=False)
print("\n95% confidence region (estimated kappa), maximum angles per axis:")
for r in reg.confidence_region_axes(res['H112'], xchi, A=res['A']):
    print(f"  axis (plate1) lat {r['lat_u']:7.2f} lon {r['lon_u']:8.2f}  "
          f"max_angle={r['max_angle_deg']:.3f} deg")

# Consistency check: kappa should be ~ (actual sd)^-2, expected ~1 if sd is set correctly
print(f"\nnote: kappa_hat ~ 1 indicates that sd={sd_km} km is consistent with the actual "
      "noise (kappa_hat>>1: sd overestimated; <<1: underestimated).")
