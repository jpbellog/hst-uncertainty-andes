"""
Example of use of chang.hotspot (method of O'Neill et al. 2005) for the uncertainty
of positions on hotspot tracks.

This example ILLUSTRATES the workflow and the equations (it does not reproduce the exact
rotations of O'Neill, which require their full dataset, not exactly legible from
the scanned table). It uses realistic numbers in the style of the Ninetyeast Ridge (India plate).

Workflow:
  1. Convert age uncertainty to along-track spatial uncertainty (eq. 2).
  2. Build error ellipses per point (along-track >> across-track).
  3. Fit an absolute rotation with the Hellinger criterion (reuse, Appendix A5).
  4. Compute the uncertainty ellipse of a reconstructed point (eq. A6 / eq 22).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from chang import geometry as geo
from chang import hotspot as hs
from chang import region as reg

print("=== 1. Age -> space conversion (O'Neill eq. 2) ===")
# typical angular velocity of the India plate in the Tertiary ~ 0.7 deg/Myr (+/- 0.05)
omega = 0.70
omega_err = 0.05
# some points in the style of the Ninetyeast Ridge (age Ma, error Ma) from Table 1
puntos = [("DSDP254", 37.8, 1.3), ("ODP756D", 43.2, 0.5),
          ("ODP757", 57.9, 1.6), ("DSDP214", 60.9, 1.3)]
print(f"omega = {omega} +/- {omega_err} deg/Myr")
print(f"{'sample':10s} {'age':>6s} {'+/-':>5s}  sigma_along(deg)")
for name, age, aerr in puntos:
    s = hs.along_track_position_sigma(age, aerr, omega, omega_err)
    print(f"{name:10s} {age:6.1f} {aerr:5.1f}  {s:8.3f}   "
          f"(= {s*111.2:.0f} km)")

print("\n=== 2. Error ellipse of a point (along >> across) ===")
# across-track of ~0.3 deg (from the full width at half maximum of the gravity profile)
ell = hs.point_error_ellipse(-30.97, 87.90,
                             along_track_sigma_deg=hs.along_track_position_sigma(
                                 37.8, 1.3, omega, omega_err),
                             across_track_sigma_deg=0.30,
                             track_azimuth_deg=10.0)
print(f"center {ell['center_latlon']}, sigma_along={ell['sigma_along']:.3f} deg, "
      f"sigma_across={ell['sigma_across']:.3f} deg, azimuth={ell['azimuth_deg']} deg")
print("cov 2x2 (East,North) deg^2:\n", np.round(ell['cov2d'], 4))

print("\n=== 3. Absolute rotation fit (Hellinger; illustrative/synthetic) ===")
# 'true' rotation of the plate to its hotspot frame
A_true = geo.rotation_matrix(20.0, 30.0, 15.0)
rng = np.random.default_rng(7)
sd_km = 40.0  # large uncertainties, typical of hotspots
sd_rad = sd_km / geo.EARTH_KM_PER_RAD
track, hotspots = [], []
for i in range(1, 5):
    eta = rng.normal(size=3); eta /= np.linalg.norm(eta)
    for _ in range(4):
        v = rng.normal(size=3); v -= (v @ eta) * eta; v /= np.linalg.norm(v)
        p = v + rng.normal(scale=sd_rad, size=3); p /= np.linalg.norm(p)
        track.append((i, *geo.xyz_to_latlon(p), sd_km))          # track geometry
        q = A_true @ v + rng.normal(scale=sd_rad, size=3); q /= np.linalg.norm(q)
        hotspots.append((i, *geo.xyz_to_latlon(q), sd_km))       # hotspot position
res = hs.fit_hotspot_rotation(track, hotspots, (22.0, 27.0, 13.0))
print(f"fitted rotation: lat {res['lat']:.2f}, lon {res['lon']:.2f}, "
      f"rho {res['rho']:.2f}  (true 20, 30, 15)")
print(f"df={res['df']}, kappa_hat={res['kappa_hat']:.3f}")

print("\n=== 4. Uncertainty of a reconstructed position (O'Neill A6 / eq 22) ===")
# covariance of the rotation on the sigma^2 scale (divide the implicit H11.2^{-1} by kappa)
cov_u = res['cov_h'] / res['kappa_hat']
punto = (-15.0, 88.0)  # track point to reconstruct
ell_rec = hs.reconstructed_hotspot_ellipse(res['A'], cov_u, punto, alpha=0.95,
                                           kappa_hat=res['kappa_hat'], df=res['df'],
                                           known_kappa=False, simultaneous=True)
plat, plon = ell_rec['point_latlon']
print(f"reconstructed point: lat {plat:.2f}, lon {plon:.2f}")
print(f"semi-major axis = {np.rad2deg(ell_rec['semi_major']):.3f} deg "
      f"({np.rad2deg(ell_rec['semi_major'])*111.2:.0f} km)")
print(f"semi-minor axis = {np.rad2deg(ell_rec['semi_minor']):.3f} deg "
      f"({np.rad2deg(ell_rec['semi_minor'])*111.2:.0f} km)")
print("\n(The error ellipse of the reconstructed position is the sought uncertainty "
      "of the hotspot-track position, propagated from the data through the covariance "
      "of the rotation.)")
