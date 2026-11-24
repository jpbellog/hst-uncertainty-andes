# chang (Python): uncertainty of plate and hotspot rotations

Python port of the Hellinger/Chang method for estimating plate rotations and their
uncertainties, and its extension to hotspot tracks (O'Neill et al. 2005). Ported from the
original Fortran code of T. Chang and checked against the papers.

## Installation and environment

Requires Python 3 with numpy, scipy and pandas. No installation is needed: it is a local
package. The scripts call `sys.path.insert` to find the `chang` package.

## Structure

```
hellinger_chang/
  chang/            main package
    geometry.py     quaternions, coordinates, exp/log maps (dtrans.f)
    stats_dist.py   chi2/F critical points (scipy; replaces XIDCH/XIDF)
    fit_uv.py       spherical regression of homologous pairs (u,v) (SPHEREREGRESS)
    hellinger.py    Hellinger fit by sections (hellinger1.f)
    triple.py       triple-junction fit, single kappa (hellinger3.f)
    region.py       confidence regions, Bingham Q, ellipse of a reconstructed point
    combine.py      inversion (eq. 20) and combination with equal kappas (eq. 19)
    tests.py        identity and difference tests (pval.f, F-test)
    hotspot.py      O'Neill extension (eq. 1, eq. 2, A6): uncertainty of hotspot tracks
  validation/       validation scripts against the papers and the Fortran
  examples/         usage example for hotspots
  reference/        Fortran reference outputs used by the validation
  PORT_NOTES.md     conventions, possible bugs of the Fortran, status
```

## Quick use

```python
from chang import fit_uv, hellinger, region, hotspot

# 1) Rotation of homologous pairs (u,v) in (lat,lon)
res = fit_uv.fit_uv(u_latlon, v_latlon)          # -> q, A, lat, lon, rho, SSE

# 2) Hellinger fit by sections (crossings per section)
#    crossings: list of (iside(1/2), isect, lat, lon, sd_km)
r = hellinger.fit_hellinger(crossings, init_latlon_rho=(30, -40, 25))
#    r['H112'], r['kappa_hat'], r['df'], r['A'] ...

# 3) Confidence region (maximum angles per axis)
xchi = region.critical_value(0.95, df=r['df'], kappa_hat=r['kappa_hat'], known_kappa=False)
axes = region.confidence_region_axes(r['H112'], xchi, A=r['A'])

# 4) Hotspots (O'Neill): age -> space, and ellipse of a reconstructed point
sig = hotspot.along_track_position_sigma(age=57.9, age_err=1.6, omega_deg_per_myr=0.7,
                                         omega_err_deg_per_myr=0.05)
ell = hotspot.reconstructed_hotspot_ellipse(r['A'], r['cov_h']/r['kappa_hat'],
                                            point_latlon=(-15, 88), alpha=0.95,
                                            kappa_hat=r['kappa_hat'], df=r['df'],
                                            known_kappa=False)
```

## Validation results

Run the scripts in `validation/` with `python`:

- `validate_geometry.py`: O'Neill A7 matrix to ~1e-4; exact round-trips.
- `validate_gulf_of_aden.py`: Gulf of Aden (Hanna 1990): axis 25.34/24.26 vs 25.31/24.29,
  angle 2.38 deg, SSE 1.284e-5 vs 1.279e-5, matrix to 4.4e-5.
- `validate_chang1988_region.py`: anomaly 13 (Chang 1988): Q to 0.05%, eigenvector of the
  smallest eigenvalue of H11.2 at 62.6 S / 116.7 W (exact).
- `validate_hellinger_synthetic.py`: recovers a known rotation to 0.03 deg with
  kappa_hat = 1.004 (consistency of the weighting).

- `validate_triple_and_tests.py`: synthetic triple junction: recovers A1 and A2 to <0.05 deg,
  df = n - 2s - 6.
- `validate_triple_vs_fortran.py`: `triple.py` against the compiled `hellinger3.exe`
  (data `reference/dataT.dat`): rotations 1->2, 1->3, 2->3 to <6e-6 deg, kappa_hat to 1e-6
  relative, identical df, H11.2 of 1->2 and 1->3 to 4e-8 relative (Frobenius norm).

Example: `examples/hotspot_oneill_example.py`.

## Scope and pending work

Ported and validated: core (fit of one rotation + region), regression of pairs (u,v),
basic inversion/combination, and the O'Neill extension for hotspots. See
`PORT_NOTES.md` for the pending parts (triple junction, combination with unequal kappas,
full 2D drawing of the Bingham region, Fstat test) and the possible defects of the
original Fortran found in the consistency check.
