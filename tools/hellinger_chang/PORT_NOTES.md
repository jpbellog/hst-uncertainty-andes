# Notes on the Python port

Summary of the notes relevant to the port (conventions, possible defects of the
original Fortran, validation status), from an equation-by-code check of the papers
against the Fortran.

## 1. Critical conventions (verified against the Fortran)

- Weights = sigma^-2 with sigma in radians. The input error is in km and is
  converted with `R = 6366.197724 km/rad` (`EARTH_KM_PER_RAD`). `kappa = 1/sigma_rad^2`.
- `iside=1` = moving plate (rotated), `iside=2` = fixed plate. The estimated rotation
  takes side 1 to side 2.
- Parametrization `q_hat = q0 * exp(M(h))` (post-multiplication; TRANS2, corrected
  in Oct-1989 from `exp(2x)` to `exp(x)`). The minimized vector `h` is the `h` of
  Kirkwood / `t` of Chang / `u` of O'Neill (A7).
- The Fortran `mul(a,b,c)` computes `C = A*B` (standard product), even though the
  comment "c=b.a" is mislabelled. This affects addrot/invrot/add2/add3/add2tj/addtj4k.
- Sign in eq. (18) of the triple junction: the code (`hellinger3`/`add*tj`) uses
  `h3 ~ A1*(h2 - h1)` (MINUS sign in the cross terms), correct by derivation and
  opposite to eq. (18) as printed. Corrected in the Fortran in March-1991.
  Port the signs of the code, not those of the paper.
- Covariance scale: hellinger1 stores `H11.2^{-1}` WITHOUT division by kappa;
  to combine rotations (eq. 19) it must be brought to the sigma^2 scale by dividing by
  kappa. In this port, `combine.combine_rotations_equal_kappa` assumes that the covariance
  is passed ALREADY on the sigma^2 scale (documented in the function).
- `hatkapc = 3/cfact` in the programs with unequal kappas is a device that lets
  `conreg` (formula for equal kappas) draw the Welch region. In a clean port it is better
  to store `(cfact, df)` or the critical value directly, and not carry `hatkapc`.

## 2. Possible defects of the original Fortran (verify numerically before porting)

> All four affect ONLY the critical value of the confidence region, NOT the covariance
> matrix or the output rotation. Found by static reading of the code.
> None of them is in the modules already ported (core + region + O'Neill); they belong to the
> family of combination with unequal kappas, not yet ported.

1. `add3.f` (L202-204 + L298-306): probable double `sqrt(kappa)` scaling. It scales each
   design block twice (its ancestors do it once).
2. `addtj4k.f` (L362-364): block B unscaled. The loop `do 309 j=1,ncolx1+1,ncolx`
   iterates only `j=1` (column zero in the rows of B), leaving B unweighted by
   `sqrt(kappab)`. Compare with `add2tj.f` L369-371.
3. `averotplus.f` (L145): `nr = 2*if + 3*nsectot` is a CONFIRMED TYPO. It must be
   `3*if + 2*nsectot`. Kirkwood and Chang (1998, Appendix B p.100): the rank of the
   design separated by fit is `r_i = 2*s_i + 3`, so `rank = sum r_i = 3*(#fits) +
   2*(#sections)`. The coefficients of L145 are swapped.
4. `trikap.f` (`abc` in `teclib.f`): a commented-out `trace1` call leaves `t`
   uninitialized in the first contribution.

## 3. What was ported and its validation status

| Python module | Fortran / paper source | Validation |
|---------------|------------------------|-----------|
| `geometry.py` | dtrans.f | O'Neill A7 (matrix to ~1e-4); round-trips to machine precision; great_circle_fit vs O'Neill Fig. A1b (Ninetyeast): eigenvalues and normal exact to the rounding of the paper, eigenvectors to ~1e-9 |
| `stats_dist.py` | numlib XIDCH/XIDF -> scipy | (replaced by scipy, more accurate) |
| `fit_uv.py` | SPHEREREGRESS (Hanna 1990) | Gulf of Aden: axis 25.34/24.26 vs 25.31/24.29, SSE 1.284e-5 vs 1.279e-5, matrix to 4.4e-5 |
| `hellinger.py` | hellinger1.f | synthetic: recovers rotation to 0.03 deg, correct df, kappa_hat=1.004 |
| `triple.py` | hellinger3.f (eq. 16-18) | synthetic: recovers A1/A2 to ~0.04 deg, correct df=n-2s-6 |
| `region.py` | conreg.f/bingham.f/ellipse.f | Chang 1988: Q to 0.05%, smallest eigenvector at 62.6 S/116.7 W (exact) |
| `combine.py` | invrot.f (eq. 20), addrot.f (eq. 19), royer2/averotkap (average), teclib/K&C1998 (Welch) | invrot/addrot/average: direct algebra; welch_df/welch_critical_value vs K&C1998 full data (v=42.38 vs 42.34, d=3.049 exact); conservative_critical_value vs K&C1998 Theorem 5.2.2 (crit=10.83 vs 10.81, factor 1.64) |
| `tests.py` | pval.f, F-test (eq. 24-27; O'Neill A7) | synthetic: rejects identity p~1e-181 for an 18 deg rotation, does not reject for 0.02 deg |
| `hotspot.py` | O'Neill 2005 (eq. 1, eq. 2, A6/A9/A10) | coherent end-to-end example (the paper gives no exact numerical output) |

The residual differences for the Gulf of Aden and Chang 1988 are explained by the rounding of
the tables/matrices printed in the papers (Table 1 to 2 decimals; H11.2 to 3 decimals).

Validation against the compiled Fortran (validation/validate_vs_fortran.py,
reference in reference/): the port agrees with hellinger1.exe and add2.exe to ~1e-8:
- hellinger.py vs hellinger1.exe: identical rotation/kappa/df, H11.2 to 1e-8 relative.
- welch_df/welch_critical_value vs add2.exe: E, Var, df=91.984, c=3.00567, crit=8.126
  EXACT (starting from the capa/capc that add2 computes from the raw data).
The Fortran was compiled with gfortran 4.5.0 (Force 2.0); the cygwin gfortran 7.4.0 gives an
internal compiler error on the files with mixed `function/subroutine + common`.

## 4. Pending (not ported yet)

- Triple junction with 3 kappas: `trikap.f` (the single-kappa one, `hellinger3.f`, is ALREADY
  in `triple.py`).
- Combination with full data and unequal kappas: `add2.f`/`add3.f`/`add2tj.f`/
  `addtj4k.f`. The CRITICAL VALUE (Welch formulas `combine.welch_df`/`welch_critical_value`)
  is ALREADY ported and validated against Kirkwood and Chang (1998). What remains is only the
  pipeline that reconstructs the restricted design matrices (C=BA) and the projections P0/P1
  from the raw data to obtain the Johansen A, B, C (the delicate indexing part; it requires
  a Fortran reference run to validate it numerically end-to-end).
- Conservative combination (summary only): `addplus.f`/`averotplus.f`: PORTED AND
  VALIDATED in `combine.conservative_critical_value` (Kirkwood and Chang 1998, Theorem 5.2.2;
  reproduces the example of the paper crit~10.81). It requires only summary statistics
  (n_i, f_i), not raw data.
- royer2 kappa correction (capa/capc/cfact terms): the covariance combination
  and the mean rotation are ALREADY in `combine.average_rotations`; the correction based on
  projections (Chang et al. 1990) is missing.
- Full graphical Bingham region (`bingham.f`: BOUNDC/GRID, which draw the 2D boundary
  of the set of admissible axes and the grids of min/max angle). The construction of Q
  (`region.bingham_Q`), the endmembers and the ellipses were ported here; the 2D drawing was not.

The conservative unequal-kappa family (`addplus`, `averotplus`) and
`addtj4k` need Kirkwood and Chang (1998), J. Multivariate Anal., 65, 71-108. The
core, the single-kappa triple junction, the basic combinations and the Welch formula do
not need it.

## 5. Numerical replacements used

`JACOBI/JACOB2` -> `numpy.linalg.eigh`; `AMOEBA` -> `scipy.optimize.minimize(method=
'Nelder-Mead')`; `XIDCH/XDCH` -> `scipy.stats.chi2`; `XIDF/XDF` -> `scipy.stats.f`;
`dpotrf/dpotri` -> `numpy.linalg.inv` (or `scipy.linalg.cho_*`).
