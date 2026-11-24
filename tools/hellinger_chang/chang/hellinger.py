"""
Rotation fit by the Hellinger criterion (by great-circle sections).

Port of hellinger1.f (Kirkwood et al. 1999 eq. 1, 12-15; Hellinger 1979 ch. III;
Chang 1988). Estimates the rotation that takes side 1 (moving) to side 2 (fixed) from
crossings organized by section, and builds the matrix H11.2 (inverse of the
covariance of the rotation).

Criterion (Kirkwood eq. 1):
    r(A, eta) = sum_i [ sum_j (u_ij^T eta_i)^2/sigma_ij^2
                       + sum_k (v_ik^T A eta_i)^2/sigma_ik^2 ]
For fixed A, r is minimized by choosing eta_i = eigenvector of the smallest eigenvalue of
    Sigma1_i + A^T Sigma2_i A,   with Sigma_i = sum kappa u u^T  (kappa = 1/sigma^2).
r(A) = sum of those smallest eigenvalues = Phi_0 (Hellinger III.5).

Conventions verified against the Fortran:
- input error in km; sigma_rad = sigma_km / 6366.197724; weight = 1/sigma_rad^2.
- parametrization A = A0 * exp(M(h)); minimization over h in R^3.
- df = n_data - 2*n_sections - 3;  kappa_hat = df / r_min.
"""
import numpy as np
from scipy.optimize import minimize
from . import geometry as geo


def _sigma_matrices(crossings):
    """
    Builds Sigma[side][isect] = sum_k axis axis^T / sigma_rad^2.
    crossings: list of (iside(1/2), isect(int), lat, lon, sd_km).
    Returns (sigma, sect_ids, ndat) with sigma[side][sect_index] a 3x3 matrix.
    """
    sects = sorted(set(c[1] for c in crossings))
    idx = {s: i for i, s in enumerate(sects)}
    ns = len(sects)
    sigma = {1: [np.zeros((3, 3)) for _ in range(ns)],
             2: [np.zeros((3, 3)) for _ in range(ns)]}
    ndat = 0
    for iside, isect, lat, lon, sd_km in crossings:
        if iside not in (1, 2):
            continue
        sd = sd_km / geo.EARTH_KM_PER_RAD
        u = geo.latlon_to_xyz(lat, lon)
        sigma[iside][idx[isect]] += np.outer(u, u) / sd ** 2
        ndat += 1
    return sigma, sects, ndat


def _O_basis(eta_vec):
    """
    Basis O_i (3x2) of the orthogonal complement of eta_i: two columns of M(eta)
    chosen by the dominant component of eta (numerical stability), replicating
    the logic of hellinger1.f (r1, lines 352-364).
    Columns of M(eta): c0=(0,z3,-z2), c1=(-z3,0,z1), c2=(z2,-z1,0).
    """
    z1, z2, z3 = eta_vec
    Meta = geo.skew(eta_vec)
    c0, c1, c2 = Meta[:, 0], Meta[:, 1], Meta[:, 2]
    if abs(z3) > 0.2:
        return np.column_stack([c0, c1])
    elif abs(z1) > 0.2:
        return np.column_stack([c1, c2])
    else:
        return np.column_stack([c0, c2])


def _fit_state(h, q0, sigma, ns):
    """Evaluates r(h) and also returns eta_i, M(eta_i), O_i and the S2_i per section."""
    q = geo.apply_perturbation(q0, h)
    A = geo.quat_to_matrix(q)
    r = 0.0
    etas, Metas, Os = [], [], []
    for i in range(ns):
        S2 = A.T @ sigma[2][i] @ A
        sig = sigma[1][i] + S2
        evals, evecs = np.linalg.eigh(sig)
        eta = evecs[:, 0]  # smallest eigenvalue
        r += evals[0]
        etas.append(eta)
        Metas.append(geo.skew(eta))
        Os.append(_O_basis(eta))
    return r, q, A, etas, Metas, Os


def fit_hellinger(crossings, init_latlon_rho, max_restarts=6, tol=1e-12):
    """
    Hellinger fit of one rotation (side 1 -> side 2).

    Parameters
    ----------
    crossings      : list of (iside, isect, lat, lon, sd_km).
    init_latlon_rho: (lat, lon, rho) initial guess of the axis and angle (degrees).
    max_restarts   : Nelder-Mead restarts (analogous to calling amoeba several times).

    Returns dict:
        q, A, lat, lon, rho : fitted rotation
        r_min               : criterion r(Ahat, eta_hat) = Phi_0
        H112                : matrix H11.2 (inverse of the covariance of h)
        cov_h               : covariance of h = H112^{-1} (WITHOUT division by kappa)
        kappa_hat, df, ndat, nsect
        XtX                 : design matrix (2*nsect+3)
    """
    sigma, sects, ndat = _sigma_matrices(crossings)
    ns = len(sects)
    lat0, lon0, rho0 = init_latlon_rho
    q0 = geo.rotation_to_quat(lat0, lon0, rho0)

    def objective(h):
        return _fit_state(h, q0, sigma, ns)[0]

    best_r = objective(np.zeros(3))
    for _ in range(max_restarts):
        res = minimize(objective, np.zeros(3), method="Nelder-Mead",
                       options=dict(xatol=1e-10, fatol=1e-14, maxiter=4000))
        r_new, q_new, _, _, _, _ = _fit_state(res.x, q0, sigma, ns)
        q0 = q_new  # move the base point (like calling amoeba again in the Fortran)
        if best_r - r_new < tol:
            best_r = r_new
            break
        best_r = r_new

    # Final state at h=0 around the converged q0
    r_min, q, A, etas, Metas, Os = _fit_state(np.zeros(3), q0, sigma, ns)

    # Design matrix X^T X (dimension 2*ns+3), by blocks (Kirkwood eq. 12-13)
    p = 2 * ns + 3
    XtX = np.zeros((p, p))
    H11 = np.zeros((3, 3))
    for i in range(ns):
        S2 = A.T @ sigma[2][i] @ A
        S1 = sigma[1][i]
        Meta, O = Metas[i], Os[i]
        H11 += Meta @ S2 @ Meta.T
        H12 = Meta @ S2 @ O            # (3x2)
        H22 = O.T @ (S1 + S2) @ O      # (2x2)
        r0, r1 = 3 + 2 * i, 3 + 2 * i + 2
        XtX[0:3, r0:r1] = H12
        XtX[r0:r1, 0:3] = H12.T
        XtX[r0:r1, r0:r1] = H22
    XtX[0:3, 0:3] = H11

    cov_full = np.linalg.inv(XtX)
    cov_h = cov_full[0:3, 0:3]     # covariance of h (= H112^{-1})
    H112 = np.linalg.inv(cov_h)
    df = ndat - 2 * ns - 3
    kappa_hat = df / r_min if r_min > 0 else np.inf
    lat, lon, rho = geo.quat_to_rotation(q)
    return dict(q=q, A=A, lat=lat, lon=lon, rho=rho, r_min=r_min,
                H112=H112, cov_h=cov_h, kappa_hat=kappa_hat, df=df,
                ndat=ndat, nsect=ns, XtX=XtX)
