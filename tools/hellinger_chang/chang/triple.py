"""
Fit of rotations at a triple junction by the Hellinger criterion,
single kappa on the three branches.

Port of hellinger3.f (Kirkwood et al. 1999 eq. 16-18; df = n - 2s - 6). Estimates two
independent rotations A1 (plate 1 -> 2) and A2 (plate 1 -> 3); the third is derived
as A3 = A2 * A1^T (closure).

Criterion (eq. 16):
    r(A1,A2,eta) = sum_i smallest_eigenvalue(Sigma1_i + A1^T Sigma2_i A1 + A2^T Sigma3_i A2)
with A1 = A10 exp(M(h1)), A2 = A20 exp(M(h2)), h = [h1, h2] in R^6.
"""
import numpy as np
from scipy.optimize import minimize
from . import geometry as geo
from .hellinger import _O_basis


def _sigma3(crossings):
    """sigma[side(1/2/3)][sect] = sum axis axis^T / sigma_rad^2."""
    sects = sorted(set(c[1] for c in crossings))
    idx = {s: i for i, s in enumerate(sects)}
    ns = len(sects)
    sigma = {s: [np.zeros((3, 3)) for _ in range(ns)] for s in (1, 2, 3)}
    ndat = 0
    for iside, isect, lat, lon, sd_km in crossings:
        if iside not in (1, 2, 3):
            continue
        sd = sd_km / geo.EARTH_KM_PER_RAD
        u = geo.latlon_to_xyz(lat, lon)
        sigma[iside][idx[isect]] += np.outer(u, u) / sd ** 2
        ndat += 1
    return sigma, sects, ndat


def _state(h, q10, q20, sigma, ns):
    q1 = geo.apply_perturbation(q10, h[0:3])
    q2 = geo.apply_perturbation(q20, h[3:6])
    A1 = geo.quat_to_matrix(q1)
    A2 = geo.quat_to_matrix(q2)
    r = 0.0
    etas, Metas, Os = [], [], []
    for i in range(ns):
        sig = sigma[1][i] + A1.T @ sigma[2][i] @ A1 + A2.T @ sigma[3][i] @ A2
        evals, evecs = np.linalg.eigh(sig)
        eta = evecs[:, 0]
        r += evals[0]
        etas.append(eta)
        Metas.append(geo.skew(eta))
        Os.append(_O_basis(eta))
    return r, q1, q2, A1, A2, etas, Metas, Os


def fit_triple(crossings, init1, init2, max_restarts=6, tol=1e-12):
    """
    Triple-junction fit.

    Parameters
    ----------
    crossings : list of (iside(1/2/3), isect, lat, lon, sd_km).
    init1     : (lat, lon, rho) guess for A1 (plate 1 -> 2), degrees.
    init2     : (lat, lon, rho) guess for A2 (plate 1 -> 3), degrees.

    Returns dict with, for each rotation (1, 2, and the combined 3 = A2 A1^T):
        A1, A2, A3, lat1/lon1/rho1, lat2/lon2/rho2, lat3/lon3/rho3
        H112_1, H112_2, H112_3, cov1, cov2, cov3, cov12 (cross)
        r_min, df, kappa_hat, ndat, nsect
    """
    sigma, sects, ndat = _sigma3(crossings)
    ns = len(sects)
    q10 = geo.rotation_to_quat(*init1)
    q20 = geo.rotation_to_quat(*init2)

    def objective(h):
        return _state(h, q10, q20, sigma, ns)[0]

    best = objective(np.zeros(6))
    for _ in range(max_restarts):
        res = minimize(objective, np.zeros(6), method="Nelder-Mead",
                       options=dict(xatol=1e-10, fatol=1e-14, maxiter=8000))
        r_new, q1n, q2n, *_ = _state(res.x, q10, q20, sigma, ns)
        q10, q20 = q1n, q2n
        if best - r_new < tol:
            best = r_new
            break
        best = r_new

    r_min, q1, q2, A1, A2, etas, Metas, Os = _state(np.zeros(6), q10, q20, sigma, ns)

    # X^T X of dimension 2*ns + 6 (blocks h1[0:3], h2[3:6], xi_i)
    p = 2 * ns + 6
    XtX = np.zeros((p, p))
    H11a = np.zeros((3, 3))  # block h1,h1
    H11b = np.zeros((3, 3))  # block h2,h2
    for i in range(ns):
        S2 = A1.T @ sigma[2][i] @ A1
        S3 = A2.T @ sigma[3][i] @ A2
        S1 = sigma[1][i]
        Meta, O = Metas[i], Os[i]
        H11a += Meta @ S2 @ Meta.T
        H11b += Meta @ S3 @ Meta.T
        H12a = Meta @ S2 @ O   # h1 x xi_i
        H12b = Meta @ S3 @ O   # h2 x xi_i
        H22 = O.T @ (S1 + S2 + S3) @ O
        c0, c1 = 6 + 2 * i, 6 + 2 * i + 2
        XtX[0:3, c0:c1] = H12a
        XtX[c0:c1, 0:3] = H12a.T
        XtX[3:6, c0:c1] = H12b
        XtX[c0:c1, 3:6] = H12b.T
        XtX[c0:c1, c0:c1] = H22
    XtX[0:3, 0:3] = H11a
    XtX[3:6, 3:6] = H11b
    # cross block h1,h2 = 0 (no datum involves A1 and A2 at the same time)

    cov_full = np.linalg.inv(XtX)
    cov1 = cov_full[0:3, 0:3]
    cov2 = cov_full[3:6, 3:6]
    cov12 = cov_full[0:3, 3:6]
    H112_1 = np.linalg.inv(cov1)
    H112_2 = np.linalg.inv(cov2)

    # Combined rotation A3 = A2 A1^T; covariance eq. (18) with the sign of the code:
    #   cov(h3) = A1 [cov1 + cov2 - cov12 - cov12^T] A1^T
    A3 = A2 @ A1.T
    cov3 = A1 @ (cov1 + cov2 - cov12 - cov12.T) @ A1.T
    H112_3 = np.linalg.inv(cov3)

    df = ndat - 2 * ns - 6
    kappa_hat = df / r_min if r_min > 0 else np.inf
    l1, o1, p1 = geo.quat_to_rotation(q1)
    l2, o2, p2 = geo.quat_to_rotation(q2)
    A3q = geo.matrix_to_axis_angle(A3)
    l3, o3 = geo.xyz_to_latlon(A3q[0])
    p3 = np.rad2deg(A3q[1])
    return dict(A1=A1, A2=A2, A3=A3,
                lat1=l1, lon1=o1, rho1=p1, lat2=l2, lon2=o2, rho2=p2,
                lat3=l3, lon3=o3, rho3=p3,
                H112_1=H112_1, H112_2=H112_2, H112_3=H112_3,
                cov1=cov1, cov2=cov2, cov3=cov3, cov12=cov12,
                r_min=r_min, df=df, kappa_hat=kappa_hat, ndat=ndat, nsect=ns)
