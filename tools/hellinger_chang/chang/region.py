"""
Confidence regions for rotations and ellipses of reconstructed points.

Implements:
- The "endmembers" of the confidence region (partial uncertainty rotations
  of Stock and Molnar): eigenvectors/eigenvalues of H11.2 -> maximum rotation angle
  per axis. (hellinger1.f / conreg.f; Kirkwood eq. 15; O'Neill A11.)
- The matrix Q = 4 M(qhat)^T H11.2 M(qhat) that feeds the graphical Bingham region
  (conreg.f / bingham.f; Hanna and Chang 1990).
- The uncertainty ellipse of a reconstructed point A w
  (ellipse.f / PTELLIPSE; Kirkwood eq. 22; O'Neill A6).

Note on kappa: H11.2 is such that the covariance of h is (kappa * H11.2)^{-1}. The
critical point is computed as:
    known kappa     : xchi = chi2_alpha(3)          (kappa=1 is assumed unless given)
    estimated kappa : xchi = 3 * F_alpha(3, df) / kappa_hat
and the maximum angle for axis j is sqrt(xchi / lambda_j), with lambda_j an eigenvalue of H11.2.
"""
import numpy as np
from . import geometry as geo
from . import stats_dist as sd


def critical_value(alpha, df=None, kappa_hat=1.0, known_kappa=True):
    """
    Critical point 'xchi' of the region (see the module docstring).
    - known_kappa=True: xchi = chi2_alpha(3) / kappa_hat (kappa_hat=1 by default).
    - known_kappa=False: xchi = 3 * F_alpha(3, df) / kappa_hat.
    """
    if known_kappa:
        return sd.chi2_crit(alpha, 3) / kappa_hat
    if df is None:
        raise ValueError("df is required when known_kappa=False")
    return 3.0 * sd.f_crit(alpha, 3, df) / kappa_hat


def confidence_region_axes(H112, xchi, A=None):
    """
    Endmembers of the confidence region from H11.2 (Kirkwood eq. 15).

    Returns a list of dicts, one per axis (eigenvector of H11.2):
        eigenvalue     : lambda_j
        axis_u         : eigenvector (axis on plate 1)
        lat_u, lon_u   : its lat/lon (degrees)
        max_angle_deg  : maximum admissible rotation angle = sqrt(xchi/lambda) (degrees)
        axis_v, lat_v, lon_v : axis on plate 2 (A @ eigenvector), if A is given
    """
    H112 = np.asarray(H112, dtype=float)
    evals, evecs = np.linalg.eigh(H112)
    out = []
    for j in range(3):
        lam = evals[j]
        z = evecs[:, j]
        ang_rad = np.sqrt(xchi / lam)
        lat_u, lon_u = geo.xyz_to_latlon(z)
        rec = dict(eigenvalue=float(lam), axis_u=z,
                   lat_u=lat_u, lon_u=lon_u,
                   max_angle_deg=float(ang_rad / geo.DEG))
        if A is not None:
            zv = np.asarray(A, dtype=float) @ z
            lat_v, lon_v = geo.xyz_to_latlon(zv)
            rec.update(axis_v=zv, lat_v=lat_v, lon_v=lon_v)
        out.append(rec)
    return out


def bingham_Q(qhat, H112):
    """
    Q = 4 M(qhat)^T H11.2 M(qhat)  (conreg.f L145-153; Hanna and Chang 1990).
    Non-negative symmetric 4x4 matrix with qhat as eigenvector of eigenvalue 0. It is
    passed to the Bingham subroutine/algorithm to draw the region in the form
    (min/max rotation angle as a function of the axis).
    """
    M = geo.M_3x4(qhat)  # 3x4
    return 4.0 * (M.T @ np.asarray(H112, dtype=float) @ M)


def reconstructed_point_ellipse(A, cov_u, w_latlon, xchi):
    """
    Uncertainty ellipse of the reconstructed point A w  (Kirkwood eq. 22; O'Neill A6;
    ellipse.f).

        cov(A w) = A M(w) cov_u M(w)^T A^T,
    with cov_u the covariance of the rotation (= H11.2^{-1}, on the appropriate kappa
    scale; see critical_value). cov(A w) has a null eigenvector in the direction A w;
    the ellipse lies in the tangent plane (A w)^perp.

    Parameters
    ----------
    A       : 3x3 rotation matrix
    cov_u   : 3x3 covariance of the rotation
    w_latlon: (lat, lon) of the point to reconstruct (degrees)
    xchi    : critical point (see critical_value); the individual ellipse uses
              2*F(2,df) or chi2(2); the simultaneous Scheffe one uses 3*F(3,df) or chi2(3).

    Returns dict:
        point_xyz, point_latlon : the reconstructed point A w
        semi_major, semi_minor  : semi-axes (radians) = sqrt(xchi * lambda)
        major_dir, minor_dir    : 3D vectors of the semi-axes
        cov_Aw                  : 3x3 matrix
    """
    A = np.asarray(A, dtype=float)
    cov_u = np.asarray(cov_u, dtype=float)
    w = geo.latlon_to_xyz(*w_latlon)
    Mw = geo.skew(w)
    cov_Aw = A @ Mw @ cov_u @ Mw.T @ A.T
    Aw = A @ w
    evals, evecs = np.linalg.eigh(cov_Aw)
    # the smallest eigenvalue (~0) corresponds to the direction A w; the other two are the ellipse
    order = np.argsort(evals)  # ascending
    lam_mid, lam_max = evals[order[1]], evals[order[2]]
    dir_mid, dir_max = evecs[:, order[1]], evecs[:, order[2]]
    semi_major = np.sqrt(xchi * lam_max)
    semi_minor = np.sqrt(xchi * lam_mid)
    return dict(point_xyz=Aw, point_latlon=geo.xyz_to_latlon(Aw),
                semi_major=float(semi_major), semi_minor=float(semi_minor),
                major_dir=dir_max, minor_dir=dir_mid, cov_Aw=cov_Aw)
