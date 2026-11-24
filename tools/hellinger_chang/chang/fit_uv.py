"""
Spherical regression for pairs of homologous points (u_i, v_i).

Port of the core of SPHEREREGRESS (Hanna and Chang 1990) / spherical regression of
Chang (1986). Unlike the Hellinger fit by sections (see hellinger.py), here the data
are "homologous" pairs: without error v_i = A0 u_i, and NO assumption is made about
pairing by segments.

Mathematics (Moran 1976; Hanna and Chang 1990 eqs. 4-6):
    SSE(A) = sum_i |v_i - A u_i|^2 = 2n - 2 sum_i v_i^T A u_i = q^T Q q,
with q the quaternion of A, Q = 2n I - 2n R, and R (4x4 symmetric) built from
X = (1/n) sum_i u_i v_i^T (the "UVT/N" matrix of the Fortran).
The fit minimizes SSE => maximizes q^T R q => q_hat = eigenvector of the LARGEST
eigenvalue of R (equivalent to the smallest eigenvalue of Q).
"""
import numpy as np
from . import geometry as geo


def _R_from_X(X):
    """Symmetric 4x4 matrix R from X (formulas of Hanna and Chang 1990)."""
    x11, x12, x13 = X[0, 0], X[0, 1], X[0, 2]
    x21, x22, x23 = X[1, 0], X[1, 1], X[1, 2]
    x31, x32, x33 = X[2, 0], X[2, 1], X[2, 2]
    r11 = x11 + x22 + x33
    r12 = x23 - x32
    r13 = x31 - x13
    r14 = x12 - x21
    r22 = x11 - x22 - x33
    r23 = x12 + x21
    r24 = x13 + x31
    r33 = -x11 + x22 - x33
    r34 = x23 + x32
    r44 = -x11 - x22 + x33
    R = np.array([[r11, r12, r13, r14],
                  [r12, r22, r23, r24],
                  [r13, r23, r33, r34],
                  [r14, r24, r34, r44]])
    return R


def fit_uv(u_latlon, v_latlon):
    """
    Fits the rotation A such that v_i ~ A u_i for pairs (u_i, v_i).

    Parameters
    ----------
    u_latlon, v_latlon : array (n, 2) of (lat, lon) in degrees.

    Returns dict with:
        q       : fitted quaternion (w,x,y,z)
        A       : 3x3 rotation matrix
        lat, lon, rho : axis (degrees) and angle (degrees)
        R, Q    : 4x4 matrices
        X       : matrix (1/n) sum u v^T
        sse     : SSE(A) = q^T Q q
        n       : number of pairs
    """
    u_latlon = np.asarray(u_latlon, dtype=float)
    v_latlon = np.asarray(v_latlon, dtype=float)
    n = len(u_latlon)
    U = np.array([geo.latlon_to_xyz(la, lo) for la, lo in u_latlon])  # (n,3)
    V = np.array([geo.latlon_to_xyz(la, lo) for la, lo in v_latlon])
    X = (U.T @ V) / n  # (1/n) sum u_i v_i^T
    R = _R_from_X(X)
    Q = 2.0 * n * np.eye(4) - 2.0 * n * R
    # q_hat = eigenvector of the largest eigenvalue of R
    evals, evecs = np.linalg.eigh(R)
    q = evecs[:, np.argmax(evals)]
    if q[0] < 0:
        q = -q  # convention: positive scalar part
    A = geo.quat_to_matrix(q)
    lat, lon, rho = geo.quat_to_rotation(q)
    sse = float(q @ Q @ q)
    return dict(q=q, A=A, lat=lat, lon=lon, rho=rho, R=R, Q=Q, X=X, sse=sse, n=n)
