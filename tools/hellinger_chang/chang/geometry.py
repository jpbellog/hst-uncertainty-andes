"""
Geometry of rotations on the sphere: coordinates, quaternions and the
exponential/logarithmic maps.

Direct port of dtrans.f (TRANS1..TRANS9, QMULT). Conventions verified against
the Fortran and against Chang (1988) / Kirkwood et al. (1999) / Hanna and Chang (1990)
/ O'Neill et al. (2005) Appendix A.

Key conventions:
- Quaternion q = (w, x, y, z), w = scalar part. q represents the right-handed
  rotation of rho radians about the unit axis u:
      q = (cos(rho/2), sin(rho/2)*u).
- The rotation matrix A acts on column vectors: v_rotated = A @ v.
- "Moving exponential" perturbation: qhat = q * exp(h), where exp(h) is the quaternion
  of the rotation of |h| rad about h/|h| (parametrization of Chang 1988 /
  Kirkwood Phi(h) = exp(M(h)) / O'Neill A7). The order is POST-multiplication
  (TRANS2 of dtrans.f, corrected in Oct-1989 to exp(x) instead of exp(2x)).
- Input/output angles of the trans_* functions in DEGREES (as in the Fortran).
"""
import numpy as np

DEG = np.pi / 180.0
# Earth radius used by the Fortran to convert km <-> radians (km per radian).
EARTH_KM_PER_RAD = 6366.197724


def latlon_to_xyz(lat_deg, lon_deg):
    """TRANS1: lat/lon (degrees) -> unit Cartesian vector (x, y, z)."""
    lat = np.asarray(lat_deg) * DEG
    lon = np.asarray(lon_deg) * DEG
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return np.array([x, y, z])


def xyz_to_latlon(u):
    """TRANS6: Cartesian vector -> (lat, lon) in degrees."""
    u = np.asarray(u, dtype=float)
    z = np.clip(u[2], -1.0, 1.0)
    lat = np.arcsin(z) / DEG
    if u[0] ** 2 + u[1] ** 2 <= 0.0:
        return lat, 0.0
    lon = np.arctan2(u[1], u[0]) / DEG
    return lat, lon


def rotation_to_quat(lat_deg, lon_deg, rho_deg):
    """TRANS3: (axis lat, lon, angle rho in degrees) -> quaternion (w,x,y,z)."""
    lat = lat_deg * DEG
    lon = lon_deg * DEG
    rho = rho_deg * DEG
    w = np.cos(rho / 2.0)
    s = np.sin(rho / 2.0)
    x = s * np.cos(lat) * np.cos(lon)
    y = s * np.cos(lat) * np.sin(lon)
    z = s * np.sin(lat)
    return np.array([w, x, y, z])


def quat_to_rotation(q):
    """TRANS5: quaternion -> (axis lat, lon, angle rho in degrees)."""
    q = np.asarray(q, dtype=float)
    w = np.clip(q[0], -1.0, 1.0)
    ang = np.arccos(w)
    rho = ang * 2.0 / DEG
    s = np.sin(ang)
    if abs(s) < 1e-30:
        return 0.0, 0.0, 0.0
    z = np.clip(q[3] / s, -1.0, 1.0)
    lat = np.arcsin(z) / DEG
    lon = np.arctan2(q[2], q[1]) / DEG
    return lat, lon, rho


def quat_mult(a, b):
    """QMULT: Hamilton quaternion product c = a*b, q=(w,x,y,z)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    w = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
    x = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2]
    y = a[0] * b[2] + a[2] * b[0] + a[3] * b[1] - a[1] * b[3]
    z = a[0] * b[3] + a[3] * b[0] + a[1] * b[2] - a[2] * b[1]
    return np.array([w, x, y, z])


def quat_conj(a):
    a = np.asarray(a, dtype=float)
    return np.array([a[0], -a[1], -a[2], -a[3]])


def exp_quat(h):
    """
    Quaternion exp(h) for h in R^3 (TRANS2): rotation of |h| rad about h/|h|.
    exp(h) = (cos(|h|/2), sin(|h|/2)/|h| * h).
    """
    h = np.asarray(h, dtype=float)
    theta = np.linalg.norm(h)
    if theta <= 0.0:
        return np.array([1.0, 0.0, 0.0, 0.0])
    fact = np.sin(theta / 2.0) / theta
    return np.array([np.cos(theta / 2.0), fact * h[0], fact * h[1], fact * h[2]])


def apply_perturbation(q, h):
    """qhat = q * exp(h)  (moving exponential parametrization, TRANS2)."""
    return quat_mult(q, exp_quat(h))


def log_quat(q):
    """
    TRANS8: solves q = exp(h) for h in R^3 (logarithmic map).
    Returns h = rho * axis (rotation vector). Safeguards for |w|>=1.
    """
    q = np.asarray(q, dtype=float)
    w = q[0]
    if abs(w) >= 1.0:
        return np.zeros(3)
    if w < 0.0:
        temp = -w
        ind = -1.0
    else:
        temp = w
        ind = 1.0
    fact = np.arccos(temp)
    if abs(fact) < 1e-30:
        fact = 0.0
    else:
        fact = 2.0 * ind * fact / np.sin(fact)
    return fact * q[1:4]


def quat_to_matrix(q):
    """TRANS4: quaternion -> 3x3 rotation matrix (v_rot = A @ v)."""
    q = np.asarray(q, dtype=float)
    a0, a1, a2, a3 = q
    A = np.empty((3, 3))
    A[0, 0] = a0 * a0 + a1 * a1 - a2 * a2 - a3 * a3
    A[1, 0] = 2.0 * (a0 * a3 + a1 * a2)
    A[2, 0] = 2.0 * (a1 * a3 - a0 * a2)
    A[0, 1] = 2.0 * (a1 * a2 - a0 * a3)
    A[1, 1] = a0 * a0 - a1 * a1 + a2 * a2 - a3 * a3
    A[2, 1] = 2.0 * (a0 * a1 + a2 * a3)
    A[0, 2] = 2.0 * (a0 * a2 + a1 * a3)
    A[1, 2] = 2.0 * (a2 * a3 - a0 * a1)
    A[2, 2] = a0 * a0 - a1 * a1 - a2 * a2 + a3 * a3
    return A


def matrix_to_axis_angle(A):
    """TRANS9: rotation matrix -> (unit axis h, rho in radians)."""
    A = np.asarray(A, dtype=float)
    c = (A[0, 0] + A[1, 1] + A[2, 2] - 1.0) / 2.0
    c = np.clip(c, -1.0, 1.0)
    rho = np.arccos(c)
    if rho == 0.0:
        return np.zeros(3), 0.0
    fact = 2.0 * np.sin(rho)
    h = np.array([A[2, 1] - A[1, 2], A[0, 2] - A[2, 0], A[1, 0] - A[0, 1]]) / fact
    h = h / np.linalg.norm(h)
    return h, rho


def rotation_matrix(lat_deg, lon_deg, rho_deg):
    """TRANS7: (axis lat, lon, rho in degrees) -> rotation matrix."""
    return quat_to_matrix(rotation_to_quat(lat_deg, lon_deg, rho_deg))


def skew(t):
    """
    Antisymmetric cross-product matrix M(t): M(t) @ v = t x v.
        M(t) = [[0, -t3, t2], [t3, 0, -t1], [-t2, t1, 0]]
    (Chang M(t), Kirkwood M(h)/M(eta), O'Neill T of eq. A5.)
    """
    t = np.asarray(t, dtype=float)
    return np.array([[0.0, -t[2], t[1]],
                     [t[2], 0.0, -t[0]],
                     [-t[1], t[0], 0.0]])


def M_3x4(q):
    """
    3x4 matrix M(q) of Hanna and Chang (1990) / conreg.f (lines 133-144), with
    q = (a0, a1, a2, a3):
        [[-a1,  a0,  a3, -a2],
         [-a2, -a3,  a0,  a1],
         [-a3,  a2, -a1,  a0]]
    It satisfies: pure quaternion part of (qhat_conj * q) ~ (1/2) M(qhat) q, and the
    construction Q = 4 M(qhat)^T H11.2 M(qhat).
    """
    q = np.asarray(q, dtype=float)
    a0, a1, a2, a3 = q
    return np.array([[-a1, a0, a3, -a2],
                     [-a2, -a3, a0, a1],
                     [-a3, a2, -a1, a0]])


def position_matrix(points_xyz):
    """
    Position matrix P = (1/n) sum_i p_i p_i^T of a set of unit points
    (O'Neill et al. 2005 eq. A2; = matrix X_i / Sigma_i without Hellinger III.4 weights).
    points_xyz: array (n, 3) of unit vectors.
    """
    P = np.asarray(points_xyz, dtype=float)
    n = len(P)
    return (P.T @ P) / n


def great_circle_fit(P):
    """
    Great-circle fit to a position matrix P (O'Neill A1).
    Returns dict:
        eigenvalues : (3,) ascending [L1<=L2<=L3]
        normal      : eigenvector of the SMALLEST eigenvalue = normal to the great circle (S1)
        along       : eigenvector of the intermediate one (S2, on the great circle)
        center      : eigenvector of the LARGEST eigenvalue = center/mean of the data (S3)
    In the limit of data exactly on a great circle, L1 -> 0.
    """
    P = np.asarray(P, dtype=float)
    evals, evecs = np.linalg.eigh(P)  # ascending
    return dict(eigenvalues=evals,
                normal=evecs[:, 0], along=evecs[:, 1], center=evecs[:, 2])


def rodrigues(u):
    """
    F(u) of O'Neill et al. (2005) eq. A7 (Rodrigues formula):
        F(u) = I + (sin|u|/|u|) T + ((1-cos|u|)/|u|^2) T^2,   T = skew(u).
    Equivalent to quat_to_matrix(exp_quat(u)); included for cross-validation.
    """
    u = np.asarray(u, dtype=float)
    n = np.linalg.norm(u)
    if n <= 0.0:
        return np.eye(3)
    T = skew(u)
    return np.eye(3) + (np.sin(n) / n) * T + ((1.0 - np.cos(n)) / n ** 2) * (T @ T)
