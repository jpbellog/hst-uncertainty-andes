"""
Inversion and combination of rotations with covariance propagation.

Port of invrot.f (Kirkwood eq. 20) and addrot.f (Kirkwood eq. 19, equal kappas).

Fortran convention (verified):
- The Fortran subroutine `mul(a,b,c)` computes c = a*b (standard product), even though
  its comment "c=b.a" is mislabelled.
- In addrot.f the product of rotations is C = B*A (second file * first), and the
  matrix transposed in the covariance is A (= A2 of the paper, the right factor).
"""
import numpy as np
from . import geometry as geo


def invert_rotation(A, cov):
    """
    Inverse of a rotation with its covariance (invrot.f; Kirkwood eq. 20).
        A_inv = A^T ;  cov_inv = A cov A^T ;  the angle is negated.
    Parameters: A (3x3), cov (3x3, covariance of the rotation h).
    Returns (A_inv, cov_inv).
    """
    A = np.asarray(A, dtype=float)
    cov = np.asarray(cov, dtype=float)
    A_inv = A.T
    cov_inv = A @ cov @ A.T
    return A_inv, cov_inv


def combine_rotations_equal_kappa(A1, cov1, kap1, df1, A2, cov2, kap2, df2):
    """
    Combines two independent rotations A3 = A1 * A2 assuming equal kappas
    (addrot.f; Kirkwood eq. 19 + case 1 of kappa combination).

    Covariance (eq. 19):  cov3 = A2^T cov1 A2 + cov2
    (with cov1, cov2 on the sigma^2 scale; see the kappa note below).
    Combined kappa (case 1):
        kap3 = (df1 + df2) / (df1/kap1 + df2/kap2),   df3 = df1 + df2.

    Scale NOTE: in the Fortran workflow the covariance stored by hellinger1 is
    H11.2^{-1} WITHOUT division by kappa; to combine, it must be brought to the sigma^2
    scale by dividing by kappa. Here cov1, cov2 are assumed to be ALREADY on the sigma^2
    scale (cov = H11.2^{-1} / kappa). Rescale outside if they are passed unscaled.

    Returns dict: A3, cov3, kap3, df3, H112_3 (=inv(cov3)).
    """
    A1 = np.asarray(A1, dtype=float)
    A2 = np.asarray(A2, dtype=float)
    cov1 = np.asarray(cov1, dtype=float)
    cov2 = np.asarray(cov2, dtype=float)
    cov3 = A2.T @ cov1 @ A2 + cov2
    A3 = A1 @ A2
    kap3 = (df1 + df2) / (df1 / kap1 + df2 / kap2)
    df3 = df1 + df2
    H112_3 = np.linalg.inv(cov3)
    return dict(A3=A3, cov3=cov3, kap3=kap3, df3=df3, H112_3=H112_3)


def average_rotations(A_list, H112_list, kappa_list, iters=50, tol=1e-12):
    """
    Average of independent estimates of the SAME rotation, combined through the
    information matrix (core of royer2.f / averotkap.f; generalizes eq. 19).

        V = sum_k kappa_k * H112_k          (combined information matrix)
        cov_comb = V^{-1}                   (combined covariance, sigma^2 scale)
    The combined rotation is the Frechet mean in SO(3) weighted by the information,
    found by iterating  delta = V^{-1} sum_k (kappa_k H112_k) log(A_mean^T A_k)  to 0.

    Parameters
    ----------
    A_list     : list of 3x3 rotation matrices.
    H112_list  : list of H11.2 matrices (inverse covariance, without kappa).
    kappa_list : list of estimated kappas.

    Returns dict: A (mean), cov (V^{-1}), H112 (=V), df_note.

    NOTE: the royer2 kappa correction (capa/capc/cfact terms) is not included; it
    comes from Chang et al. (1990) and requires the projection matrices (pending,
    see PORT_NOTES.md). Here the covariance and the mean rotation are combined, which is
    the verifiable part.
    """
    A_list = [np.asarray(A, float) for A in A_list]
    P = [k * np.asarray(H, float) for k, H in zip(kappa_list, H112_list)]  # info per fit
    V = sum(P)
    Vinv = np.linalg.inv(V)
    A_mean = A_list[0].copy()
    for _ in range(iters):
        rhs = np.zeros(3)
        for Ak, Pk in zip(A_list, P):
            # rotation vector of A_mean^T A_k (small near convergence)
            axis, ang = geo.matrix_to_axis_angle(A_mean.T @ Ak)
            rhs += Pk @ (axis * ang)
        delta = Vinv @ rhs
        A_mean = A_mean @ geo.quat_to_matrix(geo.exp_quat(delta))
        if np.linalg.norm(delta) < tol:
            break
    return dict(A=A_mean, cov=Vinv, H112=V)


def welch_df(A_john, B_john, C_john):
    """
    Welch-James-Johansen adjustment of the degrees of freedom for the critical value of
    a confidence region with unequal kappas (teclib.f `critcalc`; Kirkwood
    section "Unequal kappas").

    Given the Johansen quantities A, B, C (traces of the projections of the restricted
    and full models, subroutine `abc` of teclib), the distribution of Q is
    approximated by c*F(3, df) by matching mean and variance:
        E[Q]  = 3 + 2A + 2B
        Var[Q]= 6 + 14A + 2B + 2C
        df    = (12*Var + 2*E^2) / (3*Var - 2*E^2)
        c     = E * (1 - 2/df)
    The critical value at level alpha is c * F_alpha(3, df).

    For independent data (add2 / averotkap), B = 0.

    Returns dict: E, Var, df, cfact.
    """
    E = 3.0 + 2.0 * A_john + 2.0 * B_john
    Var = 6.0 + 14.0 * A_john + 2.0 * B_john + 2.0 * C_john
    df = (12.0 * Var + 2.0 * E ** 2) / (3.0 * Var - 2.0 * E ** 2)
    cfact = E * (1.0 - 2.0 / df)
    return dict(E=E, Var=Var, df=df, cfact=cfact)


def welch_critical_value(A_john, B_john, C_john, alpha=0.95):
    """
    Critical value of the confidence region with unequal kappas, FULL DATA
    (Kirkwood and Chang 1998, Sec. 4 + Appendix B; add2/add3/averotkap via teclib critcalc).
    Given the Johansen quantities A,B,C (from the projections P0/P1), returns the
    critical value = c * F_alpha(3, df) with (c, df) from welch_df.
    Returns dict: E, Var, df, cfact, crit.
    """
    from . import stats_dist as sd
    w = welch_df(A_john, B_john, C_john)
    crit = w["cfact"] * sd.f_crit(alpha, 3, w["df"])
    w["crit"] = float(crit)
    return w


def conservative_critical_value(nf_pairs, alpha=0.95):
    """
    CONSERVATIVE critical value of the confidence region with unequal kappas, using
    ONLY summary statistics (Kirkwood and Chang 1998, Theorem 5.2.2, pp. 94-95;
    addplus.f / averotplus.f). Raw data are not required.

    Parameters
    ----------
    nf_pairs : list of (n_i, f_i) per separate estimate, with n_i = number of points and
               f_i = degrees of freedom (= N_i - 2 s_i - 3 for a two-plate fit).
    alpha    : confidence level (0.95 by default; the bound is valid for alpha <= 0.95).

    Formulas (Theorem 5.2.2):
        n = sum n_i,   r = sum (n_i - f_i),   f_min = min f_i,   f_max = max f_i
        z = (r-3)^2/n + 6 - r
        a_min = min(z/f_min, z/f_max)
        if a_min > -3/2:  K = (6 + 54/f_min)/(3 + 2 a_min)^2,  v = (12K+2)/(3K-2)
        otherwise:        v = 4
        crit = (3 + 6/f_min) (1 - 2/v) F_alpha(3, v)

    Returns dict: n, r, f_min, f_max, z, a_min, K, v, factor, crit.
    """
    from . import stats_dist as sd
    n = sum(ni for ni, _ in nf_pairs)
    r = sum(ni - fi for ni, fi in nf_pairs)
    fs = [fi for _, fi in nf_pairs]
    f_min, f_max = min(fs), max(fs)
    z = (r - 3) ** 2 / n + 6 - r
    a_min = min(z / f_min, z / f_max)
    if a_min > -1.5:
        K = (6.0 + 54.0 / f_min) / (3.0 + 2.0 * a_min) ** 2
        v = (12.0 * K + 2.0) / (3.0 * K - 2.0)
    else:
        K = None
        v = 4.0
    factor = (3.0 + 6.0 / f_min) * (1.0 - 2.0 / v)
    crit = factor * sd.f_crit(alpha, 3, v)
    return dict(n=n, r=r, f_min=f_min, f_max=f_max, z=z, a_min=a_min,
                K=K, v=v, factor=factor, crit=float(crit))
