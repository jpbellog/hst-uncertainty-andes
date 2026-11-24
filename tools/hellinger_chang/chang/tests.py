"""
Statistical tests on rotations.

Port of pval.f (identity test) and of the F test for differences (Kirkwood eq. 24-27;
O'Neill A7). Replaces the in-house inverses with scipy.
"""
import numpy as np
from . import geometry as geo
from . import stats_dist as sd


def identity_test(A, H112, kappa_hat=1.0, df=None, known_kappa=True):
    """
    Tests whether the rotation A is the identity (pval.f).

    Statistic: T = kappa_hat * h^T H11.2 h, with h = log(A) (rotation vector).
    Under H0 (A = identity):
        known kappa     : T ~ chi2(3);  p = P(chi2(3) > T).
        estimated kappa : T/3 ~ F(3, df); p = P(F(3,df) > T/3).

    Returns dict: statistic, p_value, reject_at_95 (True if p < 0.05, i.e. A != I).
    """
    A = np.asarray(A, dtype=float)
    H112 = np.asarray(H112, dtype=float)
    h, _ = geo.matrix_to_axis_angle(A)
    _, rho = geo.matrix_to_axis_angle(A)
    hvec = h * rho  # rotation vector (equivalent to log)
    T = kappa_hat * float(hvec @ H112 @ hvec)
    if known_kappa:
        p = sd.chi2_sf(T, 3)
    else:
        if df is None:
            raise ValueError("df is required when known_kappa=False")
        p = sd.f_sf(T / 3.0, 3, df)
    return dict(statistic=T, p_value=float(p), reject_at_95=bool(p < 0.05))


def kappa_equality_test(kappa1, df1, kappa2, df2, sig=0.10):
    """
    Pre-test of kappa equality before combining (O'Neill A7; addrot.f).
    F = max(k1,k2)/min(k1,k2) compared with F_{1-sig/2}. Returns dict with F, p_value,
    equal (True if equality is not rejected at level sig).
    """
    if kappa1 >= kappa2:
        F = kappa1 / kappa2
        d1, d2 = df1, df2
    else:
        F = kappa2 / kappa1
        d1, d2 = df2, df1
    p = 2.0 * sd.f_sf(F, d1, d2)  # two-tailed
    return dict(F=float(F), p_value=float(min(p, 1.0)), equal=bool(p > sig))


def difference_test(rA1, dfA1, rA2, dfA2, rA3, dfA3, alpha=0.95):
    """
    F test of whether the combination A3 = A1 A2 differs from the identity / whether two
    plates move together (Kirkwood eq. 25; O'Neill A7):

        F = [ (rA3 - (rA1 + rA2)) / (dfA3 - (dfA1 + dfA2)) ]
            * [ (dfA1 + dfA2) / (rA1 + rA2) ]

    It is compared with F_alpha(3, dfA3). If F_alpha(3, dfA3) < F, the rotations differ
    significantly at level alpha.

    Returns dict: F, F_crit, significant (True if they differ), df_num.
    """
    dnum = dfA3 - (dfA1 + dfA2)
    F = ((rA3 - (rA1 + rA2)) / dnum) * ((dfA1 + dfA2) / (rA1 + rA2))
    F_crit = sd.f_crit(alpha, 3, dfA3)
    return dict(F=float(F), F_crit=float(F_crit),
                significant=bool(F_crit < F), df_num=dnum)
