"""
Critical points of the chi-square and F distributions.

They replace the in-house Fortran routines (XIDCH, XIDF in numlib.f) with scipy,
which is more accurate. XIDCH/XIDF solved the inverse by bisection on the incomplete
gamma/beta; scipy.stats.*.ppf does the same robustly.
"""
from scipy import stats


def chi2_crit(prob, df):
    """chi2 critical point: value x such that P(chi2_df <= x) = prob. (XIDCH)"""
    return stats.chi2.ppf(prob, df)


def f_crit(prob, df1, df2):
    """F critical point: value x such that P(F_{df1,df2} <= x) = prob. (XIDF)"""
    return stats.f.ppf(prob, df1, df2)


def chi2_sf(x, df):
    """Upper tail of chi2: P(chi2_df > x). For p-values (pval.f)."""
    return stats.chi2.sf(x, df)


def f_sf(x, df1, df2):
    """Upper tail of F: P(F_{df1,df2} > x)."""
    return stats.f.sf(x, df1, df2)
