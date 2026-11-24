"""
Extension of O'Neill, Muller and Steinberger (2005) for the uncertainty of
positions on hotspot tracks.

O'Neill et al. adapt the Hellinger/Chang method (implemented by the other
modules) to the hotspot problem. The parts SPECIFIC to O'Neill are:

1. Conversion of age uncertainty to spatial uncertainty (eq. 2 of the paper):
   the temporal uncertainty of a dated point, combined with the angular velocity
   of the plate (and its uncertainty), gives a SPATIAL uncertainty along the
   track (along-track). Together with the perpendicular uncertainty (across-track,
   estimated from the full width at half maximum of the gravity profile, sec. 2.2), it
   forms an error ellipse elongated along the track.

2. Modified Hellinger criterion (eq. 1): includes the hotspot positions
   as data (fit to their present/computed position) in addition to the track geometry.

3. Reuse of the Chang machinery (Appendix A): covariance H11.2 (via the
   Hellinger fit, A5), confidence region (A8/A10), and ellipse of a reconstructed point (A6,
   program PTELLIPSE = ellipse.f). O'Neill A9: kappa_hat = (N - 2s - 3)/r.

This module implements (1), an error-ellipse utility, the fit (2) reusing
hellinger.fit_hellinger, and the uncertainty of a reconstructed point (3).

NOTE on eq. (2): the O'Neill PDF is scanned and eq. (2) reads
ambiguously ("Dd = Dt.omega + t.Domega" vs "Dd = Dt.omega + Dt.Domega"). Here the
PHYSICALLY standard first-order propagation of d = omega * t is used:
    Delta_d = omega * Delta_t + t * Delta_omega,
which is the natural interpretation (the term t*Delta_omega usually dominates). To
reproduce the exact formula of the paper, adjust `along_track_position_sigma`.
"""
import numpy as np
from . import geometry as geo
from . import hellinger as hel
from . import region as reg


def along_track_position_sigma(age_myr, age_err_myr, omega_deg_per_myr,
                               omega_err_deg_per_myr=0.0):
    """
    Along-track spatial uncertainty (degrees) derived from the age
    (O'Neill eq. 2, propagation of d = omega * t):
        Delta_d = omega * Delta_t + t * Delta_omega.

    Parameters
    ----------
    age_myr              : age of the point (Myr).
    age_err_myr          : 1-sigma uncertainty of the age (Myr).
    omega_deg_per_myr    : angular velocity of the plate (degrees/Myr).
    omega_err_deg_per_myr: uncertainty of the angular velocity (degrees/Myr).

    Returns the along-track uncertainty in DEGREES of arc.
    """
    return (omega_deg_per_myr * age_err_myr
            + age_myr * omega_err_deg_per_myr)


def point_error_ellipse(lat, lon, along_track_sigma_deg, across_track_sigma_deg,
                        track_azimuth_deg):
    """
    Error ellipse of a track point (O'Neill sec. 2.2, Fig. 3): elongated
    along the track (along-track >> across-track).

    Parameters
    ----------
    lat, lon              : position of the point (degrees).
    along_track_sigma_deg : 1-sigma semi-axis along the track (degrees), typically
                            from along_track_position_sigma().
    across_track_sigma_deg: perpendicular 1-sigma semi-axis (degrees), from the full
                            width at half maximum of the profile.
    track_azimuth_deg     : azimuth of the track at the point (degrees from north).

    Returns dict with the ellipse in the local tangent plane:
        center_latlon, sigma_along, sigma_across, azimuth_deg,
        cov2d : 2x2 covariance in local (East, North) axes (degrees^2).
    """
    az = np.deg2rad(track_azimuth_deg)
    # along-track direction in local (East, North) axes
    t_hat = np.array([np.sin(az), np.cos(az)])
    n_hat = np.array([np.cos(az), -np.sin(az)])
    cov2d = (along_track_sigma_deg ** 2) * np.outer(t_hat, t_hat) \
        + (across_track_sigma_deg ** 2) * np.outer(n_hat, n_hat)
    return dict(center_latlon=(lat, lon),
                sigma_along=along_track_sigma_deg,
                sigma_across=across_track_sigma_deg,
                azimuth_deg=track_azimuth_deg, cov2d=cov2d)


def fit_hotspot_rotation(track_points, hotspot_points, init_latlon_rho, **kwargs):
    """
    Fit of an absolute rotation from hotspot track data, reusing the
    Hellinger criterion (O'Neill Appendix A5: "this is done by HELLINGER1").

    The problem is set up as in Hellinger: the "moving side" (iside=1) is the
    track points (track geometry, reconstructed to the hotspot position),
    and the "fixed side" (iside=2) is the points that mark the hotspot position /
    age isochrons. Each section i groups a track segment with its associated hotspot.

    Parameters
    ----------
    track_points   : list of (isect, lat, lon, sd_km)   -> moving side (track)
    hotspot_points : list of (isect, lat, lon, sd_km)   -> fixed side (hotspot/age)
    init_latlon_rho: initial guess (lat, lon, rho) in degrees.

    Returns the same dict as hellinger.fit_hellinger (includes H112, cov_h, kappa_hat).

    NOTE: O'Neill eq. (1) distinguishes hotspot position terms on both
    plates; this function covers the general case of two conjugate sets (track vs
    hotspot). For the "single hotspot" term (Sum_k h_k^T h_k), represent it
    as a section with data only on the fixed side.
    """
    crossings = []
    for isect, lat, lon, sd in track_points:
        crossings.append((1, isect, lat, lon, sd))
    for isect, lat, lon, sd in hotspot_points:
        crossings.append((2, isect, lat, lon, sd))
    return hel.fit_hellinger(crossings, init_latlon_rho, **kwargs)


def reconstructed_hotspot_ellipse(A, cov_u, point_latlon, alpha=0.95,
                                  kappa_hat=1.0, df=None, known_kappa=True,
                                  simultaneous=True):
    """
    Uncertainty of the reconstructed position of a hotspot track point (O'Neill
    A6; Kirkwood eq. 22; ellipse.f / PTELLIPSE). Given the fitted rotation A and the
    covariance of the rotation cov_u, returns the error ellipse of the point A*w.

    Parameters
    ----------
    A            : 3x3 rotation matrix (from fit_hotspot_rotation).
    cov_u        : covariance of the rotation (= H112^{-1} on the sigma^2 scale; see the
                   kappa note in region.critical_value).
    point_latlon : (lat, lon) of the point to reconstruct (degrees).
    alpha        : confidence level (e.g. 0.95).
    kappa_hat, df, known_kappa : see region.critical_value.
    simultaneous : True uses 3F(3,df)/chi2(3) (Scheffe); False uses 2F(2,df)/chi2(2)
                   (individual ellipse of eq. 22).

    Returns the dict of region.reconstructed_point_ellipse (semi-axes in radians,
    reconstructed point, directions).
    """
    from . import stats_dist as sd
    ndim = 3 if simultaneous else 2
    if known_kappa:
        xchi = sd.chi2_crit(alpha, ndim) / kappa_hat
    else:
        if df is None:
            raise ValueError("df is required when known_kappa=False")
        xchi = ndim * sd.f_crit(alpha, ndim, df) / kappa_hat
    return reg.reconstructed_point_ellipse(A, cov_u, point_latlon, xchi)
