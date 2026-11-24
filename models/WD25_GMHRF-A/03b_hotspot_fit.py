# -*- coding: utf-8 -*-
"""STEP 03b -- Hotspot position as an unknown. Replaces the "this study"
anchors of plates 13/14/15 (one rotation per dated edifice, see
models/WD25_WK08-A/refit_hotspots.py) by a parametric model of the hotspot in the
frame (plate 002), fitted by maximum likelihood to the geochronology with
additive geological dispersion.

MODEL
For a sample of age t, the predicted present position of the track is
    x(t) = B(t) . h(t),   B(t) = R_201(0)^T R_rel(0) R_rel(t)^T R_002(t)
(same chain as motion_path with plates 13/14/15 at identity relative to
002), with h(t) the position of the hotspot in the frame:
  - "fijo" (fixed):   h(t) = h0                (lat0, lon0; 2 parameters)
  - "deriva" (drift): h0 displaced along the great circle with constant
              velocity (ve, vn) in km/Myr measured at h0 (4 parameters). After
              the oldest sample of the fit h stays at h(t_max): the drift is
              not extrapolated.
The residual obs - x(t) is decomposed in the tangent plane into along (track
direction) and across. Per-sample variances (km^2):
    along:  (omega (km/Myr))^2 (sigma_t^2 + w_t^2) + sigma_pos^2
    across: sigma_pos^2 + w_x^2
sigma_t = 1 sigma analytical error, w_t (Myr) and w_x (km) additive geological
dispersion (duration of volcanism of each edifice, hotspot-eruption lag,
off-axis samples), estimated together with the parameters. Precedent:
Gaastra et al. (2022), who also use an accumulated uncertainty instead of the
analytical one.

SELECTION AND OUTLIERS
- At each iteration "fijo" and "deriva" are fitted on the same set and the one
  with the lower BIC is chosen (the AIC is also reported).
- Outlier: Mahalanobis distance d^2 = z_along^2 + z_across^2 with the
  dispersion already estimated; the worst sample is excluded if
  d^2 > chi2_2(0.9973) = 11.83 (the 2D equivalent of |z| > 3) and the fit is
  repeated, until no sample exceeds it. The prior use_fit selection of each CSV
  is kept.
- Present hotspot position: by default it does NOT enter the fit: the seeds are
  not measured active volcanism, they are the present positions of the model of
  Bello-Gonzalez et al. (2018, 2026). With PAPER_HS_P0=1 (sensitivity) the
  seed of each track enters as a pseudo-observation of age 0, without analytical error. Its age is not measured, so
  its dispersion is only spatial and isotropic, w_x along and across
  (PAPER_HS_P0_MODE=iso, base); PAPER_HS_P0_MODE=edificio treats it as a
  dated edifice (along omega^2 w_t^2, across w_x^2). PAPER_HS_P0_KM
  adds its own position error (km, default 0); the outputs carry
  "_conP0". It is never rejected as an outlier nor scaled by calibration.
- San Felix (2 samples): only "fijo", with w_t and w_x taken from JFR (same
  plate) and no rejection.
- Sensitivity PAPER_HS_LIK=t: Student t likelihood (nu = 4) without rejection.

PARAMETER UNCERTAINTY
- Statistical: Gauss-Newton with the weights of the solution, Cov = (J^T J)^-1.
- Calibration (common systematics): each systematic k is an age scale factor
  (1 + e_k), e_k ~ N(0, s_k). The sensitivity g_k = d theta / d e_k is stored
  (numerical refit with w and the sample set fixed). Step 04 samples e_k ONCE
  per realisation and shares it among tracks: theta = theta_hat + L z + sum_k g_k e_k.
- PAPER_HS_W=pred (sensitivity): adds w to the hotspot position as prediction
  dispersion (locus of an individual edifice) instead of using it only in the
  weights (confidence interval, base case).

MODELLED DRIFT (GMHRF-B frame)
PAPER_HS_PRIOR=gmhrf adds the hotspot motion predicted by the GMHRF mantle
flow model (Data Set S1 of Doubrovine et al. 2012, read with
../rotation_models/frames/gmhrf/build_rot_gmhrf.py):
    h(t) = S D(t) S^T h_lin(t),
D(t) = minimal rotation D1(0) -> D1(t) (frozen after the last age of the
D1), S = minimal rotation D1(0) -> h0, so the D1 displacement is applied at
the fitted position. The D1 drift applies from 0 Ma. "fijo" = D1 only;
"deriva" = D1 plus a residual linear drift. The D1 uncertainty is not
propagated (Data Set S1 does not provide it). With PAPER_HS_PRIOR=none
(default) the hotspot is fixed or has a linear drift in the frame, as in
WK08-A, GMHRF-A' and GWG20.

AGES
Easter: age_E and err2s_E from ages_easter_ray2012_recal.csv: plateaus
recomputed from the ERDA steps with additive dispersion, without sqrt(MSWD) or
J, and inverse isochron of DRFT-115-2 with York et al. (2004)
(analysis/ray2012_erda_ages.py), in FCT-3 = FCs 28.201 /
Min 2000. PAPER_HS_EASTER=A uses age_A (published ages, recalibrated) and B
uses age_B (factor of Daze et al. 2003). JFR: Lara et al. (2018), already in
FCs 28.201. San Felix: age_recal from ages_taltal_recal.csv.

Outputs:
  - hotspot_fit.json : parameters, covariances, g_k, w, BIC/AIC, excluded samples.
  - hotspot_residuals_{ridge}.csv : residuals and z per sample (includes excluded).
  - rot/<base>_hsb.rot : the .rot with plates 13/14/15 replaced by the
    fitted model (nodes every 1 Myr up to t_max with drift), for steps
    04, 07, 08, 10 and 11 and for GPlates (PAPER_ROT_FILE).

Usage:  py 03b_hotspot_fit.py
"""
import os
import sys
import json
import importlib
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
eng = importlib.import_module("02_rot_engine")

KM_PER_DEG = 111.19492664455873
R_KM = KM_PER_DEG * 180.0 / np.pi
REL, ANCHOR, FRAME = 902003, 201, 2   # f_PAC 0.44; identical to 0.5/0.36 in 0-50 Ma
DELTA = 0.25                          # Myr, track derivative
D2_OUT = float(chi2.ppf(0.9973, 2))
NU_T = 4.0
LIK = os.environ.get("PAPER_HS_LIK", "gauss")
W_MODE = os.environ.get("PAPER_HS_W", "conf")
USE_P0 = os.environ.get("PAPER_HS_P0", "0") == "1"   # no seed by default
P0_KM = float(os.environ.get("PAPER_HS_P0_KM", "0"))
P0_MODE = os.environ.get("PAPER_HS_P0_MODE", "iso")
P0_NAME = "HOTSPOT_0MA"
EASTER = os.environ.get("PAPER_HS_EASTER", "E")
EASTER_COL = {"E": ("age_E", "err2s_E"), "A": "age_A", "B": "age_B"}[EASTER]
EASTER_KEY = f"age_{EASTER}"           # in the output file names
PRIOR = os.environ.get("PAPER_HS_PRIOR", "none")
if PRIOR not in ("none", "gmhrf"):
    raise ValueError(f"PAPER_HS_PRIOR={PRIOR!r}: use none or gmhrf")
_D1 = {}

# seeds: present hotspot position in steps 04, 07, 08, 10 and 11. The new
# nodes take the seed to h(t), so those steps do not change.
SEEDS = {13: (-33.750097, -81.905212), 14: (-26.3138, -106.599),
         15: (-26.28791, -81.42893), 16: (-37.3764, -110.888)}
TRACKS = [("Easter", 14, "ages_easter_ray2012_recal.csv", EASTER_COL),
          ("JFR", 13, "ages_jfr_lara2018.csv", "age_ma"),
          ("TT", 15, "ages_taltal_recal.csv", "age_recal"),
          ("Foundation", 16, os.environ.get("PAPER_HS_FOUNDATION",
                                            "ages_foundation_oconnor2004_v3.csv"), "age_recal")]
# PAPER_HS_FOUNDATION: file of Foundation ages in data/. Default _v3
# (2001 recalculation without Seamounts 8, 9 and 13a). Others:
# ages_foundation_oconnor2004_v1.csv (O'Connor et al. 2001, recalculation of 1998) or _v2.csv (with
# those of 2001 and 2002 from Table 11.1 of O'Connor et al. 2004); see foundation_oconnor2004.py.
# Foundation / Iquique ridge (fourth chain): the
# hotspot builds the Foundation chain on the Pacific plate (901) since the
# age of its oldest sample (Seamount 1a, O'Connor et al. 1998) and the Iquique
# ridge on Farallon/Nazca (REL) before that age (Bello-Gonzalez et al.
# 2018). The fit uses only the dated chain, with plate 901; t_switch = t_max
# of the track. Seed: present position in the model of Bello-Gonzalez et al. (2018,
# Foundation_mod_1.gpml), the same source as the other three. No modelled drift
# in the GMHRF (Data Set S1 of Doubrovine et al. 2012 does not include Foundation).
REL_YOUNG = {"Foundation": 901}
# t_switch: age at which the hotspot crosses the Challenger FZ from north to south
# (conjugate of the Resolution FZ), criterion of Bello-Gonzalez et al. (2018, p. 9). It is the
# crossing, in present Nazca coordinates (REL), of the extended track with the FZ trace
# (data/fz_challenger_gsfml.csv); if the crossing is younger than the oldest sample of the
# Pacific chain (t_max), it contradicts the dated chain and t_max is used.
FZ_SWITCH = {"Foundation": "fz_challenger_gsfml.csv"}

# Calibration systematics, relative scale factor 1 sigma
# (recalibration with the primary sources). They affect only 40Ar/39Ar ages
# (not the K-Ar age of San Ambrosio).
SYSTEMATICS = {
    # FCs 28.201 +/- 0.046 Ma (2 sigma), Kuiper et al. (2008): common to all
    "FCs": (0.046 / 2.0) / 28.201,
    # FCT-3/FCs, Daze et al. (2003): 28.13 +/- 0.47 Ma (2 sigma); Ray 2012 only
    "FCT3": (0.47 / 2.0) / 28.13,
    # ACs/FCs = 0.04229 +/- 0.00006, Renne et al. (1998); Haase 2000 only.
    # Sigma level of the +/-0.00006 [TO BE VERIFIED]; taken as 1 sigma.
    "ACs": 0.00006 / 0.04229,
    # J factor: common to all samples of each source, as if they shared one
    # irradiation (bound: full correlation within the source).
    # Ray et al. (2012, Methods): J with 0.3-0.5 % (1 sigma); 0.5 % is used.
    # Lara et al. (2018) and Haase et al. (2000) include J in their errors without
    # giving its value: 0.5 % is our ASSUMPTION, the same as for Ray. The published
    # per-sample errors already contain J, so this counts it twice
    # (conservative). PAPER_HS_J=0 removes it.
    "J_Ray": 0.005,
    "J_Lara": 0.005,
    "J_Haase": 0.005,
    # TCs/FCs = 1.0112 +/- 0.0010, Renne et al. (1998, p. 125); sigma level not
    # stated, taken as 1 sigma (as for ACs). O'Connor et al. (1998) only.
    "TCs": 0.0010 / 1.0112,
    # O'Connor et al. (1998, p. 46): J with 0.1-0.3 % (s.e.m.); 0.3 % is used.
    "J_OConnor": 0.003,
}
SYS_BY_RIDGE = {"Easter": ("FCs", "FCT3", "J_Ray"), "JFR": ("FCs", "J_Lara"),
                "TT": ("FCs", "ACs", "J_Haase"),
                "Foundation": ("FCs", "TCs", "J_OConnor")}
# Easter with age_E (ERDA data): one J per irradiation instead of
# J_Ray. The age_E ages do not include J; the J of each ArArCALC file (J_pct1s, %1s)
# is common to all steps of the sample and is taken as common to all
# samples of the same batch (full correlation within the batch, independence
# between batches); sigma = mean J of the batch. Key "J_<lote>".
J_LOTE = {}
if EASTER == "E":
    _e = pd.read_csv(os.path.join(HERE, "data", "ages_easter_ray2012_recal.csv"), comment="#")
    _e = _e[_e.use_fit == 1]
    J_LOTE = {f"J_{k}": v / 100.0 for k, v in _e.groupby("irradiation").J_pct1s.mean().items()}
    SYSTEMATICS.pop("J_Ray")
    SYSTEMATICS.update(J_LOTE)
    SYS_BY_RIDGE["Easter"] = ("FCs", "FCT3") + tuple(sorted(J_LOTE))
# Foundation with columns irradiation and J_pct1s (foundation_oconnor2004.py): one J per batch, as
# for Easter, instead of J_OConnor.
_fn_f = [t[2] for t in TRACKS if t[0] == "Foundation"][0]
_f = pd.read_csv(os.path.join(HERE, "data", _fn_f), comment="#")
if {"irradiation", "J_pct1s"} <= set(_f.columns):
    _f = _f[_f.use_fit == 1]
    _jf = {f"J_{k}": v / 100.0 for k, v in _f.groupby("irradiation").J_pct1s.mean().items()}
    SYSTEMATICS.pop("J_OConnor")
    SYSTEMATICS.update(_jf)
    J_LOTE.update(_jf)
    SYS_BY_RIDGE["Foundation"] = ("FCs", "TCs") + tuple(sorted(_jf))
if os.environ.get("PAPER_HS_J", "1") == "0":
    for _k in ("J_Ray", "J_Lara", "J_Haase", "J_OConnor", *J_LOTE):
        SYSTEMATICS.pop(_k, None)
    SYS_BY_RIDGE = {r: tuple(k for k in v if k in SYSTEMATICS) for r, v in SYS_BY_RIDGE.items()}
# Sensitivities: PAPER_HS_EXCLUDE="DRFT-115-2,DRFT-126-1" removes
# those samples from the fit (the YORK2 isochrons of Ray et al. 2012);
# PAPER_HS_OUT_TAG adds a suffix to hotspot_fit*.json and hotspot_residuals_*;
# PAPER_HS_ROT_OUT sets the path of the written .rot (default <base>_hsb<tag>.rot).
EXCLUDE = [x.strip() for x in os.environ.get("PAPER_HS_EXCLUDE", "").split(",") if x.strip()]
OUT_TAG = os.environ.get("PAPER_HS_OUT_TAG", "")
ROT_OUT = os.environ.get("PAPER_HS_ROT_OUT", "")
# PAPER_HS_FORCE="Easter=deriva" imposes the model of a track
# (fijo and deriva are still both fitted and both BIC are reported; the forced model
# also drives the outlier rejection). Gives Easter the same parametrisation
# in all models (in WD25 with GTS2020 the BIC favours fijo by 0.7).
FORCE = dict(x.strip().split("=") for x in os.environ.get("PAPER_HS_FORCE", "").split(",") if x.strip())
# PAPER_HS_LAG = systematic lag (Myr) between the passage of the
# plate over the hotspot and the dated eruption; it is added to all ages
# (ArAr and K-Ar) and to the 0 Ma pseudo-observation (the active volcanism passed
# over the hotspot LAG Myr ago). A lag that varies between edifices goes into
# w_t. Sensitivity only (use with PAPER_HS_OUT_TAG).
LAG = float(os.environ.get("PAPER_HS_LAG", "0"))
AR_TYPES = {"plateau", "isochron", "total_fusion", "weighted_mean"}


# ------------------------------------------------------------------ geometry
def enu(p):
    la, lo = eng.xyz_to_latlon(p)
    la, lo = np.radians(la), np.radians(lo)
    e = np.array([-np.sin(lo), np.cos(lo), 0.0])
    n = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    return e, n


def min_rot(u, v):
    """Matrix of the minimal rotation that takes u to v."""
    return eng.quat_to_matrix(eng.quat_from_pole(*G_D1().minimal_rotation(u, v)))


def G_D1():
    if "mod" not in _D1:
        sys.path.insert(0, os.path.join(HERE, "..", "..", "rotation_models", "frames", "gmhrf"))
        _D1["mod"] = importlib.import_module("build_rot_gmhrf")
    return _D1["mod"]


def has_d1(pid):
    """True if plate pid carries the modelled GMHRF drift."""
    return PRIOR != "none" and pid in G_D1().D1_NAME


def prior_D(pid, t):
    """D(t) of the D1 for plate pid: (3, 3) or (n, 3, 3); None without modelled
    drift."""
    if not has_d1(pid):
        return None
    G = G_D1()
    if pid not in _D1:
        d1 = G.read_d1(G.D1_NAME[pid])
        _D1[pid] = (d1, G.d1_at(d1, 0.0), float(d1.age.max()))
    d1, u0, _ = _D1[pid]
    t = np.asarray(t, float)
    if t.ndim == 0:
        return min_rot(u0, G.d1_at(d1, float(t)))
    return np.array([min_rot(u0, G.d1_at(d1, float(a))) for a in t])


def d1_end(pid):
    prior_D(pid, 0.0)
    return _D1[pid][2]


def hotspot_xyz(theta, t, t_max, pid=None, D=None):
    """Hotspot position (unit vector, frame 002) at age t (scalar or
    array of ages; returns (3,) or (n, 3)). D: matrices D(t) of the modelled
    drift, already evaluated (otherwise computed from pid)."""
    p0 = eng.latlon_axis(theta[0], theta[1])
    t = np.asarray(t, float)
    if len(theta) == 2:
        out = p0 if t.ndim == 0 else np.tile(p0, (t.size, 1))
    else:
        e, n = enu(p0)
        v = (theta[2] * e + theta[3] * n) / R_KM          # rad/Myr
        nv = np.linalg.norm(v)
        if nv < 1e-15:
            out = p0 if t.ndim == 0 else np.tile(p0, (t.size, 1))
        else:
            a = nv * np.minimum(t, t_max)
            out = np.cos(a)[..., None] * p0 + np.sin(a)[..., None] * (v / nv)
    if pid is None or not has_d1(pid):
        return out
    if D is None:
        D = prior_D(pid, t)
    S = min_rot(_D1[pid][1], p0)
    if t.ndim == 0:
        return S @ D @ S.T @ out
    return np.einsum("ij,njk,lk,nl->ni", S, D, S, out)


def b_matrix(M, t, rel=REL):
    """B(t): hotspot position in the frame -> present position on the relative
    plate, seen from 201 at 0 Ma."""
    R = lambda p, a: eng.quat_to_matrix(M.total_rotation(p, a, 0))
    return R(ANCHOR, 0.0).T @ R(rel, 0.0) @ R(rel, t).T @ R(FRAME, t)


class Track:
    def __init__(self, M, ridge, pid, df, age_col, t_max=None):
        self.ridge, self.pid = ridge, pid
        self.rel = REL_YOUNG.get(ridge, REL)
        self.df = df.reset_index(drop=True)
        self.M = M
        age_col, err_col = age_col if isinstance(age_col, tuple) else (age_col, None)
        self.age0 = self.df[age_col].to_numpy(float)
        if err_col:
            self.sig_t = self.df[err_col].to_numpy(float) / 2.0
        else:
            a_pub = self.df.age_ma.to_numpy(float)
            scale = np.divide(self.age0, a_pub, out=np.ones_like(a_pub), where=a_pub > 0)
            self.sig_t = self.df.err2s_ma.to_numpy(float) / 2.0 * scale
        self.obs = np.array([eng.latlon_axis(a, b) for a, b in
                             zip(self.df.lat, self.df.lon)])
        pe = (pd.to_numeric(self.df["pos_err_km"], errors="coerce").fillna(0.0)
              if "pos_err_km" in self.df else pd.Series(0.0, index=self.df.index))
        self.sig_pos = pe.to_numpy(float)
        self.is_ar = self.df.age_type.isin(AR_TYPES).to_numpy()
        self.pseudo = (self.df["sample"] == P0_NAME).to_numpy()
        if LAG:   # passage over the hotspot = measured age + lag; also the active
            self.age0 = self.age0 + LAG   # volcanism (pseudo-observation): passed LAG Myr ago
        self.t_max = float(self.age0.max()) if t_max is None else t_max
        self.set_ages(self.age0)

    def set_ages(self, ages):
        self.ages = np.asarray(ages, float)
        self.tm = np.maximum(self.ages - DELTA, 0.0)
        self.tp = self.ages + DELTA
        # B at (t-d, t, t+d), arrays (n, 3, 3)
        self.Bm = np.array([b_matrix(self.M, a, self.rel) for a in self.tm])
        self.B0 = np.array([b_matrix(self.M, a, self.rel) for a in self.ages])
        self.Bp = np.array([b_matrix(self.M, a, self.rel) for a in self.tp])
        # modelled drift (PAPER_HS_PRIOR) at the same ages
        self.D = {id(a): prior_D(self.pid, a) for a in (self.tm, self.ages, self.tp)}

    def scaled(self, sysname, eps):
        """Ages with the systematic sysname scaled by (1 + eps)."""
        m = self.is_ar
        if sysname in J_LOTE:                  # only the samples of that batch
            m = m & (self.df["irradiation"] == sysname[2:]).to_numpy()
        f = np.where(m, 1.0 + eps, 1.0)
        return self.age0 * f

    def residuals(self, theta, keep):
        """(r_along, r_across, omega) in km and km/Myr for the samples keep."""
        k = np.flatnonzero(keep)
        def mul(B, t):
            D = self.D[id(t)]
            h = hotspot_xyz(theta, t[k], self.t_max, self.pid,
                            None if D is None else D[k])
            return np.einsum("nij,nj->ni", B[k], h)
        x = mul(self.B0, self.ages)
        # track towards older ages: from x(t-d) to x(t+d)
        d = (mul(self.Bp, self.tp) - mul(self.Bm, self.tm)) / \
            (self.tp - self.tm)[k, None]
        d -= np.sum(d * x, 1)[:, None] * x
        w = np.linalg.norm(d, axis=1)
        a_hat = d / w[:, None]
        c_hat = np.cross(x, a_hat)
        r = self.obs[k] - x
        return (np.sum(r * a_hat, 1) * R_KM, np.sum(r * c_hat, 1) * R_KM, w * R_KM)

    def variances(self, om, keep, wt, wx):
        sp2 = self.sig_pos[keep] ** 2
        va = om ** 2 * (self.sig_t[keep] ** 2 + wt ** 2) + sp2
        if P0_MODE == "iso":
            # present position: its age is not measured, the dispersion is only spatial
            va = np.where(self.pseudo[keep], sp2 + wx ** 2, va)
        return (va, sp2 + wx ** 2)


# ------------------------------------------------------------------ fit
def nll_terms(r, V, lik):
    if lik == "t":
        return (gammaln(NU_T / 2) - gammaln((NU_T + 1) / 2)
                + 0.5 * np.log(NU_T * np.pi * V)
                + (NU_T + 1) / 2 * np.log1p(r ** 2 / (NU_T * V)))
    return 0.5 * (r ** 2 / V + np.log(2 * np.pi * V))


def fit(tr, keep, n_theta, lik=LIK, w_fixed=None, theta0=None):
    """Maximum likelihood. Returns a dict with theta, w_t, w_x, nll, cov."""
    seed = SEEDS[tr.pid]
    x0 = np.array(list(seed) + [0.0] * (n_theta - 2)) if theta0 is None \
        else np.array(theta0, float)

    def unpack(p):
        th = p[:n_theta]
        if w_fixed is not None:
            return th, w_fixed[0], w_fixed[1]
        return th, np.exp(p[n_theta]), np.exp(p[n_theta + 1])

    def f(p):
        th, wt, wx = unpack(p)
        ra, rc, om = tr.residuals(th, keep)
        Va, Vc = tr.variances(om, keep, wt, wx)
        return float(np.sum(nll_terms(ra, Va, lik)) + np.sum(nll_terms(rc, Vc, lik)))

    # several starting points: the JFR likelihood is bimodal (small w_t and large
    # w_x, or the reverse with more drift), so the best mode is kept instead of
    # the first one found
    x0s = [x0]
    if theta0 is None and n_theta == 4:
        x0s += [np.array(list(seed) + list(v)) for v in
                ((20.0, -10.0), (-20.0, 10.0), (40.0, 0.0), (0.0, 20.0))]
    best = None
    starts_w = [(1.0, 20.0), (0.2, 5.0), (3.0, 100.0)] if w_fixed is None else [None]
    for x0, start_w in [(a, b) for a in x0s for b in starts_w]:
        p0 = x0 if start_w is None else np.concatenate([x0, np.log(start_w)])
        s = minimize(f, p0, method="BFGS", options=dict(gtol=1e-6))
        s = minimize(f, s.x, method="Nelder-Mead",
                     options=dict(xatol=1e-7, fatol=1e-9, maxiter=20000,
                                  maxfev=20000, adaptive=True))
        if best is None or s.fun < best.fun:
            best = s
    th, wt, wx = unpack(best.x)
    return dict(theta=np.array(th), w_t=float(wt), w_x=float(wx),
                nll=float(best.fun), n_theta=n_theta,
                k=n_theta + (0 if w_fixed is not None else 2),
                cov=gn_cov(tr, keep, th, wt, wx))


def gn_cov(tr, keep, th, wt, wx):
    """Cov of theta by Gauss-Newton with the weights of the solution."""
    def zres(t_):
        ra, rc, om = tr.residuals(t_, keep)
        Va, Vc = tr.variances(om, keep, wt, wx)
        return np.concatenate([ra / np.sqrt(Va), rc / np.sqrt(Vc)])
    steps = np.array([1e-4, 1e-4, 1e-2, 1e-2])[:len(th)]
    J = np.empty((2 * int(keep.sum()), len(th)))
    for j, h in enumerate(steps):
        e = np.zeros(len(th)); e[j] = h
        J[:, j] = (zres(th + e) - zres(th - e)) / (2 * h)
    return np.linalg.inv(J.T @ J)


def zscores(tr, keep_all, sol):
    ra, rc, om = tr.residuals(sol["theta"], keep_all)
    Va, Vc = tr.variances(om, keep_all, sol["w_t"], sol["w_x"])
    return ra, rc, om, ra / np.sqrt(Va), rc / np.sqrt(Vc)


def bic_aic(sol, n_data):
    return (2 * sol["nll"] + sol["k"] * np.log(n_data), 2 * sol["nll"] + 2 * sol["k"])


def select(tr, keep, lik, w_fixed=None, models=(2, 4)):
    force = {"fijo": 2, "deriva": 4}.get(FORCE.get(tr.ridge, ""))
    sols = {}
    for n in models:
        if 2 * keep.sum() <= n + (0 if w_fixed else 2):
            continue
        sols[n] = fit(tr, keep, n, lik, w_fixed)
    n_data = 2 * int(keep.sum())
    best = min(sols, key=lambda n: bic_aic(sols[n], n_data)[0])
    if force in sols:
        best = force
    return sols, best


def fit_track(tr, lik=LIK, w_fixed=None, models=(2, 4), reject=True):
    keep = np.ones(len(tr.df), bool)
    excluded = []
    while True:
        sols, best = select(tr, keep, lik, w_fixed, models)
        sol = sols[best]
        if not reject or lik == "t":
            break
        _, _, _, za, zc = zscores(tr, keep, sol)
        d2 = za ** 2 + zc ** 2
        d2_cand = np.where(tr.pseudo[keep], -1.0, d2)   # the pseudo-obs is not rejected
        j = int(np.argmax(d2_cand))
        if d2[j] <= D2_OUT or keep.sum() <= 6:
            break
        i = np.flatnonzero(keep)[j]
        excluded.append(dict(sample=tr.df.loc[i, "sample"], d2=float(d2[j]),
                             model=best))
        keep[i] = False
    return sols, best, keep, excluded


def sensitivities(tr, keep, sol, w_fixed=None):
    """g_k = d theta / d e_k by central differences (w and sample set fixed)."""
    g = {}
    for k in SYS_BY_RIDGE[tr.ridge]:
        if not tr.is_ar[keep].any():
            continue
        if k in J_LOTE and not (tr.df["irradiation"][keep] == k[2:]).any():
            continue
        h = 10 * SYSTEMATICS[k]
        th = []
        for s in (+h, -h):
            tr.set_ages(tr.scaled(k, s))
            th.append(fit(tr, keep, sol["n_theta"], "gauss",
                          (sol["w_t"], sol["w_x"]), sol["theta"])["theta"])
        tr.set_ages(tr.age0)
        g[k] = ((th[0] - th[1]) / (2 * h)).tolist()
    return g


# ------------------------------------------------------------------ nodes/.rot
def node_rows(pid, theta, t_max):
    """Nodes (age, lat, lon, angle) of plate pid relative to 002 that take
    the seed to h(t). Fixed: one node at 0 and the sentinel at 250. Drift: every
    1 Myr up to t_max (slerp between nodes ~ great circle) and sentinel. With
    modelled drift, every 1 Myr and at each age of the D1 up to its end."""
    t_end = max(t_max, d1_end(pid)) if has_d1(pid) else t_max
    ages = [0.0, 250.0] if (len(theta) == 2 and not has_d1(pid)) else \
        list(np.arange(0.0, np.floor(t_end) + 1.0)) + [t_max, t_end, 250.0]
    if has_d1(pid):
        # all ages of the D1 (spacing 0.2-1 Myr): with nodes every 1 Myr the
        # slerp departs up to ~3 km from the modelled drift at JFR
        ages += [float(a) for a in _D1[pid][0].age if a <= t_end]
    return [(a,) + node_pole(pid, theta, a, t_max)
            for a in sorted(set(round(x, 4) for x in ages))]


def node_pole(pid, theta, age, t_max):
    """(lat, lon, angle) of the minimal rotation that takes the seed to h(age)."""
    s = eng.latlon_axis(*SEEDS[pid])
    h = hotspot_xyz(theta, age, t_max, pid)
    ax = np.cross(s, h)
    sn, cs = np.linalg.norm(ax), float(s @ h)
    if sn < 1e-15:
        return (90.0, 0.0, 0.0)
    la, lo = eng.xyz_to_latlon(ax / sn)
    return (la, lo, np.degrees(np.arctan2(sn, cs)))


def write_rot(fits, out_path, base=eng.ROT_FILE):
    # surrogateescape on reading and writing: the bytes of the other rows pass
    # unchanged (the _hsb .rot files carry UTF-8 comments)
    with open(base, encoding="utf-8", errors="surrogateescape") as f:
        lines = f.readlines()
    out, done = [], set()
    for ln in lines:
        tok = ln.split()
        pid = int(tok[0]) if tok and tok[0].isdigit() else None
        if pid in fits:
            if pid not in done:
                th, t_max, model = fits[pid]
                for (a, la, lo, ang) in node_rows(pid, th, t_max):
                    out.append(f"{pid} {a:.4f} {la:.8f} {lo:.8f} {ang:.8f} {FRAME} "
                               f"!hotspot {model} (03b_hotspot_fit.py) @REF this study\n")
                done.add(pid)
            continue
        out.append(ln)
    # plates missing from the base .rot (16, Foundation, on the first run)
    for pid in sorted(set(fits) - done):
        th, t_max, model = fits[pid]
        if out and not out[-1].endswith("\n"):
            out[-1] += "\n"
        for (a, la, lo, ang) in node_rows(pid, th, t_max):
            out.append(f"{pid} {a:.4f} {la:.8f} {lo:.8f} {ang:.8f} {FRAME} "
                       f"!hotspot {model} (03b_hotspot_fit.py) @REF this study\n")
    with open(out_path, "w", encoding="utf-8", errors="surrogateescape", newline="\n") as f:
        f.writelines(out)


# ------------------------------------------------------------------ main
def main():
    M = eng.RotModel()
    res, rot_fits = {}, {}
    names = {2: "fijo", 4: "deriva"}
    w_jfr = None
    for ridge, pid, fn, col in TRACKS:
        df = pd.read_csv(os.path.join(HERE, "data", fn), comment="#")
        df = df[(df.use_fit == 1) & ~df["sample"].isin(EXCLUDE)]
        if USE_P0:
            p0 = {c: np.nan for c in df.columns}
            p0.update(sample=P0_NAME, lat=SEEDS[pid][0], lon=SEEDS[pid][1],
                      age_ma=0.0, err2s_ma=0.0, age_type="present", use_fit=1)
            for c in (col if isinstance(col, tuple) else (col,)):
                p0[c] = 0.0
            if "pos_err_km" in df.columns or P0_KM > 0:
                df = df.assign(pos_err_km=df.get("pos_err_km", 0.0))
                p0["pos_err_km"] = P0_KM
            df = pd.concat([df, pd.DataFrame([p0])], ignore_index=True)
        tr = Track(M, ridge, pid, df, col)
        if ridge == "TT":
            sols, best, keep, excl = fit_track(tr, "gauss", w_jfr, (2,), False)
        else:
            sols, best, keep, excl = fit_track(tr)
        sol = sols[best]
        if ridge == "JFR":
            w_jfr = (sol["w_t"], sol["w_x"])
        wf = w_jfr if ridge == "TT" else None
        g = sensitivities(tr, keep, sol, wf)
        cov = np.array(sol["cov"])
        if W_MODE == "pred":
            # prediction dispersion: w added to the position of h0 (east, north)
            _, _, om, _, _ = zscores(tr, keep, sol)
            s_al = np.median(om) * sol["w_t"] / KM_PER_DEG
            s_ac = sol["w_x"] / KM_PER_DEG
            s2 = 0.5 * (s_al ** 2 + s_ac ** 2)
            coslat = np.cos(np.radians(sol["theta"][0]))
            cov[0, 0] += s2
            cov[1, 1] += s2 / coslat ** 2
        n_data = 2 * int(keep.sum())
        info = {}
        for n, s in sols.items():
            b, a = bic_aic(s, n_data)
            info[names[n]] = dict(theta=s["theta"].tolist(), w_t=s["w_t"],
                                  w_x=s["w_x"], nll=s["nll"], bic=b, aic=a)
        all_keep = np.ones(len(tr.df), bool)
        ra, rc, om, za, zc = zscores(tr, all_keep, sol)
        out = tr.df[["sample", "age_ma", "err2s_ma", "lat", "lon"]].copy()
        out["age_used"] = tr.age0
        out["r_along_km"], out["r_across_km"] = ra, rc
        out["dt_myr"] = ra / om
        out["omega_km_myr"] = om
        out["z_along"], out["z_across"] = za, zc
        out["d2"] = za ** 2 + zc ** 2
        out["kept"] = keep.astype(int)
        out.to_csv(os.path.join(HERE, f"hotspot_residuals_{ridge}{OUT_TAG}.csv"),
                   index=False, float_format="%.4f")
        k = keep
        res[ridge] = dict(plate=pid, rel=tr.rel, model=names[best], theta=sol["theta"].tolist(),
                          cov=cov.tolist(), w_t=sol["w_t"], w_x=sol["w_x"],
                          w_from_jfr=ridge == "TT", t_max=tr.t_max,
                          g=g, n_used=int(k.sum()), n_input=len(k),
                          excluded=excl, models=info,
                          rms_along_km=float(np.sqrt(np.mean(ra[k] ** 2))),
                          rms_across_km=float(np.sqrt(np.mean(rc[k] ** 2))),
                          seed_shift_km=float(eng.gc_km(*SEEDS[pid], *sol["theta"][:2])))
        if tr.rel != REL:   # plate change: tr.rel up to t_switch, REL before
            res[ridge].update(t_switch=tr.t_max, rel_old=REL)
        rot_fits[pid] = (sol["theta"], tr.t_max, names[best])
        sd = np.sqrt(np.diag(cov))
        print(f"{ridge:6s} model {names[best]:6s} n={k.sum()}/{len(k)}  "
              f"h0 = ({sol['theta'][0]:.3f}, {sol['theta'][1]:.3f}) "
              f"+/- ({sd[0]*KM_PER_DEG:.1f} km, {sd[1]*KM_PER_DEG*np.cos(np.radians(sol['theta'][0])):.1f} km)"
              + (f"  drift ({sol['theta'][2]:.1f}, {sol['theta'][3]:.1f}) "
                 f"+/- ({sd[2]:.1f}, {sd[3]:.1f}) km/Myr" if len(sd) == 4 else "")
              + f"  w_t {sol['w_t']:.2f} Myr  w_x {sol['w_x']:.1f} km"
              + ("  (from JFR)" if ridge == "TT" else ""), flush=True)
        for n, s in info.items():
            print(f"         {n:6s}: BIC {s['bic']:8.2f}  AIC {s['aic']:8.2f}  "
                  f"w_t {s['w_t']:.2f}  w_x {s['w_x']:.1f}")
        for e in excl:
            print(f"         excluded {e['sample']} (d2 = {e['d2']:.1f})")
        for kk, v in g.items():
            print(f"         g[{kk}] = {np.round(v, 3).tolist()} (sigma {SYSTEMATICS[kk]:.4f})")
    meta = dict(lik=LIK, w_mode=W_MODE, easter_ages=EASTER_COL, rel=REL,
                present_obs=USE_P0, present_pos_err_km=P0_KM, present_mode=P0_MODE,
                prior=PRIOR, rot_base=eng.ROT_FILE, d2_threshold=D2_OUT, systematics=SYSTEMATICS,
                excluded_input=EXCLUDE, lag_myr=LAG, force=FORCE)
    tag = "" if (LIK, W_MODE, EASTER_KEY, USE_P0, P0_KM, P0_MODE) == \
        ("gauss", "conf", "age_E", False, 0.0, "iso") else \
        f"_{LIK}_{W_MODE}_{EASTER_KEY}" + ("_conP0" if USE_P0 else "") + \
        (f"_p0{P0_KM:g}km" if P0_KM > 0 else "") + \
        ("" if P0_MODE == "iso" else f"_p0{P0_MODE}")
    tag += OUT_TAG
    with open(os.path.join(HERE, f"hotspot_fit{tag}.json"), "w", encoding="utf-8") as f:
        json.dump(dict(meta=meta, tracks=res), f, indent=1)
    os.makedirs(os.path.join(HERE, "rot"), exist_ok=True)
    base = os.path.splitext(os.path.basename(eng.ROT_FILE))[0]
    hs = "_hsb" if PRIOR == "none" else f"_hsb{PRIOR}"
    rot_out = ROT_OUT or os.path.join(HERE, "rot", f"{base}{hs}{tag}.rot")
    write_rot(rot_fits, rot_out)
    if fz_switch(rot_out, res):
        with open(os.path.join(HERE, f"hotspot_fit{tag}.json"), "w", encoding="utf-8") as f:
            json.dump(dict(meta=meta, tracks=res), f, indent=1)
    print(f"hotspot_fit{tag}.json and {rot_out} written")
    check_rot(rot_out, res)


def fz_switch(path, res):
    """t_switch = N->S crossing of the FZ in FZ_SWITCH by the hotspot (see above), with the
    .rot just written; never younger than t_max. Returns True if any t_switch changed."""
    M2 = eng.RotModel(path)
    taus = np.round(np.arange(0.0, 60.0001, 0.02), 3)
    changed = False
    for ridge, fn in FZ_SWITCH.items():
        r = res.get(ridge)
        if r is None or "t_switch" not in r:
            continue
        F = pd.read_csv(os.path.join(HERE, "data", fn), comment="#").sort_values("lon")
        P = M2.motion_path(SEEDS[r["plate"]], r["plate"], r["rel_old"], ANCHOR, taus, 0.0)
        la, lo = P[:, 0], P[:, 1]
        inside = (lo >= F.lon.min()) & (lo <= F.lon.max())
        side = np.where(inside, la - np.interp(lo, F.lon, F.lat), np.nan)
        # N->S crossings towards the present: north (side > 0) at tau, south at tau - dtau
        k = np.flatnonzero((side[1:] > 0) & (side[:-1] <= 0))
        t_fz = [float(taus[i] + side[i] / (side[i] - side[i + 1]) * (taus[i + 1] - taus[i]))
                for i in k]
        t_fz = max(t_fz) if t_fz else None
        t_sw = r["t_max"] if t_fz is None or t_fz < r["t_max"] else t_fz
        r.update(t_switch=t_sw, t_switch_fz=t_fz,
                 t_switch_rule="fz" if t_sw == t_fz else "t_max")
        changed = True
        print(f"  t_switch {ridge}: FZ crossing {fn} = "
              + ("no crossing" if t_fz is None else f"{t_fz:.2f} Ma")
              + f"; t_max = {r['t_max']:.2f} Ma -> t_switch = {t_sw:.2f} Ma ({r['t_switch_rule']})")
    return changed


def check_rot(path, res):
    """The written .rot must reproduce x(t) = B(t) h(t) of the fit."""
    M2 = eng.RotModel(path)
    M = eng.RotModel()
    for ridge, r in res.items():
        th = np.array(r["theta"])
        err = 0.0
        for t in (0.0, 1.3, 5.0, 12.7, r["t_max"], 40.0):
            p = M2.motion_path(SEEDS[r["plate"]], r["plate"], r["rel"], ANCHOR, [t], 0.0)[0]
            x = b_matrix(M, t, r["rel"]) @ hotspot_xyz(th, t, r["t_max"], r["plate"])
            err = max(err, eng.gc_km(p[0], p[1], *eng.xyz_to_latlon(x)))
        print(f"  .rot check {ridge}: max error {err:.3f} km")


if __name__ == "__main__":
    main()
