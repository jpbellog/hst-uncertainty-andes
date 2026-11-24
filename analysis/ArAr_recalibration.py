"""Recalibration of the 40Ar/39Ar ages to FCs 28.201 Ma (Kuiper et al. 2008) with the
total decay constant of Min et al. (2000).

Does not modify the pipeline CSVs: reads `models/Q22_WK08-A/data/ages_*.csv` and writes
separate tables to `outputs/`. The pipeline CSVs change only when the pipeline is regenerated.

Equation (same irradiation, monitor s, sample u; Renne et al. 1998):
    R = (exp(lam*t_u) - 1) / (exp(lam*t_s) - 1)      (independent of the calibration)
    t_u' = ln(1 + R * (exp(lam'*t_s') - 1)) / lam'
The recalibration uses the published age (not the step-heating data).

Values and their sources (all read in the original publications):
- lam SJ = 5.543e-10 /a: Steiger and Jager (1977), p. 359 (lam_b 4.962 + lam_e 0.581).
- lam Min = 5.463e-10 /a: Min et al. (2000), Table 3 (corrected Endt and Van der Leun).
  Kuiper et al. (2008, p. 3) cite it as 5.463 +/- 0.214 at 2 sigma.
- FCs = 28.201 +/- 0.046 Ma (2 sigma, includes decay constant and astronomical age):
  Kuiper et al. (2008), p. 3 and note 6 ("all errors are stated at the 2 sigma level").
  With lam SJ the same calibration gives 28.198 +/- 0.044 Ma (p. 3).
- Easter (Ray et al. 2012, p. 7): monitor "FCT-3 biotite (28.03 +/- 0.18 Ma, 1s)
  (Renne et al., 1998)", SJ constants. Note: Renne et al. (1998) give no age for
  FCT-3 (the PDF does not contain "FCT"); 28.03 +/- 0.18 is the astronomical age of the
  FCs SANIDINE of Renne et al. (1994), according to Min et al. (2000, p. 8).
  Scenario A: FCT-3 treated as equivalent to FCs (as implied by the value used).
  Scenario B: FCT-3/FCs from Daze et al. (2003): FCT-3 = 28.13 +/- 0.47 Ma (2 sigma)
  with FCs = 27.98 Ma; the factor 1.005 +/- 0.017 (2 sigma) is taken as systematic.
- San Felix (Haase et al. 2000): "essentially following the method described by
  Sharp et al. (1996)", SJ constants, 1 sigma errors (caption of their Fig. 2). Sharp et al.
  (1996, p. 2) use Alder Creek sanidine at 1.186 Ma (Turrin et al. 1994).
  [INFERRED: Haase et al. do not name the monitor; L. Lara was asked.]
  ACs/FCs = 0.04229 +/- 0.00006 (Renne et al. 1998, abstract and p. 12).
- JFR (Lara et al. 2018): already on FCs 28.201 / Min 2000 (Methods of Lara et al.). No change.
- San Ambrosio (Gonzalez-Ferran 1987, K-Ar): no PDF, unknown constants.
  Only the SJ -> Min sensitivity is estimated.
"""
import math
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "models" / "Q22_WK08-A" / "data"
OUT = Path(__file__).resolve().parent / "outputs"
OUT.mkdir(exist_ok=True)

LAM_SJ = 5.543e-10 * 1e6     # /Myr
LAM_MIN = 5.463e-10 * 1e6
FCS_NEW, FCS_NEW_2S = 28.201, 0.046
FCS_OLD_RAY = 28.03          # used by Ray et al. as FCT-3
ACS_SHARP = 1.186            # Sharp et al. (1996)
R_ACS_FCS = 0.04229          # Renne et al. (1998)
DAZE_FCS, DAZE_FCT3, DAZE_FACTOR_2S = 27.98, 28.13, 0.017


def F(t, lam):
    return math.expm1(lam * t)


def recal(t_u, t_s_old, lam_old, t_s_new, lam_new):
    R = F(t_u, lam_old) / F(t_s_old, lam_old)
    return math.log1p(R * F(t_s_new, lam_new)) / lam_new


def age_from_R(R, t_ref, lam):
    return math.log1p(R * F(t_ref, lam)) / lam


# Monitors in the new calibration
fct3_A = FCS_NEW
R_fct3 = F(DAZE_FCT3, LAM_SJ) / F(DAZE_FCS, LAM_SJ)          # Daze et al. use SJ [assumed]
fct3_B = age_from_R(R_fct3, FCS_NEW, LAM_MIN)
acs_new = age_from_R(R_ACS_FCS, FCS_NEW, LAM_MIN)

lines = ["# Recalibracion 40Ar/39Ar a FCs 28.201 Ma / Min et al. (2000)", "",
         f"Monitor de Ray (28.03 Ma, SJ) en la calibracion nueva: A = {fct3_A:.3f} Ma; "
         f"B (Daze) = {fct3_B:.3f} Ma (R FCT-3/FCs = {R_fct3:.5f})",
         f"ACs (Sharp 1.186 Ma, SJ) en la calibracion nueva: {acs_new:.4f} Ma", ""]

# Easter
ray = pd.read_csv(DATA / "ages_easter_ray2012.csv", comment="#")
ray["age_A"] = [recal(t, FCS_OLD_RAY, LAM_SJ, fct3_A, LAM_MIN) for t in ray.age_ma]
ray["age_B"] = [recal(t, FCS_OLD_RAY, LAM_SJ, fct3_B, LAM_MIN) for t in ray.age_ma]
ray["dA_myr"] = ray.age_A - ray.age_ma
ray["dB_myr"] = ray.age_B - ray.age_ma
# common systematic (2 sigma, Myr): Kuiper monitor +/- FCT-3/FCs intercalibration
ray["sys_kuiper_2s"] = ray.age_A * FCS_NEW_2S / FCS_NEW
ray["sys_daze_2s"] = ray.age_A * DAZE_FACTOR_2S
ray.to_csv(OUT / "ages_easter_ray2012_recal.csv", index=False, float_format="%.4f")

# San Felix (40Ar/39Ar ages only)
tt = pd.read_csv(DATA / "ages_taltal.csv", comment="#")
ar = tt.age_type != "unverified"
tt["age_recal"] = tt.age_ma
tt.loc[ar, "age_recal"] = [recal(t, ACS_SHARP, LAM_SJ, acs_new, LAM_MIN) for t in tt.age_ma[ar]]
tt["d_myr"] = tt.age_recal - tt.age_ma
tt.to_csv(OUT / "ages_taltal_recal.csv", index=False, float_format="%.4f")

# K-Ar San Ambrosio: t = ln(1 + (lam/lam_e) Ar*/K)/lam; SJ -> Min sensitivity
lam_e_sj, lam_e_min = 0.581e-4, 0.580e-4   # /Myr (SJ 1977; Min 2000 Table 3)
t_kar = 2.93
ratio = F(t_kar, LAM_SJ) * lam_e_sj / LAM_SJ
t_kar_min = math.log1p(ratio * LAM_MIN / lam_e_min) / LAM_MIN

u = ray[ray.use_fit == 1]
lines += ["## Pascua (Ray et al. 2012), samples del ajuste", "",
          "| | min | mediana | max |", "|---|---|---|---|",
          f"| cambio A (Myr) | {u.dA_myr.min():+.3f} | {u.dA_myr.median():+.3f} | {u.dA_myr.max():+.3f} |",
          f"| cambio B (Myr) | {u.dB_myr.min():+.3f} | {u.dB_myr.median():+.3f} | {u.dB_myr.max():+.3f} |",
          f"| cambio A relativo (%) | {(100*u.dA_myr/u.age_ma).min():+.3f} | {(100*u.dA_myr/u.age_ma).median():+.3f} | {(100*u.dA_myr/u.age_ma).max():+.3f} |",
          f"| sistematico Kuiper 2s (Myr) | {u.sys_kuiper_2s.min():.3f} | {u.sys_kuiper_2s.median():.3f} | {u.sys_kuiper_2s.max():.3f} |",
          f"| sistematico Daze 2s (Myr) | {u.sys_daze_2s.min():.3f} | {u.sys_daze_2s.median():.3f} | {u.sys_daze_2s.max():.3f} |",
          "", "## San Felix (Haase et al. 2000)", ""]
for r in tt[ar].itertuples():
    lines.append(f"- {r.sample}: {r.age_ma:.3f} -> {r.age_recal:.4f} Ma ({1000*r.d_myr:+.1f} ka)")
lines += ["", "## San Ambrosio K-Ar (sensibilidad)", "",
          f"- 2.93 Ma con SJ -> {t_kar_min:.4f} Ma con Min ({1000*(t_kar_min-t_kar):+.1f} ka), "
          "si Gonzalez-Ferran uso SJ [sin verificar]", ""]
(OUT / "resumen_recalibracion.md").write_text("\n".join(lines), encoding="utf-8")
print("\n".join(lines))
