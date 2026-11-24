"""40Ar/39Ar ages of the Foundation chain (O'Connor et al. 1998, EPSL 164, 41-59,
Table 1, pp. 44-45, transcribed from the PDF image) recalibrated to
FCs 28.201 Ma (Kuiper et al. 2008) with the decay constant of Min et al. (2000), as those
of Easter and San Felix (ArAr_recalibration.py).

Fourth chain of the paper: the Foundation hotspot, which according to
Bello-Gonzalez et al. (2018, Gondwana Res. 64) built the Iquique Ridge
on the Farallon/Nazca plate before moving to the Pacific plate (Foundation chain).
In the pipeline the track lies on the Pacific plate (901) up to the oldest dated
age (Seamount 1a) and on Nazca/Farallon (902003) before.

Sources and assumptions:
- Plateau ages +/- 1 sigma (Table 1 = plateau of Table 2; "Age +/- 1s"). Here
  err2s_ma = 2 x 1 sigma, as in the other files.
- Position: the "on bottom" position of Table 1 (degrees and decimal minutes; S and W
  negative). For SO100 50DS-1 the "off bottom" position reads 34 51.860'S against
  35 51.492'S "on bottom" (probable erratum in the paper); "on bottom" is used.
- Monitor: TCR 85G003 = 27.92 Ma (O'Connor et al. 1998, p. 46, ref. [23]).
- TCs/FCs: R = 1.0112 +/- 0.0010 (Renne et al. 1998, Chem. Geol. 145, p. 125; the
  sigma level of the +/- is not stated in that sentence: it is taken as 1 sigma, as the
  ACs/FCs in ArAr_recalibration.py).
- Decay constant used by O'Connor et al. (1998): not stated. ASSUMPTION: Steiger and
  Jager (1977), 5.543e-10/a, the one in general use in 1998 (same assumption as for Daze).
- O'Connor et al. (1998, p. 53-54) leave Seamounts 8, 9 and 13a out of their regression
  (column oconnor_excl). Here all enter (use_fit = 1); 03b applies its
  own outlier rejection (d2 > 11.83).

Output: outputs/ages_foundation_oconnor1998_recal.csv (and a copy in <folder>/data/
of the eight models).
"""
import math
import os
import shutil
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, ".."))

LAM_SJ = 5.543e-10 * 1e6
LAM_MIN = 5.463e-10 * 1e6
FCS_NEW = 28.201
TCR_OLD = 27.92
R_TCS_FCS, R_TCS_FCS_1S = 1.0112, 0.0010


def F(t, lam):
    return math.expm1(lam * t)


def recal(t_u, t_s_old, lam_old, t_s_new, lam_new):
    R = F(t_u, lam_old) / F(t_s_old, lam_old)
    return math.log1p(R * F(t_s_new, lam_new)) / lam_new


def dm(d, m):
    return -(d + m / 60.0)


# sample, seamount, lat (deg, min), lon (deg, min) on bottom, material, age, 1s, excl. O'Connor
ROWS = [
    ("SO100-11DS-1", "1a", (32, 56.414), (130, 45.459), "whole rock", 21.2, 0.1, 0),
    ("FHDR1-3", "between 1a and 5", (32, 30.33), (127, 30.09), "whole rock", 18.2, 0.2, 0),
    ("SO100-18DS-1", "5", (32, 28.871), (126, 0.456), "whole rock", 16.3, 0.2, 0),
    ("SO100-17DS-1", "5", (32, 28.474), (126, 4.069), "whole rock", 16.1, 0.1, 0),
    ("SO100-28GTV-2", "8", (33, 41.735), (124, 54.612), "plagioclase", 8.8, 0.1, 1),
    ("SO100-25DS-1", "9", (33, 20.764), (123, 52.525), "whole rock", 13.1, 0.1, 1),
    ("SO100-26DS-1", "9", (33, 31.700), (124, 6.434), "plagioclase", 12.9, 0.3, 1),
    ("SO100-33DS-1", "10", (34, 7.384), (122, 22.038), None, 12.5, 0.1, 0),
    ("SO100-38DS-1", "11", (34, 19.119), (121, 58.729), None, 11.6, 0.1, 0),
    ("SO100-41DS-1", "12b", (34, 52.352), (121, 33.309), "plagioclase", 11.6, 0.1, 0),
    ("SO100-45DS-1", "13a", (35, 3.180), (120, 43.226), "plagioclase", 13.4, 0.1, 1),
    ("SO100-46DS-2", "16", (34, 57.336), (120, 24.490), "whole rock", 10.4, 0.1, 0),
    ("SO100-50DS-1", "18", (35, 51.492), (119, 6.724), None, 9.4, 0.1, 0),
    ("SO100-54DS-1", "19b", (35, 6.976), (118, 33.132), None, 8.4, 0.1, 0),
    ("SO100-56DS-1", "21a", (35, 22.943), (118, 5.256), None, 7.7, 0.1, 0),
    ("SO100-63DS-1", "22", (35, 48.000), (117, 26.305), None, 7.7, 0.1, 0),
    ("SO100-60DS-1", "23", (35, 27.003), (117, 11.833), None, 7.6, 0.5, 0),
    ("SO100-59DS-1", "24b", (35, 26.729), (116, 38.813), None, 7.2, 0.1, 0),
    ("SO100-67DS-4", "25", (36, 1.980), (115, 59.276), None, 6.2, 0.1, 0),
    ("SO100-66DS-1", "26", (35, 47.484), (115, 39.335), None, 4.7, 0.1, 0),
    ("SO100-69DS-1", "27b", (36, 33.804), (115, 16.627), None, 5.1, 0.1, 0),
    ("SO100-70DS-2", "28", (36, 20.907), (113, 55.696), None, 3.7, 0.1, 0),
    ("SO100-71DS-1", "30", (36, 40.812), (113, 28.350), None, 2.1, 0.1, 0),
]


def main():
    tcs_new = math.log1p(R_TCS_FCS * F(FCS_NEW, LAM_MIN)) / LAM_MIN
    rows = []
    for s, smt, la, lo, mat, age, s1, ex in ROWS:
        a = recal(age, TCR_OLD, LAM_SJ, tcs_new, LAM_MIN)
        rows.append(dict(sample=s, edifice=f"Seamount {smt}", lat=round(dm(*la), 6),
                         lon=round(dm(*lo), 6), age_ma=age, err2s_ma=2 * s1,
                         age_type="plateau", use_fit=1, segment="Foundation",
                         material=mat or "", oconnor_excl=ex, age_recal=a,
                         d_myr=a - age))
    df = pd.DataFrame(rows)
    hdr = ("# Edades 40Ar/39Ar de la cadena Foundation, O'Connor et al. (1998, EPSL 164, Tabla 1);\n"
           "# err2s_ma = 2 x el 1 sigma publicado. age_recal: TCR 27.92 Ma (SJ, supuesto) a FCs 28.201 /\n"
           f"# Min 2000 con TCs/FCs = {R_TCS_FCS} (Renne et al. 1998): TCs = {tcs_new:.4f} Ma.\n"
           "# oconnor_excl = 1: fuera de la regresion de O'Connor et al. (Seamounts 8, 9, 13a).\n"
           "# Generado por analysis/foundation_oconnor1998.py\n")
    out = os.path.join(HERE, "outputs", "ages_foundation_oconnor1998_recal.csv")
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write(hdr)
        df.to_csv(f, index=False, float_format="%.6f")
    print(f"new TCs = {tcs_new:.4f} Ma; factor = {df.age_recal.iloc[0] / df.age_ma.iloc[0]:.6f}")
    print(df[["sample", "lat", "lon", "age_ma", "age_recal", "oconnor_excl"]].to_string(index=False))
    for d in ("Q22_WK08-A", "WD25_WK08-A", "Q22_GMHRF-A", "Q22_GMHRF-B",
              "WD25_GMHRF-B", "Q22_GWG20",
              "WD25_GWG20", "SG12_WK08-A"):
        dd = os.path.join(P, "models", d, "data")
        os.makedirs(dd, exist_ok=True)
        shutil.copy2(out, dd)
    print("copied to data/ of the eight models")


if __name__ == "__main__":
    main()
