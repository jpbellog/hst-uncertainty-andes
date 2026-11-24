# -*- coding: utf-8 -*-
"""STEP 05 -- Uncertainty of the convergence velocity (Quiero et al. 2022).

Uses the Nazca-South America stage Euler vectors with covariance (Supplementary
Data 2 of Quiero et al. 2022, Nat. Commun. 13:1339; file MOESM5.txt).
Format per row:
t_ini (Ma), t_fin (Ma), pole lon (degE), pole lat (degN), omega (deg/Myr),
Cxx, Cxy, Cxz, Cyy, Cyz, Czz (rad^2/Myr^2).

The velocity at a margin point r is v = w x r (w = Euler vector).
Exact linear propagation: sigma_v = |J C J^T|^(1/2) with J = d(vc)/dw = -[r]_x
(per Cartesian component), projected onto the local frame (East, North) to
obtain the sigma of the magnitude, the azimuth and the component ORTHOGONAL to the margin.

Output: vc_sigma_quiero.csv (stage mid-age x margin latitude):
vc_mag_cm_yr, sigma_mag, azimuth, sigma_az_deg, vc_ort, sigma_ort.
Local margin azimuth: from Margin_SAM_0.1.csv at age 0 (present margin),
consistent with the use of the VC field elsewhere in the pipeline.
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SI = (__import__("os").environ["PAPER_ROOT"] + "/external/quiero2022/MOESM5.txt")
MARGIN_CSV = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "margin", "Margin_SAM_0.1.csv")  # copy of Margin_SAM (see rotation_models/README.md)
R_E_CM_YR = 6371e5 / 1e6   # cm per (rad/Myr): 6371 km = 6.371e8 cm; /Myr -> cm/yr = 6.371e8/1e6
D2R = np.pi / 180.0

def unit(lat, lon):
    la, lo = lat * D2R, lon * D2R
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])

def local_frame(lat, lon):
    la, lo = lat * D2R, lon * D2R
    e = np.array([-np.sin(lo), np.cos(lo), 0.0])
    n = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    return e, n

def margin_nodes():
    m = pd.read_csv(MARGIN_CSV)
    row = m.iloc[(m["age"] - 0).abs().idxmin()]
    nn = sum(1 for c in m.columns if c.startswith("Lon_"))
    lon = np.array([row.get(f"Lon_{k}") for k in range(1, nn + 1)], float)
    lat = np.array([row.get(f"Lat_{k}") for k in range(1, nn + 1)], float)
    k = np.isfinite(lat) & np.isfinite(lon) & (lat <= -5) & (lat >= -36)
    return lat[k], lon[k]

def gts2020(S):
    """Converts the stages (t_ini, t_fin in GTS2012) to GTS2020 and rescales rate and covariance."""
    import importlib.util
    p = os.path.join(os.path.dirname(HERE), "analysis", "timescale_GTS2020.py")
    sp = importlib.util.spec_from_file_location("esc", p)
    esc = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(esc)
    f = esc.remap("gts12")
    S = S.copy()
    t0, t1 = f(S[:, 0]), f(S[:, 1])
    k = (S[:, 1] - S[:, 0]) / (t1 - t0)
    S[:, 0], S[:, 1] = t0, t1
    S[:, 4] *= k
    S[:, 5:11] *= (k ** 2)[:, None]
    return S


if __name__ == "__main__":
    S = np.loadtxt(SI)
    # GTS2020: the stage boundaries of Quiero et al. (2022) are in GTS2012 and the
    # rotation model is in GTS2020 (analysis/rot_gts2020.py). The
    # finite rotation of each stage does not change, so the pole is kept, the rate is scaled by
    # dt_GTS2012/dt_GTS2020 and the covariance by the square of that factor.
    S = gts2020(S)
    mlat, mlon = margin_nodes()
    # local margin azimuth (southward, node to node)
    rows = []
    for srow in S:
        t0, t1, plon, plat, om = srow[:5]
        Cw = np.array([[srow[5], srow[6], srow[7]],
                       [srow[6], srow[8], srow[9]],
                       [srow[7], srow[9], srow[10]]])   # rad^2/Myr^2
        w = unit(plat, plon) * om * D2R                 # rad/Myr
        # all margin nodes (43 in the 5-36 S window), so that latitudes close to
        # 20 S / 30 S, needed by the VC uncertainty table, are sampled.
        # The margin azimuth uses nodes j -> j+10.
        for j in range(0, len(mlat), 1):
            la, lo = mlat[j], mlon[j]
            r = unit(la, lo)
            v3 = np.cross(w, r) * R_E_CM_YR             # cm/yr
            J = -np.array([[0, -r[2], r[1]],
                           [r[2], 0, -r[0]],
                           [-r[1], r[0], 0]]) * R_E_CM_YR
            Cv = J @ Cw @ J.T
            e, n = local_frame(la, lo)
            vE, vN = v3 @ e, v3 @ n
            C2 = np.array([[e @ Cv @ e, e @ Cv @ n],
                           [e @ Cv @ n, n @ Cv @ n]])
            vmag = np.hypot(vE, vN)
            u_mag = np.array([vE, vN]) / vmag
            s_mag = float(np.sqrt(u_mag @ C2 @ u_mag))
            u_az = np.array([vN, -vE]) / vmag           # transverse direction
            s_az = float(np.sqrt(u_az @ C2 @ u_az)) / vmag / D2R
            az = (np.degrees(np.arctan2(vE, vN))) % 360
            # local margin azimuth: towards the next node (south)
            j2 = min(j + 10, len(mlat) - 1)
            dE = (mlon[j2] - lo) * np.cos(la * D2R)
            dN = mlat[j2] - la
            azm = np.degrees(np.arctan2(dE, dN)) % 360
            # component orthogonal to the margin
            t_hat = np.array([np.sin(azm * D2R), np.cos(azm * D2R)])
            n_hat = np.array([t_hat[1], -t_hat[0]])
            vort = float(np.array([vE, vN]) @ n_hat)
            s_ort = float(np.sqrt(n_hat @ C2 @ n_hat))
            rows.append(dict(t_ini=t0, t_fin=t1, age_mid=0.5 * (t0 + t1),
                             lat=la, lon=lo, vc_mag=vmag, sigma_mag=s_mag,
                             az=az, sigma_az=s_az, vc_ort=abs(vort),
                             sigma_ort=s_ort))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(HERE, "vc_sigma_quiero.csv"), index=False)
    # summary at 20S and 30S
    for latq in (-20, -30):
        s = df.iloc[(df.lat - latq).abs().argsort()].groupby("age_mid").first()
        s = s.sort_index()
        print(f"--- margin ~{latq} S (VC +/- 1sigma, cm/yr; azimuth +/- deg):")
        for a, r in s.iterrows():
            print(f"  {r.t_ini:5.2f}-{r.t_fin:5.2f} Ma: "
                  f"VC {r.vc_mag:5.2f} +/- {r.sigma_mag:4.2f}  "
                  f"az {r.az:5.1f} +/- {r.sigma_az:3.1f}  "
                  f"VC_ort {r.vc_ort:5.2f} +/- {r.sigma_ort:4.2f}")
