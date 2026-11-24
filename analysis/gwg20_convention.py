# -*- coding: utf-8 -*-
"""Unit and convention of the GWG20 covariance (Table S6 of Gaastra et al. 2022).

1. Magnitudes: maximum sigma read as rad^2 and as rad^2 Ma^-2 (x t^2).
2. Side of the perturbation: compares the published covariance with a two-hotspot model
   (Hawaii and Louisville, along-track/across-track = q) expressed in the hotspot frame
   (perturbation on the left of PAC->HS) or rotated to the Pacific (on the right).
   The hotspot positions are approximate: they only serve to discriminate the convention.
3. Effect on the 1 sigma ellipse of the Easter and JFR positions on the Pacific plate: C' = C (current)
   versus C' = B C B^T, with B = A^T the rotation of row 002.

Usage:  py -I gwg20_convention.py
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
S6 = os.path.join(HERE, "..", "external", "gaastra2022",
                  "Supporting_Information", "2021tc006772-sup-0003-table si-s02.csv")


def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def rot(la, lo, ang):
    k = xyz(la, lo)
    a = np.radians(ang)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K


def cross_m(p):
    return np.array([[0, -p[2], p[1]], [p[2], 0, -p[0]], [-p[1], p[0], 0]])


def cov(r):
    return np.array([[r.a, r.b, r.c], [r.b, r.d, r.e], [r.c, r.e, r.f]])


def ellipse(Cp, p):
    e = np.cross([0, 0, 1], p)
    e /= np.linalg.norm(e)
    T = np.vstack([e, np.cross(p, e)])
    w, V = np.linalg.eigh(T @ Cp @ T.T)
    return np.sqrt(w[1]) * 6371, np.sqrt(w[0]) * 6371, np.degrees(np.arctan2(V[0, 1], V[1, 1])) % 180


def main():
    d = pd.read_csv(S6)
    d.columns = ["age", "lat", "lon", "ang", "a", "b", "c", "d", "e", "f", "k"]
    rows = {int(r.age): r for r in d.itertuples()}

    print("1. maximum sigma (deg): rad^2 | rad^2/kappa | rad^2 Ma^-2 x t^2")
    for t in (1, 10, 20, 30, 40, 50, 60, 70, 80):
        r = rows[t]
        s = np.degrees(np.sqrt(np.linalg.eigvalsh(1e-5 * cov(r)).max()))
        print(f"   {t:3d} Ma  {s:5.2f} | {s / np.sqrt(r.k):5.2f} | {s * t:6.1f}")

    hs = [xyz(19.0, -155.2), xyz(-51.0, -138.0)]
    print("\n2. normalized distance to the two-hotspot model (smaller = better): left | right")
    for t in (12, 15, 20, 25, 30, 35, 40, 45):
        r, r1 = rows[t], rows[t + 1]
        A, A1 = rot(r.lat, r.lon, r.ang), rot(r1.lat, r1.lon, r1.ang)
        Ct = cov(r) / np.trace(cov(r))
        best = {}
        for conv in ("izq", "der"):
            b = (9.0, 0.0)
            for q in np.linspace(1, 15, 57):
                J = np.zeros((3, 3))
                for h in hs:
                    dd = A1 @ A.T @ h - h
                    dd -= h * (h @ dd)
                    dd /= np.linalg.norm(dd)
                    a1, a2 = np.cross(h, dd), np.cross(h, np.cross(h, dd))
                    J += np.outer(a1, a1) / q ** 2 + np.outer(a2, a2)
                Cm = np.linalg.inv(J)
                if conv == "der":
                    Cm = A.T @ Cm @ A
                b = min(b, (np.linalg.norm(Ct - Cm / np.trace(Cm)), q))
            best[conv] = b
        print(f"   {t:3d} Ma (angle {r.ang:5.1f})  {best['izq'][0]:.3f} (q {best['izq'][1]:4.1f}) | "
              f"{best['der'][0]:.3f} (q {best['der'][1]:4.1f})")

    print("\n3. 1 sigma ellipse of the hotspot position on the Pacific plate (C/kappa), km and azimuth:")
    for nm, h in (("Easter", xyz(-27.1, -109.3)), ("JFR", xyz(-33.6, -78.8))):
        for t in (10, 20, 30, 40, 50, 60):
            r = rows[t]
            A = rot(r.lat, r.lon, r.ang)
            C = 1e-5 / r.k * cov(r)
            p = A.T @ h
            o = ellipse(cross_m(p) @ C @ cross_m(p).T, p)
            n = ellipse(cross_m(p) @ (A.T @ C @ A) @ cross_m(p).T, p)
            print(f"   {nm:6s} {t:3d} Ma  C'=C: {o[0]:4.0f} x {o[1]:3.0f} az {o[2]:4.0f} | "
                  f"C'=BCB^T: {n[0]:4.0f} x {n[1]:3.0f} az {n[2]:4.0f}")


if __name__ == "__main__":
    main()
