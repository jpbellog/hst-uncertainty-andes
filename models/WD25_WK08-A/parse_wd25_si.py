# -*- coding: utf-8 -*-
"""Extracts the finite rotations and their covariances from the Supporting Information of
Wilson and DeMets (2025), JGR Solid Earth 130:e2025JB031933.

Source: the Supporting Information PDF of that article
  (Wilson-DeMets_2025_Changes-motion-Nazca-Farallon-34Ma-flat-slab_JGR_SI.pdf).

Tables:
  S2.1 Antarctic-Nazca   (3 plates)  p10
  S2.2 Pacific-Nazca     (3 plates)  p11
  S2.3 Antarctic-Pacific (3 plates)  p12
  S2.4 Pacific-Nazca     (2 plates)  p13
  S2.5 Nazca-South America, totals  p14
  S2.6 Nazca-South America, stages   p15   (a-f in rad^2/Myr^2)

Covariance convention (footnote of Table S2.1, with `g` as a per-row column):
    C = (1/k) * [[a,b,c],[b,d,e],[c,e,f]] * 10^(-g)   [rad^2]
Tables S2.5 and S2.6 have NO k column (k = 1 is used).
The rotation restores the FIRST plate of the title to its past position relative to
the SECOND, right-hand rule (counterclockwise positive).

Output: wd25/wd25_S2_{1..6}.csv (next to this script) + console summary.
Usage:  py parse_wd25_si.py
"""
import os
import re
import numpy as np
import pandas as pd
import pypdf

SI = (__import__("os").environ["PAPER_ROOT"] + "/external/wilson2025/wilson2025_SI.pdf")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "wd25")

TABLES = {
    "S2.1": dict(page=9,  k=True,  nsnp=True,  desc="Antarctic-Nazca, 3-plate"),
    "S2.2": dict(page=10, k=True,  nsnp=True,  desc="Pacific-Nazca, 3-plate"),
    "S2.3": dict(page=11, k=True,  nsnp=True,  desc="Antarctic-Pacific, 3-plate"),
    "S2.4": dict(page=12, k=True,  nsnp=True,  desc="Pacific-Nazca, 2-plate"),
    "S2.5": dict(page=13, k=False, nsnp=False, desc="Nazca-South America, total"),
    "S2.6": dict(page=14, k=False, nsnp=False, desc="Nazca-South America, stage"),
}

NUM = r"-?\d+\.?\d*"
# The AGU PDF splits each minus sign into three characters: "-" + U+00AD + U+2010.
# The last two are soft hyphens and are REMOVED; converting them to "-" would leave
# double minus signs and break the parser.
SOFT = ("\u00ad", "\u2010", "\u2011")
DASH = ("\u2013", "\u2014", "\u2212")


def page_text(page):
    """Page text, without soft hyphens and with collapsed whitespace."""
    t = page.extract_text() or ""
    for junk in SOFT:
        t = t.replace(junk, "")
    for dash in DASH:
        t = t.replace(dash, "-")
    t = t.replace("\r", " ").replace("\u00a0", " ")
    t = re.sub(r"\s+", " ", t)
    if "--" in t:
        print("  [warning] double minus signs remain: check the normalisation")
    return t


def parse_rows(txt, has_k, has_nsnp, stage=False):
    """Data rows. Each row starts with an isochron ID (C1o, C2Ay,
    C4n.1r*, ...), except in S2.6, which starts with two ages."""
    rows = []
    if stage:
        # S2.6: age_old age_young lat lon rate a b c d e f g
        pat = re.compile(r"(?<![\w.])" + r"\s+".join(["(" + NUM + ")"] * 12)
                         + r"(?!\d)")
        for m in pat.finditer(txt):
            v = [float(x) for x in m.groups()]
            rows.append(dict(age_old=v[0], age_young=v[1], lat=v[2], lon=v[3],
                             rate=v[4], a=v[5], b=v[6], c=v[7], d=v[8], e=v[9],
                             f=v[10], g=v[11], k=1.0))
        return rows
    n = 4 + 7 + (1 if has_k else 0) + (4 if has_nsnp else 0)
    pat = re.compile(r"(C[0-9A-Za-z.]+\*?)\s+"
                     + r"\s+".join(["(" + NUM + ")"] * n) + r"(?!\d)")
    for m in pat.finditer(txt):
        v = [float(x) for x in m.groups()[1:]]
        i = 0
        age, lat, lon, ang = v[i:i + 4]; i += 4
        k = v[i] if has_k else 1.0
        if has_k:
            i += 1
        a, b, c, d, e, f, g = v[i:i + 7]; i += 7
        row = dict(isochron=m.group(1), age=age, lat=lat, lon=lon, angle=ang,
                   k=k, a=a, b=b, c=c, d=d, e=e, f=f, g=g)
        if has_nsnp:
            row.update(dict(ns_iso=v[i], np_iso=v[i + 1],
                            ns_fz=v[i + 2], np_fz=v[i + 3]))
        rows.append(row)
    return rows


def cov_matrix(r):
    """3x3 covariance in rad^2 (rad^2/Myr^2 in S2.6)."""
    s = 10.0 ** (-r["g"]) / r["k"]
    return np.array([[r["a"], r["b"], r["c"]],
                     [r["b"], r["d"], r["e"]],
                     [r["c"], r["e"], r["f"]]]) * s


def check(df, agecol):
    """Consistency checks; returns the list of problems found."""
    bad = []
    if not df[agecol].between(0, 40).all():
        bad.append("ages outside 0-40 Ma")
    if not df["lat"].between(-90, 90).all():
        bad.append("latitudes out of range")
    if not df["lon"].between(-180, 180).all():
        bad.append("longitudes out of range")
    npd = sum(1 for _, r in df.iterrows()
              if np.linalg.eigvalsh(cov_matrix(r)).min() <= 0)
    if npd:
        bad.append(f"{npd} covariances not positive definite")
    if not df[agecol].is_monotonic_increasing:
        bad.append("ages not monotonic")
    return bad


def main():
    os.makedirs(OUT, exist_ok=True)
    rd = pypdf.PdfReader(SI)
    for name, spec in TABLES.items():
        txt = page_text(rd.pages[spec["page"]])
        rows = parse_rows(txt, spec["k"], spec["nsnp"], stage=(name == "S2.6"))
        if not rows:
            print(f"[{name}] NO ROWS -- check the parser")
            continue
        df = pd.DataFrame(rows)
        agecol = "age_old" if name == "S2.6" else "age"
        if name == "S2.6":
            df = df.sort_values(agecol).reset_index(drop=True)
        problems = check(df, agecol)
        fn = os.path.join(OUT, f"wd25_{name.replace('.', '_')}.csv")
        df.to_csv(fn, index=False)
        flag = "OK" if not problems else "CHECK: " + "; ".join(problems)
        print(f"[{name}] {spec['desc']:24s} {len(df):3d} rows | "
              f"{df[agecol].min():6.3f}-{df[agecol].max():6.3f} Ma | {flag}")


if __name__ == "__main__":
    main()
