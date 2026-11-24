"""
Validation of the Python port DIRECTLY against the original compiled Fortran
(hellinger1.exe and add2.exe), using the reference in ../reference/.

Compares:
1. chang.hellinger.fit_hellinger vs hellinger1.exe (rotation, kappa, H11.2) on dataA/dataB.
2. chang.combine.welch_df / welch_critical_value vs add2.exe (capa/capc -> E,Var,df,c,crit).

The .dat files, the fit files (fitA.txt, fitB.txt) and the add2 log are generated with
../reference/gen_reference.py + the .exe programs (see ../reference/README).
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chang import hellinger, combine, geometry as geo

REF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reference")


def read_dat(path):
    """Reads a Hellinger-format .dat -> list of (iside, isect, lat, lon, sd_km)."""
    out = []
    with open(path) as f:
        f.readline()  # nsect
        for ln in f:
            p = ln.split()
            if len(p) >= 5:
                out.append((int(p[0]), int(p[1]), float(p[2]), float(p[3]), float(p[4])))
    return out


def parse_fit(path):
    """Parses the output file of hellinger1.exe."""
    with open(path) as f:
        lines = [l.rstrip() for l in f]
    d = {}
    for i, l in enumerate(lines):
        if "Fitted rotation--alat" in l:
            d["latlonrho"] = [float(x) for x in lines[i + 1].split()]
        elif "kappahat, degrees of freedom,xchi" in l:
            vals = [float(x) for x in lines[i + 1].split()]
            d["kappa"], d["df"] = vals[0], vals[1]
        elif l.strip() == "H11.2 matrix:":
            d["H112"] = np.array([[float(x) for x in lines[i + 1 + k].split()]
                                  for k in range(3)])
    return d


print("=== 1. hellinger.py vs hellinger1.exe (Fortran) ===")
for tag in ("A", "B"):
    dat = read_dat(os.path.join(REF, f"data{tag}.dat"))
    fit_f = parse_fit(os.path.join(REF, f"fit{tag}.txt"))
    lat0, lon0, rho0 = fit_f["latlonrho"]
    res = hellinger.fit_hellinger(dat, (lat0, lon0, rho0))
    print(f"\n-- dataset {tag} --")
    print(f"  Python rotation : {res['lat']:.4f}, {res['lon']:.4f}, {res['rho']:.4f}")
    print(f"  Fortran rotation: {lat0:.4f}, {lon0:.4f}, {rho0:.4f}")
    print(f"  kappa Python {res['kappa_hat']:.5f} vs Fortran {fit_f['kappa']:.5f}; "
          f"df {res['df']} vs {int(fit_f['df'])}")
    dH = np.max(np.abs(res['H112'] - fit_f['H112'])) / np.max(np.abs(fit_f['H112']))
    print(f"  max|H11.2 diff| / max|H11.2| = {dH:.2e}")

print("\n=== 2. welch_df/welch_critical_value vs add2.exe (Fortran) ===")
# Values of capa (A) and capc (C) printed by add2.exe ('A, C ...'):
capa, capc = 0.0362397676906008, 0.0718547148215389
w = combine.welch_critical_value(capa, 0.0, capc, alpha=0.95)
print("Python (a=%.6f, c=%.6f, b=0):" % (capa, capc))
print(f"  E   = {w['E']:.7f}   Fortran 3.0724795")
print(f"  Var = {w['Var']:.7f}   Fortran 6.6510662")
print(f"  df  = {w['df']:.4f}   Fortran 91.983971")
print(f"  c   = {w['cfact']:.7f}   Fortran 3.0056749")
print(f"  crit= c*F(3,df) = {w['crit']:.6f}   Fortran 8.1261717")
