# -*- coding: utf-8 -*-
"""Builds a version of the rotation model with the South Pacific branches
replaced by the rotations of Wilson and DeMets (2025).

Input:   Global_Model_R_BG_2025.rot  (base model, not modified)
Output:  models/WD25_WK08-A/rot/Global_Model_R_BG_2025_WD25.rot

WHAT CHANGES
The two circuit branches that in the base model come from Quiero et al. (2022)
and for which Wilson and DeMets publish an equivalent:

  plate 911 (NAZ-PAC, fixed 901), 0 < t <= 28.278 Ma (chron C10o)
      <- Tables S2.2 (3-plate solution, 0.781-15.140 Ma) and S2.4 (2-plate
         solution, 22.342-28.278 Ma), plus the two nodes of the 15-22 Ma gap
         (C5Cy and C6y), which appear only in Table 1 of the paper and are
         therefore WITHOUT covariance.

  plate 901 (PAC-WANT, fixed 804), 0 < t <= 15.140 Ma
      <- Table S2.3 (Antarctic-Pacific, 3-plate solution).
         Replacing 911 alone would leave the rest of the South Pacific circuit
         on Quiero et al. (2022), so the "WD25" model would be an undeclared
         hybrid.

WHAT DOES NOT CHANGE
  - The identity nodes at 0 Ma.
  - Plate 901 between 15.140 and 28.278 Ma: Table S2.3 ends at chron
    C5Bo (15.140 Ma) and there is no published equivalent further back, so those
    nodes remain from Quiero et al. (2022). The model is a HYBRID in that
    interval and this must be declared.
  - Plate 902 (FAR-PAC) between 21.767 and 28.278 Ma, also from Quiero et al.
    (2022): Wilson and DeMets do not publish Farallon-Pacific.
  - Everything older than the cut-offs: the structure of Rowan and Rowley (2014) and the
    spreading-asymmetry models f_PAC (plates 902xxx / 9020x).
  - The rest of the circuit (804, 802, 201), the hotspot frame (002) and the
    hotspot plates (13, 14, 15), which are refitted separately with
    refit_hotspots.py.

CONVENTIONS
  - Wilson and DeMets title their tables "A-B" for the rotation that restores
    plate A to its past position relative to B. The .rot needs in both cases
    the OPPOSITE rotation:
        S2.2/S2.4 "Pacific-Nazca"     -> the .rot needs Nazca rel Pacific (911)
        S2.3      "Antarctic-Pacific" -> the .rot needs Pacific rel WANT (901)
    It is the inverse rotation: same axis, negated angle, written here as the
    antipodal pole with positive angle, which is the style of the existing lines.
  - The source covariance goes UNTRANSFORMED to the line of the inverse.
    Wilson and DeMets estimate the covariances following Chang (1988), right
    convention: R_true = R exp(M(h)). The inverse is exp(M(-h)) R^T, which in the
    left convention of the engine (02_rot_engine) has covariance C.
    Writing A^T C A, with A the original rotation, would be correct only for a
    source in the left convention.
  - The source covariance is  C = (1/k) [[a,b,c],[b,d,e],[c,e,f]] 10^(-g)
    in rad^2 (footnote of Table S2.1; `g` is a column of each row).

Usage:  py build_rot_wd25.py
"""
import os
import sys
import importlib
import io
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
eng = importlib.import_module("02_rot_engine")
parse = importlib.import_module("parse_wd25_si")

SRC = (__import__("os").environ["PAPER_ROOT"] + "/inputs/rotations/Global_Model_R_BG_2025.rot")
OUTDIR = os.path.join(HERE, "rot")
OUT = os.path.join(OUTDIR, "Global_Model_R_BG_2025_WD25.rot")

# nodes of the 15-22 Ma gap of the Nazca-Pacific branch, Table 1 of the paper
# (Pacific-Nazca, without covariance in tables S2.x); checked against the main
# PDF, p. 7
GAP_NODES_911 = [("C5Cy", 15.974, 66.72, -92.19, 24.522),
                 ("C6y",  18.715, 68.74, -90.44, 29.943)]

# (plate, fixed, tag, source tables, cutover, nodes without covariance)
BRANCHES = [
    dict(plate=911, fixed=901, tag="NAZ-PAC",
         tables=["wd25_S2_2.csv", "wd25_S2_4.csv"], cutover=33.705,
         gap=GAP_NODES_911),
    dict(plate=901, fixed=804, tag="PAC-WANT",
         tables=["wd25_S2_3.csv"], cutover=15.140, gap=[]),
]

# CROSSOVER Nazca/Farallon.
# Wilson and DeMets (2025) treat Nazca and Farallon as ONE single plate over
# the last 34 Myr (this is the title of their paper). The base model instead
# separates the two plates at 21.767 Ma and from there to 28.278 Ma uses the
# Farallon-Pacific branch of Quiero et al. (2022). For the WD25 model the crossover
# must be where their solution ends, 33.705 Ma (chron C13o).
#
# Procedure (GPlates, Tutorial 2.2 "Changing Rotations, Equivalent Finite
# Rotations, and Cross-Overs", Shephard et al., EarthByte): at the crossover age
# the file carries TWO rows for the moving plate, one with the young fixed
# plate and one with the new one, and the second contains the EQUIVALENT rotation
# of the moving plate relative to the new fixed plate (in GPlates: that plate is
# set as anchored plate and the equivalent rotation is read). Here the moving
# plate 902 is attached to 911 by identity up to the crossover, so its
# equivalent rotation relative to 901 is exactly the row of 911
# relative to 901 at that age. The TOTAL rotation of 902 is thus continuous.
CROSSOVER = dict(plate=902, young_fixed=911, old_fixed=901, age=33.705,
                 tag="FAR-PAC")


def inverse_pole(lat, lon, ang):
    """Antipodal pole and positive angle = inverse rotation."""
    return -lat, (lon + 180.0 + 180.0) % 360.0 - 180.0, ang


def branch_lines(br):
    """New .rot lines of a branch, sorted by age, and their source table."""
    src = pd.concat([pd.read_csv(os.path.join(HERE, "wd25", t))
                     for t in br["tables"]], ignore_index=True)
    src = src[src.age <= br["cutover"] + 1e-9].sort_values("age")
    src = src.reset_index(drop=True)
    lines = []
    for _, r in src.iterrows():
        C = parse.cov_matrix(r)        # Chang (1988), right convention: the inverse carries C
        la, lo, an = inverse_pole(r.lat, r.lon, r.angle)
        cov = "   ".join(f"{v:.4E}" for v in
                         (C[0, 0], C[0, 1], C[0, 2], C[1, 1], C[1, 2], C[2, 2]))
        lines.append((float(r.age),
                      f"{br['plate']}   {r.age:.3f}   {la:.6f}   {lo:.6f}   "
                      f"{an:.7f}   {br['fixed']}   !{br['tag']} "
                      f"@REF Wilson & DeMets 2025 "
                      f'@DOI"10.1029/2025JB031933" @CHRONID"{r.isochron}"   {cov}\n'))
    for iso, age, lat, lon, ang in br["gap"]:
        la, lo, an = inverse_pole(lat, lon, ang)
        lines.append((age,
                      f"{br['plate']}   {age:.3f}   {la:.6f}   {lo:.6f}   "
                      f"{an:.7f}   {br['fixed']}   !{br['tag']} "
                      f"@REF Wilson & DeMets 2025 Table 1 "
                      f'@DOI"10.1029/2025JB031933" @CHRONID"{iso}" @NOCOV\n'))
    lines.sort(key=lambda x: x[0])
    return lines, src


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    built = {br["plate"]: branch_lines(br) for br in BRANCHES}
    by_plate = {br["plate"]: br for br in BRANCHES}

    # crossover row: that of 911 relative to 901 at the crossover age
    X = CROSSOVER
    xline = [l for a, l in built[911][0] if abs(a - X["age"]) < 1e-9]
    assert len(xline) == 1, f"node 911 at {X['age']} Ma is missing for the crossover"
    xf = xline[0].split()
    x_pole = (float(xf[2]), float(xf[3]), float(xf[4]))

    out = []
    dropped = {p: 0 for p in by_plate}
    inserted = {p: False for p in by_plate}
    x_dropped, x_done = 0, False
    for line in io.open(SRC, encoding="latin-1"):
        parts = line.split()
        # ---- crossover plate (902): its young segment is rebuilt ----
        if len(parts) >= 6 and parts[0].isdigit() and int(parts[0]) == X["plate"]:
            age = float(parts[1]); fixed = int(parts[5].split("!")[0])
            drop = ((fixed == X["young_fixed"] and age > 0.0)
                    or (fixed == X["old_fixed"] and age <= X["age"] + 1e-9))
            if drop:
                x_dropped += 1
                if not x_done:
                    x_done = True
                    # identity relative to the young plate up to the crossover
                    out.append(f"{X['plate']}   {X['age']:.3f}   90.000000   "
                               f"0.000000   0.0000000   {X['young_fixed']}   "
                               f"!{X['tag']} crossover Nazca/Farallon movido a "
                               f"{X['age']} Ma (Wilson & DeMets 2025 tratan las "
                               f"dos placas como una en 0-34 Ma)\n")
                    # second row at the SAME age, now relative to the old plate:
                    # equivalent rotation = that of 911 rel 901 (GPlates Tut. 2.2)
                    out.append(f"{X['plate']}   {X['age']:.3f}   {x_pole[0]:.6f}   "
                               f"{x_pole[1]:.6f}   {x_pole[2]:.7f}   "
                               f"{X['old_fixed']}   !{X['tag']} crossover, rotacion "
                               f"equivalente rel {X['old_fixed']} = fila de la 911 "
                               f"(Wilson & DeMets 2025, Tabla S2.4 C13o)\n")
                continue
        if len(parts) >= 6 and parts[0].isdigit() and int(parts[0]) in by_plate:
            plate = int(parts[0])
            br = by_plate[plate]
            age = float(parts[1])
            # EVERY node of the plate within the cutover is deleted, regardless of
            # the fixed plate. This is deliberate: 911 in the base model has a
            # HANDOVER node at 21.767 Ma (identity relative to 902) that marked the
            # end of the Nazca-Pacific branch of Quiero. Since the Wilson and
            # DeMets branch reaches 28.278 Ma, that node falls INSIDE the new
            # segment, and finite_rotation() sorts by age without separating by fixed
            # plate: when crossing a change of fixed plate it snaps to the nearest
            # node, so at t ~ 21 Ma it would take the identity and Nazca
            # would equal Farallon. That would produce a 5-7 degree jump in
            # the collision latitude of the THREE tracks at that age only.
            if 0.0 < age <= br["cutover"] + 1e-9:
                dropped[plate] += 1
                if not inserted[plate]:
                    out.extend(l for _, l in built[plate][0])
                    inserted[plate] = True
                continue
        out.append(line)
    for p in by_plate:
        assert inserted[p], f"no node {p} found in the range"

    assert x_done, "no row of the crossover plate found"
    io.open(OUT, "w", encoding="latin-1", newline="\n").writelines(out)
    print(f"crossover plate {X['plate']}: {x_dropped} old rows removed, "
          f"2 written at {X['age']} Ma (rel {X['young_fixed']} and rel "
          f"{X['old_fixed']})")
    for br in BRANCHES:
        p = br["plate"]
        lines, src = built[p]
        print(f"plate {p} ({br['tag']}, fixed {br['fixed']}): "
              f"{dropped[p]} nodes of Quiero et al. (2022) removed "
              f"(0 < t <= {br['cutover']} Ma), {len(lines)} inserted "
              f"({len(src)} with covariance + {len(br['gap'])} without)")
    print(f"wrote: {OUT}")

    # --- verification: re-read the new model and check node by node ---
    M = eng.RotModel(OUT)
    for br in BRANCHES:
        p, lines, src = br["plate"], *built[br["plate"]]
        got = [r for r in M.table[p]
               if 0 < r[0] <= br["cutover"] + 1e-9 and r[4] == br["fixed"]]
        assert len(got) == len(lines), (p, len(got), len(lines))
        by_age = {round(float(a), 3): (la, lo, an)
                  for (a, la, lo, an, fx, C) in got}
        worst = 0.0
        for _, r in src.iterrows():                # only those with covariance
            la, lo, an = by_age[round(float(r.age), 3)]
            R_new = eng.quat_to_matrix(eng.quat_from_pole(la, lo, an))
            R_wd = eng.quat_to_matrix(eng.quat_from_pole(r.lat, r.lon, r.angle))
            # the written node must be the INVERSE of the published rotation
            worst = max(worst, np.abs(R_new - R_wd.T).max())
        ncov = sum(1 for r in got if r[5] is not None)
        print(f"plate {p}: {len(src)} nodes verified | "
              f"max |R_new - R_WD^T| = {worst:.2e} | "
              f"with covariance read by the parser: {ncov}")

        # GUARD: within the replaced segment the fixed plate must be
        # UNIQUE. finite_rotation() sorts by age without separating by fixed plate and,
        # when crossing a change, snaps to the nearest node; a handover node
        # left inside the segment produces a discontinuity of several degrees
        # at a single age (see the comment of the deletion filter).
        inside = [r for r in M.table[p] if 0.0 < r[0] <= br["cutover"] + 1e-9]
        fixed_set = sorted({r[4] for r in inside})
        assert fixed_set == [br["fixed"]], (
            f"plate {p}: within 0-{br['cutover']} Ma there are nodes with fixed "
            f"plates {fixed_set}, expected only [{br['fixed']}]")
        print(f"plate {p}: guard OK, single fixed plate {fixed_set} in the segment")

    # Crossover GUARD: the TOTAL rotation of the moving plate must be
    # continuous across the crossover age (this is the condition GPlates
    # imposes). It is measured as the displacement of a point on the plate.
    # It must be measured with epsilon -> 0. With a finite epsilon the difference
    # is the plate MOTION over 2*eps (here ~101 km/Myr); it scales
    # linearly with eps and says nothing about continuity.
    pt = eng.latlon_axis(-25.0, -95.0)
    eps = 1e-7
    rows = []
    for t in (X["age"] - eps, X["age"] + eps):
        q = M.total_rotation(X["plate"], t, 0)
        rows.append(eng.xyz_to_latlon(eng.quat_to_matrix(q) @ pt))
    jump = float(eng.gc_km(rows[0][0], rows[0][1], rows[1][0], rows[1][1]))
    print(f"crossover plate {X['plate']} at {X['age']} Ma: jump of the total "
          f"rotation (eps={eps:.0e}) = {jump:.8f} km")
    assert jump < 1e-4, f"the crossover is not continuous: jump of {jump:.5f} km"

    # and the algebraic check, which is exact: the crossover row rel the
    # old plate must be identical to that of 911 rel 901 at that age
    got = [r[1:4] for r in M.table[X["plate"]]
           if abs(r[0] - X["age"]) < 1e-9 and r[4] == X["old_fixed"]]
    ref = [r[1:4] for r in M.table[911]
           if abs(r[0] - X["age"]) < 1e-9 and r[4] == 901]
    assert len(got) == 1 and len(ref) == 1 and np.allclose(got[0], ref[0]), (
        f"the crossover row does not match that of 911: {got} vs {ref}")
    print(f"crossover: row rel {X['old_fixed']} identical to that of 911 "
          f"{tuple(round(v, 6) for v in ref[0])}")

    # and the moving plate must follow 911 just before the crossover
    q1 = eng.quat_to_matrix(M.total_rotation(X["plate"], X["age"] - eps, 0)) @ pt
    q2 = eng.quat_to_matrix(M.total_rotation(911, X["age"] - eps, 0)) @ pt
    la1, lo1 = eng.xyz_to_latlon(q1); la2, lo2 = eng.xyz_to_latlon(q2)
    d = float(eng.gc_km(la1, lo1, la2, lo2))
    print(f"plates {X['plate']} and 911 just before the crossover: {d:.6f} km "
          f"(they must coincide, they are the same plate)")
    assert d < 0.01, f"902 and 911 do not coincide before the crossover: {d:.3f} km"


if __name__ == "__main__":
    main()
