# -*- coding: utf-8 -*-
"""Common timescale. Converts to GTS2012 the ages of the plate 002 nodes taken
from Table S3 "Pacific" of Doubrovine et al. (2012), which are in Cande and Kent
(1995) for the Cenozoic (Methods of the GMHRF).

REMAPPING FUNCTION
Piecewise linear CK95 -> GTS2012 through the boundaries of the normal chrons that
appear in both tables (Cande and Kent 1995, Table 2; Ogg 2012, Table 5.2),
read by ../../analysis/timescale_CK95_GTS2012.py into
outputs/chrons_CK95_vs_GTS2012_full.csv, plus the point (0, 0). After the
last boundary read (base of C33n, 79.075 Ma CK95 = 79.900 Ma GTS2012) that
difference is kept constant. C25n was not read from the PDF; that segment is
interpolated between C24n.3n and C26n.

WHAT IS REMAPPED
Only the rows of plate 002 with "Table S3 Pacific (" in the comment (10-80 Ma,
with covariance; the covariance does not change). These stay as they are:
  - the 90-120 Ma rows "Pacific via WK08 stage rotations" and the WK08-A rows
    older than 120 Ma: the WK08 ages are radiometric, the paper names no
    polarity timescale;
  - the drift of Data Set S1 (version B and PAPER_HS_PRIOR=gmhrf): the timescale
    of the flow model is not verified [TO BE VERIFIED].
The circuit (Quiero et al. 2022 and Wilson and DeMets 2025) is already in GTS2012.

Outputs: rot/Global_Model_R_BG_2025_{TAG}GMHRF_gts12.rot for TAG "" and "WD25_",
and the table of ages in rot/remap_gmhrf_gts2012.log.

Usage:  py -I remap_gmhrf_gts2012.py
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CRONES = os.path.join(ROOT, "analysis", "outputs",
                      "chrons_CK95_vs_GTS2012_full.csv")
TAGS = ("", "WD25_")
KEY = "Table S3 Pacific ("


def remap_function():
    df = pd.read_csv(CRONES)
    ck = np.concatenate([[0.0], df.ck95_joven, df.ck95_vieja])
    gt = np.concatenate([[0.0], df.gts12_joven, df.gts12_vieja])
    o = np.argsort(ck)
    ck, gt = ck[o], gt[o]
    if np.any(np.diff(gt) <= 0):
        raise ValueError("the chron table is not monotonic")

    def f(t):
        t = np.asarray(t, float)
        return np.where(t <= ck[-1], np.interp(t, ck, gt), t + (gt[-1] - ck[-1]))
    return f, (ck[-1], gt[-1])


def main():
    f, last = remap_function()
    log = [f"CK95 -> GTS2012 function: {CRONES}",
           f"last boundary {last[0]:.3f} -> {last[1]:.3f} Ma, constant difference after it"]
    for tag in TAGS:
        src = os.path.join(HERE, "rot", f"Global_Model_R_BG_2025_{tag}GMHRF.rot")
        out = os.path.join(HERE, "rot", f"Global_Model_R_BG_2025_{tag}GMHRF_gts12.rot")
        lines, n = [], 0
        with open(src, encoding="utf-8", errors="replace") as fh:
            for ln in fh:
                tok = ln.split()
                if tok and tok[0] in ("002", "2") and KEY in ln:
                    a = float(tok[1])
                    b = float(f(a))
                    # only the age field changes; the rest of the line stays the same
                    i = ln.index(tok[1])
                    ln = ln[:i] + f"{b:.3f}" + ln[i + len(tok[1]):]
                    ln = ln.replace(KEY, "Table S3 Pacific, age CK95->GTS2012 (", 1)
                    log.append(f"  {tag or 'Quiero_'}GMHRF  002  {a:6.2f} Ma CK95 -> {b:7.3f} Ma GTS2012 "
                               f"({b - a:+.3f})")
                    n += 1
                lines.append(ln)
        with open(out, "w", encoding="utf-8", newline="\n") as fh:
            fh.writelines(lines)
        log.append(f"{out}: {n} nodes remapped")
    with open(os.path.join(HERE, "rot", "remap_gmhrf_gts2012.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
