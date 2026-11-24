"""Common GTS2020 timescale (Ogg 2020, Table 5.3 of chrons).
Collects, per normal chron, the young/old boundaries of CK95 (Cande and Kent 1995),
GTS2012 (Ogg 2012) and GTS2020 (Ogg 2020), and provides the piecewise linear
remapping functions CK95 -> GTS2020 and GTS2012 -> GTS2020.

GTS2020: each row of the PDF gives, after the chron name, the age of the
BASE (Ma, 3 decimals) and the duration (myr, 3 decimals); the top of a normal
chron is base - duration. GTS2012 is reread here (top, base, duration), including
the chrons with a name in parentheses (Brunhes, Olduvai...); CK95 comes
from timescale_CK95_GTS2012.py (chrons_CK95_vs_GTS2012_full.csv).

Output: outputs/chrons_CK95_GTS2012_GTS2020.csv
Usage:  py timescale_GTS2020.py
"""
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pdfplumber

LIB = Path(__import__("os").environ["PAPER_ROOT"] + "/external/timescale_pdfs")
HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
TABLA = OUT / "chrons_CK95_GTS2012_GTS2020.csv"


CRON = r"\b(C\d+[A-Z]*n(?:\.\dn)?)(?:\([^)]*\))?(?=\s)"
NUM = r"(?<![\d.])(\d+\.\d{3})(?![\d])"


def _leer(patron, cols):
    """(young, old) boundaries per normal chron; cols selects the table format."""
    pdf = next(LIB.glob(patron))
    rows = {}
    with pdfplumber.open(pdf) as p:
        for pg in p.pages:
            for line in (pg.extract_text() or "").split("\n"):
                for m in re.finditer(CRON, line):
                    nums = [float(x) for x in re.findall(NUM, line[m.end():])]
                    if m.group(1) in rows:
                        continue
                    if cols == "gts20" and len(nums) >= 2 and nums[1] <= nums[0]:
                        rows[m.group(1)] = (round(nums[0] - nums[1], 3), nums[0])
                    elif cols == "gts12" and len(nums) >= 3:
                        rows[m.group(1)] = (nums[-3], nums[-2])
                    elif cols == "gts12" and m.group(1) == "C1n" and len(nums) == 2:
                        rows["C1n"] = (0.0, nums[0])        # top "0" without decimals
    return rows


def gts20():
    """GTS2020: the chron is followed by base (Ma) and duration (myr)."""
    return _leer("Ogg_2020_*.pdf", "gts20")


def gts12():
    """GTS2012: the last three numbers with 3 decimals are top, base, duration."""
    return _leer("Ogg_2012_*.pdf", "gts12")


def tabla():
    """One row per normal chron present in GTS2020 and in at least one other timescale.
    CK95 from chrons_CK95_vs_GTS2012_full.csv (reading of Cande and Kent 1995)."""
    if TABLA.exists():
        return pd.read_csv(TABLA)
    ck = pd.read_csv(OUT / "chrons_CK95_vs_GTS2012_full.csv").set_index("cron")
    g12, g20 = gts12(), gts20()
    rows = []
    for c, (j20, v20) in g20.items():
        r = dict(cron=c, gts20_joven=j20, gts20_vieja=v20)
        r["gts12_joven"], r["gts12_vieja"] = g12.get(c, (np.nan, np.nan))
        if c in ck.index:
            r["ck95_joven"], r["ck95_vieja"] = ck.loc[c, "ck95_joven"], ck.loc[c, "ck95_vieja"]
        rows.append(r)
    df = pd.DataFrame(rows).sort_values("gts20_joven")
    df = df[["cron", "ck95_joven", "ck95_vieja", "gts12_joven", "gts12_vieja",
             "gts20_joven", "gts20_vieja"]]
    for s in ("joven", "vieja"):
        df[f"d20_12_{s}"] = df[f"gts20_{s}"] - df[f"gts12_{s}"]
        df[f"d20_ck_{s}"] = df[f"gts20_{s}"] - df[f"ck95_{s}"]
    df.to_csv(TABLA, index=False, float_format="%.3f")
    return df


# Top of C34n (end of the Cretaceous superchron), outside the normal-chron
# tables: CK95 83.0 (Cande and Kent 1995); GTS2012 83.64 (Ogg 2012, row C33r);
# GTS2020 83.650 (Ogg 2020, "C33r extension ... assigned basal age").
C34_TECHO = {"ck95": 83.0, "gts12": 83.64, "gts20": 83.65}


def remap(origen):
    """Function t_origen -> t_GTS2020, piecewise linear through the boundaries of the
    common normal chrons, with (0, 0); beyond the last boundary, constant
    difference. origen: "ck95" or "gts12"."""
    df = tabla().dropna(subset=[f"{origen}_joven"])
    x = np.concatenate([[0.0], df[f"{origen}_joven"], df[f"{origen}_vieja"],
                        [C34_TECHO[origen]]])
    y = np.concatenate([[0.0], df.gts20_joven, df.gts20_vieja, [C34_TECHO["gts20"]]])
    k = x > 0
    x, y = np.r_[0.0, x[k]], np.r_[0.0, y[k]]
    o = np.argsort(x)
    x, y = x[o], y[o]
    if np.any(np.diff(x) <= 0) or np.any(np.diff(y) <= 0):
        raise ValueError(f"chron table not monotonic ({origen})")

    def f(t):
        t = np.asarray(t, float)
        return np.where(t <= x[-1], np.interp(t, x, y), t + (y[-1] - x[-1]))
    return f


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    if TABLA.exists():
        TABLA.unlink()
    df = tabla()
    pd.set_option("display.width", 200)
    print(df.round(3).to_string(index=False))
    print("n chrons", len(df))
