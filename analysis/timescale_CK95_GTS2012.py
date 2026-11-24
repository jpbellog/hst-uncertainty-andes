"""Compares the ages of the normal chrons in CK95 (Cande and Kent 1995, Table 2)
and GTS2012 (Ogg 2012, Table 5.2), 0-40 Ma. Reads the PDFs in
external/timescale_pdfs; does not touch the pipeline.

CK95: Table 2 gives normal intervals "young age - old age | chron".
GTS2012: in each row of Table 5.2 the numbers with 3 decimals are
Top(Ma), Base(Ma) and Duration(myr); the distances have 2 decimals.
"""
import re
import sys
from pathlib import Path
import pdfplumber
import pandas as pd

LIB = Path(__import__("os").environ["PAPER_ROOT"] + "/external/timescale_pdfs")
OUT = Path(__file__).resolve().parent / "outputs"
OUT.mkdir(exist_ok=True)
sys.stdout.reconfigure(encoding="utf-8")


def norm(c):
    return c.replace(" ", "").replace("l", "1").rstrip(".")


def ck95():
    pdf = next(LIB.glob("Cande-Kent_1995_*.pdf"))
    rows = {}
    with pdfplumber.open(pdf) as p:
        for pg in p.pages:
            t = pg.extract_text() or ""
            for m in re.finditer(r"(\d+\.\d{3})\s*-\s*(\d+\.\d{3})\s+(C\s?\d+[A-Z]*n(?:\.\s?[l1\d]n)?)(?!\s*-)", t):
                rows.setdefault(norm(m.group(3)), (float(m.group(1)), float(m.group(2))))
    return rows


def gts12():
    pdf = next(LIB.glob("Ogg_2012_*.pdf"))
    rows = {}
    with pdfplumber.open(pdf) as p:
        for pg in p.pages:
            for line in (pg.extract_text() or "").split("\n"):
                m = re.search(r"\b(C\d+[A-Z]*n(?:\.\dn)?)\s", line)
                if not m:
                    continue
                nums = re.findall(r"(?<![\d.])(\d+\.\d{3})(?![\d])", line[m.end():])
                if len(nums) >= 3:
                    rows[m.group(1)] = (float(nums[-3]), float(nums[-2]))
    return rows


a, b = ck95(), gts12()
df = pd.DataFrame([dict(cron=c, ck95_joven=a[c][0], ck95_vieja=a[c][1],
                        gts12_joven=b[c][0], gts12_vieja=b[c][1])
                   for c in a if c in b])
df = df.sort_values("ck95_joven")
# full table (up to C34n, for remapping the GMHRF nodes
# of 10-80 Ma) and the 0-40 Ma table
full = df.assign(d_joven=df.gts12_joven - df.ck95_joven,
                 d_vieja=df.gts12_vieja - df.ck95_vieja)
full.to_csv(OUT / "chrons_CK95_vs_GTS2012_full.csv", index=False, float_format="%.3f")
df = df[df.ck95_joven <= 40]
df["d_joven"] = df.gts12_joven - df.ck95_joven
df["d_vieja"] = df.gts12_vieja - df.ck95_vieja
df.to_csv(OUT / "chrons_CK95_vs_GTS2012.csv", index=False, float_format="%.3f")
print(f"chrons in both tables (0-40 Ma): {len(df)}  (CK95 read {len(a)}, GTS12 read {len(b)})")
print(df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
for lo, hi in ((0, 10), (10, 20), (20, 30), (30, 40)):
    s = df[(df.ck95_joven >= lo) & (df.ck95_joven < hi)]
    d = pd.concat([s.d_joven, s.d_vieja])
    print(f"{lo:2d}-{hi:2d} Ma: n={len(s):2d}  difference GTS12-CK95 min {d.min():+.3f} max {d.max():+.3f} |mean| {d.abs().mean():.3f} Myr")
