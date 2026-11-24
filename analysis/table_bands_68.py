# -*- coding: utf-8 -*-
"""68 % rows of Table 2 of the text: at each age
of the table, union of the 68 % intervals (p16--p84) of the nine models with f_PAC = 0.44,
clipped by the frame_coverage of each frame, together with the median of the medians and the 95 % envelope of
all models and f_PAC (the same as in manuscript_numbers.sec_tabla1).
Usage:  py table_bands_68.py   (output: outputs/table_bands_68.csv)"""
import os
import pandas as pd
import manuscript_numbers as C

d = C.bandas().dropna(subset=["lat_p50"])
a = C.bandas(pf=None).dropna(subset=["lat_p50"])
rows = []
for r in ("Easter", "TT", "JFR", "Foundation"):
    for age in range(0, 61, 10):
        s, t = d[(d.ridge == r) & (d.age == age)], a[(a.ridge == r) & (a.age == age)]
        if not len(s):
            continue
        rows.append(dict(cadena=C.CADENAS[r], edad=age, n=s.modelo.nunique(), mediana=s.lat_p50.median(),
                         lo68=s.lat_p16.min(), hi68=s.lat_p84.max(), lo95=t.lat_p2_5.min(), hi95=t.lat_p97_5.max()))
T = pd.DataFrame(rows)
T.to_csv(os.path.join(C.SAL, "table_bands_68.csv"), index=False, float_format="%.1f")
print(T.round(1).to_string(index=False))
