# -*- coding: utf-8 -*-
"""Merges the step 04 parts of a pipeline folder
(one process per variant x f_PAC, PAPER_SUFFIX = <tag>_<var3><pf>)
into mc_collision_bands_{variant}{tag}.csv. Same procedure as
merge_mc_parts.py, with the suffix of option (b) ("_hsb" or "_hsb_pred").

Usage:  py merge_mc_parts_tag.py <folder> <tag> [comma-separated variants]
"""
import os
import sys
import pandas as pd

d, tag = sys.argv[1], sys.argv[2]
variants = sys.argv[3].split(",") if len(sys.argv) > 3 else \
    ["total", "circuit", "wk08", "hotspot"]
ok = True
for v in variants:
    fns = [os.path.join(d, f"mc_collision_bands_{v}{tag}_{v[:3]}{pf.replace('.', '')}.csv")
           for pf in ("0.5", "0.44", "0.36")]
    missing = [os.path.basename(f) for f in fns if not os.path.exists(f)]
    if missing:
        print(f"[{v}] MISSING: {missing}")
        ok = False
        continue
    df = pd.concat([pd.read_csv(f, dtype={"pf": str}) for f in fns], ignore_index=True)
    df = df.sort_values(["ridge", "pf", "age"]).reset_index(drop=True)
    out = os.path.join(d, f"mc_collision_bands_{v}{tag}.csv")
    df.to_csv(out, index=False)
    print(f"[{v}] written {out}: {len(df)} rows, pf = {sorted(df.pf.unique())}, "
          f"median N = {int(df.n.median())}")
sys.exit(0 if ok else 1)
