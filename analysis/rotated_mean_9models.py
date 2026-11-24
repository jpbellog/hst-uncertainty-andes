# -*- coding: utf-8 -*-
"""Runs rotated_mean.py in the nine models of the text with the .rot of each folder
(sens9_summary.ROT) and collects the outputs in outputs/rotated_mean_9models.txt (read by manuscript_numbers.py).
Usage:  py rotated_mean_9models.py"""
import os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "figures"))
import paper_models as MP  # noqa: E402
import sens9_summary as S9  # noqa: E402
out = []
for c in MP.CARPETAS:
    r = subprocess.run([sys.executable, os.path.join(HERE, "rotated_mean.py"), os.path.join(ROOT, "models", c),
                        os.path.join(ROOT, "models", c, "rot", S9.ROT[c])], check=True, capture_output=True, text=True,
                       env=dict(os.environ, PYTHONIOENCODING="utf-8"), encoding="utf-8")
    out.append(f"=== {c}\n{r.stdout.rstrip()}\n")
with open(os.path.join(HERE, "outputs", "rotated_mean_9models.txt"), "w", encoding="utf-8") as f:
    f.write("".join(out))
print("".join(out))
