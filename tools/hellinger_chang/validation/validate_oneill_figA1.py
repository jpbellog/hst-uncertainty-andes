"""
Validation of geometry.great_circle_fit against Figure A1 (panel b) of
O'Neill et al. (2005): position matrix of a digitized section of the Ninetyeast
Ridge, with its eigenvectors (S) and eigenvalues (Lambda).

This is the CENTRAL operation of the Hellinger / O'Neill criterion: the normal to the great
circle of a section is the eigenvector of the smallest eigenvalue of the position matrix.

Values supplied by the author (read from the figure of the paper; not checked
directly against the PDF).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from chang import geometry as geo

P = np.array([
    [0.000580234, 0.014465508, -0.00475319991],
    [0.014465508, 0.99106697, -0.0879223599],
    [-0.00475319991, -0.0879223599, 0.0431457396],
])

# from the paper (v1, v2, v3) and their eigenvalues
evecs_paper = {
    "v1": np.array([0.995191281, -0.0058517261, 0.097775623]),
    "v2": np.array([-0.096817352, 0.092635994, 0.990981823]),
    "v3": np.array([0.0148564926, 0.995682846, -0.0916239853]),
}
latlon_paper = {"v1": (5.611, -0.337), "v2": (82.300, 136.264), "v3": (-5.257, 89.145)}
evals_paper = [2.818e-05, 0.035, 0.999]  # "rounded"

print("=== Validation of great_circle_fit vs O'Neill Fig. A1b (Ninetyeast Ridge) ===")
print("Symmetric position matrix:", np.allclose(P, P.T))
print("trace(P) =", np.round(np.trace(P), 5),
      " sum of paper eigenvalues =", round(sum(evals_paper), 5), "(must agree)")

fit = geo.great_circle_fit(P)
ev = fit["eigenvalues"]
print(f"\neigenvalues (asc): {ev[0]:.4e}, {ev[1]:.4f}, {ev[2]:.4f}")
print(f"paper (asc)      : {evals_paper[0]:.4e}, {evals_paper[1]:.4f}, {evals_paper[2]:.4f}")


def align(a, b):
    """aligns the sign of a with that of b and returns max|diff|."""
    a = a if (a @ b) >= 0 else -a
    return a, np.max(np.abs(a - b))


print("\nEigenvectors (sign-aligned):")
# the normal (smallest eigenvalue) must correspond to S1; the paper lists v1 with lat 5.6 (almost
# equatorial) -> each computed eigenvector is compared with the closest one of the paper.
computed = {"normal (min)": fit["normal"], "along (mid)": fit["along"],
            "center (max)": fit["center"]}
for name, vec in computed.items():
    lat, lon = geo.xyz_to_latlon(vec)
    # find the closest v of the paper
    best, bestd = None, 1e9
    for k, vp in evecs_paper.items():
        _, d = align(vec.copy(), vp)
        if d < bestd:
            bestd, best = d, k
    va, d = align(vec.copy(), evecs_paper[best])
    lap, lop = latlon_paper[best]
    print(f"  {name:14s} ~ paper {best}: max|vec diff|={d:.2e}  "
          f"lat/lon calc ({lat:.3f},{lon:.3f}) vs paper ({lap},{lop})")

print("\n=> The NORMAL to the great circle (eigenvector of the smallest eigenvalue) is what the")
print("   Hellinger/O'Neill criterion uses as eta_i. Agreement with the paper validates")
print("   the central operation of the fit (used in hellinger.py and triple.py).")
