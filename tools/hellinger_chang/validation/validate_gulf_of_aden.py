"""
Validation of fit_uv.py against the Gulf of Aden example of Hanna and Chang (1990),
Table 1 and output of the SPHEREREGRESS program (pp. 165-167).

Data: 11 pairs (u on the Somali plate, v on the Arabian plate). Expected result from the paper:
    axis lat 25.31, lon 24.29, angle 2.38 deg, SSE = 0.12787467e-04
    fitted rotation matrix:
        [[ 0.99972, -0.01747,  0.01574],
         [ 0.01800,  0.99926, -0.03405],
         [-0.01513,  0.03433,  0.99930]]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from chang import fit_uv

# Table 1 of Hanna and Chang (1990): (lat, lon) in degrees
u = [(13.05, 57.56), (13.34, 57.07), (13.89, 56.50), (14.19, 55.96),
     (14.10, 55.92), (14.21, 55.38), (12.68, 50.95), (11.97, 47.56),
     (12.06, 47.35), (11.63, 45.80), (11.73, 45.36)]
v = [(14.28, 58.12), (14.54, 57.67), (15.00, 57.16), (15.33, 56.51),
     (15.25, 56.48), (15.37, 55.93), (13.59, 51.51), (12.78, 48.11),
     (12.86, 47.89), (12.44, 46.39), (12.58, 45.87)]

A_paper = np.array([[0.99972, -0.01747, 0.01574],
                    [0.01800, 0.99926, -0.03405],
                    [-0.01513, 0.03433, 0.99930]])

print("=== Validation of fit_uv vs Gulf of Aden (Hanna and Chang 1990) ===")
res = fit_uv.fit_uv(u, v)
print(f"fitted axis  : lat {res['lat']:.2f}, lon {res['lon']:.2f}, "
      f"angle {res['rho']:.2f} deg")
print("expected     : lat 25.31, lon 24.29, angle 2.38 deg")
print(f"fitted SSE   : {res['sse']:.8e}   expected 1.2787467e-05")
print("\nfitted matrix A:\n", np.round(res['A'], 5))
print("\nmax|A - A_paper|     =", np.max(np.abs(res['A'] - A_paper)))
print("max|A^T - A_paper|   =", np.max(np.abs(res['A'].T - A_paper)),
      "(in case the direction convention is reversed)")
print("\nmatrix Q (x1e4) to compare structure:\n", np.round(res['Q'] / 1e4, 3))
