"""
chang: Python port of the plate-rotation uncertainty method of
Chang (1988) / Hellinger (1979, 1981) / Kirkwood et al. (1999) / Hanna and Chang
(1990), extended to hotspot tracks following O'Neill et al. (2005).

Modules:
    geometry    - quaternions, coordinates, exp/log maps (port of dtrans.f)
    stats_dist  - chi2 / F critical points (replace XIDCH/XIDF)
    fit_uv      - spherical regression for homologous pairs (u,v) (SPHEREREGRESS)
    hellinger   - Hellinger fit by sections (port of hellinger1.f)
    triple      - triple-junction fit, single kappa (port of hellinger3.f)
    region      - confidence regions and ellipses of reconstructed points
    combine     - inversion and combination of rotations (invrot, addrot, average, Welch)
    tests       - identity and difference tests (pval, F-test)
    hotspot     - O'Neill extension for hotspot-track uncertainty
"""
from . import (geometry, stats_dist, fit_uv, hellinger, triple, region,  # noqa: F401
               combine, tests, hotspot)
