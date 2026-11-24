# -*- coding: utf-8 -*-
"""Crop of ETOPO1 (Bedrock, grid registered, 1 arc-minute) used as the background of Fig. 1.

Reads the uncompressed global GeoTIFF with a memmap (int16, one row per strip, origin
-180.00833/90.00833, step 1') and saves the box of the figure at a 2' step in
figures/data/etopo1_F1.npz (lon, lat, z in m), so the figure does not need the 466 MB file.

Input: external/etopo1/ETOPO1_Bed_g_geotiff.tif (see THIRD_PARTY.md).
Usage:  py -I extract_etopo1_F1.py
"""
import os
import numpy as np

ETOPO = os.path.join(os.environ["PAPER_ROOT"], "external", "etopo1", "ETOPO1_Bed_g_geotiff.tif")
NX, NY, OFF0 = 21601, 10801, 87470
LON0, LAT0, D = -180.0, 90.0, 1 / 60.0
BOX = (-116.0, -51.0, -46.0, 6.0)    # extent of Fig. 1 plus a margin
STEP = 2
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "etopo1_F1.npz")


def main():
    with open(ETOPO, "rb") as f:
        assert f.read(2) == b"II"
    z = np.memmap(ETOPO, dtype="<i2", mode="r", offset=OFF0, shape=(NY, NX))
    i0, i1 = int(round((BOX[0] - LON0) / D)), int(round((BOX[1] - LON0) / D))
    j0, j1 = int(round((LAT0 - BOX[3]) / D)), int(round((LAT0 - BOX[2]) / D))
    sub = np.array(z[j0:j1 + 1:STEP, i0:i1 + 1:STEP], dtype=np.int16)
    lon = LON0 + D * np.arange(i0, i1 + 1, STEP)
    lat = LAT0 - D * np.arange(j0, j1 + 1, STEP)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez_compressed(OUT, lon=lon, lat=lat, z=sub)
    print(f"{OUT}: {sub.shape}, z {sub.min()} to {sub.max()} m")
    # check: summit of Aconcagua and trench off Antofagasta, of the expected order
    for nm, la, lo in (("Aconcagua", -32.65, -70.01), ("trench 23.5S", -23.5, -71.4)):
        print(f"  {nm}: {sub[int(round((lat[0] - la) / (D * STEP))), int(round((lo - lon[0]) / (D * STEP)))]} m")


if __name__ == "__main__":
    main()
