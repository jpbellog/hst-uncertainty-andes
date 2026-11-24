# Third-party inputs

These published datasets are read by the code but are not redistributed in this repository. Download
each one from its source and place it at the path given (relative to the repository root, `PAPER_ROOT`).
All of them are listed in `.gitignore`.

"Required" inputs are read by the main chain (`run_all.sh`). "Optional" inputs are only needed to
regenerate an intermediate file that is already included in this repository.

## Required

| Dataset | Source | Files | Place at | Read by |
|---|---|---|---|---|
| Nazca–Pacific–Antarctica rotations and covariances, Tables S2.2–S2.4 | Wilson and DeMets (2025), Supporting Information, doi:10.1029/2025JB031933 | Supporting Information PDF, saved as `wilson2025_SI.pdf` | `external/wilson2025/` | `models/WD25_WK08-A/parse_wd25_si.py` (writes `models/WD25_WK08-A/wd25/wd25_S2_*.csv`), then `rotation_models/build_rotation_models.py` |
| GMHRF finite rotations (Table S3) and hotspot motion (Data Set S1) | Doubrovine et al. (2012), Supporting Information, doi:10.1029/2011JB009072 | `jgrb17279-sup-0004-ts03.txt`, `jgrb17279-sup-0008-ds01.txt` | `external/doubrovine2012/Supporting_Information/` | `build_rotation_models.py`, `03b_hotspot_fit.py` (GMHRF-B) |
| GWG20 Pacific hotspot frame (their Table S6) | Gaastra et al. (2022), Supporting Information, doi:10.1029/2021TC006772 | `2021tc006772-sup-0003-table si-s02.csv` | `external/gaastra2022/Supporting_Information/` | `build_rotation_models.py` |
| Nazca–South America stage-vector covariances | Quiero et al. (2022), Supplementary Data, doi:10.1038/s41467-022-29055-4 | `MOESM5.txt` | `external/quiero2022/` | `models/*/05_vc_uncertainty.py` |
| Present-day coastlines and static plate polygons | Müller et al. (2019), plate motion model v2.0, doi:10.1029/2018TC005462 | `Coastlines_USGS_MedRes_AREPS_polyline.shp` (+ `.shx`, `.dbf`), `Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons_2019_v1.shp` (+ companions) | `external/muller2019/` | `figures/make_F1_track_ellipses.py`, `analysis/hellinger_circuit.py` |
| Utrecht paleogeography model | van Hinsbergen et al. (2026), Supporting Information, doi:10.1371/journal.pone.0346817 | `Utrecht_Paleogeography_15Dec25.rot`, `Continental_polygons.gpml`, `gAPWP25.rot` | `external/vanhinsbergen2026/` | `analysis/andean_shortening_utrecht.py` (Andean shortening test) |
| Minimum-continent-motion mantle frame | Wagenaar et al. (2025), Zenodo doi:10.5281/zenodo.17536784 (version 3.0) | `Rotation_files.zip` (unzipped to `Rotation_files/`) | `external/wagenaar2025/` | `analysis/build_rot_wagenaar.py` |
| ETOPO1 bedrock relief | Amante and Eakins (2009), doi:10.7289/V5C8276M | `ETOPO1_Bed_g_geotiff.tif` | `external/etopo1/` | `analysis/bathymetry_axes.py`, `analysis/track_vs_ridge_axis.py` |

## Optional

| Dataset | Source | Place at | Regenerates |
|---|---|---|---|
| Magnetic anomaly picks (GSFML) | Seton et al. (2014), doi:10.1002/2013GC005176 | `external/gsfml/GSFML.global.picks.gmt` | covariances fitted in this work at C32y (`analysis/hellinger_circuit.py`; needs GDAL `ogr2ogr`, set `OSGEO4W_ROOT`) |
| Fracture-zone traces (GSFML) | Matthews et al. (2011), doi:10.1029/2011JB008413; Wessel et al. (2015), doi:10.1002/2015GC005853 | `external/gsfml/SF_shp/` | control plot of the Challenger FZ (`analysis/foundation_t_switch_fz.py`); the trace used is already in `models/*/data/fz_challenger_gsfml.csv` |
| Step-heating data of the Easter chain, EarthRef Digital Archive (ERDA) entries 1305–1342 | Ray et al. (2012), doi:10.1093/petrology/egs021; ERDA at earthref.org [download address to be verified] | `external/erda_ray2012/<entry>.xls` | Easter plateau ages (`analysis/ray2012_erda.py`, `ray2012_erda_ages.py`); the result is `models/*/data/ages_easter_ray2012_recal.csv` |
| Polarity time scales (article PDFs) | Cande and Kent (1995); Ogg (2012); Ogg (2020), doi:10.1016/B978-0-12-824360-2.00005-X | `external/timescale_pdfs/` | chron tables used to move node ages to GTS2020 (`analysis/timescale_*.py`); the tables are included in `analysis/outputs/` |
| Foundation chain ages (article PDFs) | O'Connor et al. (2001), doi:10.1016/S0012-821X(01)00461-7; O'Connor et al. (2004), doi:10.1007/978-3-642-18782-7_12 | `external/oconnor_pdfs/` | `analysis/foundation_oconnor2004.py`; the result is `models/*/data/ages_foundation_oconnor2004_v3.csv` |
