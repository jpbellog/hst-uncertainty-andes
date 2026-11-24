# Reproduction of the article results

Each number, table and figure of the article and its supplementary material comes from the
published inputs and the saved scripts, with no ad hoc commands.

The final check is `analysis/manuscript_numbers.py`: it recomputes each number of the text from the
saved outputs and compares it with the value printed in the article. A number that is not there is
not backed by the repository.

## Models

Nine in the main text (3 circuits × 3 frames) and three in the supplementary material (GMHRF-A').
Each one has its folder under `models/`; the list with the `.rot` and the prior of each one is in
`MODELOS=` of `analysis/run_models.sh` and in `analysis/sens9_summary.py` (`ROT`).

| Circuit \ frame | WK08-A | GMHRF-B | GWG20 | GMHRF-A' (suppl.) |
|---|---|---|---|---|
| Q22 (Quiero et al. 2022) | `Q22_WK08-A` | `Q22_GMHRF-B` | `Q22_GWG20` | `Q22_GMHRF-A` |
| WD25 (Wilson and DeMets 2025) | `WD25_WK08-A` | `WD25_GMHRF-B` | `WD25_GWG20` | `WD25_GMHRF-A` |
| SG12 (Somoza and Ghidella 2012) | `SG12_WK08-A` | `SG12_GMHRF-B` | `SG12_GWG20` | `SG12_GMHRF-A` |

GMHRF-A' (folders `*_GMHRF-A`) and GMHRF-B share the `.rot`; only the prior of step 03b
(`PAPER_HS_PRIOR`) distinguishes them. The scripts of `models/Q22_WK08-A` (`02`, `03b`, `04` and
`data/`) are identical in the twelve folders.

## Inputs

**Inside the repository**

- `inputs/rotations/Global_Model_R_BG_2025.rot`: starting rotation model, with the source of each
  row in its comment and in its original timescale (declared input; sha256 `d861d537...`).
- Published tables used by `build_rotation_models.py`: Wilson and DeMets 2025, S2.2–S2.4 (extracted
  by `models/WD25_WK08-A/parse_wd25_si.py` into `models/WD25_WK08-A/wd25/`); Doubrovine et al. 2012,
  Table S3 (read by `rotation_models/frames/gmhrf/`); Gaastra et al. 2022, Table S6 (read by
  `rotation_models/frames/gwg20/` from `external/gaastra2022/`); Somoza and Ghidella (2012, Table 1)
  is written in `analysis/build_rot_sg12.py`.
- Geochronology: `models/*/data/ages_*.csv` (source of each age in its comment; Easter ages
  recomputed from ERDA with `analysis/ray2012_erda.py` and `ray2012_erda_ages.py`).
- `inputs/deposits/`: `Porphyry_Plate_ID_geodyn_v7.csv` and
  `Epithermal_Plate_Id_geodyn_v6_thesis.csv`, deposit database of this work, an update of the
  compilation of Bello-González et al. (2026) (sha256 `6999d373...` and `d21b4f4d...`). Rosario and
  Ujina are separate entries (in earlier compilations both were listed as Collahuasi). It is read by
  step 07 (`07_deposit_tracks.py`) and figures F1, F3, F4 and FigS1_decomp.
- `inputs/geochemistry/Full_Data_v11_tt_with_Clusters_v12_geodyn2.xlsx`: dated igneous samples
  (magmatic background), database of this work, an update of the compilation of Bello-González et al. (2026)
  (sha256 `b4e7fc54...`); same readers.
- `inputs/margin/`: `Margin_SAM_0.1.csv` (trench nodes every 0.1 Myr: discretisation of the
  restoration model of Arriagada et al. (2008) made by the author for Bello-González et al. (2018)) and
  `Margin_SAM.shp/.shx/.dbf/.prj/.cpg` (present-day trench, restored by
  `analysis/andean_shortening_utrecht.py`) (sha256 of the csv `2fb8eded...`). They are read by steps
  04, 05, 07 (via 04) and 09, `tables/margin_map.py`, figures F1, F3, FigS1c, FigS1d and FigS4, and
  the Andean shortening test.
- `inputs/validation/` (only the validations of Section S5; identical byte by byte to the
  originals): `motion_paths/` (27 `.xy` files, motion paths exported from GPlates
  with the same starting rotation model, 10, 36 and 56 Ma, f_PAC 0.36, 0.44 and 0.5), read by the `__main__` of
  `02_rot_engine.py`; and `HST_Margin/` (15 csv files, ages 0, 5, ..., 70 Ma). The rotation file is
  `inputs/rotations/Global_Model_R_BG_2025.rot`, which `02_rot_engine.py` uses by default without
  `PAPER_ROT_FILE`.
- `tools/hellinger_chang/`: own Python port of the Hellinger–Chang toolkit (package `chang`), with
  `validation/validate_vs_fortran.py` and `reference/` (dataA/dataB, outputs of `hellinger1.exe`
  fitA/fitB, and `gen_reference.py`). The `.exe` files are not included: regenerating the references
  requires compiling Chang's Fortran sources (T. Chang, personal communication 2019), which are not
  distributed. It is imported by `analysis/hellinger_circuit.py`.
- Utrecht model (Andean shortening test): `external/vanhinsbergen2026/`
  (`Utrecht_Paleogeography_15Dec25.rot`, `Continental_polygons.gpml`, `gAPWP25.rot`), from the
  supplementary material of van Hinsbergen et al. (2026); the Andean blocks 25xx are those of Schepers et al. (2017),
  updated by Boschman (2021).
- Wagenaar frame: `external/wagenaar2025/` (`Continent_frame701.rot` inside `Rotation_files.zip`),
  from Wagenaar et al. (2025), Zenodo doi:10.5281/zenodo.17536784 (version 3.0).

**Third-party inputs not included** (see `THIRD_PARTY.md`; place them under `external/`)

| Input | Source | Location | Read by |
|---|---|---|---|
| ETOPO1 | Amante and Eakins (2009) | `external/etopo1/` | `analysis/track_vs_ridge_axis.py`; F1 |
| Müller plate model | Müller et al. (2019) | `external/muller2019/` | F1; `analysis/hellinger_circuit.py` |
| GSFML, magnetic anomaly picks (isochron C32y) | Seton et al. (2014) | `external/gsfml/GSFML.global.picks.gmt` | `analysis/hellinger_circuit.py` (own covariances of 801 and 802 at 71.1 Ma, supplementary material) |
| GSFML, fracture zones (trace of the Challenger FZ) | Matthews et al. (2011), Wessel et al. (2015) | `external/gsfml/SF_shp/`; the trace used is already copied in `models/*/data/fz_challenger_gsfml.csv` | 03b (t_switch of Foundation, supplementary material); `analysis/foundation_t_switch_fz.py` (control plot) |

The trace `FZ_RM_0002-0189` has no name in GSFML; it is identified as the Challenger FZ because it is
the continuous trace that reaches the margin at ~33°S (its conjugate in the Pacific is the
Resolution FZ).

## Input sources

The citations are copied from the bibliography of the article.

**Circuit, frames and timescale** (`build_rotation_models.py`, `circuit_covariance.py`): Q22
Quiero et al. (2022); WD25 Wilson and DeMets (2025) (`models/WD25_WK08-A/wd25/`, Tables S2.2–S2.4); SG12 Somoza and Ghidella (2012)
(Table 1, written in `analysis/build_rot_sg12.py`); WK08-A Wessel and Kroenke (2008); GMHRF Doubrovine et al. (2012)
(Table S3, read by `rotation_models/frames/gmhrf/`); GWG20 Gaastra et al. (2022) (Table S6, in
`external/gaastra2022/` and `rotation_models/frames/gwg20/`); West–East Antarctica Granot et al. (2013),
Cande et al. (2000) and Granot and Dyment (2018); covariances added to the circuit Croon et al. (2008), Cande et al. (2010),
DeMets et al. (2021); GTS2020 timescale Ogg (2020); C32y picks Seton et al. (2014). The rotation models update the
dataset published with Bello-González et al. (2026).

**Geochronology** (`models/*/data/ages_*.csv`): Easter Ray et al. (2012) (ERDA heating steps); Juan
Fernández Lara et al. (2018); San Félix Haase et al. (2000); San Ambrosio González-Ferrán (1987); Foundation
O'Connor et al. (1998), O'Connor et al. (2001), O'Connor et al. (2004), with the plate change at the Challenger FZ (Bello-González et al. (2018);
trace from Matthews et al. (2011) and Wessel et al. (2015)); recalibration to Fish Canyon Kuiper et al. (2008) and decay
constant Min et al. (2000).

**Deposits, magmatic background and margin**: databases of this work, which update Bello-González et al. (2026);
trench restoration Arriagada et al. (2008), discretised for Bello-González et al. (2018).

**Sensitivity tests**: Utrecht model van Hinsbergen et al. (2026), with the blocks of Schepers et al. (2017) and
Boschman (2021); frame of Wagenaar et al. (2025); relief Amante and Eakins (2009); Müller et al. (2019).

References:

- Quiero et al. (2022): Quiero, F. et al. (2022). Growth of Neogene Andes linked to changes in plate convergence using high-resolution kinematic models. Nature Communications, 13, 1339. doi:10.1038/s41467-022-29055-4
- Wilson and DeMets (2025): Wilson, D. S. and DeMets, C. (2025). Changes in motion of the Nazca/Farallon plate over the last 34 million years: Implications for flat-slab subduction and the propagation of plate kinematic changes. Journal of Geophysical Research: Solid Earth, 130, e2025JB031933. doi:10.1029/2025JB031933
- Somoza and Ghidella (2012): Somoza, R. and Ghidella, M. E. (2012). Late Cretaceous to recent plate motions in western South America revisited. Earth and Planetary Science Letters, 331–332, 152–163. doi:10.1016/j.epsl.2012.03.003
- Wessel and Kroenke (2008): Wessel, P. and Kroenke, L. W. (2008). Pacific absolute plate motion since 145 Ma: An assessment of the fixed hot spot hypothesis. Journal of Geophysical Research, 113, B06101. doi:10.1029/2007JB005499
- Doubrovine et al. (2012): Doubrovine, P. V., Steinberger, B. and Torsvik, T. H. (2012). Absolute plate motions in a reference frame defined by moving hot spots in the Pacific, Atlantic, and Indian oceans. Journal of Geophysical Research, 117, B09101. doi:10.1029/2011JB009072
- Gaastra et al. (2022): Gaastra, K. M., Gordon, R. G. and Woodworth, D. T. (2022). Quantification of Pacific Plate hotspot tracks since 80 Ma. Tectonics, 41, e2021TC006772. doi:10.1029/2021TC006772
- Granot et al. (2013): Granot, R. et al. (2013). Revised Eocene–Oligocene kinematics for the West Antarctic rift system. Geophysical Research Letters, 40, 279–284. doi:10.1029/2012GL054181
- Cande et al. (2000): Cande, S. C. et al. (2000). Cenozoic motion between East and West Antarctica. Nature, 404, 145–150. doi:10.1038/35004501
- Granot and Dyment (2018): Granot, R. and Dyment, J. (2018). Late Cenozoic unification of East and West Antarctica. Nature Communications, 9, 3189. doi:10.1038/s41467-018-05270-w
- Croon et al. (2008): Croon, M. B., Cande, S. C. and Stock, J. M. (2008). Revised Pacific-Antarctic plate motions and geophysics of the Menard Fracture Zone. Geochemistry, Geophysics, Geosystems, 9, Q07001. doi:10.1029/2008GC002019
- Cande et al. (2010): Cande, S. C., Patriat, P. and Dyment, J. (2010). Motion between the Indian, Antarctic and African plates in the early Cenozoic. Geophysical Journal International, 183, 127–149. doi:10.1111/j.1365-246X.2010.04737.x
- DeMets et al. (2021): DeMets, C., Merkouriev, S. and Sauter, D. (2021). High resolution reconstructions of the Southwest Indian Ridge, 52 Ma to present: implications for the breakup and absolute motion of the Africa plate. Geophysical Journal International, 226, 1461–1497. doi:10.1093/gji/ggab107
- Ogg (2020): Ogg, J. G. (2020). Geomagnetic Polarity Time Scale. Geologic Time Scale 2020, 159–192. doi:10.1016/B978-0-12-824360-2.00005-X
- Bello-González et al. (2018): Bello-González, J. P., Contreras-Reyes, E. and Arriagada, C. (2018). Predicted path for hotspot tracks off South America since Paleocene times: Tectonic implications of ridge–trench collision along the Andean margin. Gondwana Research, 64, 216–234. doi:10.1016/j.gr.2018.07.008
- Bello-González et al. (2026): Bello-González, J. P. et al. (2026). Cenozoic geodynamic and tectonomagmatic evolution of the Andean margin: Implications for metallogenic endowment of orogenic belts. Physics and Chemistry of the Earth, 142, 104202. doi:10.1016/j.pce.2025.104202
- Ray et al. (2012): Ray, J. S. et al. (2012). Chronology and geochemistry of lavas from the Nazca Ridge and Easter Seamount Chain: A ~30 Myr hotspot record. Journal of Petrology, 53, 1417–1448. doi:10.1093/petrology/egs021
- Lara et al. (2018): Lara, L. E. et al. (2018). 40Ar/39Ar geochronological constraints on the age progression along the Juan Fernández Ridge, SE Pacific. Frontiers in Earth Science, 6, 194. doi:10.3389/feart.2018.00194
- Haase et al. (2000): Haase, K. M. et al. (2000). Sr-Nd-Pb isotope ratios, geochemical compositions, and 40Ar/39Ar data of lavas from San Félix Island (southeast Pacific): Implications for magma genesis and sources. Terra Nova, 12, 90–96. doi:10.1111/j.1365-3121.2000.00278.x
- González-Ferrán (1987): González-Ferrán, O. (1987). Evolución geológica de las islas chilenas en el Océano Pacífico. Islas oceánicas chilenas: conocimiento científico y necesidades de investigaciones, 37–54.
- O'Connor et al. (1998): O'Connor, J. M., Stoffers, P. and Wijbrans, J. R. (1998). Migration rate of volcanism along the Foundation Chain, SE Pacific. Earth and Planetary Science Letters, 164, 41–59. doi:10.1016/S0012-821X(98)00165-4
- O'Connor et al. (2001): O'Connor, J. M., Stoffers, P. and Wijbrans, J. R. (2001). En echelon volcanic elongate ridges connecting intraplate Foundation Chain volcanism to the Pacific–Antarctic spreading center. Earth and Planetary Science Letters, 192, 633–648. doi:10.1016/S0012-821X(01)00461-7
- O'Connor et al. (2004): O'Connor, J. M., Stoffers, P. and Wijbrans, J. R. (2004). The Foundation Chain: Inferring hotspot-plate interaction from a weak seamount trail. Oceanic Hotspots: Intraplate Submarine Magmatism and Tectonism, 349–374. doi:10.1007/978-3-642-18782-7_12
- Kuiper et al. (2008): Kuiper, K. F. et al. (2008). Synchronizing rock clocks of Earth history. Science, 320, 500–504. doi:10.1126/science.1154339
- Min et al. (2000): Min, K. et al. (2000). A test for systematic errors in 40Ar/39Ar geochronology through comparison with U/Pb analysis of a 1.1-Ga rhyolite. Geochimica et Cosmochimica Acta, 64, 73–98. doi:10.1016/S0016-7037(99)00204-5
- Seton et al. (2014): Seton, M. et al. (2014). Community infrastructure and repository for marine magnetic identifications. Geochemistry, Geophysics, Geosystems, 15, 1629–1641. doi:10.1002/2013GC005176
- Matthews et al. (2011): Matthews, K. J. et al. (2011). The tectonic fabric of the ocean basins. Journal of Geophysical Research, 116, B12109. doi:10.1029/2011JB008413
- Wessel et al. (2015): Wessel, P. et al. (2015). Semiautomatic fracture zone tracking. Geochemistry, Geophysics, Geosystems, 16, 2462–2472. doi:10.1002/2015GC005853
- Arriagada et al. (2008): Arriagada, C. et al. (2008). Paleogene building of the Bolivian Orocline: Tectonic restoration of the central Andes in 2-D map view. Tectonics, 27, TC6014. doi:10.1029/2008TC002269
- van Hinsbergen et al. (2026): van Hinsbergen, D. J. J. et al. (2026). Paleolatitude.org 3.0: a calculator for paleoclimate and paleobiology studies based on a new global paleogeography model. PLoS One, 21, e0346817. doi:10.1371/journal.pone.0346817
- Schepers et al. (2017): Schepers, G. et al. (2017). South-American plate advance and forced Andean trench retreat as drivers for transient flat subduction episodes. Nature Communications, 8, 15249. doi:10.1038/ncomms15249
- Boschman (2021): Boschman, L. M. (2021). Andean mountain building since the Late Cretaceous: a paleoelevation reconstruction. Earth-Science Reviews, 220, 103640. doi:10.1016/j.earscirev.2021.103640
- Wagenaar et al. (2025): Wagenaar, S. D. M., Vaes, B. and van Hinsbergen, D. J. J. (2025). Toward reconstructing mantle convection using a minimum-continent-motion reference frame. Journal of Geophysical Research: Solid Earth, 130, e2024JB030430. doi:10.1029/2024JB030430
- Amante and Eakins (2009): Amante, C. and Eakins, B. W. (2009). ETOPO1 1 Arc-Minute Global Relief Model: Procedures, Data Sources and Analysis. NOAA National Geophysical Data Center. doi:10.7289/V5C8276M
- Müller et al. (2019): Müller, R. D. et al. (2019). A global plate model including lithospheric deformation along major rifts and orogens since the Triassic. Tectonics, 38, 1884–1907. doi:10.1029/2018TC005462

`02_rot_engine.py` uses `inputs/rotations/Global_Model_R_BG_2025.rot` by default if
`PAPER_ROT_FILE` is not set; all the orchestrators set it.

## Processing chain, in order

All the orchestrators are in `analysis/` except the first two.

Steps 3 to 10 run together with `bash run_chain.sh --run` (in `analysis/`). Without arguments it
is a test mode that runs nothing: it checks scripts, `.rot` files and fits, lists the commands and
flags the outputs older than the 03b fit of their folder. With `--run` it first moves the old
step 04 outputs to `<folder>/backup_chain_<date>/` (otherwise `run_models.sh` skips them as
complete) and stops at the first error. It adds steps 05 and 09 (convergence velocity), which no
other orchestrator runs.

1. **Rotations.** `py rotation_models/build_rotation_models.py`: builds the twelve `.rot` files from
   the starting one (Granot et al. 2013 for 804-802, node of Cande et al. 2000 at 43.8 Ma, notes on
   the Farallon–Nazca change, ages to GTS2020 with `rot_gts2020.py`, which is not idempotent) and
   compares them with the `.rot` files of the runs.
2. **Hotspot fit (03b).** `py rotation_models/fit_hotspots.py`: the twelve fits with the environment
   of the final run (without the present-day position as an age-0 observation, the BIC chooses fixed
   or drifting, Foundation ages `_v3`, w as confidence, gmhrf prior in GMHRF-B). It writes with tag
   `_rep` and compares with `hotspot_fit.json` of each folder (difference 0 in the twelve).
3. **Monte Carlo and following steps.** `bash run_models.sh`: step 04 (N = 1000; variants total,
   circuit, frame and hotspot; three f_PAC values), merge, 07 (deposits), 08 (ellipses), 10, 11 and
   12 (deposit–track test). About 25 min per step 04 job.
4. **Deposits with lag and Utrecht.** Functions `lag` (step 07 with `PAPER_DEP_LAG` = 1, 2, 5 Myr)
   and `utrecht` of `run_chain.sh`, which calls `run_shortening_step07.sh` (Andean shortening test).
5. **Bands at the present-day margin.** Function `presentes` of `run_chain.sh`
   (`tables/make_present_bands.py`, outputs in `results/`).
6. **Supplementary sensitivities (S4).** `run_sens9.sh` (lag, without YORK2, interpolation),
   `run_bic_noforce.sh` (Easter BIC), function `wagenaar` of `run_chain.sh` (frame of Wagenaar et
   al. 2025), `run_sens_selkirk.sh` (Selkirk microplate); summaries with
   `sens9_summary.py`, `sens9_deposits.py`, `sens9_wagenaar.py`, `sens_selkirk.py`.
7. **Summaries read by `manuscript_numbers.py`.** `bash run_summaries.sh` (deposit offset, Nazca
   Ridge, rotated mean, numbers of the deposit distribution) and finally `manuscript_numbers.py`.
8. **Figures.** `bash run_figures.sh` (Figs. 1–6 and the supplementary figures that depend on the
   fit; copies to `figures/output/`). These do not depend on the fit: FigS_circuit_nodes,
   FigS1b_vcsigma*, FigS1c, FigS4_workflow.
9. **Supplementary tables.** `supplementary/make_supplementary_tables.py` (`tableS1.tex`) and
   `make_supplementary_excel.py` (`TablesS_data.xlsx`).
10. **Verification.** `py manuscript_numbers.py` (all sections) or `py manuscript_numbers.py <section>`.

Step 04 uses a fixed seed (20260714): with the same `.rot` and the same fit it gives identical bands
(checked in `Q22_WK08-A`, total f_PAC 0.44, 258 rows with difference 0).

## Validations of Section S5

They run with the inputs of `inputs/validation/` and `tools/hellinger_chang/`:

- Hellinger–Chang vs Fortran: `py validation/validate_vs_fortran.py` in `tools/hellinger_chang/`.
  It compares with the saved outputs of the `.exe` files (not included): identical rotations and
  degrees of freedom, H11.2 with relative difference 1.6e-8 (A) and 1.2e-8 (B); Welch df 91.9840 vs
  91.983971.
- Rotation engine vs GPlates motion paths: `py 02_rot_engine.py` in `models/Q22_WK08-A/` without
  `PAPER_ROT_FILE`; maximum 0.014 km (article: 0.02 km).
