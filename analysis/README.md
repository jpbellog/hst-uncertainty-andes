# Analysis

Orchestrators, tools that build parts of the rotation models, sensitivity tests and the summaries
read by `manuscript_numbers.py`. The order in which everything runs is set out in
`rotation_models/README.md` (section "Processing chain, in order").

## Orchestrators

| Script | Role |
|---|---|
| `run_chain.sh` | Whole chain after the hotspot fit. Without arguments it only checks inputs and lists the commands; `--run` runs it |
| `run_models.sh` | Steps 04 to 12 in the twelve models (holds the `MODELOS` table with the rotation file, prior and frame covariance of each model) |
| `run_sens9.sh` | Sensitivity tests of Supplementary Section S4 in the nine main-text models (hotspot-eruption lag, without the YORK2 isochrons, covariance interpolation) |
| `run_bic_noforce.sh` | Fits for the BIC of the Easter fixed/drift choice |
| `run_sens_selkirk.sh` | Foundation chain and the Selkirk microplate |
| `run_shortening_step07.sh` | Step 07 with the deposit positions and trench restored for Andean shortening |
| `run_summaries.sh` | Summaries read by `manuscript_numbers.py`, then `manuscript_numbers.py` |
| `run_figures.sh` | Figures of the article and the Supplementary Material (see `figures/README.md`) |

## Verification

`manuscript_numbers.py` recomputes the numbers of the main text and the Supplementary Material from the
saved outputs and compares each one with the value printed in the article (`py manuscript_numbers.py`
for all sections, `py manuscript_numbers.py <section>` for one). It writes
`outputs/manuscript_numbers.csv`. `table_bands_68.py` builds the 68 % rows of Table 2 and imports it.

## Rotation models and covariances

| Script | Role |
|---|---|
| `build_rot_sg12.py` | SG12 circuit: Nazca (Farallon)–South America motion of Somoza and Ghidella (2012, Table 1) |
| `build_rot_wagenaar.py` | Model in the minimum-continent-motion frame of Wagenaar et al. (2025) (sensitivity test) |
| `circuit_covariance.py` | Covariances of the circuit links that the starting rotation file carries without covariance |
| `hellinger_circuit.py` | Hellinger fits of circuit nodes without a usable published covariance (GSFML picks; uses `tools/hellinger_chang/`) |
| `interp_node_covariances.py` | Effect of drawing the perturbations at the nodes and interpolating between them (plate 002) |
| `frame_coverage.py` | Age limit of each model given its oldest frame covariance |
| `gwg20_convention.py`, `gwg20_rotate_covariance.py` | Unit and convention of the GWG20 covariances (Gaastra et al., 2022, Table S6) |
| `rot_gts2020.py`, `timescale_GTS2020.py`, `timescale_CK95_GTS2012.py` | Node ages moved to GTS2020 and the chron tables used for it |
| `rot_sens_granot2013.py`, `g13_sign.py` | West–East Antarctica nodes of Granot et al. (2013) and the check of their sign convention |

## Geochronology

| Script | Role |
|---|---|
| `ArAr_recalibration.py` | 40Ar/39Ar ages recalibrated to FCs 28.201 Ma with the decay constant of Min et al. (2000) |
| `ray2012_erda.py`, `ray2012_erda_ages.py`, `ray2012_isochrons.py` | Easter ages rebuilt from the step-heating data of Ray et al. (2012) in ERDA |
| `foundation_oconnor1998.py`, `foundation_oconnor2004.py`, `foundation_t_switch_fz.py` | Foundation chain ages and the age of its plate change at the Challenger fracture zone |
| `decay_constant_test.py` | Effect of the uncertainty of the 40K decay constant |
| `rotated_mean.py`, `rotated_mean_9models.py` | Mean of the dated samples rotated back to the hotspot frame, compared with the fit |

## Sensitivity tests and summaries

| Script | Role |
|---|---|
| `sens9_summary.py`, `sens9_deposits.py`, `sens9_wagenaar.py`, `sens_selkirk.py`, `sens_noP0.py` | Summaries of the sensitivity tests |
| `andean_shortening_utrecht.py`, `step07_utrecht.py` | Restoration of Andean shortening with the Utrecht block model and step 07 on it |
| `deposit_lag_9models.py` | Collision–deposit lag in the nine models |
| `assign_nearest_track.py` | Assignment of deposits to the nearest of the four tracks |
| `nominal_collision.py`, `nominal_sensitivity.py`, `merge_mc_parts_tag.py` | Unperturbed collision points, and the merge of the parallel parts of step 04 |
| `bathymetry_axes.py`, `track_vs_ridge_axis.py`, `ridge_axes_figure.py` | Tracks against the ridge axes in ETOPO1 (Fig. S12) |
| `flat_slab_tracks.py` | Arrival of the tracks beneath the flat-slab segments of Wilson and DeMets (2025, Table 3) (Fig. 5) |

## Data and outputs

`data/wd25_table3_flat_slabs.csv` transcribes Table 3 of Wilson and DeMets (2025). `outputs/` holds the
outputs read by `manuscript_numbers.py` and the figures: `frame_coverage/` (with the Hellinger fits in
`hellinger/`), `sens9/`, `selkirk/`, `andean_shortening/`, the chron tables, `deposit_lag_9models.csv`,
`rotated_mean_9models.txt`, `track_vs_ridge_axis_g13.csv` and `ray2012_erda_samples.csv`.

Some column names and some text lines in these outputs are in Spanish. They are parsed by
`manuscript_numbers.py` and are left as written.
