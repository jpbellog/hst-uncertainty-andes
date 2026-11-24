# Figures

One script per figure. `analysis/run_figures.sh` runs them all and copies the PDFs to
`figures/output/main/` and `figures/output/supplementary/` under the article's numbering.

| Article | Script | Main inputs |
|---|---|---|
| Fig. 1 | `make_F1_track_ellipses.py` | `models/*/track_ellipses_0.44_hsb_cc_g13.csv`, `track_nominal_*`, `data/etopo1_F1.npz` |
| Fig. 2 | `make_F2_collision_bands.py` | `results/mc_collision_bands_*_present_*.csv` |
| Fig. 3 | `make_F5_model_comparison.py` | `results/mc_collision_bands_total_present_*.csv`; writes `results/F5_model_comparison_numbers.csv` |
| Fig. 4 | `make_F3_reconstruction.py` | `models/*/recon_ellipses_0.44_hsb_cc_g13.csv`, `recon_model_tracks_*`, `inputs/` |
| Fig. 5 | `analysis/flat_slab_tracks.py` | `analysis/data/wd25_table3_flat_slabs.csv`, collision bands |
| Fig. 6 | `make_F4b_deposit_summary.py` | `models/*/deposit_tracks_*.csv`, `models/Q22_WK08-A/deposit_track_test_by_model.csv` |
| Fig. S1 | `make_FigS_circuit_nodes.py` | `analysis/outputs/frame_coverage/hellinger/`; writes `analysis/outputs/circuit_nodes_g13.csv` |
| Figs. S2, S3 | `make_FigS1_decomp_vcsigma.py` (`FigS1_decomp`, `FigS1b_vcsigma*`) | collision bands by variant (`total`, `circuit`, `hotspot`, `wk08`), `models/*/vc_sigma_*.csv` |
| Fig. S4 | `make_FigS1c_velocity_field.py` | collision bands, the rotation file in `PAPER_ROT_FILE`, `inputs/margin/` |
| Figs. S5, S6 | `make_F2_collision_bands.py`, `make_F5_model_comparison.py` with the GMHRF-A' models | `results/` |
| Fig. S7 | `make_FigS_wagenaar.py` | `analysis/outputs/sens9/wagenaar_9models.csv` |
| Fig. S8 | `make_FigS_shortening.py` | `analysis/outputs/andean_shortening/` |
| Fig. S9 | `make_FigS_lag.py` | `analysis/outputs/deposit_lag_9models.csv` |
| Fig. S10 | `make_FigS_deposit_cdf_chain.py` | `models/*/deposit_tracks_*.csv` |
| Fig. S11 | `make_FigS4_workflow.py` | `models/Q22_WK08-A/track_*_0.44_hsb_cc_g13.csv` |
| Fig. S12 | `analysis/ridge_axes_figure.py` | outputs of `analysis/bathymetry_axes.py` (needs ETOPO1), `models/*/track_nominal_0.44_hsb_cc_g13.csv` |

`make_F4_deposit_distance.py` has no figure in the article; it writes the distances quoted in the text
to `results/F4_deposit_distance_numbers.csv`. `paper_models.py` defines the model sets (nine in the main
text, three GMHRF-A' in the Supplementary Material) shared by the scripts.

## Environment

| Variable | Default | Role |
|---|---|---|
| `PAPER_ROOT` | required | repository root |
| `PAPER_FIG_TAG` | empty (`run_figures.sh` sets `_hsb_cc_g13`) | tag of the outputs read |
| `PAPER_TABLES_DIR` | `tables/` (`run_figures.sh` sets `results/`) | folder with the present-day collision bands |
| `PAPER_FIG_OUT` | this folder | where the PDFs are written (the supplementary figures S7 to S9 go to `results/`) |
| `PAPER_PIPE`, `PAPER_PIPE_WD` | `Q22_WK08-A`, `WD25_WK08-A` | model folders drawn by the single-model figures |

`data/etopo1_F1.npz` is a crop of ETOPO1 (Amante and Eakins, 2009) used as the background of Fig. 1, so
that the figure can be drawn without the full grid.
