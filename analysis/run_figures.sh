#!/bin/bash
# Article figures that depend on the 03b fit and on the g13 bands, from the
# outputs of run_models.sh and the bands at the present margin (results/), copied to
# figures/output/main and figures/output/supplementary.
# Mapping: Fig1 F1_track_ellipses, Fig2 F2_collision_bands, Fig3 F5_model_comparison,
# Fig4 F3_reconstruction, Fig5 flat_slab_tracks_test (analysis/flat_slab_tracks.py), Fig6
# F4b_deposit_summary; supplement: FigS1_decomp, F2/F5 _supl_gmhrfA, FigS_deposit_cdf_chain.
# Independent of the fit (not remade): FigS_circuit_nodes, FigS1b_vcsigma*, FigS1c, FigS4_workflow.
# Usage:  bash run_figures.sh > run_figures.log 2>&1
P="${PAPER_ROOT:?define PAPER_ROOT}"; A="$P/analysis"; F="$P/figures"
G="$P/results"; E="$P/figures/output"
mkdir -p "$E/main" "$E/supplementary"
export PAPER_FIG_TAG=_hsb_cc_g13 PAPER_TABLES_DIR="$G" PAPER_FIG_OUT="$G" PYTHONIOENCODING=utf-8
cd "$F"
for s in make_F1_track_ellipses make_F2_collision_bands make_F5_model_comparison make_F3_reconstruction \
         make_F4b_deposit_summary make_FigS1_decomp_vcsigma make_FigS_deposit_cdf_chain \
         make_FigS_wagenaar make_FigS_shortening make_FigS_lag; do
  py $s.py > "$G/log_$s.txt" 2>&1; echo "$s exit $?"
done
for s in make_F2_collision_bands make_F5_model_comparison; do
  PAPER_MODELOS=supl py $s.py > "$G/log_${s}_supl.txt" 2>&1; echo "$s supl exit $?"
done
( cd "$A" && py flat_slab_tracks.py > "$G/log_flat_slab_tracks.txt" 2>&1 ); echo "flat_slab exit $?"
cp "$G/F1_track_ellipses.pdf" "$E/main/Fig1.pdf"; cp "$G/F2_collision_bands.pdf" "$E/main/Fig2.pdf"
cp "$G/F5_model_comparison.pdf" "$E/main/Fig3.pdf"; cp "$G/F3_reconstruction.pdf" "$E/main/Fig4.pdf"
cp "$G/flat_slab_tracks_test.pdf" "$E/main/Fig5.pdf"; cp "$F/F4b_deposit_summary.pdf" "$E/main/Fig6.pdf"
for f in FigS1_decomp F2_collision_bands_supl_gmhrfA F5_model_comparison_supl_gmhrfA FigS_deposit_cdf_chain \
         FigS_wagenaar FigS_shortening FigS_lag; do
  cp "$G/$f.pdf" "$E/supplementary/"
done
echo "copied $(date +%T)"
