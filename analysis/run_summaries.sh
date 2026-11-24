#!/bin/bash
# Summaries read by manuscript_numbers.py that depend on the 03b fit, on the deposit tracks
# (step 07, with and without lag) and on the g13 bands. Run after run_models.sh, after
# step 07 with lag (run_chain.sh, step 3) and after the bands at the present margin.
#   deposit_lag_9models.py   -> outputs/deposit_lag_9models.csv
#   track_vs_ridge_axis.py   -> outputs/track_vs_ridge_axis_g13.csv
#   rotated_mean_9models.py  -> outputs/rotated_mean_9models.txt
#   make_F4_deposit_distance -> results/F4_deposit_distance_numbers.csv
#   manuscript_numbers.py    -> outputs/manuscript_numbers.csv (each number of the text against the computed one)
# Usage:  bash run_summaries.sh > run_summaries.log 2>&1
P="${PAPER_ROOT:?define PAPER_ROOT}"; A="$P/analysis"; G="$P/results"
export PYTHONIOENCODING=utf-8
cd "$A"
py deposit_lag_9models.py > /dev/null; echo "deposit_lag exit $?"
py -I track_vs_ridge_axis.py > /dev/null; echo "track_vs_ridge_axis exit $?"
py rotated_mean_9models.py > /dev/null; echo "rotated_mean exit $?"
( cd "$P/figures" && PAPER_FIG_TAG=_hsb_cc_g13 PAPER_TABLES_DIR="$G" PAPER_FIG_OUT="$G" \
    py make_F4_deposit_distance.py > "$G/log_make_F4_deposit_distance.txt" 2>&1 ); echo "F4 numbers exit $?"
py manuscript_numbers.py | tail -1
