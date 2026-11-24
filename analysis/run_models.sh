#!/bin/bash
# Twelve-model grid (3 circuits x 4 frames) with Foundation _v3
# (foundation_oconnor2004.py). Uses Granot et al. (2013) for the 804-802 nodes
# (12o, 13o, 16y and 18o, and the Cande 2000 node at 43.462 in place of
# "18o constant at 45.03"); .rot in GTS2020; circuit covariances (PAPER_CIRCUIT_COV=b).
# Tags _hsb_cc_g13 and _hsb_pred_cc_g13. Includes SG12_WK08-A (circuit of
# Somoza and Ghidella 2012, build_rot_sg12.py; uses the 03b fit of the base model) and the
# 03b fit with the ERDA ages (pred: hotspot_fit_gauss_pred_age_E.json).
# Each folder uses its 03b fit (PAPER_HS_FIT) and its _hsb .rot.
#   step 04: N = 1000, variants total, circuit, wk08 (frame) and hotspot, three
#            f_PAC; with w as confidence (tag _hsb) and, for total and hotspot,
#            with w as prediction (tag _hsb_pred). NP processes at a time (8 by default).
#   steps 08 (N 500) and 10 (N 400) with both tags; 11 with _hsb.
#   step 07 (07_deposit_tracks.py, assignment of each deposit to the nearest of the four
#   chains) in all folders with the confidence fit, then step 12 (12_deposit_track_test.py,
#   test across the models).
# The outputs carry the tag.
# Usage:  bash run_models.sh > run_models.log 2>&1
P="${PAPER_ROOT:?define PAPER_ROOT}"
A="$P/analysis"
N04=${N04:-1000}; N08=${N08:-500}; N10=${N10:-400}   # N04=2 N08=3 N10=3 for a test run
NP=${NP:-8}   # simultaneous processes (NP=6 leaves cores free for other work)
export PAPER_CIRCUIT_COV=b PAPER_TAG_EXTRA=_g13 PAPER_INTERP_COV=norm
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# folder | conf rot | pred rot | prior | frame covariance
MODELOS="Q22_WK08-A|Global_Model_R_BG_2025_hsb.rot|Global_Model_R_BG_2025_hsb_gauss_pred_age_A.rot|none|wk08
WD25_WK08-A|Global_Model_R_BG_2025_WD25_hsb.rot|Global_Model_R_BG_2025_WD25_hsb_gauss_pred_age_A.rot|none|wk08
Q22_GMHRF-A|Global_Model_R_BG_2025_GMHRF_gts12_hsb.rot|Global_Model_R_BG_2025_GMHRF_gts12_hsb_gauss_pred_age_A.rot|none|rot
Q22_GMHRF-B|Global_Model_R_BG_2025_GMHRF_gts12_hsbgmhrf.rot|Global_Model_R_BG_2025_GMHRF_gts12_hsbgmhrf_gauss_pred_age_A.rot|gmhrf|rot
WD25_GMHRF-B|Global_Model_R_BG_2025_WD25_GMHRF_gts12_hsbgmhrf.rot|Global_Model_R_BG_2025_WD25_GMHRF_gts12_hsbgmhrf_gauss_pred_age_A.rot|gmhrf|rot
Q22_GWG20|Global_Model_R_BG_2025_GWG20_hsb.rot|Global_Model_R_BG_2025_GWG20_hsb_gauss_pred_age_A.rot|none|rot
WD25_GWG20|Global_Model_R_BG_2025_WD25_GWG20_hsb.rot|Global_Model_R_BG_2025_WD25_GWG20_hsb_gauss_pred_age_A.rot|none|rot
SG12_WK08-A|Global_Model_R_BG_2025_SG12_hsb.rot|Global_Model_R_BG_2025_SG12_hsb_gauss_pred_age_A.rot|none|wk08
WD25_GMHRF-A|Global_Model_R_BG_2025_WD25_GMHRF_gts12_hsb.rot|Global_Model_R_BG_2025_WD25_GMHRF_gts12_hsb_gauss_pred_age_A.rot|none|rot
SG12_GMHRF-A|Global_Model_R_BG_2025_SG12_GMHRF_gts12_hsb.rot|Global_Model_R_BG_2025_SG12_GMHRF_gts12_hsb_gauss_pred_age_A.rot|none|rot
SG12_GMHRF-B|Global_Model_R_BG_2025_SG12_GMHRF_gts12_hsbgmhrf.rot|Global_Model_R_BG_2025_SG12_GMHRF_gts12_hsbgmhrf_gauss_pred_age_A.rot|gmhrf|rot
SG12_GWG20|Global_Model_R_BG_2025_SG12_GWG20_hsb.rot|Global_Model_R_BG_2025_SG12_GWG20_hsb_gauss_pred_age_A.rot|none|rot"

# runs one step with the environment of (folder, w): w = conf | pred
step() {  # folder w log command...
  local d="$1" w="$2" log="$3"; shift 3
  local line; line=$(echo "$MODELOS" | grep "^$d|")
  IFS="|" read -r _ rc rp prior cov <<< "$line"
  local rot="$rc" fit="hotspot_fit.json"
  if [ "$w" = pred ]; then rot="$rp"; fit="hotspot_fit_gauss_pred_age_E.json"; fi
  ( cd "$P/models/$d" && PAPER_ROT_FILE="$P/models/$d/rot/$rot" PAPER_HS_FIT="$P/models/$d/$fit" \
      PAPER_HS_PRIOR="$prior" PAPER_FRAME_COV="$cov" "$@" > "$log" 2>&1 )
  echo "$(date +%T) $d [$w] $* -> exit $?"
}
export -f step; export P MODELOS N04

echo "== start $(date '+%F %T')"
# list of step 04 jobs (longest first: total)
JOBS=""
for v in total hotspot circuit wk08; do
  for line in $MODELOS; do d="${line%%|*}"
    for pf in 0.5 0.44 0.36; do p="${pf//./}"
      JOBS+="$d conf _hsb_cc_g13_${v:0:3}$p $v $pf"$'\n'
      if [ $v = total ] || [ $v = hotspot ]; then
        JOBS+="$d pred _hsb_pred_cc_g13_${v:0:3}$p $v $pf"$'\n'
      fi
    done
  done
done
# a part that is already complete is skipped (its log records N04 realizations and the CSV
# exists), so an interrupted run resumes by launching the same script again
printf "%s" "$JOBS" | xargs -P "$NP" -L 1 bash -c \
  'if grep -q ": $N04 realizaciones" "$P/$0/mc$2.log" 2>/dev/null &&
      [ -f "$P/$0/mc_collision_bands_$3$2.csv" ]; then echo "$0 $2 already complete"; exit 0; fi
   PAPER_SUFFIX="$2" step "$0" "$1" "mc$2.log" py 04_mc_collisions.py "$N04" "$3" "$4"'
echo "== 04 end $(date '+%F %T')"

for line in $MODELOS; do d="${line%%|*}"
  py "$A/merge_mc_parts_tag.py" "$P/models/$d" _hsb_cc_g13 > "$P/models/$d/merge_g13_hsb.log" 2>&1; echo "merge $d _hsb_cc_g13 exit $?"
  py "$A/merge_mc_parts_tag.py" "$P/models/$d" _hsb_pred_cc_g13 total,hotspot > "$P/models/$d/merge_g13_hsb_pred.log" 2>&1; echo "merge $d _hsb_pred_cc_g13 exit $?"
done

# steps 07, 08, 10, 11, NP at a time
J2=""
for line in $MODELOS; do d="${line%%|*}"
  for w in conf pred; do
    J2+="$d $w ds08_g13_$w.log py 08_track_ellipses.py $N08"$'\n'
    J2+="$d $w ds10_g13_$w.log py 10_recon_ellipses.py $N10"$'\n'
  done
  J2+="$d conf ds11_g13_conf.log py 11_recon_model_tracks.py"$'\n'
  J2+="$d conf ds07_g13_conf.log py 07_deposit_tracks.py"$'\n'
done
printf "%s" "$J2" | xargs -P "$NP" -L 1 bash -c 'step "$@"' _
( cd "$P/models/Q22_WK08-A" && py 12_deposit_track_test.py > ds12_g13.log 2>&1; echo "12 exit $?" )
echo "== end $(date '+%F %T')"
