#!/bin/bash
# Sensitivity tests of the supplement with the NINE models of the text (g13 final run).
#   (a) hotspot-eruption lag: 03b without the 0 Ma observation (PAPER_HS_P0=0) with
#       PAPER_HS_LAG = 0, 1, 2 Myr (tags _s9_lag0/1/2).
#   (b) without the two YORK2 isochrons of Ray et al. 2012 (PAPER_HS_EXCLUDE, tag _s9_noYORK2).
#   (c) interpolation: step 04 (total, f_PAC 0.44, N 1000, same seed) WITHOUT
#       PAPER_INTERP_COV=norm; compared with the g13 bands (which use norm).
# Nothing overwrites the final outputs: the .rot files go to rot/sens9/ and the json/band files carry the tag _s9.
# Usage:  bash run_sens9.sh > run_sens9.log 2>&1
P="${PAPER_ROOT:?define PAPER_ROOT}"; A="$P/analysis"
NP=${NP:-7}; N04=${N04:-1000}
export PAPER_CIRCUIT_COV=b PAPER_TAG_EXTRA=_g13 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# no PAPER_HS_FORCE (the BIC prefers the Easter drift without the 0 Ma observation)
MODELOS="Q22_WK08-A|Global_Model_R_BG_2025_hsb.rot|none|wk08
Q22_GMHRF-B|Global_Model_R_BG_2025_GMHRF_gts12_hsbgmhrf.rot|gmhrf|rot
Q22_GWG20|Global_Model_R_BG_2025_GWG20_hsb.rot|none|rot
WD25_WK08-A|Global_Model_R_BG_2025_WD25_hsb.rot|none|wk08
WD25_GMHRF-B|Global_Model_R_BG_2025_WD25_GMHRF_gts12_hsbgmhrf.rot|gmhrf|rot
WD25_GWG20|Global_Model_R_BG_2025_WD25_GWG20_hsb.rot|none|rot
SG12_WK08-A|Global_Model_R_BG_2025_SG12_hsb.rot|none|wk08
SG12_GMHRF-B|Global_Model_R_BG_2025_SG12_GMHRF_gts12_hsbgmhrf.rot|gmhrf|rot
SG12_GWG20|Global_Model_R_BG_2025_SG12_GWG20_hsb.rot|none|rot"

uno() {  # folder type tag
  local d="$1" tipo="$2" tag="$3" line rc prior cov
  line=$(echo "$MODELOS" | grep "^$d|"); IFS="|" read -r _ rc prior cov <<< "$line"
  mkdir -p "$P/models/$d/rot/sens9"
  cd "$P/models/$d" || return
  case "$tipo" in
    lag) env PAPER_ROT_FILE="$P/models/$d/rot/$rc" PAPER_HS_PRIOR="$prior" PAPER_HS_P0=0 \
           PAPER_HS_LAG="${tag#_s9_lag}" PAPER_HS_OUT_TAG="$tag" PAPER_HS_ROT_OUT="$P/models/$d/rot/sens9/$tag.rot" \
           py 03b_hotspot_fit.py > "hotspot_fit$tag.log" 2>&1 ;;
    york) env PAPER_ROT_FILE="$P/models/$d/rot/$rc" PAPER_HS_PRIOR="$prior" PAPER_HS_EXCLUDE="DRFT-115-2,DRFT-126-1" \
           PAPER_HS_OUT_TAG="$tag" PAPER_HS_ROT_OUT="$P/models/$d/rot/sens9/$tag.rot" \
           py 03b_hotspot_fit.py > "hotspot_fit$tag.log" 2>&1 ;;
    interp) env PAPER_ROT_FILE="$P/models/$d/rot/$rc" PAPER_HS_FIT="$P/models/$d/hotspot_fit.json" PAPER_HS_PRIOR="$prior" \
           PAPER_FRAME_COV="$cov" PAPER_SUFFIX="$tag" py 04_mc_collisions.py "$N04" total 0.44 > "mc$tag.log" 2>&1 ;;
  esac
  echo "$(date +%T) $d $tipo $tag exit $?"
}
export -f uno; export P MODELOS N04
echo "== start $(date '+%F %T')"
J=""
for line in $MODELOS; do d="${line%%|*}"; J+="$d interp _s9_interp0_tot044"$'\n'; done
for line in $MODELOS; do d="${line%%|*}"
  for L in 0 1 2; do J+="$d lag _s9_lag$L"$'\n'; done
  J+="$d york _s9_noYORK2"$'\n'
done
printf "%s" "$J" | xargs -P "$NP" -L 1 bash -c 'uno "$@"' _
echo "== end $(date '+%F %T')"
