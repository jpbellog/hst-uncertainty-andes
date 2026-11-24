#!/bin/bash
# Single orchestrator of steps 3 to 10 of rotation_models/README.md (the ordered chain):
# everything after the 03b fit (rotation_models/fit_hotspots.py) up to manuscript_numbers.py.
#
#   bash run_chain.sh            test mode (default): runs and writes nothing. Checks that
#                                 the scripts and inputs exist, that the 03b fit of the twelve folders
#                                 is the base case (no 0 Ma observation and no forced drift), lists
#                                 the commands in order and flags as OLD each output older than the
#                                 hotspot_fit.json of its folder (or MISSING).
#   bash run_chain.sh --run   runs the chain. It first MOVES to <folder>/backup_chain_<date>/ the
#                                 step 04 outputs older than the fit, so that
#                                 run_models.sh does not skip them as complete. Stops at the
#                                 first step with an error. Log in outputs/chain_<date>/.
#
# Steps (in parentheses, the script that is called or replaced here):
#   1. convergence velocity: step 05 (pipeline) and step 09 (Q22, WD25, SG12); read by
#      manuscript_numbers.py (sec_velocidad)
#   2. step 04, merge, 07, 08, 10, 11 and 12 in the twelve models (run_models.sh)
#   3. step 07 with deposit lag 1, 2 and 5 Myr in the nine models
#   4. Andean shortening, Utrecht (run_shortening_step07.sh) and its summary (andean_shortening_utrecht.py)
#   5. bands at the present margin
#   6. supplementary sensitivities: run_sens9.sh, run_bic_noforce.sh, Wagenaar frame (with a backup
#      of the previous outputs), run_sens_selkirk.sh and their summaries
#   7. summaries read by manuscript_numbers.py (run_summaries.sh)
#   8. figures (run_figures.sh)
#   9. supplementary tables (make_supplementary_tables.py, make_supplementary_excel.py)
#  10. manuscript_numbers.py
# Does not run the 03b (step 2 of the README) or the seeded variant (sens_conP0.py), which is not in the text.
P="${PAPER_ROOT:?define PAPER_ROOT}"
A="$P/analysis"; G="$P/results"
S="$P/supplementary"; M="$A/outputs/wagenaar_frame"
MODO=prueba; [ "$1" = "--run" ] && MODO=correr
FECHA=$(date +%Y%m%d_%H%M); L="$A/outputs/chain_$FECHA"
NP=${NP:-8}
export P OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONIOENCODING=utf-8
MODELOS=$(sed -n '/^MODELOS="/,/"$/p' "$A/run_models.sh" | sed 's/^MODELOS="//; s/"$//')
NUEVE=$(sed -n '/^MODELOS="/,/"$/p' "$A/run_sens9.sh" | sed 's/^MODELOS="//; s/"$//' | cut -d"|" -f1)
CIRCUITOS="Q22_WK08-A WD25_WK08-A SG12_WK08-A"
export MODELOS
ERR=0

# ---------------------------------------------------------------- utilities
viejo() {  # file reference label: OLD if older than the reference, MISSING if absent
  if [ ! -e "$1" ]; then echo "    MISSING $3"; return 0; fi
  if [ "$1" -ot "$2" ]; then echo "    OLD     $3"; return 0; fi
  return 1
}
paso() {  # number description command...: test mode shows it; run mode executes it and stops on failure
  local n="$1" desc="$2"; shift 2
  echo "-- $n. $desc"
  if [ "$MODO" = prueba ]; then echo "    \$ $*"; return; fi
  local log="$L/paso$n.log"
  echo "   $(date +%T) start -> $log"
  "$@" > "$log" 2>&1; local rc=$?
  if [ $rc -ne 0 ] || grep -q "exit [1-9]" "$log"; then
    echo "   ERROR in step $n (exit $rc); lines with errors:"; grep "exit [1-9]" "$log" | head
    echo "== chain stopped $(date '+%F %T')"; exit 1
  fi
  echo "   $(date +%T) end, exit 0"
}
rot_de() { echo "$MODELOS" | grep "^$1|" | cut -d"|" -f2; }   # confidence .rot of a folder

# ---------------------------------------------------------------- checks (both modes)
echo "== run_chain.sh, mode $MODO, $(date '+%F %T')"
echo "-- checks"
for s in run_models.sh run_shortening_step07.sh run_sens9.sh run_bic_noforce.sh run_sens_selkirk.sh \
         run_summaries.sh run_figures.sh andean_shortening_utrecht.py sens9_summary.py \
         sens9_deposits.py sens9_wagenaar.py sens_selkirk.py build_rot_wagenaar.py nominal_collision.py \
         merge_mc_parts_tag.py manuscript_numbers.py; do
  [ -f "$A/$s" ] || { echo "    MISSING script $s"; ERR=1; }
done
for s in "$S/make_supplementary_tables.py" "$S/make_supplementary_excel.py" "$P/tables/make_present_bands.py" \
         "$P/models/Q22_WK08-A/05_vc_uncertainty.py"; do
  [ -f "$s" ] || { echo "    MISSING script $s"; ERR=1; }
done
for line in $MODELOS; do IFS="|" read -r d rc rp prior cov <<< "$line"
  for f in "rot/$rc" "rot/$rp" hotspot_fit.json hotspot_fit_gauss_pred_age_E.json; do
    [ -f "$P/models/$d/$f" ] || { echo "    MISSING $d/$f"; ERR=1; }
  done
  # base case: without the present position as a 0 Ma observation (mandatory)
  py -c "import json,sys; sys.exit(1 if json.load(open(sys.argv[1]))['meta']['present_obs'] else 0)" \
    "$P/models/$d/hotspot_fit.json" || { echo "    $d: hotspot_fit.json uses the present position (not the base case)"; ERR=1; }
  # meta.force: some fits were written with PAPER_HS_FORCE=Easter=deriva; without FORCE the
  # BIC chooses the drift in all twelve and fit_hotspots.py reproduces them with difference 0. Reported only.
  F=$(py -c "import json,sys; print(json.load(open(sys.argv[1]))['meta'].get('force') or '')" "$P/models/$d/hotspot_fit.json")
  [ -n "$F" ] && NF=$((${NF:-0}+1))
done
[ $ERR -eq 0 ] && echo "    scripts, .rot and fits of the twelve folders: ok (no present position in any)"
[ -n "$NF" ] && echo "    note: $NF fits with meta.force (written with the Easter drift forced; fit_hotspots.py without FORCE reproduces them with difference 0)"
[ $ERR -ne 0 ] && [ "$MODO" = correr ] && { echo "== chain stopped: fix the above"; exit 1; }

# ---------------------------------------------------------------- old outputs
echo "-- outputs older than the 03b fit of their folder (or missing)"
NV=0
for line in $MODELOS; do d="${line%%|*}"; R="$P/models/$d/hotspot_fit.json"
  for v in total_hsb_cc_g13 total_hsb_pred_cc_g13 hotspot_hsb_cc_g13 circuit_hsb_cc_g13; do
    viejo "$P/models/$d/mc_collision_bands_$v.csv" "$R" "$d step 04 $v" && NV=$((NV+1))
  done
  viejo "$P/models/$d/deposit_tracks_porphyries.csv" "$R" "$d step 07" && NV=$((NV+1))
done
# step 09 is not compared: it uses only the circuit plates of the _hsb .rot (the 03b rewrites 13-16) and
# the chain always runs it (about 7 min per circuit)
for d in $NUEVE; do R="$P/models/$d/hotspot_fit.json"
  viejo "$P/models/$d/deposit_tracks_porphyries_lag1.csv" "$R" "$d step 07 lag" && NV=$((NV+1))
  for t in _s9_lag1 _s9_noYORK2 _s9_noFORCE _s9_sinSelkirk; do
    viejo "$P/models/$d/hotspot_fit$t.json" "$R" "$d sensitivity$t" && NV=$((NV+1)); done
  viejo "$P/models/$d/mc_collision_bands_total_s9_interp0_tot044.csv" "$R" "$d sensitivity interp" && NV=$((NV+1))
done
R=$(ls -t "$P"/models/*/hotspot_fit.json | head -1)   # most recent fit of all
echo "    (summaries, figures and tables against the most recent fit: ${R#$P/})"
for f in "$A/outputs/andean_shortening/position_summary.csv" "$A/outputs/sens9/wagenaar_9models.csv" \
         "$A/outputs/sens9/nominal_lag_york_9models.csv" "$A/outputs/sens9/test_porphyries_variants.csv" \
         "$A/outputs/selkirk/fits_selkirk.csv" "$A/outputs/deposit_lag_9models.csv" \
         "$A/outputs/rotated_mean_9models.txt" "$G/F4_deposit_distance_numbers.csv" \
         "$P/figures/output/main/Fig2.pdf" "$S/tableS1.tex" "$S/TablesS_data.xlsx"; do
  viejo "$f" "$R" "${f#$P/}" && NV=$((NV+1))
done
echo "    $NV old or missing outputs"

# ---------------------------------------------------------------- chain
if [ "$MODO" = correr ]; then
  mkdir -p "$L"
  # step 04 outputs older than the fit: moved (not deleted) so that they are not skipped
  for line in $MODELOS; do d="${line%%|*}"
    if viejo "$P/models/$d/mc_collision_bands_total_hsb_cc_g13.csv" "$P/models/$d/hotspot_fit.json" "$d" > /dev/null; then
      B="$P/models/$d/backup_chain_$FECHA"; mkdir -p "$B"
      ( cd "$P/models/$d" && mv mc_hsb*_cc_g13_*.log mc_collision_bands_*_hsb*_cc_g13*.csv "$B/" 2>/dev/null )
      echo "   step 04 of $d moved to ${B#$P/}"
    fi
  done
fi

vel() {  # step 05 in Q22_WK08-A and 09 in the three circuits, with the confidence .rot
  ( cd "$P/models/Q22_WK08-A" && py 05_vc_uncertainty.py > 05_vc_cadena.log 2>&1 ); echo "05 Q22_WK08-A exit $?"
  for d in $CIRCUITOS; do
    ( cd "$P/models/$d" && PAPER_ROT_FILE="$P/models/$d/rot/$(rot_de $d)" PAPER_CIRCUIT_COV=b PAPER_INTERP_COV=norm \
        py 09_vc_uncertainty_circuit.py 400 > 09_vc_cadena.log 2>&1 ); echo "09 $d exit $?"
  done
}
lag() {  # folder lag: step 07 with the deposit lag
  local d="$1" L="$2" line rc rp prior cov
  line=$(echo "$MODELOS" | grep "^$d|"); IFS="|" read -r _ rc rp prior cov <<< "$line"
  ( cd "$P/models/$d" && PAPER_CIRCUIT_COV=b PAPER_TAG_EXTRA=_g13 PAPER_INTERP_COV=norm PAPER_DEP_LAG="$L" \
      PAPER_ROT_FILE="$P/models/$d/rot/$rc" PAPER_HS_FIT="$P/models/$d/hotspot_fit.json" PAPER_HS_PRIOR="$prior" \
      PAPER_FRAME_COV="$cov" py 07_deposit_tracks.py > "ds07_g13_lag$L.log" 2>&1 )
  echo "$(date +%T) 07 lag$L $d exit $?"
}
retraso() {
  local J=""; for d in $NUEVE; do for L in 1 2 5; do J+="$d $L"$'\n'; done; done
  printf "%s" "$J" | xargs -P "$NP" -L 1 bash -c 'lag "$@"' _
}
utrecht() {  # first the Utrecht positions and trench (they depend on the deposit database), then step 07
  ( cd "$A" && py andean_shortening_utrecht.py ); echo "shortening positions exit $?"; NP=3 bash "$A/run_shortening_step07.sh"; }
presentes() {  # bands at the present-day margin
  for line in $MODELOS; do d="${line%%|*}"
    for t in _hsb_cc_g13 _hsb_pred_cc_g13; do
      PAPER_PIPE="$d" PAPER_BANDS_SUFFIX="_${d}" PAPER_FIG_TAG="$t" PAPER_TABLES_DIR="$G" \
        py "$P/tables/make_present_bands.py" > /dev/null 2>&1; echo "present $d $t exit $?"
    done
  done
}
wagenaar() {  # Wagenaar et al. (2025) frame, with a backup of the previous outputs
  local W="$M/Q22_WVH25" B="$M/backup_chain_$FECHA"
  local ROT="$W/rot/Global_Model_R_BG_2025_WVH25_hsb.rot"
  mkdir -p "$B" && cp -p "$M"/*.csv "$M"/*.log "$B/" && cp -rp "$W" "$B/" || { echo "backup exit 1"; return; }
  ( cd "$A" && py build_rot_wagenaar.py > "$M/build_g13.log" 2>&1 ); echo "build exit $?"
  ( cd "$W" && env -u PAPER_HS_P0 -u PAPER_HS_FORCE PAPER_ROT_FILE="$ROT" PAPER_HS_ROT_OUT="$ROT" \
      py 03b_hotspot_fit.py > hotspot_fit_g13.log 2>&1 ); echo "03b exit $?"
  ( cd "$W" && PAPER_ROT_FILE="$ROT" PAPER_HS_FIT="$W/hotspot_fit.json" py 07_deposit_tracks.py > ds07_g13.log 2>&1 )
  echo "07 exit $?"
  ( cd "$W" && PAPER_ROT_FILE="$ROT" PAPER_HS_FIT="$W/hotspot_fit.json" \
      py "$A/nominal_collision.py" "$W" "$M/coll_WVH25_g13.csv" ); echo "coll exit $?"
}
sensibilidades() {
  NP=$NP bash "$A/run_sens9.sh"
  NP=2 bash "$A/run_bic_noforce.sh"
  wagenaar
  NP=2 bash "$A/run_sens_selkirk.sh"
  cd "$A" || return
  for s in sens9_summary sens9_wagenaar sens9_deposits sens_selkirk; do py $s.py > /dev/null; echo "$s exit $?"; done
}
tablas() { cd "$S" || return; py make_supplementary_tables.py; echo "tableS1 exit $?"; py make_supplementary_excel.py; echo "excel exit $?"; }
manuscript_numbers() {
  cd "$A" || return; py manuscript_numbers.py > "$L/manuscript_numbers.txt" 2>&1; tail -1 "$L/manuscript_numbers.txt"
  if grep -q "^DIF\|Traceback" "$L/manuscript_numbers.txt"; then grep "^DIF" "$L/manuscript_numbers.txt"; echo "manuscript_numbers exit 1"
  else echo "manuscript_numbers exit 0"; fi
}
export -f lag rot_de; export NUEVE CIRCUITOS NP

[ "$MODO" = prueba ] && echo "-- commands that --run would run (in order; nothing runs in this mode)"
paso 1  "convergence velocity (05 and 09)"               vel
paso 2  "step 04 and following, twelve models"          env NP=$NP bash "$A/run_models.sh"
paso 3  "step 07 with lag 1, 2, 5 Myr, nine models"     retraso
paso 4  "Andean shortening (Utrecht) and summary"       utrecht
paso 5  "bands at the present margin"                   presentes
paso 6  "supplementary sensitivities and summaries"    sensibilidades
paso 7  "summaries read by manuscript_numbers.py"       bash "$A/run_summaries.sh"
paso 8  "figures"                                      bash "$A/run_figures.sh"
paso 9  "supplementary tables"                         tablas
paso 10 "manuscript_numbers.py (stops if a number does not match)"   manuscript_numbers
if [ "$MODO" = prueba ]; then
  echo "   (step 1: ${CIRCUITOS// /, } with PAPER_ROT_FILE=<rot _hsb>, PAPER_CIRCUIT_COV=b, PAPER_INTERP_COV=norm, N 400)"
  echo "   (step 6, Wagenaar: backup in ${M#$P/}/backup_chain_<date>/)"
  echo "== end of test mode; to run: bash run_chain.sh --run > run_chain.log 2>&1"
else
  echo "== chain complete $(date '+%F %T')"
fi
