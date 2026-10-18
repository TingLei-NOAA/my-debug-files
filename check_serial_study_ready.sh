#!/bin/bash
# Pre-flight check for the serial vs 10-group timing study (perf_pair_job_serial.sh). Read-only.
#   usage: check_serial_study_ready.sh     (run it from the directory that holds the study scripts)
# Every path, YAML name and rank count is read from perf_pair_job_serial.sh, so this checks exactly
# what the jobs will use. Ends with a count of errors (must be 0) and warnings (read them).

set -u
cd "$(dirname "$0")" || exit 1
JOB=perf_pair_job_serial.sh
NERR=0; NWARN=0
ok()   { echo "  ok    $*"; }
warn() { echo "  WARN  $*"; NWARN=$((NWARN + 1)); }
err()  { echo "  ERROR $*"; NERR=$((NERR + 1)); }
need_file() { if [ -e "$1" ]; then ok "$1"; else err "missing: $1"; fi; }

# value(s) of "key:" in a YAML, quotes and trailing comments removed (also matches "- key: value")
yaml_values() {
  sed -n "s/^[[:space:]]*\(-[[:space:]]*\)\{0,1\}$1:[[:space:]]*//p" "$2" \
    | sed 's/[[:space:]]#.*$//; s/^"//; s/"[[:space:]]*$//; s/[[:space:]]*$//' | sort -u
}

# ------------------------------------------------------------------ 1. study scripts
echo "=== 1. Study scripts in $(pwd)"
for s in perf_pair_job_serial.sh submit_perf_series_serial.sh collect_perf_stats_serial.sh; do
  if [ ! -f "$s" ]; then err "missing: $s"; continue; fi
  if grep -q $'\r' "$s"; then err "$s has Windows line endings (run: sed -i 's/\r$//' $s)"; fi
  if bash -n "$s" 2> /dev/null; then ok "$s: bash syntax"; else err "$s: bash syntax error (bash -n $s)"; fi
  [ -x "$s" ] || warn "$s is not executable (chmod +x $s)"
done
[ -f "$JOB" ] || { echo "cannot continue without $JOB"; exit 1; }

eval "$(grep -E '^(BASE_RUNDIR|YAML_A|YAML_B|EXEC|moduledir)=' "$JOB")"
eval "$(grep -E '^NPROC_A=' "$JOB")"
eval "$(grep -E '^NPROC_B=' "$JOB")"
NODES=$(sed -n 's/^#PBS -l select=\([0-9]*\):.*/\1/p' "$JOB")
echo "  run dir $BASE_RUNDIR"
echo "  A = $YAML_A ($NPROC_A ranks, $PPN_A per node)"
echo "  B = $YAML_B ($NPROC_B ranks, $PPN_B per node)"
for s in A B; do
  eval "np=\$NPROC_$s; pp=\$PPN_$s"
  nn=$(( (np + pp - 1) / pp ))
  if [ "$nn" -le "$NODES" ]; then ok "setup $s needs $nn nodes, job requests $NODES"
  else err "setup $s needs $nn nodes, job requests only $NODES"; fi
done

# ------------------------------------------------------------------ 2. executable and environment
echo; echo "=== 2. Executable and environment"
need_file "$EXEC"
[ -d "$moduledir" ] && ok "module dir $moduledir" || err "missing module dir $moduledir"
command -v qsub > /dev/null && ok "qsub found" || err "qsub not in PATH"
LIBDIR=$(dirname "$(dirname "$EXEC")")
SABERLIB=$(ls "$LIBDIR"/lib*/libsaber.so 2> /dev/null | head -1)
if [ -z "$SABERLIB" ]; then warn "libsaber.so not found under $LIBDIR/lib*"
elif strings "$SABERLIB" | grep -qi l_loc_filter_once; then ok "$SABERLIB contains l_loc_filter_once"
else err "$SABERLIB has no l_loc_filter_once: the build does not have the filter-once switch"; fi

# ------------------------------------------------------------------ 3. per YAML
check_yaml() {  # $1 = setup (A|B), $2 = yaml, $3 = ranks
  local s=$1 y=$BASE_RUNDIR/$2 np=$3 f v lay prod init
  echo; echo "=== 3$s. Setup $s: $2"
  [ -f "$y" ] || { err "missing: $y"; return; }
  if grep -q $'\r' "$y"; then warn "$2 has Windows line endings"; fi
  cd "$BASE_RUNDIR" || return

  echo "  -- fixed files and namelists"
  for k in "namelist filename" "field table filename" akbk fieldset; do
    for f in $(yaml_values "$k" "$y"); do need_file "$f"; done
  done

  echo "  -- FV3 layout of every geometry namelist (product must be $np, or $np/groups for group geometries)"
  for f in $(yaml_values "namelist filename" "$y" | grep -v fmsmpp); do
    [ -f "$f" ] || continue
    lay=$(grep -iE '^[[:space:]]*layout[[:space:]]*=' "$f" | head -1 | sed 's/.*=//; s/[[:space:]]//g; s/,$//')
    prod=$(echo "$lay" | awk -F, '{print $1 * $2}')
    if [ -z "$lay" ]; then warn "$f: no layout line"
    elif [ "$prod" = "$np" ]; then ok "$f: layout $lay = $prod (full communicator)"
    elif [ "$s" = B ] && [ "$prod" = $((np / 10)) ]; then ok "$f: layout $lay = $prod (one of 10 groups)"
    else err "$f: layout $lay = $prod does not fit $np ranks"; fi
  done
  v=$(sed -n 's/^_layout: &fv3_layout \[\(.*\)\].*/\1/p' "$y" | tr -d ' ')
  prod=$(echo "$v" | awk -F, '{print $1 * $2}')
  [ "$prod" = "$np" ] && ok "YAML fv3_layout [$v] = $prod" || err "YAML fv3_layout [$v] = $prod, job runs $np ranks"

  echo "  -- background and ensemble members"
  for f in $(yaml_values datapath "$y" | grep -v '%mem%' | grep -v '^\./*$'); do
    for v in $(yaml_values "filename_[a-z]*" "$y"); do
      [ -e "$f/$v" ] || err "missing background file $f/$v"
    done
  done
  v=$(yaml_values datapath "$y" | grep '%mem%' | head -1)
  local nmem=$(yaml_values nmembers "$y" | head -1) miss=0 m dir
  for m in $(seq -f %03g 1 "${nmem:-0}"); do
    dir=${v//%mem%/$m}
    for f in $(yaml_values "filename_\(core\|trcr\|cplr\)" "$y"); do
      [ -e "$dir/$f" ] || { miss=$((miss + 1)); [ $miss -le 5 ] && err "missing member file $dir/$f"; }
    done
  done
  [ $miss -eq 0 ] && ok "$nmem members complete under ${v}" || err "$miss member files missing"

  echo "  -- obs input files"
  miss=0
  for f in $(yaml_values obsfile "$y" | grep -v '^jdiag'); do [ -e "$f" ] || { miss=$((miss + 1)); warn "missing obs file $f"; }; done
  [ $miss -eq 0 ] && ok "all obs input files present"

  echo "  -- MGBF namelists"
  for init in $(yaml_values "mgbf sdl and vdl init namelist file" "$y"); do
    if [ ! -f "$init" ]; then err "missing MGBF init namelist $init"; continue; fi
    ok "init $init"
    grep -qi l_loc_filter_once "$init" && err "$init sets l_loc_filter_once: it belongs only in the group files (forrtl 19)"
    local ngrp=0 nfo=0 g gp
    # group files = quoted names in the init namelist; tried relative to the run dir and to the init's dir
    for g in $(grep -oE "['\"][^'\"]+\.nml['\"]" "$init" | tr -d "'\"" | sort -u); do
      gp=$g; [ -f "$gp" ] || gp=$(dirname "$init")/$g
      if [ ! -f "$gp" ]; then err "group namelist $g (named in $init) not found"; continue; fi
      ngrp=$((ngrp + 1))
      if grep -qiE '^[[:space:]]*l_loc_filter_once[[:space:]]*=[[:space:]]*\.true\.' "$gp"; then nfo=$((nfo + 1)); fi
      printf "  group %-55s %s\n" "$gp" "$(grep -ioE '^[[:space:]]*(km2|km3|l_for_localization|l_loc_filter_once)[[:space:]]*=[^,]*' "$gp" | tr -d ' ' | tr '\n' ' ')"
    done
    if [ $ngrp -eq 0 ]; then err "no group namelists found in $init"
    elif [ $nfo -eq $ngrp ]; then ok "filter-once ON in all $ngrp group namelists"
    else err "filter-once ON in only $nfo of $ngrp group namelists (study needs it on in all, in both setups)"; fi
  done
  cd - > /dev/null
}
check_yaml A "$YAML_A" "$NPROC_A"
check_yaml B "$YAML_B" "$NPROC_B"

# ------------------------------------------------------------------ 4. results directory
echo; echo "=== 4. Results directory"
RES=$(dirname "$BASE_RUNDIR")/perf_results_serial
if ls "$RES"/*.meta > /dev/null 2>&1; then
  warn "$RES already holds $(ls "$RES"/*.meta | wc -l) runs; the collector would mix them with the new ones"
else ok "$RES is empty or does not exist yet"; fi

echo; echo "=== Result: $NERR errors, $NWARN warnings"
[ $NERR -eq 0 ]
