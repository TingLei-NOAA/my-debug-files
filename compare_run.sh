#!/bin/bash
# Compare the netCDF output files of several test runs with those of a control run.
#
# Edit the settings below, then run:  ./compare_runs.sh
#
# Each run's output prefix is read from the "prefix:" key of the top-level "output:" block
# of its YAML.  All output files are assumed to be in the current directory, named
# <prefix>.<anything>.nc.  For every file a test produced, the control file with the same
# name after the prefix is compared with it.  Files a test did not produce are skipped.
# Each comparison's output goes to compare_logs/<test prefix>.<file>.log; a summary is
# printed at the end.  The script exits 1 if any comparison differs or a control file is missing.

# ------------------------------- settings ------------------------------------------------
CTRL_YAML=deb2-hyb-dual-norm-vdl_v1-p1936.yaml                 # the control run

TEST_YAMLS=(                           # the test runs, one per line
  deb2-parallel_io-hyb-dual-norm-vdl_v1-p1936.yaml
  deb2-parallel_io_inplace-hyb-dual-norm-vdl_v1-p1936.yaml
  deb2-parallel_io_inplace_only_anavar-hyb-dual-norm-vdl_v1-p1936.yaml
)

# Comparison command, run as:  $CMP CONTROL_FILE TEST_FILE
# exit status 0 = identical, non-zero = different.
# Leave empty to use "python compare-nc.py" if compare-nc.py is here, else "nccmp -dfs".
CMP=""
# ------------------------------------------------------------------------------------------

set -u
if [ -z "$CMP" ]; then
  if [ -f compare-nc.py ]; then CMP="python compare-nc.py"
  elif command -v nccmp > /dev/null; then CMP="nccmp -dfs"
  else echo "no comparison tool: set CMP, or put compare-nc.py here"; exit 1
  fi
fi

# Value of KEY in the top-level "output:" block of a YAML file (quotes and comments removed)
output_key() {
  awk -v key="$2" '
    /^[^[:space:]#-][^:]*:/ { inout = ($0 ~ /^output:/); next }      # a top-level key
    inout && $0 ~ "^[[:space:]]+" key ":" {
      sub("^[[:space:]]+" key ":[[:space:]]*", "")
      sub(/[[:space:]]+#.*$/, "")
      gsub(/^["\047]|["\047]$/, "")
      print; exit
    }' "$1"
}

[ -f "$CTRL_YAML" ] || { echo "control YAML not found: $CTRL_YAML"; exit 1; }
CTRL_PREFIX=$(output_key "$CTRL_YAML" prefix)
[ -n "$CTRL_PREFIX" ] || { echo "no output prefix in $CTRL_YAML"; exit 1; }
mkdir -p compare_logs

echo "Control: $CTRL_YAML  (prefix $CTRL_PREFIX)"
echo "Compare: $CMP"
printf "\n%-45s %-40s %s\n" "test" "file" "result"

status=0
for TEST_YAML in "${TEST_YAMLS[@]}"; do
  if [ ! -f "$TEST_YAML" ]; then
    printf "%-45s %-40s %s\n" "$TEST_YAML" "-" "YAML NOT FOUND"; status=1; continue
  fi
  TEST_PREFIX=$(output_key "$TEST_YAML" prefix)
  if [ -z "$TEST_PREFIX" ]; then
    printf "%-45s %-40s %s\n" "$TEST_YAML" "-" "NO OUTPUT PREFIX IN YAML"; status=1; continue
  fi
  shopt -s nullglob
  files=( "$TEST_PREFIX".*.nc )
  shopt -u nullglob
  if [ ${#files[@]} -eq 0 ]; then
    printf "%-45s %-40s %s\n" "$TEST_YAML" "-" "no output files ($TEST_PREFIX.*.nc)"; continue
  fi
  for test_file in "${files[@]}"; do
    suffix=${test_file#"$TEST_PREFIX"}            # e.g. .20240527.010000.fv_core.res.nc
    ctrl_file=$CTRL_PREFIX$suffix
    log=compare_logs/$TEST_PREFIX$suffix.log
    if [ ! -f "$ctrl_file" ]; then
      result="CONTROL MISSING ($ctrl_file)"; status=1
    elif $CMP "$ctrl_file" "$test_file" > "$log" 2>&1; then
      result="identical"
    else
      result="DIFFERENT (see $log)"; status=1
    fi
    printf "%-45s %-40s %s\n" "$TEST_YAML" "${suffix#.}" "$result"
  done
done
exit $status

