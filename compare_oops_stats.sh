#!/bin/bash
# Compare the OOPS_STATS timing of two JEDI runs.
#
# Usage:  compare_oops_stats.sh OLD_LOG NEW_LOG [THRESHOLD_MS]
#
#   OLD_LOG, NEW_LOG  files holding the rank-0 OOPS_STATS lines of each run (labelled
#                     lines such as "nid003104... 0: OOPS_STATS ..." are fine)
#   THRESHOLD_MS      list timers whose total changed by more than this (default 200)
#
# Prints:
#   1. Milestones ("Run start", "Run end", ... - Runtime: <s> sec), seconds since start
#   2. Where the time went: whole process, the OOPS run, and the rest (start-up etc.)
#   3. Timers whose total changed by more than THRESHOLD_MS, largest change first
#   4. Timers present in only one of the runs
#
# Timer lines are "<name> : <total ms> <count> <average ms>".  Only lines with exactly
# these three numbers are used, and only the first occurrence of each name, so the
# object-count table and the second (parallel) statistics table are not mixed in.

set -u
if [ $# -lt 2 ]; then
  echo "usage: $0 OLD_LOG NEW_LOG [THRESHOLD_MS (default 200)]"
  exit 1
fi
OLD=$1
NEW=$2
THR=${3:-200}
for f in "$OLD" "$NEW"; do
  [ -r "$f" ] || { echo "cannot read $f"; exit 1; }
  grep -q "OOPS_STATS" "$f" || { echo "no OOPS_STATS lines in $f (use the rank-0 output file)"; exit 1; }
done

extract() { grep "OOPS_STATS" "$1" | sed 's/^.*OOPS_STATS //'; }

awk -v thr="$THR" -v oldname="$OLD" -v newname="$NEW" '
function trim(s) { sub(/^[ \t]+/, "", s); sub(/[ \t]+$/, "", s); return s }

# "<name> - Runtime: <seconds> sec, ..."  ->  mname, mval
function milestone(line,   p, w) {
  p = index(line, " - Runtime:")
  if (!p) return 0
  mname = trim(substr(line, 1, p - 1))
  split(substr(line, p + 11), w, " ")
  if (w[1] !~ /^[0-9.]+$/) return 0
  mval = w[1] + 0
  return 1
}

# "<name> : <total ms> <count> <average ms>"  ->  tname, ttot, tcnt
function timer(line,   n, i, w) {
  if (!match(line, / : /)) return 0
  tname = trim(substr(line, 1, RSTART - 1))
  n = split(substr(line, RSTART + 3), w, " ")
  if (n != 3) return 0
  for (i = 1; i <= 3; i++) if (w[i] !~ /^[0-9.eE+-]+$/) return 0
  ttot = w[1] + 0
  tcnt = w[2] + 0
  return 1
}

FNR == 1 { run++ }

milestone($0) {
  if (!((run, mname) in M)) {
    M[run, mname] = mval
    if (run == 1) mord[++nm] = mname
    else if (!((1, mname) in M)) mnew[++nmnew] = mname
  }
  next
}

timer($0) {
  if ((run, tname) in T) { dup[run]++; next }
  T[run, tname] = ttot
  C[run, tname] = tcnt
  if (run == 1) tord[++nt] = tname
  else if (!((1, tname) in T)) tnew[++ntnew] = tname
  next
}

END {
  printf "OLD: %s\nNEW: %s\n", oldname, newname

  # 1. Milestones
  printf "\n=== 1. Milestones (seconds since process start)\n"
  printf "%-45s %10s %10s %9s\n", "milestone", "old", "new", "change"
  for (i = 1; i <= nm; i++) {
    m = mord[i]
    if ((2, m) in M) printf "%-45s %10.2f %10.2f %+9.2f\n", m, M[1, m], M[2, m], M[2, m] - M[1, m]
    else             printf "%-45s %10.2f %10s\n", m, M[1, m], "(none)"
  }
  for (i = 1; i <= nmnew; i++) printf "%-45s %10s %10.2f\n", mnew[i], "(none)", M[2, mnew[i]]

  # 2. Where the time went
  printf "\n=== 2. Where the time went (seconds)\n"
  have_end = ((1, "Run end") in M) && ((2, "Run end") in M)
  have_run = ((1, "oops::Run::execute") in T) && ((2, "oops::Run::execute") in T)
  if (have_end) {
    dend = M[2, "Run end"] - M[1, "Run end"]
    printf "%-45s %+9.2f\n", "whole process (Run end)", dend
  }
  if (have_run) {
    drun = (T[2, "oops::Run::execute"] - T[1, "oops::Run::execute"]) / 1000
    printf "%-45s %+9.2f\n", "the OOPS run (oops::Run::execute)", drun
  }
  if (have_end && have_run)
    printf "%-45s %+9.2f\n", "outside the OOPS run (start-up, shutdown)", dend - drun
  if (!have_end || !have_run) printf "(Run end or oops::Run::execute missing in a run)\n"

  # 3. Timers that changed
  printf "\n=== 3. Timers whose total changed by more than %s ms (largest first)\n", thr
  printf "%10s %12s %12s %8s  %s\n", "change ms", "old ms", "new ms", "change", "timer  (* = call count differs)"
  cmd = "sort -t \"|\" -k1,1gr | cut -d \"|\" -f2-"
  for (i = 1; i <= nt; i++) {
    t = tord[i]
    if (!((2, t) in T)) continue
    d = T[2, t] - T[1, t]
    if (d <= thr && d >= -thr) continue
    pct = (T[1, t] > 0) ? sprintf("%+7.1f%%", 100 * d / T[1, t]) : "     n/a"
    flag = (C[1, t] != C[2, t]) ? sprintf(" *(%d -> %d calls)", C[1, t], C[2, t]) : ""
    printf "%.3f|%+10.0f %12.0f %12.0f %8s  %s%s\n", d, d, T[1, t], T[2, t], pct, t, flag | cmd
  }
  close(cmd)

  # 4. Timers in only one run
  printf "\n=== 4. Timers present in only one run\n"
  none = 1
  for (i = 1; i <= nt; i++) if (!((2, tord[i]) in T)) { printf "  old only: %s (%.0f ms)\n", tord[i], T[1, tord[i]]; none = 0 }
  for (i = 1; i <= ntnew; i++) { printf "  new only: %s (%.0f ms)\n", tnew[i], T[2, tnew[i]]; none = 0 }
  if (none) print "  (none)"

  printf "\n(%d timers compared; repeated timer names skipped: old %d, new %d)\n", nt, dup[1] + 0, dup[2] + 0
}
' <(extract "$OLD") <(extract "$NEW")

