#!/bin/bash
# Collect the timing study written by perf_pair_job.sh: every timer of every run, and their statistics.
#
#   usage: collect_perf_stats_serial.sh [RESULTS_DIR] [THRESHOLD_S]
#          RESULTS_DIR  default /lfs/h2/emc/da/noscrub/Ting.Lei/dr-NA-Ens3Dvar-MGBF/perf_results_serial
#          THRESHOLD_S  screen summary lists timers whose mean paired difference exceeds this (default 1 s)
#
# Setups: A = 1936 ranks / serial ensemble (no groups),  B = 1960 ranks / 10 ensemble groups.
# The two B paths differ (A: oops::EnsembleCovariance, B: SABER parallel block chain), so many
# covariance timers exist in one setup only; compare the timers both setups report.
#
# Everything is read from rank 0's OOPS_STATS output of each run. Four kinds of values ("table"):
#   r0       "Timing Statistics": rank 0's total time and call count of every timer
#   par      "Parallel Timing Statistics": min/max/avg over ranks (only the timers OOPS lists there);
#            the value used is the MAX over ranks, i.e. the critical path
#   mile     milestones "<name> - Runtime: <s> sec": seconds since process start
#   derived  NonB = r0 oops::Run::execute - HybridCovariance - hybrid::multiply
#            (everything outside B setup and B application)
# Timer names are cut to 51 characters by OOPS; a name repeated within a table gets "#2", "#3", ...
#
# Files written to RESULTS_DIR (all tab-separated, for any later, narrower analysis):
#   perf_runs.tsv          one row per run: job, pair, setup, position, rc, wall, final J, Run end
#   perf_timers_long.tsv   one row per run and timer: every value read
#   perf_stats_all.tsv     per table, timer and setup: n, mean, sd, cv, median, min, max (+ call counts)
#   perf_paired_all.tsv    per table and timer: paired B vs A on the same nodes (diff, ratio, 95% interval)
# Screen summary: sanity checks, milestones, timers with a paired difference above THRESHOLD_S,
# per-job (node set) means of the main timers, outliers.

set -u
RESULTS=${1:-/lfs/h2/emc/da/noscrub/Ting.Lei/dr-NA-Ens3Dvar-MGBF/perf_results_serial}
THR=${2:-1}
RUNS=$RESULTS/perf_runs.tsv
LONG=$RESULTS/perf_timers_long.tsv
STATS=$RESULTS/perf_stats_all.tsv
PAIRED=$RESULTS/perf_paired_all.tsv
ls "$RESULTS"/*.meta > /dev/null 2>&1 || { echo "no .meta files in $RESULTS"; exit 1; }

# ---------------------------------------------------------------- stage 1: read every run
printf "job\tpair\tsetup\tpos\trc\twall_s\tstart\tJ\tRunEnd_s\n" > "$RUNS"
printf "job\tpair\tsetup\tpos\ttable\ttimer\tvalue_s\tcount\tmin_s\tmax_s\tavg_s\tpct_total\timbal_pct\n" > "$LONG"
for meta in "$RESULTS"/*.meta; do
  awk -v metafile="$meta" -v runsfile="$RUNS" -v longfile="$LONG" '
  function trim(s) { sub(/^[ \t]+/, "", s); sub(/[ \t]+$/, "", s); return s }
  function isnum(s) { return s ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/ }
  function emit(table, name, value, count, mn, mx, av, pc, im) {
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", M["job"], M["pair"], M["setup"], M["pos"], \
      table, name, value, count, mn, mx, av, pc, im >> longfile
  }
  BEGIN {
    getline line < metafile
    n = split(line, kv, "\t")
    for (i = 1; i <= n; i++) { split(kv[i], q, "="); M[q[1]] = q[2] }
    runend = "NA"; J = "NA"; sect = ""
  }
  {
    # keep rank-0 lines only ("<host> 0: ..."); unlabelled lines are accepted as they are
    if ($0 ~ /^[^ ]+ [0-9]+: /) { if ($2 != "0:") next; sub(/^[^ ]+ 0: /, "") }
  }
  /Nonlinear J *=/ { s = $0; sub(/^.*Nonlinear J *= */, "", s); split(s, w, /[ ,]/); if (isnum(w[1])) J = w[1] }
  !/OOPS_STATS/ { next }
  { sub(/^.*OOPS_STATS ?/, "") }
  # section banners appear at the start and at the end of each section
  /-+ Parallel Timing Statistics/ { sect = (sect == "par") ? "" : "par"; next }
  /-+ Timing Statistics -+/       { sect = (sect == "r0")  ? "" : "r0";  next }
  /-+ Object counts -+/           { sect = (sect == "obj") ? "" : "obj"; next }
  / - Runtime:/ {
    name = trim(substr($0, 1, index($0, " - Runtime:") - 1))
    split(substr($0, index($0, " - Runtime:") + 11), w, " ")
    if (isnum(w[1])) {
      if (++mc[name] > 1) name = name "#" mc[name]
      emit("mile", name, w[1], "NA", "NA", "NA", "NA", "NA", "NA")
      if (name == "Run end") runend = w[1]
    }
    next
  }
  (sect == "r0" || sect == "par") && match($0, / : /) {
    name = trim(substr($0, 1, RSTART - 1))
    k = split(substr($0, RSTART + 3), w, " ")
    for (i = 1; i <= k; i++) if (!isnum(w[i])) next
    if (sect == "r0" && k == 3) {
      if (++rc[name] > 1) name = name "#" rc[name]
      R0[name] = w[1] / 1000
      emit("r0", name, sprintf("%.4f", w[1] / 1000), w[2], "NA", "NA", "NA", "NA", "NA")
    } else if (sect == "par" && k == 5) {
      if (++pc[name] > 1) name = name "#" pc[name]
      emit("par", name, sprintf("%.4f", w[2] / 1000), "NA", sprintf("%.4f", w[1] / 1000), \
           sprintf("%.4f", w[2] / 1000), sprintf("%.4f", w[3] / 1000), w[4], w[5])
    }
  }
  END {
    # both setups wrap their ensemble B in the same top-level hybrid, so use its multiply
    ex = "oops::Run::execute"; hc = "oops::Covariance::HybridCovariance"; bm = "oops::Covariance::hybrid::multiply"
    if ((ex in R0) && (hc in R0) && (bm in R0))
      emit("derived", "NonB", sprintf("%.4f", R0[ex] - R0[hc] - R0[bm]), "NA", "NA", "NA", "NA", "NA", "NA")
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", M["job"], M["pair"], M["setup"], M["pos"], M["rc"], \
      M["wall_s"], M["start"], J, runend >> runsfile
  }' "${meta%.meta}.log"
done
echo "runs:   $RUNS ($(($(wc -l < "$RUNS") - 1)) runs)"
echo "values: $LONG ($(($(wc -l < "$LONG") - 1)) rows)"

# ---------------------------------------------------------------- stage 2: statistics
awk -F'\t' -v thr="$THR" -v statsfile="$STATS" -v pairedfile="$PAIRED" '
function tq(df,   t) {   # two-sided 95% Student t quantile
  split("12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 2.228 2.201 2.179 2.160 2.145 2.131 2.120 2.110 2.101 2.093 2.086 2.080 2.074 2.069 2.064 2.060 2.056 2.052 2.048 2.045 2.042", t, " ")
  return (df >= 1 && df <= 30) ? t[df] : 1.96
}
function sortv(a, n,   i, j, x) { for (i = 2; i <= n; i++) { x = a[i]; for (j = i - 1; j >= 1 && a[j] > x; j--) a[j + 1] = a[j]; a[j + 1] = x } }
function median(a, n) { sortv(a, n); return (n % 2) ? a[(n + 1) / 2] : (a[n / 2] + a[n / 2 + 1]) / 2 }
function stats(a, n,   i, s, ss, b) {   # sets MEAN SD MED MIN MAX
  s = 0; for (i = 1; i <= n; i++) s += a[i]; MEAN = s / n
  ss = 0; for (i = 1; i <= n; i++) ss += (a[i] - MEAN)^2; SD = (n > 1) ? sqrt(ss / (n - 1)) : 0
  for (i = 1; i <= n; i++) b[i] = a[i]; MED = median(b, n); MIN = b[1]; MAX = b[n]
}
FNR == 1 { next }
NR == FNR {   # perf_runs.tsv
  id = $1 SUBSEP $2 SUBSEP $3
  nrun++; rid[nrun] = id; rjob[nrun] = $1; rpair[nrun] = $2; rset[nrun] = $3; rpos[nrun] = $4
  rrc[nrun] = $5; rJ[nrun] = $8; rend[nrun] = $9
  ok[id] = ($5 == "0" && $9 != "NA")
  if (!($1 in seenjob)) { seenjob[$1] = 1; jobs[++njob] = $1 }
  next
}
{             # perf_timers_long.tsv
  id = $1 SUBSEP $2 SUBSEP $3
  if (!ok[id]) next
  k = $5 SUBSEP $6
  if (!(k in seenk)) { seenk[k] = 1; keys[++nk] = k }
  V[k, id] = $7 + 0
  if ($8 != "NA") C[k, id] = $8 + 0
}
END {
  # ---- 1. Sanity
  print "\n=== 1. Sanity   (A = 1936 ranks/serial, B = 1960 ranks/10 groups)"
  for (r = 1; r <= nrun; r++) {
    s = rset[r]; total[s]++
    if (ok[rid[r]]) {
      good[s]++
      if (!((s, rJ[r]) in Jc)) Jl[s] = Jl[s] " " rJ[r]
      Jc[s, rJ[r]]++
    } else printf "  NOT USED: %s pair %s setup %s (rc=%s, Run end %s)\n", rjob[r], rpair[r], s, rrc[r], rend[r]
  }
  for (si = 1; si <= 2; si++) {
    S = (si == 1) ? "A" : "B"
    printf "  setup %s: %d runs, %d complete; final J:", S, total[S], good[S]
    n = split(Jl[S], jj, " "); for (i = 1; i <= n; i++) printf " %s (x%d)", jj[i], Jc[S, jj[i]]
    printf "%s\n", (n > 1) ? "   <-- DIFFERENT: these runs did not do the same work" : ""
  }
  # call counts must not change between runs of one setup
  ndiff = 0
  for (ki = 1; ki <= nk; ki++) for (si = 1; si <= 2; si++) {
    S = (si == 1) ? "A" : "B"; cmin = ""; cmax = ""
    for (r = 1; r <= nrun; r++) if (rset[r] == S && ((keys[ki], rid[r]) in C)) {
      c = C[keys[ki], rid[r]]; if (cmin == "" || c < cmin) cmin = c; if (cmax == "" || c > cmax) cmax = c
    }
    CMIN[ki, S] = cmin; CMAX[ki, S] = cmax
    if (cmin != cmax) { split(keys[ki], kk, SUBSEP); if (++ndiff <= 10) printf "  call count varies within setup %s: %s (%s..%s)\n", S, kk[2], cmin, cmax }
  }
  printf "  timers whose call count varies between runs of one setup: %d\n", ndiff

  # ---- per setup statistics of every value -> perf_stats_all.tsv
  printf "table\ttimer\tsetup\tn\tmean_s\tsd_s\tcv_pct\tmedian_s\tmin_s\tmax_s\tcount_min\tcount_max\n" > statsfile
  for (ki = 1; ki <= nk; ki++) {
    split(keys[ki], kk, SUBSEP)
    for (si = 1; si <= 2; si++) {
      S = (si == 1) ? "A" : "B"; n = 0
      for (r = 1; r <= nrun; r++) if (rset[r] == S && ((keys[ki], rid[r]) in V)) a[++n] = V[keys[ki], rid[r]]
      if (n == 0) continue
      stats(a, n); delete a
      SMEAN[ki, S] = MEAN
      printf "%s\t%s\t%s\t%d\t%.4f\t%.4f\t%.2f\t%.4f\t%.4f\t%.4f\t%s\t%s\n", kk[1], kk[2], S, n, MEAN, SD, \
        (MEAN > 0 ? 100 * SD / MEAN : 0), MED, MIN, MAX, (CMIN[ki, S] == "" ? "NA" : CMIN[ki, S]), \
        (CMAX[ki, S] == "" ? "NA" : CMAX[ki, S]) > statsfile
    }
  }
  close(statsfile)

  # ---- paired B vs A within each pair (same nodes) -> perf_paired_all.tsv
  printf "table\ttimer\tn_pairs\tmeanA_s\tmeanB_s\tdiff_mean_s\tdiff_ci_lo\tdiff_ci_hi\tratio_mean\tratio_sd\tratio_ci_lo\tratio_ci_hi\tsignificant\n" > pairedfile
  for (r = 1; r <= nrun; r++) { pk = rjob[r] SUBSEP rpair[r]; if (rset[r] == "A") PA[pk] = rid[r]; else PB[pk] = rid[r] }
  for (ki = 1; ki <= nk; ki++) {
    split(keys[ki], kk, SUBSEP); n = 0; m = 0; sa = 0; sb = 0
    for (pk in PA) if (pk in PB) {
      ia = PA[pk]; ib = PB[pk]
      if (!(ok[ia] && ok[ib] && ((keys[ki], ia) in V) && ((keys[ki], ib) in V))) continue
      va = V[keys[ki], ia]; vb = V[keys[ki], ib]
      d[++n] = vb - va; sa += va; sb += vb
      if (va > 0) q[++m] = vb / va
    }
    if (n == 0) continue
    stats(d, n); dm = MEAN; dh = (n > 1) ? tq(n - 1) * SD / sqrt(n) : 0
    if (m > 0) { stats(q, m); qm = MEAN; qs = SD; qh = (m > 1) ? tq(m - 1) * SD / sqrt(m) : 0 } else { qm = qs = qh = 0 }
    sig = (n > 1 && (dm - dh > 0 || dm + dh < 0)) ? "yes" : "no"
    line = sprintf("%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\t%s\t%s\t%s", kk[1], kk[2], n, sa / n, sb / n, dm, dm - dh, dm + dh, \
      (m ? sprintf("%.4f", qm) : "NA"), (m ? sprintf("%.4f", qs) : "NA"), (m ? sprintf("%.4f", qm - qh) : "NA"), (m ? sprintf("%.4f", qm + qh) : "NA"), sig)
    print line > pairedfile
    PL[ki] = line; PD[ki] = dm; PN[ki] = n
    delete d; delete q
  }
  close(pairedfile)

  # ---- 2. Milestones
  print "\n=== 2. Milestones: seconds since start, paired over same-node pairs (B - A < 0: B earlier)"
  printf "%-40s %3s %10s %10s %9s %21s %4s\n", "milestone", "n", "mean A", "mean B", "B - A", "95% interval", "sig"
  for (ki = 1; ki <= nk; ki++) {
    split(keys[ki], kk, SUBSEP)
    if (kk[1] != "mile" && kk[1] != "derived") continue
    if (!(ki in PL)) continue
    split(PL[ki], f, "\t")
    printf "%-40s %3d %10.2f %10.2f %+9.2f   [%+8.2f, %+8.2f] %4s\n", (kk[1] == "derived" ? kk[2] " (derived)" : kk[2]), f[3], f[4], f[5], f[6], f[7], f[8], f[13]
  }

  # ---- 3. Timers with a mean paired difference above thr
  printf "\n=== 3. Timers whose mean paired difference |B - A| > %s s, largest first (table: par = max over ranks, r0 = rank 0)\n", thr
  printf "%-4s %-51s %3s %10s %10s %9s %21s %7s %4s\n", "tab", "timer", "n", "mean A", "mean B", "B - A", "95% interval", "B/A", "sig"
  cmd = "sort -t \"|\" -k1,1gr | cut -d \"|\" -f2-"
  fflush()   # print the header before sort writes its output
  for (ki = 1; ki <= nk; ki++) {
    split(keys[ki], kk, SUBSEP)
    if ((kk[1] != "r0" && kk[1] != "par") || !(ki in PL)) continue
    if (PD[ki] <= thr && PD[ki] >= -thr) continue
    split(PL[ki], f, "\t")
    printf "%.4f|%-4s %-51s %3d %10.2f %10.2f %+9.2f   [%+8.2f, %+8.2f] %7s %4s\n", (PD[ki] < 0 ? -PD[ki] : PD[ki]), \
      kk[1], kk[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[13] | cmd
  }
  close(cmd)

  # ---- 4. Per job (node set)
  print "\n=== 4. Per job (each job = one node set): mean over its complete runs [s]"
  ns = split("mile|Run end;par|oops::Run::execute;par|oops::Covariance::hybrid::multiply;par|saber::mgbf::Covariance::multiply;par|saber::SaberHybridBlockChain::parallelScatter;par|saber::SaberHybridBlockChain::parallelGatherAndSum;par|oops::State::State;derived|NonB", sel, ";")
  split("RunEnd RunExec Bmult MGBFmult Scatter Gather StateIO NonB", sh, " ")
  printf "%-8s %-3s %3s", "job", "set", "n"; for (k = 1; k <= ns; k++) printf " %10s", sh[k]; printf "\n"
  for (j = 1; j <= njob; j++) for (si = 1; si <= 2; si++) {
    S = (si == 1) ? "A" : "B"; nr = 0
    for (r = 1; r <= nrun; r++) if (rjob[r] == jobs[j] && rset[r] == S && ok[rid[r]]) nr++
    printf "%-8s %-3s %3d", jobs[j], S, nr
    for (k = 1; k <= ns; k++) {
      split(sel[k], p, "|"); key = p[1] SUBSEP p[2]; n = 0; sum = 0
      for (r = 1; r <= nrun; r++) if (rjob[r] == jobs[j] && rset[r] == S && ((key, rid[r]) in V)) { n++; sum += V[key, rid[r]] }
      printf " %10s", (n ? sprintf("%.2f", sum / n) : "NA")
    }
    printf "\n"
  }

  # ---- 5. Outliers on Run end
  print "\n=== 5. Outliers: Run end more than 3 scaled MADs from the setup median"
  none = 1; key = "mile" SUBSEP "Run end"
  for (si = 1; si <= 2; si++) {
    S = (si == 1) ? "A" : "B"; n = 0
    for (r = 1; r <= nrun; r++) if (rset[r] == S && ((key, rid[r]) in V)) { a[++n] = V[key, rid[r]]; idx[n] = r }
    if (n >= 3) {
      for (i = 1; i <= n; i++) b[i] = a[i]; med = median(b, n)
      for (i = 1; i <= n; i++) b[i] = (a[i] > med) ? a[i] - med : med - a[i]; mad = 1.4826 * median(b, n)
      for (i = 1; i <= n; i++) if (mad > 0 && (a[i] - med > 3 * mad || med - a[i] > 3 * mad)) {
        r = idx[i]; none = 0
        printf "  %s pair %s setup %s: %.1f s (median %.1f, MAD %.1f), position %s in its pair\n", rjob[r], rpair[r], S, a[i], med, mad, rpos[r]
      }
    }
    delete a; delete b; delete idx
  }
  if (none) print "  (none)"
  printf "\nAll values: %s\nPer setup statistics: %s\nPaired comparison:    %s\n", "perf_timers_long.tsv", "perf_stats_all.tsv", "perf_paired_all.tsv"
}' "$RUNS" "$LONG"
