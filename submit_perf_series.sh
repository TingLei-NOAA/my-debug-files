#!/bin/bash
# Submit the timing study: NJOBS jobs x NPAIRS pairs. Each pair is one 1960-rank run (A) and one
# 980-rank run (B) on the same nodes; jobs alternate which setup goes first.
# The jobs may run at the same time: each one works in its own run directory.
#   usage: submit_perf_series.sh            (3 jobs x 4 pairs = 12 runs of each setup)
#          NJOBS=2 NPAIRS=2 submit_perf_series.sh
set -eu
cd "$(dirname "$0")"
NJOBS=${NJOBS:-3}
NPAIRS=${NPAIRS:-4}
for j in $(seq 1 "$NJOBS"); do
  if [ $((j % 2)) -eq 1 ]; then first=A; else first=B; fi
  id=$(qsub -N perf_job$j -o perf_job$j.log -v JOBTAG=job$j,FIRST=$first,NPAIRS=$NPAIRS perf_pair_job.sh)
  echo "job$j  first=$first  pairs=$NPAIRS  $id"
done
