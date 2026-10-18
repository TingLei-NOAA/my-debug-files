#!/bin/bash
# Submit the timing study: NJOBS jobs x NPAIRS pairs. Each pair is one serial 1936-rank run (A) and
# one 10-group 1960-rank run (B) on the same nodes; jobs alternate which setup goes first.
# The jobs may run at the same time: each one works in its own run directory.
#   usage: submit_perf_series_serial.sh     (3 jobs x 4 pairs = 12 runs of each setup)
#          NJOBS=2 NPAIRS=2 submit_perf_series_serial.sh
set -eu
cd "$(dirname "$0")"
NJOBS=${NJOBS:-3}
NPAIRS=${NPAIRS:-4}
for j in $(seq 1 "$NJOBS"); do
  if [ $((j % 2)) -eq 1 ]; then first=A; else first=B; fi
  id=$(qsub -N perf_ser$j -o perf_ser$j.log -v JOBTAG=ser$j,FIRST=$first,NPAIRS=$NPAIRS perf_pair_job_serial.sh)
  echo "job$j  first=$first  pairs=$NPAIRS  $id"
done
