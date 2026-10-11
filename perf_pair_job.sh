#!/bin/bash
#PBS -A RRFS-DEV
#PBS -q dev
#PBS -l walltime=05:00:00
#PBS -l select=62:mpiprocs=128:ncpus=128:ompthreads=4
#PBS -l place=vscatter:exclhost
#PBS -N perf_pair
#PBS -j oe

# Timing study, one job: run the two setups back to back NPAIRS times on THIS job's nodes.
#   A = 1960 ranks, 32 per node, 10 ensemble groups (dev_filter_once_10groups ... p1960)
#   B =  980 ranks, 16 per node,  5 ensemble groups (dev_filter_once_5groups  ... p980)
# Pair 1 runs FIRST then the other setup; the order flips every pair (AB, BA, AB, ...).
# Submit with submit_perf_series.sh, which sets JOBTAG, FIRST and NPAIRS.
#
# Outputs, in $RESULTS (shared by all jobs of the study):
#   <JOBTAG>_pair<i>_<A|B>.log   stdout/stderr of that run (rank-labelled, holds OOPS_STATS)
#   <JOBTAG>_pair<i>_<A|B>.meta  one tab-separated line: job, pair, setup, position, times, rc, nodes
#   <JOBTAG>.nodes               the node list of this job

set -u
JOBTAG=${JOBTAG:?JOBTAG not set - submit with submit_perf_series.sh}
FIRST=${FIRST:-A}
NPAIRS=${NPAIRS:-4}
WALL_SEC=${WALL_SEC:-18000}   # keep equal to the walltime requested above (5 h)

BASE_RUNDIR=/lfs/h2/emc/da/noscrub/Ting.Lei/dr-NA-Ens3Dvar-MGBF/fv3_2024052700_nadomain
YAML_A=dev_filter_once_10groups-ens3dvar-online_dual-vdl_dbz-p1960.yaml
YAML_B=dev_filter_once_5groups-ens3dvar-online_dual-vdl_dbz-p980.yaml
NPROC_A=1960; PPN_A=32
NPROC_B=980;  PPN_B=16

# Environment: identical to the single-run scripts
moduledir="/lfs/h2/emc/da/noscrub/Ting.Lei/dr-rdasapp/RDASApp/modulefiles"
module use $moduledir
module load RDAS/wcoss2.intel
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/opt/cray/pe/mpich/8.1.19/ofi/intel/19.0/lib"
ulimit -s unlimited
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=512M
export MPICH_ALLTOALLV_THROTTLE=484
EXEC=/lfs/h2/emc/da/noscrub/Ting.Lei/dr-jedi-bundle-dev/jedi-bundle/build/bin/fv3jedi_var.x
export JEDI_LIBS="/lfs/h2/emc/da/noscrub/Ting.Lei/dr-jedi-bundle-dev/jedi-bundle/build/lib64:/lfs/h2/emc/da/noscrub/Ting.Lei/dr-jedi-bundle-dev/jedi-bundle/build/lib"
module load intel-oneapi-mkl intel-oneapi-compiler
export LD_LIBRARY_PATH="$JEDI_LIBS:$MKLROOT/lib/intel64:$LD_LIBRARY_PATH"
export MPICH_OFI_STARTUP_CONNECT=1

# Private run directory for this job, a SIBLING of BASE_RUNDIR, so the YAMLs' relative paths
# ("../dr-convert2state/...") resolve exactly as in BASE_RUNDIR. Inputs are symlinked; outputs
# (obs diagnostics, analysis, .out, FMS logfiles) are written here, so concurrent jobs never
# overwrite each other or the files in BASE_RUNDIR.
STUDY_DIR=$(dirname "$BASE_RUNDIR")
RESULTS=$STUDY_DIR/perf_results
JOBDIR=$STUDY_DIR/perf_${JOBTAG}_${PBS_JOBID%%.*}
mkdir -p "$RESULTS" "$JOBDIR"
for f in "$BASE_RUNDIR"/*; do
  name=$(basename "$f")
  if [ ! -d "$f" ]; then
    # Never link a file a run could write: writing through the link would change BASE_RUNDIR.
    case "$name" in
      *.nc|*.out|*.log|*.res|logfile.*|core*|*.o[0-9]*|perf_*) continue ;;
    esac
  fi
  ln -s "$f" "$JOBDIR/$name"
done
cd "$JOBDIR" || exit 1

NODEFILE=$RESULTS/${JOBTAG}.nodes
sort -u "$PBS_NODEFILE" > "$NODEFILE"
NNODES=$(wc -l < "$NODEFILE")
echo "job $JOBTAG ($PBS_JOBID): $NNODES nodes, $NPAIRS pairs, first=$FIRST, run dir $JOBDIR"

declare -A DUR=([A]=0 [B]=0)   # longest observed duration of each setup, for the walltime guard

run_one() {  # $1 = setup (A|B), $2 = pair index, $3 = position in the pair (1|2)
  local s=$1 i=$2 pos=$3 yaml nproc ppn tag t0 t1 rc
  if [ "$s" = A ]; then yaml=$YAML_A; nproc=$NPROC_A; ppn=$PPN_A
  else                  yaml=$YAML_B; nproc=$NPROC_B; ppn=$PPN_B; fi
  tag=${JOBTAG}_pair${i}_${s}
  t0=$(date +%s)
  echo "$(date -u +%FT%TZ) start $tag ($nproc ranks, $ppn per node)"
  mpiexec -l --line-buffer -n $nproc -ppn $ppn --cpu-bind core --depth 4 --label \
    -env LD_LIBRARY_PATH "$LD_LIBRARY_PATH" \
    $EXEC $yaml ${tag}.out > "$RESULTS/${tag}.log" 2>&1
  rc=$?
  t1=$(date +%s)
  printf "job=%s\tpbs_jobid=%s\tpair=%d\tsetup=%s\tpos=%d\tnproc=%d\tppn=%d\tstart=%d\tend=%d\twall_s=%d\trc=%d\tnnodes=%d\tnodes=%s\n" \
    "$JOBTAG" "$PBS_JOBID" "$i" "$s" "$pos" "$nproc" "$ppn" "$t0" "$t1" $((t1 - t0)) "$rc" "$NNODES" "$NODEFILE" \
    > "$RESULTS/${tag}.meta"
  echo "$(date -u +%FT%TZ) end   $tag rc=$rc wall=$((t1 - t0)) s"
  [ $((t1 - t0)) -gt ${DUR[$s]} ] && DUR[$s]=$((t1 - t0))
}

other() { [ "$1" = A ] && echo B || echo A; }

for i in $(seq 1 "$NPAIRS"); do
  if [ $((i % 2)) -eq 1 ]; then s1=$FIRST; else s1=$(other "$FIRST"); fi
  s2=$(other "$s1")
  # Walltime guard: start a pair only if it can finish (10% margin + 5 min), once both
  # setups have been timed in this job.
  need=$(( (${DUR[A]} + ${DUR[B]}) * 11 / 10 + 300 ))
  if [ ${DUR[A]} -gt 0 ] && [ ${DUR[B]} -gt 0 ] && [ $((SECONDS + need)) -gt $WALL_SEC ]; then
    echo "stopping before pair $i: ${SECONDS} s used, a pair needs ~${need} s, walltime ${WALL_SEC} s"
    break
  fi
  run_one "$s1" "$i" 1
  run_one "$s2" "$i" 2
done
echo "job $JOBTAG done after $SECONDS s"
