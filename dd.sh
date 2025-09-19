#!/bin/bash
#SBATCH --account=fv3-cam
#SBATCH --qos=debug
#SBATCH -t 00:30:00
#SBATCH --job-name=RRFS_mpas_jedi
#SBATCH -o deb.log
#SBATCH --open-mode=truncate

#SBATCH --nodes=9
#SBATCH --ntasks=196
#SBATCH --ntasks-per-node=22      # leave CPU slack (192 cores/node  22*8=176 used)
#SBATCH --cpus-per-task=8
#SBATCH --exclusive

source /apps/lmod/lmod/init/bash

# Safer shell opts for bash
set -euo pipefail

ulimit -s unlimited
ulimit -v unlimited

# --- your env/toolchain ---
srcfile=/scratch3/NCEPDEV/fv3-cam/Ting.Lei/dr-jedi-bundle-current/src-jedi-bundle.sh
. "$srcfile"

export OMP_NUM_THREADS=1
export OOPS_TRACE=1
export OOPS_DEBUG=1
export I_MPI_DEBUG=10

# Paths  adjust if needed
jexec="/scratch3/NCEPDEV/fv3-cam/Ting.Lei/dr-debug/jedi-bundle/build/bin/mpasjedi_error_covariance_toolbox.x"
libdir="/scratch3/NCEPDEV/fv3-cam/Ting.Lei/dr-debug/jedi-bundle/build/lib64"
yaml="./dev-example-no_vert_flt-ps-dirac.yaml"
out="./deb.out"

# 1) Launch main MPI step (keeps per-node slack)
srun -l -n 196 --ntasks-per-node=22 --cpus-per-task=8 --exclusive \
  "$jexec" "$yaml" "$out" &
main_step=$!

# 2) Wait briefly, then find rank-0 PID and host
sleep 8
pid=""
host=""
i=0
while [ $i -lt 30 ]; do
  host=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | sed -n '1p')
  if [ -n "$host" ]; then
    pid=$(srun --overlap -N1 -n1 -w "$host" --ntasks-per-node=1 --cpus-per-task=1 \
      pgrep -fn "$jexec" 2>/dev/null || true)
  fi
  if [ -n "$pid" ] && [ -n "$host" ]; then
    break
  fi
  sleep 2
  i=$((i+1))
done

if [ -n "$pid" ] && [ -n "$host" ]; then
  echo "Attaching to host=$host pid=$pid" 1>&2

  cat > gdb.cmds <<'EOF'
set pagination off
set print pretty on
set logging file gdb_rank0.log
set logging on
set breakpoint pending on
set backtrace limit 1000
python
import gdb

def attach_commands(message):
    bps = gdb.breakpoints()
    if not bps:
        return
    bp = bps[-1]
    gdb.execute("commands {}\n  silent\n  echo {}\n\n  bt full\n  continue\nend".format(bp.number, message), to_string=True)

attach_commands('=== Hit MPI_Finalize (rank0) ===')
gdb.execute('rbreak MPI_Finalize', to_string=True)
attach_commands('=== Hit MPI_Finalize (rank0) ===')

gdb.execute('rbreak PMPI_Comm_dup', to_string=True)
attach_commands('=== Hit PMPI_Comm_dup (rank0) ===')

gdb.execute('handle SIGABRT stop print pass', to_string=True)

gdb.execute('break abort', to_string=True)
attach_commands('=== Hit abort() (rank0) ===')
end
continue
EOF

  srun --overlap -N1 -n1 -w "$host" --ntasks-per-node=1 --cpus-per-task=1 \
    gdb -q -p "$pid" \
      -ex "set solib-search-path $libdir" \
      -ex "sharedlibrary" \
      -x gdb.cmds
else
  echo "Skipping gdb attach; missing PID or host" 1>&2
fi

# 5) Wait for the main step to finish
wait "$main_step"