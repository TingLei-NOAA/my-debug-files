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

. /apps/lmod/lmod/init/sh

# Safer shell opts for /bin/sh (no pipefail here)

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

# 2) Wait briefly, then find rank-0 PID and host (pure /bin/sh)
sleep 8
pid=""
host=""
i=0
while [ $i -lt 20 ]; do
  set -- $(scontrol listpids "$SLURM_JOB_ID" | awk 'NR==2{print $2, $(NF-2)}')
  pid="${1:-}"
  host="${2:-}"
  [ -n "$pid" ] && [ -n "$host" ] && break
  sleep 1
  i=$((i+1))
done
echo "Attaching to host=$host pid=$pid" 1>&2

# 3) gdb command file (no variable expansion inside)
cat > gdb.cmds <<'EOF'
set pagination off
set print pretty on
set logging file gdb_rank0.log
set logging on
set backtrace limit 1000
rbreak MPI_Finalize
commands
  silent
  echo \n=== Hit MPI_Finalize (rank0) ===\n
  bt full
  continue
end
rbreak PMPI_Comm_dup
commands
  silent
  echo \n=== Hit PMPI_Comm_dup (rank0) ===\n
  bt full
  continue
end
handle SIGABRT stop print pass
break abort
commands
  silent
  echo \n=== Hit abort() (rank0) ===\n
  bt full
  continue
end
continue
EOF

# 4) Attach gdb to rank 0 on its node
#    IMPORTANT: override per-step resources so it doesnt inherit 228.
#    If your site disallows --overlap, you can remove it because we left slack.
srun --overlap -N1 -n1 -w "$host" --ntasks-per-node=1 --cpus-per-task=1 \
  gdb -q -p "$pid" \
    -ex "set solib-search-path $libdir" \
    -ex "sharedlibrary" \
    -x gdb.cmds

# 5) Wait for the main step to finish
wait "$main_step"
