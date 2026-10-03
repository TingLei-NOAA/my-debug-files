#!/bin/bash
#PBS -A RRFS-DEV
#PBS -q dev
#PBS -l walltime=01:45:00
#PBS -l select=484:mpiprocs=128:ncpus=128
#PBS -l place=vscatter:exclhost
#PBS -N deb-ens3vardual_12kmS1-qsubv1-30mem-halo_ob-1936-mgbf 
##PBS -l debug=true
#PBS -j oe -o deb-ens3dvardual_12kmS1-qsubv1-30mem-halo_ob-p1936.log


#/lfs/h2/emc/da/noscrub/Shun.Liu/rrfs/testD/ufs-srweather-app/regional_workflow/ush/load_modules_run_task.sh run_anal_gsi /lfs/h2/emc/da/noscrub/Shun.Liu/rrfs/testD/ufs-srweather-app/regional_workflow/jobs/JREGIONAL_RUN_ANAL

#source /lfs/h2/emc/lam/noscrub/emc.lam/rrfs/v0.6.9/ufs-srweather-app/env/build_wcoss2_intel.env
moduledir="/lfs/h2/emc/da/noscrub/Ting.Lei/dr-sams-rdasapp/RDASApp/modulefiles"
module use $moduledir
module load RDAS/wcoss2.intel
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/opt/cray/pe/mpich/8.1.19/ofi/intel/19.0/lib"


ulimit -s unlimited
ulimit -a

echo "start time:"
date


#rundir=/lfs/h2/emc/ptmp/emc.lam/Shun.Liu/anal_conv_gsi_spinup
rundir=/lfs/h2/emc/da/noscrub/Ting.Lei/dr-NA-Ens3Dvar-MGBF/fv3_2024052700_nadomain

cd $rundir

#last_run=$(ls nodiagwrt_nosat_stdout.* 2>/dev/null | sort -n -t. -k 2 | tail -1)

EXEC=/lfs/h2/emc/da/noscrub/Ting.Lei/dr-jedi-bundle/jedi-bundle/build/bin/fv3jedi_var.x
export JEDI_LIBS="/lfs/h2/emc/da/noscrub/Ting.Lei/dr-jedi-bundle/jedi-bundle/build/lib64:/lfs/h2/emc/da/noscrub/Ting.Lei/dr-jedi-bundle/jedi-bundle/build/lib"
# add subproject libs if present:
# JEDI_LIBS="$JEDI_LIBS:/lfs/.../oops/build/lib64:/lfs/.../saber/build/lib64:/lfs/.../ufo/build/lib64"

# include MKL/compiler runtimes if you linked them:
module load intel-oneapi-mkl intel-oneapi-compiler  # or whatever sets MKLROOT
export LD_LIBRARY_PATH="$JEDI_LIBS:$MKLROOT/lib/intel64:$LD_LIBRARY_PATH"
export MPICH_OFI_STARTUP_CONNECT=1

# make sure mpiexec propagates env (Hydra/PALS understands -env or -genv):
mpiexec -n 1936 -ppn 4 --cpu-bind core --depth 16 --label \
  -env LD_LIBRARY_PATH "$LD_LIBRARY_PATH" \
  $EXEC \
   compare-p1936-3dvar.yaml   compare-p1936-3dvar.out
exit

