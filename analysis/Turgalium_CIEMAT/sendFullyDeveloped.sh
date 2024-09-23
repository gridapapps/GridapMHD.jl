#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -t 10:00:00
#SBATCH --partition=volta

#SBATCH -o outputFD_test
#SBATCH -e errorFD_test
#SBATCH --job-name=FD_test
#SBATCH --mem=0


SLURM_NPROCS=`expr $SLURM_JOB_NUM_NODES \* $SLURM_NTASKS_PER_NODE`

srun hostname -s > hosts.$SLURM_JOB_ID
echo "================================================================"
hostname
echo "Using: ${SLURM_NPROCS} procs in ${SLURM_JOB_NUM_NODES} nodes"
echo "================================================================"
echo ""


SECONDS=0

source env.sh
#export OMPI_MCA_btl_openib_allow_ib=1
#export OMPI_MCA_btl_openib_if_include="mlx5_0:1"

#mpiexec -n ${SLURM_NPROCS} julia --project=$GRIDAPMHD -J $GRIDAPMHD/compile/Turgalium_CIEMAT/GridapMHD36c.so -O3 --check-bounds=no -e\
mpiexec -n ${SLURM_NPROCS} julia --project=$GRIDAPMHD -O3 --check-bounds=no -e\
'
using GridapMHD:FullyDeveloped 
FullyDeveloped(
  nc=(60,60),
  np=(2,2),
  backend=:mpi,
  Ha = 1000.0,
  b = 1.5,
  dir_B = (0.0,1.0,0.0),
  cw_s = 0.01,
  τ_Ha = 1e5,
  cw_Ha = 0.01,
  τ_s = 1e5,
  nsums = 100,
  debug = false,
  vtk = true,
  title="FD_test",
  mesh = false,
  solver=:petsc,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_28 1 -mat_mumps_icntl_29 2 -mat_mumps_icntl_4 3 -mat_mumps_cntl_1 0.001"
 )'

duration=$SECONDS
rm -f hosts.$SLURM_JOB_ID
#rm -f $MACHINE_FILE

STATUS=$?
echo "================================================================"
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "================================================================"
echo ""
echo "STATUS = $STATUS"
echo ""

