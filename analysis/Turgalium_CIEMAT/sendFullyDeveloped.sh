#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 20:00:00
#SBATCH --partition=cpu36c

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
julia --project=$GRIDAPMHD -e\
'
using GridapMHD:FullyDeveloped 
FullyDeveloped(
  nc=(6,6),
#  np=(4,5),
#  backend=:mpi,
  Ha = 10.0,
  b = 1.5,
  dir_B=(0.0,1.0,0.0),
  cw_s = 1.0,
  cw_Ha = 0.0,
  debug=true,
  vtk=true,
  title="FD_test",
  mesh = true,
  solver=:julia,
#  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_14 1000",
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

