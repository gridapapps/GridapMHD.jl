#!/bin/bash
#SBATCH -N 2
#SBATCH --ntasks-per-node=1
#SBATCH -t 24:00:00
#SBATCH --partition=volta

#SBATCH -o outputExp_Ha50_ser
#SBATCH -e errorExp_Ha50_ser
###SBATCH --mail-user=fernando.roca@ciemat.es
#SBATCH --job-name=Exp_Ha50_ser
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

#mpiexec -n ${SLURM_NPROCS} julia --project=.. -J ../GridapMHD.so -O3 --check-bounds=no -e\
mpiexec -n ${SLURM_NPROCS} julia --project=$GRIDAPMHD -J $GRIDAPMHD/compile/Turgalium_CIEMAT/GridapMHD36c.so -O3 --check-bounds=no -e\
'
using GridapMHD: expansion
expansion(;
  mesh="68k", 
  np=2,
  backend=:mpi,
  Ha = 50.0,
  N = 3740.0,
  cw = 0.01,
  debug=false,
  vtk=true,
  title="Expansion_Ha50_serial",
  solver=:julia,
#  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_7 0",
 )'



duration=$SECONDS
rm -f hosts.$SLURM_JOB_ID

STATUS=$?
echo "================================================================"
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "================================================================"
echo ""
echo "STATUS = $STATUS"
echo ""

