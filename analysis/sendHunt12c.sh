#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 03:00:00
#SBATCH --partition=cpu12c

#SBATCH -o outputHunt_Ha1000_seq
#SBATCH -e errorHunt_Ha1000_seq
#SBATCH --job-name=Hunt_1000_seq
####SBATCH --chdir=/home/froca/blankets/GridapMHD_2/GridapMHD.jl/analysis
#SBATCH --mem=0


SLURM_NPROCS=`expr $SLURM_JOB_NUM_NODES \* $SLURM_NTASKS_PER_NODE`

srun hostname -s > hosts.$SLURM_JOB_ID
echo "================================================================"
hostname
echo "Using: ${SLURM_NPROCS} procs"
echo "================================================================"
echo ""


SECONDS=0

source env.sh
#export OMPI_MCA_btl_openib_allow_ib=1
#export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

#mpiexec -n ${SLURM_NPROCS} julia --project=.. -e\
julia --project=.. -e\
'
using GridapMHD: hunt
hunt(
  nc=(40,40),
#  np=(2,4),
#  backend=:mpi,
  L=1.0,
  B=(0.,1000.,0.),
  nsums = 800,
  debug=false,
  vtk=true,
  title="Hunt_Ha1000_seq",
  mesh = false,
  BL_adapted = true,
#  kmap_x = 2,
#  kmap_y = 3,
  solver=:julia,
#  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps",
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

