#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 20:00:00
#SBATCH --partition=cpu36c

#SBATCH -o outputHunt_Ha500_mumps_n1
#SBATCH -e errorHunt_Ha500_mumps_n1
###SBATCH --mail-user=fernando.roca@ciemat.es
#SBATCH --job-name=Hunt_mumps_500
#SBATCH --mem=0


SLURM_NPROCS=`expr $SLURM_JOB_NUM_NODES \* $SLURM_NTASKS_PER_NODE`

srun hostname -s > hosts.$SLURM_JOB_ID
echo "================================================================"
hostname
echo "Using: ${SLURM_NPROCS} procs in ${SLURM_JOB_NUM_NODES} nodes"
echo "================================================================"
echo ""


SECONDS=0

#MACHINE_FILE="machinefile-${SLURM_JOB_ID}"

# Generate machine file
#echo $SLURM_NODELIST | scontrol show hostnames | sed 's/$/:'$SLURM_NTASKS_PER_NODE'/' > $MACHINE_FILE


source env.sh
#export OMPI_MCA_btl_openib_allow_ib=1
#export OMPI_MCA_btl_openib_if_include="mlx5_0:1"

#mpiexec -n ${SLURM_NPROCS} --mca btl_openib_allow_ib 1 --mca btl_openib_if_include mlx5_0 julia --project=$GRIDAPMHD -J $GRIDAPMHD/GridapMHD.so -O3 --check-bounds=no -e\
julia --project=$GRIDAPMHD -J $GRIDAPMHD/compile/Turgalium_CIEMAT/GridapMHD36c.so -O3 --check-bounds=no -e\
'
using GridapMHD: hunt
hunt(
  nc=(50,50),
#  np=(2,2),
#  backend=:mpi,
  L=1.0,
  B=(0.,500.,0.),
  nsums = 1000,
  debug=false,
  vtk=true,
  title="hunt_500_petsc_n1",
  mesh = false,
  BL_adapted = true,
  solver=:petsc,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps",
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

