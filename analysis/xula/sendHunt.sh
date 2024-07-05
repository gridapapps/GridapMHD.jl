#!/bin/bash
#SBATCH -N 1 
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:10:00
#SBATCH --partition=cpu36c

#SBATCH -o outputHunt_Ha20
#SBATCH -e errorHunt_Ha20
###SBATCH --mail-user=fernando.roca@ciemat.es
#SBATCH --job-name=Hunt_Ha20
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

mpiexec -n ${SLURM_NPROCS} julia --project=$GRIDAPMHD -J $GRIDAPMHD/compile/Turgalium_CIEMAT/GridapMHD36c.so -O3 --check-bounds=no -e\
'
using GridapMHD: hunt
hunt(
  nc=(20,20),
  np=(2,2),
  backend=:mpi,
  L=1.0,
  B=(0.,20.,0.),
  nsums = 100,
  debug=false,
  vtk=true,
  title="hunt_Ha20",
  mesh = false,
  BL_adapted = true,
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

