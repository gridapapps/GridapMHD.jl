#!/bin/bash
#SBATCH -N 2
#SBATCH --ntasks-per-node=20
#SBATCH -t 20:00:00
#SBATCH --partition=cpu36c

#SBATCH -o outputHunt_Ha1000_N2n20
#SBATCH -e errorHunt_Ha1000_N2n20
###SBATCH --mail-user=fernando.roca@ciemat.es
#SBATCH --job-name=Hunt_mpi_1000
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

mpiexec --np ${SLURM_NPROCS} julia --project=$GRIDAPMHD -e\
'
using GridapMHD: hunt
hunt(
  nc=(60,60),
  np=(5,8),
  backend=:mpi,
  L=1.0,
  B=(0.,1000.,0.),
  nsums = 1000,
  debug=false,
  vtk=true,
  title="hunt_1000_mumps",
  mesh = false,
  BL_adapted = true,
#  kmap_x = 2,
#  kmap_y = 3,
  solver=:petsc,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps",
#  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type gmres -pc_type hypre -pc_hypre_type pilut"
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

