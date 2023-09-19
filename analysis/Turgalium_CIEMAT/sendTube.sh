#!/bin/bash
#SBATCH -N 2 
#SBATCH --ntasks-per-node=4
#SBATCH -t 48:00:00
#SBATCH --partition=cpu36c

#SBATCH -o outputTube_Ha500_Ins_N2n4
#SBATCH -e errorTube_Ha500_Ins_N2n4
###SBATCH --mail-user=fernando.roca@ciemat.es
#SBATCH --job-name=Tube_Ha500_Ins_N2n4
#SBATCH --mem=0

SLURM_NPROCS=`expr $SLURM_JOB_NUM_NODES \* $SLURM_NTASKS_PER_NODE`
#SLURM_NPROCS = '16'

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

#Generate a file for passing some SLURM parameters to julia

PASS_FILE="pass_params.jl"

echo NPROCS=$SLURM_NPROCS > $PASS_FILE
echo JOB_NAME=\"$SLURM_JOB_NAME\" >> $PASS_FILE
Ha=${SLURM_JOB_NAME##*Ha}
#Ha=${Ha%_tau*}
Ha=${Ha%_Ins*}
echo Hartmann=${Ha}.0 >> $PASS_FILE
#tau=${SLURM_JOB_NAME##*tau}
#tau=${tau%_*}
#echo Tau=${tau}.0 >> $PASS_FILE


#mpiexec -n ${SLURM_NPROCS}  julia --project=$GRIDAPMHD -J $GRIDAPMHD/compile/Turgalium_CIEMAT/GridapMHD36c.so -O3 --check-bounds=no -e\
mpiexec -n ${SLURM_NPROCS}  julia --project=$GRIDAPMHD -O3 --check-bounds=no -e\
'
include("pass_params.jl")
using GridapMHD: tube
tube(;
  mesh="64k", 
  np=NPROCS,
  backend=:mpi,
  Ha = Hartmann,
  N = 10000.0,
  inlet = :plane,
  cw = 0.0,
#  τ = Tau,
  debug=false,
  vtk=true,
  title=JOB_NAME,
  solver=:petsc,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_7 0 -mat_mumps_icntl_28 1 -mat_mumps_icntl_29 2 -mat_mumps_icntl_4 3 -mat_mumps_cntl_1 0.001"
 )'


duration=$SECONDS
rm -f hosts.$SLURM_JOB_ID
rm -f pass_params.jl

STATUS=$?
echo "================================================================"
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "================================================================"
echo ""
echo "STATUS = $STATUS"
echo ""
