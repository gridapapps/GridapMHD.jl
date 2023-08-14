#!/bin/bash
#SBATCH -N 4 
#SBATCH --ntasks-per-node=4
#SBATCH -t 48:00:00
#SBATCH --partition=cpu36c

#SBATCH -o outputExp_Ha100_tau1000_N4n4
#SBATCH -e errorExp_Ha100_tau1000_N4n4
###SBATCH --mail-user=fernando.roca@ciemat.es
#SBATCH --job-name=Exp_Ha100_tau1000_N4n4
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
Ha=${Ha%_tau*}
echo Hartmann=${Ha}.0 >> $PASS_FILE
tau=${SLURM_JOB_NAME##*tau}
tau=${tau%_*}
echo Tau=${tau}.0 >> $PASS_FILE


mpiexec -n ${SLURM_NPROCS}  julia --project=$GRIDAPMHD -J $GRIDAPMHD/compile/Turgalium_CIEMAT/GridapMHD36c.so -O3 --check-bounds=no -e\
'
include("pass_parms.jl")
using GridapMHD: expansion
expansion(;
  mesh="68k", 
  np=NPROCS,
  backend=:mpi,
  Ha = Hartmann,
  N = 3740.0,
  cw = 0.028,
  τ = tau,
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

