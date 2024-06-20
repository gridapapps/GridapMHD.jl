#!/bin/bash
#SBATCH -N 2
#SBATCH --ntasks-per-node=8
#SBATCH -t 01:00:00
#SBATCH --partition=xula3

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

srun --mpi=pmix  julia --project=$GRIDAPMHD -O3 --check-bounds=no -e\
'
# make a local temp directory for the cache
d = strip(String(read(`mktemp -d`)))
mkdir(joinpath(d, "compiled"))
# copy the cache to the temp directory
# set the temp directory as the first depot so any (re)compilation will happen there and not interfere with other jobs
run(`rsync -au $(DEPOT_PATH[1])/compiled $d/compiled`)
pushfirst!(DEPOT_PATH, d)

using GridapPETSc
using SparseMatricesCSR
using GridapMHD:FullyDeveloped 

#Monolithic MUMPS
solver = Dict(
    :solver => :petsc,
    :matrix_type    => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type    => Vector{PetscScalar},
    :petsc_options  => "-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_28 1 -mat_mumps_icntl_29 2 -mat_mumps_icntl_4 3 -mat_mumps_cntl_1 0.001",
    :niter          => 100,
    :rtol           => 1e-5,
  )

FullyDeveloped(
  nc=(60,60),
  np=(4,4),
  backend=:mpi,
  Ha = 1000.0,
  b = 1.5,
  dir_B = (0.0,1.0,0.0),
  cw_s = 0.01,
  τ_Ha = 1e6,
  cw_Ha = 0.01,
  τ_s = 1e6,
  nsums = 100,
  debug = false,
  vtk = true,
  title="FD_test",
  mesh = false,
  solver=solver
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

