#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 3-00:00:00
#SBATCH -p xula3

#SBATCH -o outputExpWall_test
#SBATCH -e errorExpWall_test
#SBATCH --job-name=ExpWall_test
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

#Generate a file for passing some SLURM parameters to julia

PASS_FILE="input_params.jl"

echo NPROCS=$SLURM_NPROCS > $PASS_FILE
echo JOB_NAME=\"$SLURM_JOB_NAME\" >> $PASS_FILE

#Definition of the problem input paramenters
echo Ha=1000.0 >>		$PASS_FILE 	#Hartmann number defined on the outlet channel
echo Re=100.0 >> 		$PASS_FILE 	#Reynolds number defined on the outlet channesl
echo r=4.0 >>	 		$PASS_FILE	#Expansion ratio
echo betta = 1.0 >> 		$PASS_FILE	#Outlet channel aspect ratio
echo cw = 0.028 >>		$PASS_FILE	#Wall conductance ratio
echo L_in=8.0 >>		$PASS_FILE	#Inlet channel normalized lenght
echo L_out=8.0 >> 		$PASS_FILE      #Outlet channel normalized lenght
echo t = 0.028 >>		$PASS_FILE	#Thickness of the wall

echo N_Ha = 6 >>      $PASS_FILE      #Number of cells in the inlet channel along the B-direction
echo N_s = 6 >>       $PASS_FILE      #Number of cells in the inlet channel along the direction perpendicular to B
echo exp_Ha = 1.0 >>  $PASS_FILE      #Exponent of the Smolentsev bunching law in the B-direction
echo exp_s = 0.7 >>   $PASS_FILE      #Exponent of the Smolentsev bunching law in the direction perpendicular to B
echo dx = 0.1 >>      $PASS_FILE      #Distante of the first node in the sudden expansion
echo N_1 = 10 >>      $PASS_FILE      #Number of cells in the outlet channel along the flow direction
echo N_2 = 10 >>      $PASS_FILE      #Number of cells in the inlet channel along the flow direction
echo N_w = 1 >>	      $PASS_FILE      #Number of cells along the wall direction

#Mesh computation using gmsh

source $GRIDAPMHD/meshes/expansion/MeshGenerator-wall.sh
: <<'END_COMENT'
#Parallel run using GridapMHD

julia --project=$GRIDAPMHD -O3 --check-bounds=no -e\
'
# make a local temp directory for the cache
d = strip(String(read(`mktemp -d`)))
mkdir(joinpath(d, "compiled"))
# copy the cache to the temp directory
# set the temp directory as the first depot so any (re)compilation will happen there and not interfere with other jobs
run(`rsync -au $(DEPOT_PATH[1])/compiled $d/compiled`)
pushfirst!(DEPOT_PATH, d)

include("input_params.jl")
using GridapPETSc
using SparseMatricesCSR
using GridapMHD: expansion

#Monolithic MUMPS
solver = Dict(
    :solver => :petsc,
    :matrix_type    => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type    => Vector{PetscScalar},
    :petsc_options  => "-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_28 1 -mat_mumps_icntl_29 2 -mat_mumps_icntl_4 3 -mat_mumps_cntl_1 0.001",
    :niter          => 100,
    :rtol           => 1e-5,
  )
"""
expansion(;
  mesh="compile",
  Ha = 10.0,
  N = 100.0,
  cw = 0.0,
  debug=false,
  Z = 4.0,
  b = 0.2,
  inlet=:shercliff,
  vtk=false,
  title="warmup_gmsh_petsc",
  solver=solver,
)
"""
expansion(;
  mesh="wall_computed", 
  np=NPROCS,
  backend=:mpi,
  path =".",
  Ha = Ha,
  N = Ha^2/Re, 
  cw = cw,
  Z = r,
  b = betta,
  t = t,
  inlet=:parabolic,
  debug=true,
  vtk=true,
  title=JOB_NAME,
  solver=solver
 )'

END_COMENT

duration=$SECONDS
rm -f hosts.$SLURM_JOB_ID
rm -f input_params.jl
rm -f mesh_params

STATUS=$?
echo "================================================================"
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "================================================================"
echo ""
echo "STATUS = $STATUS"
echo ""

