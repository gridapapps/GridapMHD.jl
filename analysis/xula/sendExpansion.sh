#!/bin/bash
#SBATCH -N 4 
#SBATCH --ntasks-per-node=48 
#SBATCH -t 56:00:00
#SBATCH -p xula3

#SBATCH -o outputExp_r4Ha1000Re50_N4n4
#SBATCH -e errorExp_r4Ha1000Re50_N4n4
###SBATCH --mail-user=fernando.roca@ciemat.es
#SBATCH --job-name=Exp_r4Ha1000Re50_N4n4
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
echo Ha=1000.0 >>	$PASS_FILE 	#Hartmann number defined on the outlet channel
echo Re=50.0 >>	 	$PASS_FILE 	#Reynolds number defined on the outlet channesl
echo r=4.0 >>	 	$PASS_FILE	#Expansion ratio
echo betta = 0.2 >> 	$PASS_FILE	#Outlet channel aspect ratio
echo L_in=1.0 >>	$PASS_FILE	#Inlet channel normalized lenght
echo L_out=2.0 >> 	$PASS_FILE      #Outlet channel normalized lenght

#Definition of the mesh parameters
echo N_Ha=28 >>	 	$PASS_FILE	#Number of cells in the inlet channel along the B-direction
echo n_Ha=4 >> 		$PASS_FILE	#Number of cells in the Hartmann BL
echo N_side=20 >> 	$PASS_FILE	#Number of cells in along the direction perpendicular to B
echo n_side=5 >> 	$PASS_FILE	#Number of cells in side BL
echo n_inlet=20 >>	$PASS_FILE	#Number of axial cells in the inlet channel
echo n_outlet=40 >> 	$PASS_FILE	#Number of outlet cells in the outlet channel
echo R_exp=1.075 >>	$PASS_FILE	#Ratio of the geometric clustering towards the expansion point (x=0)

#Mesh computation using gmsh

source $GRIDAPMHD/meshes/expansion/MeshGenerator.sh

#Parallel run using GridapMHD

srun --mpi=pmix -n ${SLURM_NPROCS} julia --project=$GRIDAPMHD -O3 --check-bounds=no -e\
'
# make a local temp directory for the cache
d = strip(String(read(`mktemp -d`)))
mkdir(joinpath(d, "compiled"))
# copy the cache to the temp directory
# set the temp directory as the first depot so any (re)compilation will happen there and not interfere with other jobs
run(`rsync -au $(DEPOT_PATH[1])/compiled $d/compiled`)
pushfirst!(DEPOT_PATH, d)

include("input_params.jl")
using GridapMHD: expansion
expansion(;
  mesh="computed", 
  np=NPROCS,
  backend=:mpi,
  Ha = Ha,
  N = Ha^2/Re, 
  cw = 0.0,
  Z = r,
  b = betta,
#  τ = 1e6,
  inlet=:shercliff,
  debug=false,
  vtk=true,
  title=JOB_NAME,
  solver=:petsc,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_7 0 -mat_mumps_icntl_28 1 -mat_mumps_icntl_29 2 -mat_mumps_icntl_4 3 -mat_mumps_cntl_1 0.001"
 )'


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

