#!/bin/bash
#SBATCH -N 1 
#SBATCH --ntasks-per-node=12
#SBATCH -t 90:00:00
#SBATCH --partition=volta

#SBATCH -o outputMaPLEHa100Re1-block
#SBATCH -e errorMaPLEHa100Re1-block
#SBATCH --job-name=MaPLEHa100Re1-block
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

#Load the enviroment
source env.sh

#Generate a file for passing some SLURM parameters to julia 

PASS_FILE="pass_params.jl"

echo NPROCS=$SLURM_NPROCS > $PASS_FILE
echo JOB_NAME=\"$SLURM_JOB_NAME\" >> $PASS_FILE

#Operational parameters defining the problem 

echo Ha= 100.0   >>     $PASS_FILE
echo Re= 1.0	 >>	$PASS_FILE
echo cw = 0.0684 >>	$PASS_FILE
echo τ= 1e6 	 >>	$PASS_FILE

#Parameters for the mesh construction

echo R = 1.0    >>      $PASS_FILE	#Radious of the pipe (used for normalization)
echo p = 0.5    >>      $PASS_FILE	#Edge of the BL region
echo q = 0.5    >>      $PASS_FILE	#Core radious
echo L = 4.0    >>      $PASS_FILE	#Length of the pipe

echo N_r = 12 >>        $PASS_FILE	#Radial nodes in the BL region
echo n =4 >>            $PASS_FILE	#Radial nodes in the Hartmann BL (1/Ha)
echo N_a= 40 >>         $PASS_FILE	#Azimutal nodes 
echo N_L= 12 >>         $PASS_FILE	#Nodes in the axial direction
echo n_c = 0.2 >>       $PASS_FILE	#Maximum size of the cells edges in the core region

#Generate the mesh

#For a hybrid (block-unstructured) mesh
#export mesh_path=/ws/blankets/GridapMHD.jl/meshes/tube_hybrid

#For a block mesh (stil tetra)

export mesh_path=/ws/blankets/GridapMHD.jl/meshes/tube_block

source $mesh_path/MeshGenerator.sh

#GridapMHD computation
mpiexec -n ${SLURM_NPROCS} --mca btl_openib_if_include mlx5_0 julia --project=$GRIDAPMHD -O3 --check-bounds=no -e\
'
include("pass_params.jl")
using GridapMHD: tube
tube(;
  mesh="MaPLE", 
  np=NPROCS,
  backend=:mpi,
  Ha = Ha,
  N = Ha^2/Re,
  inlet = :parabolic,
  cw = cw,
  τ = τ,
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
