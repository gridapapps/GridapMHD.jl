#!/bin/bash
#PBS -q normal 
#PBS -l walltime=00:15:00
#PBS -l ncpus=4
#PBS -l mem=190gb
#PBS -N jobmpi
#PBS -l wd

source env.sh

mpiexec -n 4 julia --project=$GRIDAPMHD -J $GRIDAPMHD/GridapMHD.so -O3 --check-bounds=no -e\
'
using GridapMHD: hunt
hunt(
  nc=(48,48),
  np=(2,2),
  backend=:mpi,
  L=1.0,
  B=(0.,50.,0.),
  nruns=2,
  nsums=200,
  vtk=true,
  title="hunt_48_petsc",
  solver=:petsc,
  BL_adapted = true,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_4 3  -mat_mumps_icntl_28 2  -mat_mumps_icntl_29 2"
 )'

