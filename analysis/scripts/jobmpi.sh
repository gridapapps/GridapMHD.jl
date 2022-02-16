#!/bin/bash
#PBS -q normal 
#PBS -l walltime=00:10:00
#PBS -l ncpus=4
#PBS -l mem=160gb
#PBS -N jobmpi
#PBS -l wd

source env.sh

#mpiexec -n 4 julia --project=$GRIDAPMHD -J $GRIDAPMHD/GridapMHD.so -O3 --check-bounds=no -e\
mpiexec -n 4 julia --project=$GRIDAPMHD  -O0 --check-bounds=yes -e\
'
using Gridap
using GridapMHD
using PartitionedArrays
GridapMHD.hunt(
  nc=(4,4),
  np=(2,2),
  backend=mpi,
  L=1.0,
  B=VectorValue(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt",
  solver=:petsc,
 )'

