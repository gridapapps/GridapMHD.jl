#!/bin/bash
#PBS -q normal 
#PBS -l walltime=00:10:00
#PBS -l ncpus=1
#PBS -l mem=160gb
#PBS -N jobmpi
#PBS -l wd

source env.sh

julia --project=$GRIDAPMHD -J $GRIDAPMHD/GridapMHD.so -O3 --check-bounds=no -e\
'
using GridapMHD: hunt
hunt(
  nc=(4,4),
  np=(1,1),
  backend=:mpi,
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt",
  solver=:petsc,
 )'

