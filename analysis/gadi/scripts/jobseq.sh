#!/bin/bash
#PBS -q normal 
#PBS -l walltime=04:00:00
#PBS -l ncpus=1
#PBS -l mem=190gb
#PBS -N jobseq
#PBS -l wd

source env.sh

julia --project=$GRIDAPMHD -J $GRIDAPMHD/GridapMHD.so -O3 --check-bounds=no -e\
'
using GridapMHD: hunt
hunt(
  nc=(50,50),
  L=1.0,
  B=(0.,50.,0.),
  vtk=true,
  nsums=500,
  title="hunt_50_julia",
  solver=:julia,
  BL_adapted = false
 )'
