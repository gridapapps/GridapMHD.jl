#!/bin/bash
#PBS -q {{q}}
#PBS -l walltime={{walltime}}
#PBS -l ncpus={{ncpus}}
#PBS -l mem={{mem}}
#PBS -l jobfs={{jobfs}}
#PBS -N {{{name}}}
#PBS -l wd
#PBS -o {{{o}}}
#PBS -e {{{e}}} 

source ${HOME}/GridapMHD.jl/analysis/gadi/scripts/env.sh

mpiexec -n {{n}} julia --project=$GRIDAPMHD -J $GRIDAPMHD/GridapMHD.so -O3 --check-bounds=no -e\
'
# using Profile
# using FileIO
# using MPI
using GridapMHD: hunt
hunt(
  nc={{nc}},
  np={{np}},
  backend=:mpi,
  L=1.0,
  B={{B}},
  vtk={{vtk}},
  title="{{{title}}}",
  path="{{{path}}}",
  nruns={{nruns}},
  nsums={{nsums}},
  kmap_x=1,
  kmap_y=1,
  BL_adapted = true,
  solver=:petsc,
  petsc_options="{{{petsc_options}}}",
  res_assemble=false,
  jac_assemble=false,
  solve=true
 )
#  Profile.clear()
#  @profile hunt(
#   nc={{nc}},
#   np={{np}},
#   backend=:mpi,
#   L=1.0,
#   B={{B}},
#   debug={{debug}},
#   vtk={{vtk}},
#   title="{{{title}}}",
#   path="{{{path}}}",
#   nruns={{nruns}},
#   nsums={{nsums}},
#   kmap={{kmap}},
#   solver=:petsc,
#   petsc_options="{{{petsc_options}}}",
#   only_assemble=true,
#   solve=false)
#   save(joinpath("{{{path}}}","{{{name}}}_$(MPI.Comm_rank(MPI.COMM_WORLD)).jlprof"), Profile.retrieve()...)
 '