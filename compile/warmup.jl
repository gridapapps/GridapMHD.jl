using GridapMHD
using Gridap
using PartitionedArrays

info = GridapMHD.hunt(
  nc=(4,4),
  L=1.0,
  B=VectorValue(0.,50.,0.),
  debug=false,
  vtk=true,
  title="warmup_hunt",
  solver=:petsc,
 )

info = GridapMHD.hunt(
  nc=(4,4),
  np=(1,1),
  backend=mpi,
  L=1.0,
  B=VectorValue(0.,50.,0.),
  debug=false,
  vtk=true,
  title="warmup_hunt_mpi",
  solver=:petsc,
 )
