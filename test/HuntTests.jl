module HuntTests

using GridapMHD
using Gridap
using PartitionedArrays

info = GridapMHD.hunt(
  nc=(4,4),
  np=(1,1),
  backend=mpi,
  L=2.0,
  B=VectorValue(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt",
 )

info = GridapMHD.hunt(
  nc=(4,4),
  np=(2,2),
  backend=sequential,
  L=2.0,
  B=VectorValue(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt",
 )

info = GridapMHD.hunt(
  nc=(4,4),
  L=1.0,
  B=VectorValue(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt",
 )
#display(info)

end # module

