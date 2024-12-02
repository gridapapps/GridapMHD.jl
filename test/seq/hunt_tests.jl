module HuntTestsSequential

using GridapMHD: hunt
using GridapPETSc, SparseMatricesCSR

hunt(
  nc=(4,4),
  L=1.0,
  B=(0.,50.,0.),
  debug=true,
  vtk=true,
  title="hunt",
  solver=:julia,
)

hunt(
  nc=(10,10),
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt",
  solver=:julia,
  fluid_disc = :RT,
)

hunt(
  nc=(10,10),
  L=1.0,
  B=(0.,50.,0.),
  order=3,
  title="hunt",
  solver=:julia,
  simplexify=true,
)

hunt(
  nc=(4,4),
  np=(1,1),
  backend=:mpi,
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt",
  solver=:julia,
)

hunt(
  nc=(4,4),
  np=(2,2),
  backend=:sequential,
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt",
  solver=:julia,
)

hunt(
  nc=(12,12),
  backend=:sequential,
  np=(2,2),
  L=1.0,
  B=(0.,50.,0.),
  tw=0.2,
  debug=false,
  vtk=true,
  title="hunt-solid",
  solver=:julia,
  BL_adapted=false,
  kmap_x = 3,
  kmap_y = 3
)

end # module
