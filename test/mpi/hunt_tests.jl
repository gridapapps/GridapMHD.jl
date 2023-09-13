module HuntTestsMPI

using GridapPETSc
using SparseMatricesCSR

using GridapMHD: hunt

hunt(
  nc=(4,4),
  np=(2,2),
  backend=:mpi,
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt",
  solver=:petsc,
)

end