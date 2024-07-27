module ChannelTestsSequential

using GridapMHD: channel
using GridapPETSc, SparseMatricesCSR


channel(nc=(4,4,4),vtk=true)

channel(backend=:mpi,np=(1,1,1),solver=:petsc,
  sizes=(4,2,2),nc=(4,4,4),vtk=true)

channel(backend=:mpi,np=(1,1,1),solver=:petsc,
  sizes=(8,2,2),nc=(4,8,8),B0=20,vtk=true,inlet=:shercliff)

channel(sizes=(4,1,1),nc=(16,4,4),vtk=true)

channel(sizes=(4,1,1),nc=(16,8,8,),debug=false,vtk=true,
  bl_orders=(1,2,2),nonuniform_B=true,B0=20,inlet=:shercliff)


end
