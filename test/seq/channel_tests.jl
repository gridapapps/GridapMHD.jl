module ChannelTestsSequential

using GridapMHD: channel
using GridapPETSc, SparseMatricesCSR, SparseArrays

solver = Dict(
  :solver        => :badia2024,
  :matrix_type   => SparseMatrixCSC{Float64,Int},
  :vector_type   => Vector{Float64},
  :block_solvers => [:julia,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason",
  :rtol => 1.e-10
)

channel(
  nc=(4,2,2),
  order = 3,
  order_j = 4,
  L = 2.0,
  w = 1.0,
  ζ = 50.0,
  B = (0.0,50.0,0.0),
  bl_orders=(1,2,2),
  rt_scaling = true,
  convection = :none,
  simplexify = true,
  fluid_disc = :Pk_dPkm1,
  current_disc = :BDM,
  solver = solver,
  vtk=true
)

solver = Dict(
  :solver        => :badia2024,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:petsc_mumps,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason",
  :rtol => 1.e-10
)

channel(
  backend=:mpi,
  np=(1,1,1),
  solver=solver,
  sizes=(6,3,3),
  nc=(4,4,4),
  order=3,
  order_j=3,
  simplexify=true,
  bl_orders=(1,2,2),
  convection=:none,
  inlet=:shercliff,
  rt_scaling = true,
  B = (0.0,50.0,0.0),
  ζ = 100.0,
)

#channel(backend=:mpi,np=(1,1,1),solver=:petsc,
#  sizes=(8,2,2),nc=(4,8,8),B0=20,vtk=true,inlet=:shercliff)
#
#channel(sizes=(4,1,1),nc=(16,4,4),vtk=true)
#
#channel(sizes=(4,1,1),nc=(16,8,8,),debug=false,vtk=true,
#  bl_orders=(1,2,2),nonuniform_B=true,B0=20,inlet=:shercliff)


end
