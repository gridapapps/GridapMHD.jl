module HuntTestsSequential

using SparseMatricesCSR, SparseArrays
using GridapPETSc
using GridapMHD: hunt

hunt(
  nc=(4,4),
  L=1.0,
  B=(0.,50.,0.),
  debug=true,
  vtk=true,
  title="hunt-debug",
  solver=:julia
)

hunt(
  nc=(10,10),
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt-H1HDiv", # default
  solver=:julia
)

hunt(
  nc=(10,10),
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt-H1H1",
  solver=:julia,
  fluid_disc = :Qk_dPkm1,
  current_disc = :H1
)

hunt(
  nc=(10,10),
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt-HDivHDiv",
  solver=:julia,
  fluid_disc = :RT
)

hunt(
  nc=(4,4),
  np=(1,1),
  backend=:mpi,
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt-mpi-backend",
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
  title="hunt-partitioned",
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

hunt(
  nc=(8,8),
  fluid_disc = :Qk_dPkm1,
  current_disc = :H1,
  solver = Dict(
    :solver => :h1h1blocks,
    :matrix_type    => SparseMatrixCSC{Float64,Int},
    :vector_type    => Vector{Float64},
    :block_solvers  => [:julia,:julia,:julia],
    :petsc_options  => "-ksp_error_if_not_converged true -ksp_converged_reason",
  ),
  ζᵤ = 10.0,
  ζⱼ = 10.0,
)

hunt(
  nc=(8,8),
  fluid_disc = :RT,
  current_disc = :H1,
  solver = Dict(
    :solver => :h1h1blocks,
    :matrix_type    => SparseMatrixCSC{Float64,Int},
    :vector_type    => Vector{Float64},
    :block_solvers  => [:julia,:julia,:julia],
    :petsc_options  => "-ksp_error_if_not_converged true -ksp_converged_reason",
  ),
  ζᵤ = 10.0,
  ζⱼ = 10.0,
)

end # module
