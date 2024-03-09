module CavityBadia2024TestsSequential

using SparseArrays
using GridapMHD: cavity

solver = Dict(
  :solver         => :badia2024,
  :matrix_type    => SparseMatrixCSC{Float64,Int64},
  :vector_type    => Vector{Float64},
  :petsc_options  => "-ksp_error_if_not_converged true -ksp_converged_reason",
  :block_solvers  => [:julia,:cg_jacobi,:cg_jacobi],
)
cavity(np=1,backend=:mpi,solver=solver,ζ=10.0)

solver = Dict(
  :solver         => :badia2024,
  :matrix_type    => SparseMatrixCSC{Float64,Int64},
  :vector_type    => Vector{Float64},
  :petsc_options  => "-ksp_error_if_not_converged true -ksp_converged_reason",
  :block_solvers  => [:gmg,:cg_jacobi,:cg_jacobi],
)
cavity(np=1,backend=:mpi,solver=solver,ζ=10.0,ranks_per_level=[1,1])

end