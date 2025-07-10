module CavityTestsSequential

using GridapPETSc
using SparseMatricesCSR, SparseArrays
using GridapMHD: cavity

# Serial, LUSolver
cavity(
  fluid_disc = :Qk_dPkm1,
  current_disc = :RT,
)

cavity(
  fluid_disc = :Qk_dPkm1,
  current_disc = :H1,
)

cavity(
  fluid_disc = :RT,
  current_disc = :RT,
)

cavity(
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

cavity(
  nc = (2,2,2),
  fluid_disc = :Qk_dPkm1,
  current_disc = :H1,
  solver = Dict(
    :solver => :h1h1blocks,
    :matrix_type    => SparseMatrixCSC{Float64,Int},
    :vector_type    => Vector{Float64},
    :block_solvers  => [:gmg,:julia,:julia],
    :petsc_options  => "-ksp_error_if_not_converged true -ksp_converged_reason",
  ),
  ranks_per_level = [(1,1,1),(1,1,1)],
  ζᵤ = 10.0,
  ζⱼ = 10.0,   
)

end # module
