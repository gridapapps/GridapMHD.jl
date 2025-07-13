module CavityTestsMPI

using GridapPETSc
using SparseMatricesCSR
using SparseArrays

using GridapMHD: cavity

# PETSc - SNES + MUMPS
# cavity(np=4,nc=(4,4,4),backend=:mpi,solver=:petsc)

# Cavity - H1H1 - Block solver + MUMPS

# np = (2,2,1)
# cavity(
#   nc = (4,4,4),
#   np = np,
#   backend = :mpi,
#   fluid_disc = :Qk_dPkm1,
#   current_disc = :H1,
#   solver = Dict(
#     :solver => :h1h1blocks,
#     :matrix_type    => SparseMatrixCSR{0,Float64,Int},
#     :vector_type    => Vector{Float64},
#     :block_solvers  => [:petsc_mumps,:petsc_cg_jacobi,:petsc_gmres_amg],
#     :petsc_options  => "-ksp_error_if_not_converged true -ksp_converged_reason",
#   ),
#   ζᵤ = 10.0,
#   solid = true
# )

np = (2,2,1)
cavity(
  nc = (8,8,8),
  np = np,
  backend = :mpi,
  fluid_disc = :Qk_dPkm1,
  current_disc = :H1,
  solver = Dict(
    :solver => :h1h1blocks,
    :matrix_type    => SparseMatrixCSC{Float64,Int},
    :vector_type    => Vector{Float64},
    :block_solvers  => [:gmg,:petsc_cg_jacobi,:petsc_gmres_amg],
    :petsc_options  => "-ksp_error_if_not_converged true -ksp_converged_reason",
  ),
  ranks_per_level = [np,np],
  ζᵤ = 10.0,
  solid = true
)

# np = (2,2,1)
# cavity(
#   nc = (8,8,8),
#   np = np,
#   backend = :mpi,
#   fluid_disc = :RT,
#   current_disc = :RT,
#   solver = Dict(
#     :solver => :badia2024,
#     :matrix_type    => SparseMatrixCSC{Float64,Int},
#     :vector_type    => Vector{Float64},
#     :block_solvers  => [:gmg,:petsc_cg_jacobi,:petsc_cg_jacobi],
#     :petsc_options  => "-ksp_error_if_not_converged true -ksp_converged_reason",
#   ),
#   ranks_per_level = [np,np],
#   ζᵤ = 10.0,
#   ζⱼ = 10.0,
#   solid = true
# )

end # module