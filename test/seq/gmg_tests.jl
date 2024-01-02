
using GridapMHD: expansion
using GridapPETSc, SparseMatricesCSR

mesh = Dict(
  :num_refs_coarse => 0,
  :ranks_per_level => [1,1],
)
solver = Dict(
  :solver        => :badia2024,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:gmg,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
)
expansion(np=1,backend=:mpi,mesh=mesh,solver=solver,ζ=0.0,cw=0.0)
