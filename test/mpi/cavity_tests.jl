module CavityTestsMPI

using GridapPETSc
using SparseMatricesCSR

using GridapMHD: cavity

# PETSc - SNES + MUMPS
cavity(np=4,nc=(4,4,4),backend=:mpi,solver=:petsc)

# GMRES + block LU solvers
cavity(np=4,nc=(4,4,4),backend=:mpi,solver=:li2019)

# GMRES + block preconditioners
solver = Dict(
  :solver        => :li2019,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:amg,:gmres_swartz,:amg,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
)
cavity(np=4,nc=(4,4,4),backend=:mpi,solver=solver)

end # module