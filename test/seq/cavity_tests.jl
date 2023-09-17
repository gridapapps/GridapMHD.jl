module CavityTestsSequential

using GridapPETSc
using SparseMatricesCSR

using GridapMHD: cavity

# Serial, LUSolver
cavity()

# Serial, GMRES + block LU solvers
cavity(solver=:block_gmres_li2019)

# Sequential, PETSc - SNES + MUMPS
cavity(np=2,backend=:sequential,solver=:petsc)

# Sequential, GMRES + block LU solvers
cavity(np=2,backend=:sequential,solver=:block_gmres_li2019)

# Sequential, GMRES + block preconditioners
solver = Dict(
  :solver        => :block_gmres_li2019,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:amg,:gmres_swartz,:amg,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
)
cavity(np=2,backend=:sequential,solver=solver)

# Distributed MPI, LUSolver
# cavity(np=2,backend=:mpi,solver=:petsc)

end # module
