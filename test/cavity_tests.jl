module cavity_tests

using GridapMHD
using GridapPETSc
using SparseMatricesCSR

# Serial, LUSolver
GridapMHD.cavity()

# Serial, GMRES + block preconditioner
GridapMHD.cavity(solver=:block_gmres)

# Sequential, PETSc - SNES + MUMPS
GridapMHD.cavity(np=2,backend=:sequential,solver=:petsc)

# Sequential, GMRES + block LU solvers
GridapMHD.cavity(np=2,backend=:sequential,solver=:block_gmres)

# Sequential, GMRES + block preconditioners
solver_options = Dict(
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:mumps,:gmres_swartz,:amg,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
)
GridapMHD.cavity(np=2,backend=:sequential,solver=:block_gmres,solver_params=solver_options)

# Distributed MPI, LUSolver
#GridapMHD.cavity(np=2,backend=:mpi,solver=:petsc)

end # module
