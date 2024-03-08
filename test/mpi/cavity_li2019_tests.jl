module CavityLi2019TestsMPI

using GridapPETSc
using SparseMatricesCSR

using GridapMHD: cavity

# GMRES + Li2019 with MUMPS for Dj
#cavity(np=4,nc=(4,4,4),backend=:mpi,solver=:li2019)

# GMRES + Li2019 with Overlapping Additive Schwarz preconditioner for Dj 
solver = Dict(
  :solver        => :li2019,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:petsc_gmres_schwarz,:petsc_gmres_schwarz,:petsc_cg_jacobi,:petsc_cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
)
#cavity(np=4,nc=(4,4,4),backend=:mpi,solver=solver)

end # module