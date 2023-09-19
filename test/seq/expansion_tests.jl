module ExpansionTestsSequential

using GridapMHD: expansion
using GridapPETSc, SparseMatricesCSR

expansion()
expansion(np=2,backend=:sequential)
expansion(np=1,backend=:mpi)

solver = Dict(
  :solver        => :block_gmres_li2019,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:mumps,:gmres_schwarz,:gmres_schwarz,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
)
expansion(np=2,backend=:sequential,solver=solver)

end # module