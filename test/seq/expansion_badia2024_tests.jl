module ExpansionBadia2024TestsSequential

using GridapMHD: expansion
using SparseArrays
using SparseMatricesCSR
using GridapPETSc

Re = 1.0
Ha = 10.0
N = Ha^2/Re

# solver = Dict(
#   :solver        => :badia2024,
#   :matrix_type   => SparseMatrixCSC{Float64,Int64},
#   :vector_type   => Vector{Float64},
#   :block_solvers => [:julia,:cg_jacobi,:cg_jacobi],
#   :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
# )

# solver = Dict(
#   :solver        => :badia2024,
#   :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
#   :vector_type   => Vector{PetscScalar},
#   :block_solvers => [:petsc_from_options,:cg_jacobi,:cg_jacobi],
#   :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps",
# )

solver = Dict(
  :solver        => :badia2024,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:petsc_mumps,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason",
)

expansion(np=4,backend=:mpi,mesh="710",solver=solver,order=2,ζ=10.0,N=N,Ha=Ha,title="Expansion")

end # module