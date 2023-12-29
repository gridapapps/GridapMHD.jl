module ExpansionTestsSequential

using GridapMHD: expansion
using GridapPETSc, SparseMatricesCSR

expansion()
expansion(np=2,backend=:sequential)
expansion(np=1,backend=:mpi)

# Li2019, with MUMPS solver for Dj
solver = Dict(
  :solver        => :li2019,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:mumps,:gmres_schwarz,:gmres_amg,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
)
expansion(np=2,backend=:sequential,solver=solver)

# Li2019, with Schwarz solver for Dj
petsc_options = """
  -ksp_type gmres
  -ksp_rtol 1.0e-5
  -ksp_atol 1.0e-14
  -ksp_converged_reason
  -pc_type asm
  -pc_asm_overlap 10
  -pc_asm_type restrict
  -pc_asm_blocks 32
  -sub_ksp_type preonly
  -sub_pc_type lu
"""
solver = Dict(
  :solver        => :li2019,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:from_options,:gmres_schwarz,:cg_jacobi,:cg_jacobi],
  :petsc_options => petsc_options
)
expansion(np=2,backend=:sequential,solver=solver)

end # module