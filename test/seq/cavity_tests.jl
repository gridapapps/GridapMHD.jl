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
  :solver        => :block_gmres_li2019,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:from_options,:gmres_schwarz,:gmres_amg,:cg_jacobi,:cg_jacobi],
  :petsc_options => petsc_options
)
cavity(np=2,backend=:sequential,solver=solver)

# Distributed MPI, LUSolver
# cavity(np=2,backend=:mpi,solver=:petsc)

end # module
