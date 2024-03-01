module CavityTestsSequential

using GridapPETSc
using SparseMatricesCSR, SparseArrays

using GridapMHD: cavity

# Serial, LUSolver
cavity(np=1,backend=:mpi,solver=:julia)

solver = Dict(
  :solver         => :badia2024,
  :matrix_type    => SparseMatrixCSC{Float64,Int64},
  :vector_type    => Vector{Float64},
  :petsc_options  => "-ksp_error_if_not_converged true -ksp_converged_reason",
  :block_solvers  => [:gmg,:cg_jacobi,:cg_jacobi],
)
cavity(np=1,backend=:mpi,solver=solver,ζ=10.0,ranks_per_level=[1,1])

# Serial, GMRES + block LU solvers
cavity(solver=:li2019)

# Sequential, PETSc - SNES + MUMPS
cavity(np=2,backend=:sequential,solver=:petsc)

# Sequential, GMRES + block LU solvers
cavity(np=2,backend=:sequential,solver=:li2019)

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
  :solver        => :li2019,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:from_options,:gmres_schwarz,:cg_jacobi,:cg_jacobi],
  :petsc_options => petsc_options
)
cavity(np=2,backend=:sequential,solver=solver)

end # module
