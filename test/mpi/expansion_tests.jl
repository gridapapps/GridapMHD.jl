module ExpansionTestsMPI

using GridapPETSc
using SparseMatricesCSR
using GridapMHD: expansion

#expansion(np=4,backend=:mpi,savelines=true)

mesh = Dict{Symbol,Any}(
  :mesher    => :p4est_SG,
  :base_mesh => "meshes/Expansion_710.msh",
  :num_refs  => 1,
)
solver = Dict{Symbol,Any}(
  :solver        => :badia2024,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:petsc_mumps,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason",
)
expansion(
  np=4,
  backend=:mpi,
  mesh=mesh,
  solver=solver,
  ζ  = 10.0,
  savelines=true
)

end