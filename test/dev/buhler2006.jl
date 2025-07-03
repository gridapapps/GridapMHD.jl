using GridapMHD: expansion;
using SparseMatricesCSR;
using GridapPETSc;

solver = Dict(
  :solver        => :badia2024,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:petsc_mumps,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason",
)

expansion(
  backend = :mpi,
  np = 1,
  mesh = "/scratch/bt62/jm3247/mhd-validation/experiments/buhler2006/meshes/expansion_Ha_1000_N_1000000_np_16.msh",
  Ha = 1000,
  N = 1000000,
  ζ = 10.0,
  cw = 0.028,
  Z = 4.0,
  b = 1.0,
  inlet = :parabolic,
  solid_coupling = :solid,
  solver = solver,
  savelines = true,
  title = "expansion_Ha_1000_N_1000000_np_16",
  path = "/scratch/bt62/jm3247/mhd-validation/experiments/buhler2006/data"
)


