mpiexec -n 4 julia --project=. -e'
  using GridapMHD: hunt;
  using SparseMatricesCSR;
  using GridapPETSc;

  solver = Dict(
    :solver        => :badia2024,
    :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type   => Vector{PetscScalar},
    :block_solvers => [:petsc_mumps,:cg_jacobi,:cg_jacobi],
    :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason",
    :niter         => 10,
  )

  hunt(backend=:mpi,np=(2,2),solver=solver,title="hunt",ζ=10.0)
  '
