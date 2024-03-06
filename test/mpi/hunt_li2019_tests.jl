module HuntLi2019TestsMPI

using GridapPETSc
using SparseMatricesCSR

using GridapMHD: hunt

# Li2019, MUMPS for Dj
solver = Dict(
  :solver        => :li2019,
  :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type   => Vector{PetscScalar},
  :block_solvers => [:mumps,:gmres_schwarz,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
)
hunt(
  nc=(6,6),
  np=(2,2),
  backend=:mpi,
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt",
  solver=solver,
)

end # module