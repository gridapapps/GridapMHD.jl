module HuntLi2019TestsMPI

using GridapPETSc
using SparseMatricesCSR

using GridapMHD: hunt

function main(parts)

  # Default li2019
  hunt(
    nc=(4,4),
    np=parts,
    backend=:mpi,
    L=1.0,
    B=(0.0, 20.0, 0.0),
    nsums=100,
    debug=false,
    vtk=true,
    title="hunt",
    solver=:li2019,
  )

  # Li2019, MUMPS for Dj
  solver = Dict(
    :solver => :li2019,
    :matrix_type => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type => Vector{PetscScalar},
    :block_solvers => [:petsc_mumps, :petsc_gmres_schwarz, :petsc_cg_jacobi, :petsc_cg_jacobi],
    :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
  )
  hunt(
    nc=(6, 6),
    np=parts,
    backend=:mpi,
    L=1.0,
    B=(0.0, 50.0, 0.0),
    debug=false,
    vtk=true,
    title="hunt",
    solver=solver,
  )

  # Li2019, GMRES + Overlapping Additive Schwarz preconditioner for Dj
  solver = Dict(
    :solver => :li2019,
    :matrix_type => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type => Vector{PetscScalar},
    :block_solvers => [:petsc_gmres_schwarz, :petsc_gmres_schwarz, :petsc_cg_jacobi, :petsc_cg_jacobi],
    :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
  )
  hunt(
    nc=(4, 4),
    np=parts,
    backend=:sequential,
    L=1.0,
    B=(0.0, 50.0, 0.0),
    debug=false,
    vtk=true,
    title="hunt",
    solver=solver,
  )

  # Li2019, GMRES + Overlapping Additive Schwarz preconditioner for Dj (tunned by options)
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
    :solver => :li2019,
    :matrix_type => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type => Vector{PetscScalar},
    :block_solvers => [:petsc_from_options, :petsc_gmres_schwarz, :petsc_cg_jacobi, :petsc_cg_jacobi],
    :petsc_options => petsc_options
  )
  hunt(
    nc=(4, 4),
    np=parts,
    backend=:sequential,
    L=1.0,
    B=(0.0, 50.0, 0.0),
    debug=false,
    vtk=true,
    title="hunt",
    solver=solver,
  )

end # main

end # module