module HuntGMGTestsMPI

using GridapPETSc
using SparseMatricesCSR, SparseArrays
using GridapMHD: hunt

function main(parts)

np = parts
hunt(
  nc=(4,4),
  np=np,
  backend=:mpi,
  fluid_disc = :Qk_dPkm1,
  current_disc = :H1,  
  order = 2,
  solver= Dict(
      :solver => :h1h1blocks,
      :matrix_type    => SparseMatrixCSC{Float64,Int},
      :vector_type    => Vector{Float64},
      :block_solvers  => [:gmg,:petsc_cg_jacobi,:petsc_gmres_amg],
      :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
      ),
  title="hunt-H1H1-gmg",
  ranks_per_level = [np,np],
  ζᵤ = 10.0
)

end

end