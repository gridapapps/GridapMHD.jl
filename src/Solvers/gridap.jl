
get_block_solver(::Val{:julia},params)     = LUSolver()
get_block_solver(::Val{:cg_jacobi},params) = CG_jacobi_solver(params)

function CG_jacobi_solver(params)
  return CGSolver(JacobiLinearSolver();maxiter=10)
end
