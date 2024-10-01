
function Badia2024Solver(op::FEOperator,params)

  # Preconditioner
  model = params[:model]
  k  = params[:fespaces][:k]
  Ω = params[:Ω]
  dΩ = Measure(Ω,2*k)
  Ωf = params[:Ωf]
  dΩf = Measure(Ωf,2*k)
  α_p = -1.0/(params[:fluid][:β] + params[:fluid][:ζ])
  α_φ = -1.0/(1.0 + params[:fluid][:ζ])
  a_Ip(p,v_p) = ∫(α_p*p*v_p)*dΩf
  a_Iφ(φ,v_φ) = ∫(α_φ*φ*v_φ)*dΩ

  U_u, U_p, U_j, U_φ = get_trial(op)
  V_u, V_p, V_j, V_φ = get_test(op)

  diag_solvers = map(s -> get_block_solver(Val(s),params), params[:solver][:block_solvers])

  uj_block = NonlinearSystemBlock(1)
  p_block  = BiformBlock(a_Ip,U_p,V_p)
  φ_block  = BiformBlock(a_Iφ,U_φ,V_φ)
  blocks = [    uj_block        LinearSystemBlock() LinearSystemBlock();
            LinearSystemBlock()     p_block         LinearSystemBlock();
            LinearSystemBlock() LinearSystemBlock()      φ_block       ]
  coeffs = [1.0 1.0 1.0;
            0.0 1.0 0.0; 
            0.0 0.0 1.0]  
  P = BlockTriangularSolver(blocks,diag_solvers,coeffs,:upper)

  # Linear Solver
  verbose = i_am_main(get_parts(model))
  nl_rtol = params[:solver][:rtol]
  l_rtol  = nl_rtol/10.0
  
  m = params[:solver][:niter]
  l_solver = FGMRESSolver(m,P;rtol=l_rtol,atol=1e-14,verbose=verbose,name="Global System - FGMRES + Badia2024")
  #SolverInterfaces.set_depth!(l_solver,2)
  l_solver.log.depth = 2

  # Nonlinear Solver
  nl_solver = GridapSolvers.NewtonSolver(l_solver,maxiter=10,atol=1e-14,rtol=nl_rtol,verbose=verbose)
  return nl_solver
end
