
function Badia2024Solver(op::FEOperator,params)

  # Preconditioner
  model = params[:model]
  k  = params[:fespaces][:k]
  Ωf = _interior(model,params[:fluid][:domain])
  dΩ = Measure(Ωf,2*k)
  a_Ip(p,v_p) = ∫(p*v_p)*dΩ
  a_Iφ(φ,v_φ) = ∫(-φ*v_φ)*dΩ

  U_u, U_p, U_j, U_φ = get_trial(op)
  V_u, V_p, V_j, V_φ = get_test(op)

  NB = 3
  diag_solvers = map(s -> get_block_solver(Val(s)),params[:solver][:block_solvers])
  diag_blocks  = [NonlinearSystemBlock(),BiformBlock(a_Ip,U_p,V_p),BiformBlock(a_Iφ,U_φ,V_φ)]
  blocks = map(CartesianIndices((NB,NB))) do I
    (I[1] == I[2]) ? diag_blocks[I[1]] : LinearSystemBlock()
  end
  coeffs = [1.0 1.0 1.0;
            0.0 1.0 0.0; 
            0.0 0.0 1.0]  
  P = BlockTriangularSolver(blocks,diag_solvers,coeffs,:upper)

  # Linear Solver
  verbose = i_am_main(get_parts(model))
  nl_rtol = params[:solver][:rtol]
  l_rtol  = nl_rtol/10.0
  
  m = params[:solver][:niter]
  l_solver = FGMRESSolver(m,P;rtol=l_rtol,atol=1e-14,verbose=verbose)

  # Nonlinear Solver
  nl_solver = GridapSolvers.NewtonSolver(l_solver,maxiter=10,atol=1e-14,rtol=nl_rtol,verbose=verbose)
  return nl_solver
end

function badia2024_gmg(params)
  mh = params[:multigrid][:mh]
  trials = MultiFieldFESpace(params[:multigrid][:trials][[1,3]])
  tests  = MultiFieldFESpace(params[:multigrid][:tests][[1,3]])
  
  nlevs = num_levels(mh)
  k = params[:fespaces][:k]
  qdegree = map(lev -> 2*k,1:nlevs)
  
  function jac_uj(x,y,dΩ)
    u,j = x
    v_u,v_j = y
    # TODO: Add thin_wall bcs... how we deal with more than one triangulation?
    r = a_mhd_u_u(u,v_u,β,dΩ) + a_mhd_u_j(j,v_u,γ,B,dΩ) + a_mhd_j_u(u,v_j,σ,B,dΩ) + a_mhd_j_j(j,v_j,dΩ)
    return r
  end

  return gmg_solver(mh,trials,tests,jac_uj,qdegree)
end
