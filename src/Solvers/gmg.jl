
function get_block_solver(::Val{:gmg},params)
  vars = Tuple(params[:multigrid][:variables])
  gmg_solver(Val(vars),params)
end

function gmg_solver(::Val{(1,3)},params)
  mh = params[:multigrid][:mh]

  trials = MultiFieldFESpace(map(s -> s, params[:multigrid][:trials][[1,3]]))
  tests  = MultiFieldFESpace(map(s -> s, params[:multigrid][:tests][[1,3]]))
  
  nlevs = num_levels(mh)
  k = params[:fespaces][:k]
  qdegree = map(lev -> 2*k,1:nlevs)
  _, _, α, β, γ, σf, f, B, ζ = retrieve_fluid_params(params,k)
  
  function jac_uj(x,y,dΩ)
    u, j = x
    v_u, v_j = y
    # TODO: Add thin_wall bcs... how we deal with more than one triangulation?
    # TODO: Add nonlinear terms to u-u
    r = a_mhd_u_u(u,v_u,β,dΩ) + a_mhd_u_j(j,v_u,γ,B,dΩ) + a_mhd_j_u(u,v_j,σf,B,dΩ) + a_mhd_j_j(j,v_j,dΩ)
    if abs(ζ) > eps(typeof(ζ))
      r = r + a_al_u_u(u,v_u,ζ,dΩ) + a_al_j_j(j,v_j,ζ,dΩ)
    end
    return r
  end

  return gmg_solver(mh,trials,tests,jac_uj,qdegree)
end

function gmg_solver(mh,trials,tests,biform,qdegree)
  ranks = get_level_parts(mh,1)
  smatrices = compute_gmg_matrices(mh,trials,tests,biform,qdegree)
  restrictions, prolongations = setup_transfer_operators(trials,
                                                         qdegree;
                                                         mode=:residual,
                                                         solver=IS_ConjugateGradientSolver(;reltol=1.e-6))

  smoothers = gmg_patch_smoothers(mh,tests,biform,qdegree)

  # TODO: Reuse top level matrix in gmg. In fact, all matrices should inputed (and modified) by numerical_setup
  gmg = GMGLinearSolver(mh,
                        smatrices,
                        prolongations,
                        restrictions,
                        pre_smoothers=smoothers,
                        post_smoothers=smoothers,
                        coarsest_solver=PETScLinearSolver(petsc_mumps_setup),
                        maxiter=1,
                        rtol=1.0e-8,
                        verbose=false,
                        mode=:preconditioner)
  solver = GMRESSolver(5;Pr=gmg,m_add=3,maxiter=15,rtol=1.0e-8,verbose=i_am_main(ranks))
  return solver
end

function compute_gmg_matrices(mh,trials,tests,biform,qdegree)
  nlevs = num_levels(trials)

  mats = Vector{PSparseMatrix}(undef,nlevs)
  for lev in 2:nlevs
    parts = get_level_parts(mh,lev)
    if i_am_in(parts)
      model = GridapSolvers.get_model(mh,lev)
      U = GridapSolvers.get_fe_space(trials,lev)
      V = GridapSolvers.get_fe_space(tests,lev)
      Ω = Triangulation(model)
      dΩ = Measure(Ω,qdegree[lev])
      a(u,v) = biform(u,v,dΩ)
      mats[lev] = assemble_matrix(a,U,V)
    end
  end
  return mats
end

function gmg_patch_smoothers(mh,tests,biform,qdegree)
  patch_decompositions = PatchDecomposition(mh)
  patch_spaces = PatchFESpace(tests,patch_decompositions)

  nlevs = num_levels(mh)
  smoothers = Vector{RichardsonSmoother}(undef,nlevs-1)
  for lev in 1:nlevs-1
    parts = get_level_parts(mh,lev)
    if i_am_in(parts)
      PD = patch_decompositions[lev]
      Ph = GridapSolvers.get_fe_space(patch_spaces,lev)
      Vh = GridapSolvers.get_fe_space(tests,lev)
      Ω  = Triangulation(PD)
      dΩ = Measure(Ω,qdegree[lev])
      local_solver   = LUSolver()
      patch_smoother = PatchBasedLinearSolver(biform,Ph,Vh,dΩ,local_solver)
      smoothers[lev] = RichardsonSmoother(patch_smoother,5,0.2)
    end
  end
  return smoothers
end
