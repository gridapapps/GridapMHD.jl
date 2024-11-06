
function get_block_solver(::Val{:gmg},params)
  vars = Tuple(params[:multigrid][:variables])
  gmg_solver(Val(vars),params)
end

function gmg_solver(::Val{(:u,:j)},params)
  mh = params[:multigrid][:mh]

  trials_u = params[:multigrid][:trials][:u]
  tests_u  = params[:multigrid][:tests][:u]
  trials_j = params[:multigrid][:trials][:j]
  tests_j  = params[:multigrid][:tests][:j]
  trials = MultiFieldFESpace([trials_u, trials_j])
  tests  = MultiFieldFESpace([tests_u, tests_j])
  
  nlevs = num_levels(mh)
  k = params[:fespaces][:k]
  qdegree = map(lev -> 2*k+1,1:nlevs)
  is_nonlinear = has_convection(params)
  
  reffe_p = params[:fespaces][:reffe_p]
  Πp = MultilevelTools.LocalProjectionMap(divergence,reffe_p,2*k)
  
  function jacobian_uj(model)
    Ωf = Triangulation(model)
    dΩf = Measure(Ωf,2*k)
    _, _, α, β, γ, σf, f, B, ζ = retrieve_fluid_params(model,params,k)
    _, params_tw, _, _, params_Λ = retrieve_bcs_params(model,params,k)
    function a(x,dx,y)
      u, j = x
      du, dj = dx
      v_u, v_j = y
      r = a_mhd_u_u(du,v_u,β,dΩf) + a_mhd_u_j(dj,v_u,γ,B,dΩf) + a_mhd_j_u(du,v_j,σf,B,dΩf) + a_mhd_j_j(dj,v_j,dΩf)
      if is_nonlinear
        r = r + dc_mhd_u_u(u,du,v_u,α,dΩf)
      end
      for p in params_tw
        r = r + a_thin_wall_j_j(dj,v_j,p...)
      end
      if abs(ζ) > eps(typeof(ζ))
        r = r + a_al_u_u(du,v_u,ζ,Πp,dΩf) + a_al_j_j(dj,v_j,ζ,dΩf)
      end
      return r
    end
    return a
  end

  function biform_uj(model)
    Ωf = Triangulation(model)
    dΩf = Measure(Ωf,2*k)
    _, _, α, β, γ, σf, f, B, ζ = retrieve_fluid_params(model,params,k)
    _, params_tw, _, _, params_Λ = retrieve_bcs_params(model,params,k)
    function a(dx,y)
      du, dj = dx
      v_u, v_j = y
      r = a_mhd_u_u(du,v_u,β,dΩf) + a_mhd_u_j(dj,v_u,γ,B,dΩf) + a_mhd_j_u(du,v_j,σf,B,dΩf) + a_mhd_j_j(dj,v_j,dΩf)
      for p in params_tw
        r = r + a_thin_wall_j_j(dj,v_j,p...)
      end
      if abs(ζ) > eps(typeof(ζ))
        r = r + a_al_u_u(du,v_u,ζ,Πp,dΩf) + a_al_j_j(dj,v_j,ζ,dΩf)
      end
      return r
    end
    return a
  end

  function jacobian_u(model)
    _, _, α, β, γ, σf, f, B, ζ = retrieve_fluid_params(model,params,k)
    Ωf = Triangulation(model)
    dΩf = Measure(Ωf,2*k)
    function a(u,du,v_u)
      r = a_mhd_u_u(du,v_u,β,dΩf) 
      if is_nonlinear
        dc_mhd_u_u(u,du,v_u,α,dΩf)
      end
      if abs(ζ) > eps(typeof(ζ))
        r = r + a_al_u_u(du,v_u,ζ,Πp,dΩf)
      end
      return r
    end
    return a
  end

  function rhs(model)
    _, _, α, β, γ, σf, f, B, ζ = retrieve_fluid_params(model,params,k)
    Ωf = Triangulation(model)
    dΩf = Measure(Ωf,2*k)
    l((u,j),(v_u,v_j)) = a_al_u_u(u,v_u,ζ,Πp,dΩf) + a_al_j_j(j,v_j,ζ,dΩf)
    return l
  end
  function rhs_u(model)
    _, _, α, β, γ, σf, f, B, ζ = retrieve_fluid_params(model,params,k)
    Ωf = Triangulation(model)
    dΩf = Measure(Ωf,2*k)
    l(u,v_u) = a_al_u_u(u,v_u,ζ,Πp,dΩf)
    return l
  end

  weakforms = map(mhl -> jacobian_uj(GridapSolvers.get_model(mhl)),mh)
  smoothers = gmg_patch_smoothers(mh,trials,jacobian_uj;is_nonlinear,w=-0.2)

  mf_prolongation = false
  if mf_prolongation
    prolongations_u = gmg_patch_prolongations(tests_u,jacobian_u,rhs_u;is_nonlinear)
    prolongations_j = setup_prolongation_operators(tests_j,qdegree)
    prolongations = MultiFieldTransferOperator(tests,[prolongations_u,prolongations_j];op_type=:prolongation)
    restrictions_u = gmg_patch_restrictions(tests_u,prolongations_u,rhs_u,qdegree;solver=LUSolver())
    restrictions_j = setup_restriction_operators(tests_j,qdegree)
    restrictions = MultiFieldTransferOperator(tests,[restrictions_u,restrictions_j];op_type=:restriction)
  else
    prolongations = gmg_patch_prolongations(tests,biform_uj,biform_uj;is_nonlinear=false)
    restrictions = gmg_patch_restrictions(tests,prolongations,biform_uj,qdegree;solver=LUSolver())
  end

  return gmg_solver(mh,trials,tests,weakforms,restrictions,prolongations,smoothers)
end

function gmg_solver(mh,trials,tests,weakforms,restrictions,prolongations,smoothers)
  ranks = get_level_parts(mh,1)

  cranks = get_level_parts(mh,num_levels(mh))
  coarsest_solver = (length(cranks) == 1) ? LUSolver() : PETScLinearSolver(petsc_mumps_setup)

  gmg = GMGLinearSolver(
    mh, trials, tests, weakforms, prolongations, restrictions,
    pre_smoothers=smoothers, post_smoothers=smoothers, coarsest_solver=coarsest_solver,
    maxiter=3, verbose=i_am_main(ranks), mode=:preconditioner, is_nonlinear=true
  )
  gmg.log.depth = 6
  solver = FGMRESSolver(10,gmg;m_add=5,maxiter=30,rtol=1.0e-6,verbose=i_am_main(ranks),name="UJ Block - FGMRES+GMG")
  solver.log.depth = 4
  return solver
end

function gmg_patch_smoothers(
  mh,tests,weakform;
  niter = 5,
  w = 0.2,
  is_nonlinear = true,
  patch_decompositions = PatchDecomposition(mh)
)
  spaces = view(map(GridapSolvers.get_fe_space,tests),1:num_levels(tests)-1)
  patch_spaces = PatchFESpace(tests,patch_decompositions)
  smoothers = map(patch_decompositions,patch_spaces,spaces) do PD, Ph, Vh
    psolver = PatchBasedLinearSolver(weakform(PD),Ph,Vh;is_nonlinear)
    if w < 0
      solver = GMRESSolver(niter;Pr=psolver,maxiter=niter)
      patch_smoother = RichardsonSmoother(solver,1,1.0)
    else
      patch_smoother = RichardsonSmoother(psolver,niter,w)
    end
  end
  return smoothers
end

function gmg_patch_prolongations(sh,lhs,rhs;is_nonlinear=true)
  map(view(linear_indices(sh),1:num_levels(sh)-1)) do lev
    cparts = get_level_parts(sh,lev+1)
    if i_am_in(cparts)
      model = get_model_before_redist(sh,lev)
      PD = GridapSolvers.PatchBasedSmoothers.CoarsePatchDecomposition(model)
      lhs_i = lhs(PD)
      rhs_i = rhs(PD)
    else
      PD, lhs_i, rhs_i = nothing, nothing, nothing
    end
    # TODO: We should give it both tests and trials in the nonlinear case
    PatchProlongationOperator(lev,sh,PD,lhs_i,rhs_i;is_nonlinear)
  end
end

function gmg_patch_restrictions(sh,patch_prolongations,rhs,qdegrees;kwargs...)
  map(view(linear_indices(sh),1:num_levels(sh)-1)) do lev
    qdegree = qdegrees[lev]
    cparts = get_level_parts(sh,lev+1)
    if i_am_in(cparts)
      model = get_model_before_redist(sh,lev)
      rhs_i = rhs(model)
    else
      rhs_i = nothing
    end
    Ip = patch_prolongations[lev]
    PatchRestrictionOperator(lev,sh,Ip,rhs_i,qdegree;kwargs...)
  end
end
