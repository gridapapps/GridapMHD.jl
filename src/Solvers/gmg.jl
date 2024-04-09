
function get_block_solver(::Val{:gmg},params)
  vars = Tuple(params[:multigrid][:variables])
  gmg_solver(Val(vars),params)
end

function gmg_solver(::Val{(1,3)},params)
  mh = params[:multigrid][:mh]

  trials = MultiFieldFESpace(map(s -> s, params[:multigrid][:trials][[1,3]]))
  tests  = MultiFieldFESpace(map(s -> s, params[:multigrid][:tests][[1,3]]))
  tests_u, tests_j = params[:multigrid][:tests][1], params[:multigrid][:tests][3]
  
  nlevs = num_levels(mh)
  k = params[:fespaces][:k]
  qdegree = map(lev -> 2*k+1,1:nlevs)
  
  function jacobian_uj(model)
    _, dΩf, α, β, γ, σf, f, B, ζ = retrieve_fluid_params(model,params,k)
    _, dΩs, σs = retrieve_solid_params(model,params,k)
    _, params_thin_wall, _, _ = retrieve_bcs_params(model,params,k)
    Πp = params[:fespaces][:Πp]
    function a(x,dx,y)
      u, j = x
      du, dj = dx
      v_u, v_j = y
      r = a_mhd_u_u(du,v_u,β,dΩf) + a_mhd_u_j(dj,v_u,γ,B,dΩf) + a_mhd_j_u(du,v_j,σf,B,dΩf) + a_mhd_j_j(dj,v_j,dΩf)
      r = r + dc_mhd_u_u(u,du,v_u,α,dΩf)
      for p in params_thin_wall
        r = r + a_thin_wall_j_j(dj,v_j,p...)
      end
      if !isnothing(dΩs)
        r = r + a_solid_j_j(dj,v_j,dΩs)
      end
      if abs(ζ) > eps(typeof(ζ))
        r = r + a_al_u_u(du,v_u,ζ,Πp,dΩf) + a_al_j_j(dj,v_j,ζ,dΩf)
      end
      return r
    end
    return a
  end

  function jacobian_u(model)
    _, dΩf, α, β, γ, σf, f, B, ζ = retrieve_fluid_params(model,params,k)
    Πp = params[:fespaces][:Πp]
    function a(u,du,v_u)
      r = a_mhd_u_u(du,v_u,β,dΩf) + dc_mhd_u_u(u,du,v_u,α,dΩf)
      if abs(ζ) > eps(typeof(ζ))
        r = r + a_al_u_u(du,v_u,ζ,Πp,dΩf)
      end
      return r
    end
    return a
  end

  function projection_rhs(model)
    _, dΩf, α, β, γ, σf, f, B, ζ = retrieve_fluid_params(model,params,k)
    Πp = params[:fespaces][:Πp]
    l(du,v_u) = a_al_u_u(du,v_u,ζ,Πp,dΩf)
    return l
  end

  weakforms = map(mhl -> jacobian_uj(GridapSolvers.get_model(mhl)),mh)
  smoothers = gmg_patch_smoothers(mh,tests,jacobian_uj)
  restrictions = setup_restriction_operators(tests,qdegree;mode=:residual,solver=LUSolver())
  prolongations_u = gmg_patch_prolongations(tests_u,jacobian_u,projection_rhs)
  prolongations_j = setup_prolongation_operators(tests_j,qdegree)
  prolongations = MultiFieldTransferOperator(tests,[prolongations_u,prolongations_j])

  return gmg_solver(mh,trials,tests,weakforms,restrictions, prolongations,smoothers)
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
  gmg.log.depth += 4
  #solver = FGMRESSolver(10,gmg;m_add=5,maxiter=30,rtol=1.0e-6,verbose=i_am_main(ranks),name="UJ Block - FGMRES+GMG")
  #solver.log.depth += 3 # For printing purposes
  return gmg
end

function gmg_patch_smoothers(mh,tests,weakform)
  patch_decompositions = PatchDecomposition(mh)
  spaces = view(map(GridapSolvers.get_fe_space,tests),1:num_levels(tests)-1)
  patch_spaces = PatchFESpace(tests,patch_decompositions)
  smoothers = map(patch_decompositions,patch_spaces,spaces) do PD, Ph, Vh
    psolver = PatchBasedLinearSolver(weakform(PD),Ph,Vh,is_nonlinear=true)
    RichardsonSmoother(psolver,10,0.2)
  end
  return smoothers
end

function gmg_patch_prolongations(sh,lhs,rhs)
  map(view(linear_indices(sh),1:num_levels(sh)-1)) do lev
    cparts = get_level_parts(sh,lev+1)
    if i_am_in(cparts)
      model = get_model_before_redist(sh,lev)
      PD = PatchDecomposition(model)
      lhs_i = lhs(PD)
      rhs_i = rhs(PD)
    else
      PD, lhs_i, rhs_i = nothing, nothing, nothing
    end
    PatchProlongationOperator(lev,sh,PD,lhs_i,rhs_i;is_nonlinear=true)
  end
end
