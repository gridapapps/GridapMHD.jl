
function get_block_solver(::Val{:gmg},params)
  formulation = params[:fespaces][:formulation]
  solver = params[:solver][:solver]
  gmg_solver(Val(formulation),Val(solver),params)
end

function gmg_solver(
  trials,tests,weakforms,restrictions,prolongations,smoothers; 
  name = "GMG Solver", gmg_maxiter = 3, maxiter = 10, atol = 1.e-6, rtol = 1.e-10
)
  ranks = get_level_parts(trials,1)
  coarsest_solver = PETScLinearSolver(petsc_mumps_setup)
  gmg = GMGLinearSolver(
    trials, tests, weakforms, prolongations, restrictions,
    pre_smoothers=smoothers, post_smoothers=smoothers, coarsest_solver=coarsest_solver,
    maxiter=gmg_maxiter, verbose=i_am_main(ranks), mode=:preconditioner, is_nonlinear=true
  )
  gmg.log.depth = 6
  solver = FGMRESSolver(2,gmg;m_add=1,maxiter=maxiter,rtol=rtol,atol=atol,verbose=i_am_main(ranks),name=name)
  solver.log.depth = 4
  return solver
end

function gmg_coarse_solver(tests)
  nlevs = num_levels(tests)
  cparts = get_level_parts(tests,nlevs)
  if isa(cparts,MPIArray)
    return PETScLinearSolver(petsc_mumps_setup)
  else
    return LUSolver
  end
end

function gmg_patch_smoothers(tests,weakform;w=0.2)
  nlevs = num_levels(tests)
  smoothers = map(view(tests,1:nlevs-1)) do test
    model = get_model(test)
    if isa(model,AdaptedDiscreteModelTypes)
      model = get_model(model)
    end
    space = get_fe_space(test)
    ptopo = Geometry.PatchTopology(ReferenceFE{0},model)
    jac = weakform(PatchModel(model,ptopo))
    solver = PatchBasedSmoothers.PatchSolver(
      ptopo, space, space, jac;
      assembly = :star,
      collect_factorizations = true,
      is_nonlinear = true
    )
    if w > 0.0
      return RichardsonSmoother(solver,10,w)
    else
      return FGMRESSolver(10,solver;maxiter=10)
    end
  end
  return smoothers
end

function gmg_patch_prolongations(tests, weakform)
  nlevs = num_levels(tests)
  map(view(linear_indices(tests),1:nlevs-1)) do lev
    cparts = get_level_parts(tests,lev+1)
    if i_am_in(cparts)
      model = get_model_before_redist(tests,lev)
      ptopo = CoarsePatchTopology(model)
      jac = weakform(PatchModel(model,ptopo))
    else
      ptopo, jac = nothing, nothing
    end
    PatchProlongationOperator(
      lev, tests, ptopo, jac, jac;
      is_nonlinear = true,
      collect_factorizations = false
    )
  end
end

########################################################################################
# H1-H1 discretization, u-block solver

function gmg_solver(::Val{:H1H1},::Val{:h1h1blocks},params)
  # GMG with patch solvers and prolongation,
  # for a Navier-Stokes Augmented Lagrangian formulation

  mh = params[:multigrid][:mh]
  trials = params[:multigrid][:trials][:u]
  tests  = params[:multigrid][:tests][:u]
  
  weakform(model) = weak_form_h1_h1_u(model,params)
  jacs = map(mhl -> weakform(get_model(mhl)), mh)

  smoothers = gmg_patch_smoothers(tests, weakform)
  prolongations = gmg_patch_prolongations(tests, weakform)
  restrictions = setup_restriction_operators(
    tests,params[:fespaces][:q];mode=:residual
  )

  return gmg_solver(
    trials, tests, jacs, restrictions, prolongations, smoothers;
    name = "GMG H1-H1 u-block Solver"
  )
end

function weak_form_h1_h1_u(model,params)
  weakform_params = (
    retrieve_fluid_params(model,params), 
    nothing, # No solid params needed for u 
    nothing, nothing, [] # No boundary conditions for now
  )
  jac(x,dx,dy) = jac_h1_h1_u(x,dx,dy,weakform_params)
  return jac
end

function jac_h1_h1_u(_x,_dx,_dy, params)
  fluid_params, solid_params, params_φ, params_thin_wall, params_Λ = params

  x = setup_variable_u(_x)
  dx = setup_variable_u(_dx)
  dy = setup_variable_u(_dy)

  r = jac_fluid_h1_h1_u(x,dx,dy,fluid_params...)
  for p in params_Λ
    r = r + a_Λ(x,dx,p...)
  end

  return r
end

function jac_fluid_h1_h1_u(x,dx,dy,α,β,γ,B,σ,f,g,divg,ζᵤ,ζⱼ,Πp,convection,dΩ)
  u, ∇u = x[:u], x[:∇u]
  du, v = dx[:u], dy[:u]
  ∇du, ∇v = dx[:∇u], dy[:∇u]
  div_du, div_v = dx[:divu], dy[:divu]

  duB, vB = du×B, v×B
  u_block = β*(∇du⊙∇v) + γ*duB⋅vB
  if !iszero(ζᵤ)
    u_block += ζᵤ*(Πp(du)*div_v)
  end
  if convection == :picard
    u_block += α*v⋅(conv∘(u,∇du))
  elseif convection == :newton
    u_block += α*v⋅(conv∘(u,∇du) + conv∘(du,∇u))
  end

  return ∫(u_block) * dΩ
end

########################################################################################
# HDiv-HDiv discretization, uj-block solver

function gmg_solver(::Val{:HDivHDiv},::Val{:badia2024},params)
  mh = params[:multigrid][:mh]
  trials_u = params[:multigrid][:trials][:u]
  trials_j = params[:multigrid][:trials][:j]
  tests_u  = params[:multigrid][:tests][:u]
  tests_j  = params[:multigrid][:tests][:j]
  trials = MultiFieldFESpace([trials_u, trials_j])
  tests = MultiFieldFESpace([tests_u, tests_j])
  
  weakform(model) = weak_form_hdiv_hdiv_uj(model,params)
  jacs = map(mhl -> weakform(get_model(mhl)), mh)

  smoothers = gmg_patch_smoothers(tests, weakform)
  prolongations = setup_prolongation_operators(
    tests,params[:fespaces][:q];mode=:residual
  )
  restrictions = setup_restriction_operators(
    tests,params[:fespaces][:q];mode=:residual
  )

  return gmg_solver(
    trials, tests, jacs, restrictions, prolongations, smoothers;
    name = "GMG H1-H1 u-block Solver", gmg_maxiter = 1
  )
end

function weak_form_hdiv_hdiv_uj(model,params)
  weakform_params = (
    retrieve_fluid_params(model,params), 
    nothing, # No solid params needed for u 
    nothing, nothing, [], # No boundary conditions for now
    retrieve_hdiv_fluid_params(model,params)
  )
  jac(x,dx,dy) = jac_h1_h1_u(x,dx,dy,weakform_params)
  return jac
end

function jac_hdiv_hdiv_uj(_x,_dx,_dy, params)
  fluid_params, solid_params, params_φ, params_thin_wall, params_Λ, hdiv_params = params

  x = setup_variable_u(_x)
  dx = setup_variable_u(_dx)
  dy = setup_variable_u(_dy)

  r = jac_fluid_h1_hdiv_uj(x,dx,dy,fluid_params...)
  r = r + jac_fluid_hdiv_stab(x,dx,dy,hdiv_params...)
  for p in params_Λ
    r = r + a_Λ(x,dx,p...)
  end

  return r
end

function jac_fluid_h1_hdiv_uj(x,dx,dy,α,β,γ,B,σ,f,g,divg,ζᵤ,ζⱼ,Πp,convection,dΩ)
  u, ∇u = x[:u], x[:∇u]
  du, v = dx[:u], dy[:u]
  dj, s = dx[:j], dy[:j]
  ∇du, ∇v = dx[:∇u], dy[:∇u]
  div_du, div_v = dx[:divu], dy[:divu]
  div_dj, div_s = dx[:divj], dy[:divj]

  u_block = β*(∇du⊙∇v)
  j_block = dj⋅s 

  # Augmented Lagrangian terms
  if !iszero(ζᵤ)
    u_block += ζᵤ*(div_du*div_v)
  end
  if !iszero(ζⱼ)
    j_block += ζⱼ*(div_dj*div_s)
  end

  # Convection term
  if convection == :picard
    u_block += α*v⋅(conv∘(u,∇du))
  elseif convection == :newton
    u_block += α*v⋅(conv∘(u,∇du) + conv∘(du,∇u))
  end

  return ∫(u_block - γ*(dj×B)⋅v + j_block - σ*(du×B)⋅s)dΩ
end
