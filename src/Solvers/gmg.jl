
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
  _, params_thin_wall, _, _ = retrieve_bcs_params(params,k)
  
  function jacobian_uj(x,dx,dy,dΩ)
    dΩf, dΩs, nΓ_tw, dΓ_tw = dΩ
    du, dj = dx
    v_u, v_j = dy

    # TODO: Again, this is VERY BAD
    _u, _j = x
    if isa(dΩf,Gridap.Measure)
      u = getany(local_views(_u))
    else
      u = _u
    end

    r = a_mhd_u_u(du,v_u,β,dΩf) + a_mhd_u_j(dj,v_u,γ,B,dΩf) + a_mhd_j_u(du,v_j,σf,B,dΩf) + a_mhd_j_j(dj,v_j,dΩf)
    #r = r + dc_mhd_u_u(u,du,v_u,α,dΩf)
    for (i,p) in enumerate(params_thin_wall)
      τ,cw,jw,_,_ = p
      r = r + a_thin_wall_j_j(dj,v_j,τ,cw,jw,nΓ_tw[i],dΓ_tw[i])
    end
    if !isnothing(dΩs)
      r = r + a_solid_j_j(dj,v_j,dΩs)
    end
    if abs(ζ) > eps(typeof(ζ))
      r = r + a_al_u_u(du,v_u,ζ,dΩf) + a_al_j_j(dj,v_j,ζ,dΩf)
    end
    return r
  end

  function build_measures_uj(model)
    _, dΩf, _, _, _, _, _, _, _ = retrieve_fluid_params(model,params,k)
    _, dΩs, _ = retrieve_solid_params(model,params,k)

    nΓ_tw = GridapDistributed.DistributedCellField[]
    dΓ_tw = GridapDistributed.DistributedMeasure[]
    for i in 1:length(params[:bcs][:thin_wall])
      Γ = _boundary(model,params[:bcs][:thin_wall][i][:domain])
      push!(dΓ_tw,Measure(Γ,2*k))
      push!(nΓ_tw,get_normal_vector(Γ))
    end

    return dΩf, dΩs, nΓ_tw, dΓ_tw
  end

  return gmg_solver(mh,trials,tests,jacobian_uj,build_measures_uj,qdegree)
end

# TODO: This is VERY BAD, we need to find a better way to do this
function GridapDistributed.local_views(a::Tuple{<:GridapDistributed.DistributedMeasure,Nothing,Vector{<:GridapDistributed.DistributedCellField},Vector{<:GridapDistributed.DistributedMeasure}})
  _dΩf, _dΩs, _nΓ_tw, _dΓ_tw = a
  dΩf = local_views(_dΩf)
  #dΩs = local_views(_dΩs)
  if isempty(_nΓ_tw)
    map(dΩf) do dΩf
      (dΩf,nothing,[],[])
    end
  else
    nΓ_tw = map(local_views,_nΓ_tw) |> GridapDistributed.to_parray_of_arrays
    dΓ_tw = map(local_views,_dΓ_tw) |> GridapDistributed.to_parray_of_arrays
    map(dΩf,nΓ_tw,dΓ_tw) do dΩf,nΓ_tw,dΓ_tw
      (dΩf,nothing,nΓ_tw,dΓ_tw)
    end
  end
end

function gmg_solver(mh,trials,tests,biform,measures,qdegree)
  ranks = get_level_parts(mh,1)
  smatrices = compute_gmg_matrices(mh,trials,tests,biform,measures,qdegree)
  restrictions, prolongations = setup_transfer_operators(trials,
                                                         qdegree;
                                                         mode=:residual,
                                                         solver=CGSolver(JacobiLinearSolver();rtol=1.e-6))

  smoothers = gmg_patch_smoothers(mh,trials,biform,measures,qdegree)

  cranks = get_level_parts(mh,num_levels(mh))
  coarsest_solver = (length(cranks) == 1) ? LUSolver() : PETScLinearSolver(petsc_mumps_setup)

  gmg = GMGLinearSolver(mh,smatrices,prolongations,restrictions,
                        pre_smoothers=smoothers,
                        post_smoothers=smoothers,
                        coarsest_solver=coarsest_solver,
                        maxiter=1,verbose=false,mode=:preconditioner)
  solver = GMRESSolver(10;Pr=gmg,m_add=5,maxiter=30,rtol=1.0e-6,verbose=i_am_main(ranks))
  solver.log.depth += 1 # For printing purposes
  return solver
end

function compute_gmg_matrices(mh,trials,tests,biform,measures,qdegree)
  nlevs = num_levels(trials)

  mats = Vector{PSparseMatrix}(undef,nlevs)
  for lev in 2:nlevs
    parts = get_level_parts(mh,lev)
    if i_am_in(parts)
      model = GridapSolvers.get_model(mh,lev)
      U = GridapSolvers.get_fe_space(trials,lev)
      V = GridapSolvers.get_fe_space(tests,lev)
      u0 = zero(U)
      dΩ = measures(model)
      a(u,v) = biform(u0,u,v,dΩ)
      mats[lev] = assemble_matrix(a,U,V)
    end
  end
  return mats
end

function gmg_patch_smoothers(mh,tests,biform,measures,qdegree)
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
      u0 = zero(Vh)
      dΩ = measures(PD)
      local_solver = LUSolver()
      a(u,v,dΩ) = biform(u0,u,v,dΩ)
      patch_smoother = PatchBasedLinearSolver(a,Ph,Vh,dΩ,local_solver)
      smoothers[lev] = RichardsonSmoother(patch_smoother,10,0.1)
    end
  end
  return smoothers
end
