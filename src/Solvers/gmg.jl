
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
  qdegree = map(lev -> 2*k+1,1:nlevs)
  
  function jacobian_uj(model)
    _, dΩf, α, β, γ, σf, f, B, ζ = retrieve_fluid_params(model,params,k)
    _, dΩs, σs = retrieve_solid_params(model,params,k)
    _, params_thin_wall, _, _ = retrieve_bcs_params(model,params,k)
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
        r = r + a_al_u_u(du,v_u,ζ,dΩf) + a_al_j_j(dj,v_j,ζ,dΩf)
      end
      return r
    end
    return a 
  end

  return gmg_solver(mh,trials,tests,jacobian_uj,qdegree)
end

function gmg_solver(mh,trials,tests,weakform,qdegree)
  ranks = get_level_parts(mh,1)
  smatrices = compute_gmg_matrices(mh,trials,tests,weakform)
  projection_solver = LUSolver()#CGSolver(JacobiLinearSolver();rtol=1.e-6)
  restrictions, prolongations = setup_transfer_operators(trials,
                                                         qdegree;
                                                         mode=:residual,
                                                         solver=projection_solver)

  smoothers = gmg_patch_smoothers(mh,trials,weakform)

  cranks = get_level_parts(mh,num_levels(mh))
  coarsest_solver = (length(cranks) == 1) ? LUSolver() : PETScLinearSolver(petsc_mumps_setup)

  gmg = GMGLinearSolver(mh,smatrices,prolongations,restrictions,
                        pre_smoothers=smoothers,
                        post_smoothers=smoothers,
                        coarsest_solver=coarsest_solver,
                        maxiter=1,verbose=true,mode=:preconditioner)
  gmg.log.depth += 3
  solver = FGMRESSolver(10,gmg;m_add=5,maxiter=30,rtol=1.0e-6,verbose=i_am_main(ranks),name="UJ Block - FGMRES+GMG")
  solver.log.depth += 3 # For printing purposes
  return solver
end

function compute_gmg_matrices(mh,trials,tests,weakform)
  nlevs = num_levels(trials)

  mats = Vector{PSparseMatrix}(undef,nlevs)
  for lev in 1:nlevs
    parts = get_level_parts(mh,lev)
    if i_am_in(parts)
      model = GridapSolvers.get_model(mh,lev)
      U = GridapSolvers.get_fe_space(trials,lev)
      V = GridapSolvers.get_fe_space(tests,lev)
      u0 = zero(U)
      jac = weakform(model)
      a(u,v) = jac(u0,u,v)
      mats[lev] = assemble_matrix(a,U,V)
    end
  end
  return mats
end

function gmg_patch_smoothers(mh,trials,weakform)
  patch_decompositions = PatchDecomposition(mh)
  patch_spaces = PatchFESpace(trials,patch_decompositions)

  nlevs = num_levels(mh)
  smoothers = Vector{RichardsonSmoother}(undef,nlevs-1)
  for lev in 1:nlevs-1
    parts = get_level_parts(mh,lev)
    if i_am_in(parts)
      PD = patch_decompositions[lev]
      Ph = GridapSolvers.get_fe_space(patch_spaces,lev)
      Vh = GridapSolvers.get_fe_space(trials,lev)
      a = weakform(PD)
      patch_smoother = PatchBasedLinearSolver(a,Ph,Vh,is_nonlinear=true)
      smoothers[lev] = RichardsonSmoother(patch_smoother,10,0.2)
    end
  end
  return smoothers
end

"""
function Gridap.Algebra.numerical_setup(ss::GridapSolvers.LinearSolvers.GMGSymbolicSetup,mat::AbstractMatrix)
  mh              = ss.solver.mh
  pre_smoothers   = ss.solver.pre_smoothers
  post_smoothers  = ss.solver.post_smoothers
  smatrices       = ss.solver.smatrices
  coarsest_solver = ss.solver.coarsest_solver

  A = all(PartitionedArrays.getany(own_values(mat)) .== PartitionedArrays.getany(own_values(smatrices[1])))
  B = all(PartitionedArrays.getany(partition(mat)) .== PartitionedArrays.getany(partition(smatrices[1])))
  C = PartitionedArrays.matching_local_indices(axes(mat,1),axes(smatrices[1],1))
  D = PartitionedArrays.matching_local_indices(axes(mat,2),axes(smatrices[1],2))
  println(" >>>> GMG checks: ",A,B,C,D)

  finest_level_cache = GridapSolvers.LinearSolvers.gmg_finest_level_cache(mh,smatrices)
  work_vectors = GridapSolvers.LinearSolvers.gmg_work_vectors(mh,smatrices)
  pre_smoothers_caches = GridapSolvers.LinearSolvers.gmg_smoothers_caches(mh,pre_smoothers,smatrices)
  if !(pre_smoothers === post_smoothers)
    post_smoothers_caches = GridapSolvers.LinearSolvers.gmg_smoothers_caches(mh,post_smoothers,smatrices)
  else
    post_smoothers_caches = pre_smoothers_caches
  end
  coarsest_solver_cache = GridapSolvers.LinearSolvers.gmg_coarse_solver_caches(mh,coarsest_solver,smatrices,work_vectors)

  return GridapSolvers.LinearSolvers.GMGNumericalSetup(ss.solver,finest_level_cache,pre_smoothers_caches,post_smoothers_caches,coarsest_solver_cache,work_vectors)
end

function GridapSolvers.LinearSolvers.apply_GMG_level!(lev::Integer,xh::Union{PVector,Nothing},rh::Union{PVector,Nothing},ns::GridapSolvers.LinearSolvers.GMGNumericalSetup)
  mh = ns.solver.mh
  parts = get_level_parts(mh,lev)
  if i_am_in(parts)
    if (lev == num_levels(mh)) 
      ## Coarsest level
      solve!(xh, ns.coarsest_solver_cache, rh)
    else 
      ## General case
      Ah = ns.solver.smatrices[lev]
      restrict, interp = ns.solver.restrict[lev], ns.solver.interp[lev]
      dxh, Adxh, dxH, rH = ns.work_vectors[lev]

      println("       >>>> Initial norm : ",norm(rh))
      # Pre-smooth current solution
      solve!(xh, ns.pre_smoothers_caches[lev], rh)
      println("       >>>> Norm after pre-smoothing : ",norm(rh))

      # Restrict the residual
      mul!(rH,restrict,rh)

      # Apply next_level
      !isa(dxH,Nothing) && fill!(dxH,0.0)
      GridapSolvers.LinearSolvers.apply_GMG_level!(lev+1,dxH,rH,ns)

      # Interpolate dxH in finer space
      mul!(dxh,interp,dxH)

      # Update solution & residual
      xh .= xh .+ dxh
      mul!(Adxh, Ah, dxh)
      rh .= rh .- Adxh
      println("       >>>> Norm after coarse correction : ",norm(rh))

      # Post-smooth current solution
      solve!(xh, ns.post_smoothers_caches[lev], rh)
      println("       >>>> Norm after post-smoothing : ",norm(rh))
    end
  end
end

function Gridap.Algebra.solve!(y::AbstractVector,ns::GridapSolvers.LinearSolvers.GMGNumericalSetup,b::AbstractVector)
  mode = ns.solver.mode
  log  = ns.solver.log

  x = similar(y)
  rh = ns.finest_level_cache
  if (mode == :preconditioner)
    fill!(x,0.0)
    copy!(rh,b)
  else
    Ah = ns.solver.smatrices[1]
    mul!(rh,Ah,x)
    rh .= b .- rh
  end

  res  = norm(rh)
  done = GridapSolvers.SolverInterfaces.init!(log,res)
  while !done
    GridapSolvers.LinearSolvers.apply_GMG_level!(1,x,rh,ns)
    res  = norm(rh)
    done = GridapSolvers.SolverInterfaces.update!(log,res)
  end

  GridapSolvers.SolverInterfaces.finalize!(log,res)
  copy!(y,x)
  return y
end
"""