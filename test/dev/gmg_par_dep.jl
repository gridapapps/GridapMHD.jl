# Testing GMG for the velocity block
using FileIO
using BSON
using PartitionedArrays
using Gridap
using Gridap.Helpers
using Gridap.Geometry
using Gridap.MultiField

using GridapDistributed
using GridapSolvers
using GridapSolvers.BlockSolvers
using GridapSolvers.MultilevelTools
using GridapSolvers.PatchBasedSmoothers

using GridapMHD
using GridapMHD.Meshers

# Smoothers (see GridapSolvers GMG test)
function get_patch_smoothers_hdiv(sh,biform,qdegree)
  nlevs = num_levels(sh)
  smoothers = map(view(sh,1:nlevs-1)) do shl
    model = get_model(shl)
    ptopo = Geometry.PatchTopology(ReferenceFE{0},model)
    space = get_fe_space(shl)
    Ω = Geometry.PatchTriangulation(model,ptopo)
    Λ = Skeleton(Ω)
    Γ = Boundary(Ω)
    ap = (u,v) -> biform(u,v,Ω,Λ,Γ,qdegree)
    solver = PatchBasedSmoothers.PatchSolver(
      ptopo, space, space, ap;
      assembly = :star,
      collect_factorizations = true,
      is_nonlinear = false
    )
    return RichardsonSmoother(solver,10,0.2)
  end
  return smoothers
end
function get_block_jacobi_smoothers(sh)
  nlevs = num_levels(sh)
  smoothers = map(view(sh,1:nlevs-1)) do shl
    model = get_model(shl)
    ptopo = Geometry.PatchTopology(ReferenceFE{0},model)
    space = get_fe_space(shl)
    solver = PatchBasedSmoothers.BlockJacobiSolver(space, ptopo; assembly=:star)
    return RichardsonSmoother(solver,10,0.2)
  end
  return smoothers
end

# Formulation
function get_scaling(info)
  ν = info[:ν]
  γ = info[:γ]
  ζ = info[:ζ]
  if info[:scaling] == :small
    ν = ν/ζ
    γ = γ/ζ
    ζ = 1.0
  end
  return ν,γ,ζ
end

function get_bilinear_form(model,info)
  Ω = Triangulation(model)
  if info[:FESpace] == :RT
    Λ = Skeleton(model)
    Γ = Boundary(Ω)  # Γ = Boundary(model) gives an error with #polytope branches
    return (u0,u,v) -> biform_dg(u,v,Ω,Λ,Γ,info)
  elseif info[:FESpace] == :Qk
    return (u0,u,v) -> biform(u,v,Ω,info)
  else
    error("Only RaviartThomas and Lagrangian FE spaces are implemented")
  end
  return a
end
# DG formulation
function biform_dg(u,v,Ω,Λ,Γ,info) 
  ν,γ,ζ = get_scaling(info)
  B = info[:B]
  order = info[:FE_order]
  μ = order*(order+1)*ν

  dΩ = Measure(Ω,2*order)
  dΛ = Measure(Λ,2*order)
  dΓ = Measure(Γ,2*order)

  n_Λ = get_normal_vector(Λ)
  n_Γ = get_normal_vector(Γ)

  if num_cell_dims(Γ) == 1
    h_Λ = CellField(map(get_array,local_views(∫(1)dΛ)),Λ)
    h_Γ = CellField(map(get_array,local_views(∫(1)dΓ)),Γ)
  elseif num_cell_dims(Γ) == 2
    h_Λ = CellField(map(x->x.^(1/2),map(get_array,local_views(∫(1)dΛ))),Λ)
    h_Γ = CellField(map(x->x.^(1/2),map(get_array,local_views(∫(1)dΓ))),Γ)
  end

  a = ∫( ν*(∇(v)⊙∇(u)) + ζ*(∇⋅v)*(∇⋅u) + γ*((v×B)⋅(u×B)) )*dΩ + # + s*(v⋅u)
      ∫((μ/h_Γ)*(v⋅u) - ν*(v⋅(n_Γ⋅∇(u))+(n_Γ⋅∇(v))⋅u) )*dΓ +
      ∫(
        (μ/h_Λ)*(jump(v⊗n_Λ)⊙jump(u⊗n_Λ)) -
          ν*(jump(v⊗n_Λ)⊙mean(∇(u)) + mean(∇(v))⊙jump(u⊗n_Λ))
      )*dΛ
  return a
end
# H1 formulation
function biform(u,v,Ω,info)
  ν,γ,ζ = get_scaling(info)
  B = info[:B]
  order = info[:FE_order]
  dΩ = Measure(Ω,2*order)
  return ∫( ν*(∇(v)⊙∇(u)) + ζ*(∇⋅v)*(∇⋅u) + γ*((v×B)⋅(u×B)))*dΩ 
end

function get_force(info)
  if info[:solution] == :hunt
    f = x -> VectorValue(0.,0.,1.)
  elseif info[:solution] == :hunt2
    γ = info[:γ]
    f = x -> VectorValue(0.,0.,1.) + γ*VectorValue(1.,0.,0.)
  elseif info[:solution] == :exact
    uexBxB = info[:uexBxB]
    ue = info[:ue]
    ν,γ,_ = get_scaling(info)
    f = x -> (-ν*Δ(ue)(x)-γ*uexBxB(x))
  end
  return f
end

function liform_dg(v,Ω,Γ,info)
  
  ν,γ,_ = get_scaling(info)
  order = info[:FE_order]
  μ = order*(order+1)*ν
  ue = info[:ue]
  f = get_force(info)

  dΩ = Measure(Ω,2*order)
  dΓ = Measure(Γ,2*order)
  n_Γ = get_normal_vector(Γ)

  if num_cell_dims(Γ) == 1
    h_Γ = CellField(map(get_array,local_views(∫(1)dΓ)),Γ)
  elseif num_cell_dims(Γ) == 2
    h_Γ = CellField(map(x->x.^(1/2),map(get_array,local_views(∫(1)dΓ))),Γ)
  end

  return ∫(v⋅f)dΩ + ∫( (μ/h_Γ)*(v⋅ue) - ν*((n_Γ⋅∇(v))⋅ue))*dΓ
end

function liform(v,Ω,info)
  f = get_force(info)
  order = info[:FE_order]
  dΩ = Measure(Ω,2*order)
  return ∫(v⋅f)dΩ
end

# Geometry
function create_model(parts,np_per_level,info)
    if info[:is_periodic]  # 3D periodic
      nc = info[:nc]
      cells = (nc,nc,3)

      if info[:domain] == :unit
        domain = (0,1,0,1,0,0.1)
      elseif info[:domain] == :std
        domain = (-1,1,-1,1,0,0.1)
      elseif info[:domain] == :std2
        domain = (-1,1,-1,1,0,0.2)
      elseif info[:domain] == :lz
        num_levels = length(np_per_level)
        nx = nc*2^(num_levels-1)
        Lz = 3.0/nx
        domain = (0,1,0,1,0,Lz)
      end

      _np_per_level = map(x->(x[1],x[2],1),np_per_level)
      if info[:is_stretched]
        Ha = sqrt(info[:γ])
        mh = CartesianModelHierarchy(parts,_np_per_level,domain,cells;nrefs = (2,2,1),isperiodic = (false,false,true),map = GridapMHD.Meshers.hunt_stretch_map(1.0,Ha,1,1,info[:is_stretched]))
      else
        mh = CartesianModelHierarchy(parts,_np_per_level,domain,cells;nrefs = (2,2,1),isperiodic = (false,false,true))
      end
    else
      D = length(np_per_level[1])
      domain = (D == 2) ? (0,1,0,1) : (0,1,0,1,0,1)
      cells = (D == 2) ? (4,4) : (4,4,4)
      mh = CartesianModelHierarchy(parts,np_per_level,domain,cells)
      # println("Created $(D)D model hierarchy")
    end
    return mh
end

function gmg_hdiv(parts,t,info,mh)

    tic!(t;barrier=true)
    k = info[:FE_order]
    if info[:FESpace] == :RT
      RTk_FE = ReferenceFE(raviart_thomas,Float64,k-1)
      # tests  = TestFESpace(mh,RTk_FE,dirichlet_tags="boundary")
      Ωh = Triangulation(mh)
      tests  = TestFESpace(Ωh,RTk_FE,dirichlet_tags="boundary")
    elseif info[:FESpace] == :Qk
      Qk_FE = ReferenceFE(lagrangian,VectorValue{info[:D],Float64},k)
      tests  = TestFESpace(mh,Qk_FE,dirichlet_tags="boundary")
      Pkm1_FE = ReferenceFE(lagrangian,Float64,k-1,space=:P)
      Πp = MultilevelTools.LocalProjectionMap(divergence,Pkm1_FE,2*(k-1)) # Is this the correct order?
      info[:Πp] = Πp
    else
      error("Only RaviartThomas and Lagrangian FE spaces are implemented")
    end
    ue = info[:ue]
    trials = TrialFESpace(tests,ue)
    
    # Bilinear forms at all levels
    biforms = map(mhl -> get_bilinear_form(get_model(mhl),info),mh)

    # FE operator (finest level)
    model = get_model(mh,1)
    Ω = Triangulation(model)
    dΩ = Measure(Ω,2*k)
    a = get_bilinear_form(model,info)
    U = get_fe_space(trials,1)
    V = get_fe_space(tests,1)

    X = MultiFieldFESpace([U];style=BlockMultiFieldStyle(1,(1,),(1,)))
    Y = MultiFieldFESpace([V];style=BlockMultiFieldStyle(1,(1,),(1,)))

    if info[:FESpace] == :RT
      Γ = Boundary(Ω)  # Γ = Boundary(model) gives an error with #polytope branches

      # Linear operator
      # op = AffineFEOperator(a,v->liform_dg(v,Ω,Γ,info),U,V)
      # Nonlinear operator

      # res(u,v) = a(u,u,v) - liform_dg(v,Ω,Γ,info)
      # jac(u0,u,v) = a(u0,u,v)
      # op = FEOperator(res,jac,U,V)

      # Nonlinear block operator
      res(u,v) = a(u...,u...,v...) - liform_dg(v...,Ω,Γ,info)
      jac(u0,u,v) = a(u0...,u...,v...)
      op = FEOperator(res,jac,X,Y)

    elseif info[:FESpace] == :Qk
      op = AffineFEOperator(a,v->liform(v,Ω,info),U,V)
    else
      error("Only RaviartThomas and Lagrangian FE spaces are implemented")
    end

    # Multigrid components
    if info[:projection_solver] == :CG_Jacobi
      projection_solver = CGSolver(
        JacobiLinearSolver();
        maxiter=10,
        atol=info[:projection_solver_atol],
        rtol=info[:projection_solver_rtol],
        # verbose=i_am_main(parts),
        name = "Projection solver (CG_Jacobi)"
      )
      projection_solver.log.depth = 12
    else
      projection_solver = LUSolver()
    end

    qdegree = info[:qdegree]
    if info[:prolongation] == :default
      prolongations = setup_prolongation_operators(tests,qdegree;mode=:residual,solver=projection_solver)
    elseif info[:prolongation] == :patch_block_jacobi
      prolongations = PatchBasedSmoothers.setup_block_jacobi_prolongation_operators(tests)
    else
      error("Unknown prolongation operator type")
    end

    if info[:restriction] == :default
      restrictions = setup_restriction_operators(tests,qdegree;mode=:residual,solver=projection_solver)
    end

    if info[:smoother] == :block_jacobi
      smoothers = get_block_jacobi_smoothers(trials)
    end

    # coarsest_solver = PETScLinearSolver(GridapMHD.petsc_mumps_setup)
    coarsest_solver = Gridap.Algebra.LUSolver()  # this is the default

    # GMG solver
    gmg = GMGLinearSolver(
        trials,tests,biforms,
        prolongations,restrictions,
        pre_smoothers=smoothers,
        post_smoothers=smoothers,
        coarsest_solver=coarsest_solver,
        maxiter=1,
        cycle_type=info[:cycle_type],
        verbose=i_am_main(parts)
        ,is_nonlinear=true # Nonlinear solver
    )
    gmg.log.depth = 8

    atol = info[:solver_atol]
    rtol = info[:solver_rtol]
    if info[:solver] == :CG
      solver = CGSolver(gmg;maxiter=30,atol=atol,rtol=rtol,verbose=i_am_main(parts),name = "System (CG_GMG)")
    elseif info[:solver] == :FGMRES
      solver = FGMRESSolver(2,gmg;m_add=1,maxiter=30,atol=atol,rtol=rtol,verbose=i_am_main(parts),name="System (FGMRES_GMG)")
    end
    solver.log.depth = 4
    toc!(t,"setup")
    
    # Block preconditioner
    blocks = [NonlinearSystemBlock(1);;]
    P = BlockTriangularSolver(blocks,[solver])
    l_solver = FGMRESSolver(1,P;maxiter=1,rtol=rtol,atol=atol,verbose=i_am_main(parts),name="Global System - FGMRES + Block")
    l_solver.log.depth = 2

    # Nonlinear solver
    # nl_solver = GridapSolvers.NewtonSolver(solver,maxiter=1,atol=atol,rtol=rtol,verbose=i_am_main(parts))
    nl_solver = GridapSolvers.NewtonSolver(l_solver,maxiter=1,atol=atol,rtol=rtol,verbose=i_am_main(parts))

    # Nonlinear block solve
    tic!(t;barrier=true)
    xh = zero(X)
    xh, cache = solve!(xh, nl_solver,op)
    uh = xh[1]

    # Nonlinear solve
    # uh = zero(U)
    # uh, cache = solve!(uh, nl_solver,op)

    # Linear solver
    # uh = solve(solver,op)
    toc!(t,"solve")

    eh = ue-uh
    e_l2  = sum(∫(eh⋅eh)dΩ)

    if i_am_main(parts)
      println("L2 error = ", e_l2)
    end

    info[:vtk_out] && writevtk(Ω,info[:vtk_name],cellfields=["uh"=>uh])

    e_l2, solver.log.num_iters
end

const ζs = [1.0,1.0e2,1.0e4,1.0e6,1.0e8,1.0e10,1.0e12]
const npar = 7 # length(ζs)
const nruns = 2 # To measure CPU times
const path = @__DIR__

function gmg_par_dep(;D=3,
  is_periodic=true,
  is_stretched=false,
  domain=:std,
  nlev=5,
  nc=4,
  ν=1.0,
  γ0=1.0,
  γ_law=:constant, # or :variable,  in this case γ=γ0*ζ
  solution=:hunt,  # or :exact
  scaling=:large,
  fe_space=:RT,
  fe_order=1,
  qdegree=2*fe_order,
  prolongation=:default,
  restriction=:default,
  projection_solver=:CG_Jacobi,
  smoother=:block_jacobi,
  cycle_type=:v_cycle,
  solver=:FGMRES,
  tolerances=[1e-11,1e-5,1e-12,1e-6],
  vtk_out=false
  )

  if D==2 
    np=(2,2)
    name = "2D"
    ue = x -> VectorValue(x[1],-x[2]) 
  elseif D==3
    if is_periodic
      np = (2,2) 
      name = "3DP"
    else
      np=(2,2,2)
      name = "3D"
    end
    if solution == :hunt || solution == :hunt2
      ue = x -> VectorValue(0.0,0.0,0.0)
      B = VectorValue(0.0,1.0,0.0)
      uexBxB = x->VectorValue(0.0,0.0,0.0)
    elseif solution == :exact
      ue = x -> VectorValue(x[1],-x[2],0.0)
      B = VectorValue(0.0,0.0,1.0)
      uexBxB = x->(VectorValue(x[1],-x[2],0.)×B)×B
    end
  end
  is_stretched ? mesh_type = "stretched" : mesh_type = "uniform"
  γ_law == :constant ? γ_law_str = "γ$(γ0)" : γ_law_str = "γ$(γ0)xζ"
  title = "NL_BP_$(name)_$(domain)_domain_$(mesh_type)_nc$(nc)_$(solution)_$(γ_law_str)_$(fe_space)$(fe_order)_scal_$(scaling)_qdeg_$(qdegree)_$(solver)_$(cycle_type)_S_10$(smoother)_P_$(prolongation)_R_$(restriction)_Ps_$(projection_solver)_C_LU"

  info = Dict{Symbol,Any}()
  info[:title]      = title
  info[:max_levels] = nlev
  info[:max_params] = npar
  info[:D]          = D
  info[:nc]         = nc
  info[:ζs]         = ζs
  info[:ν]          = ν
  info[:γ_law]      = γ_law      # :constant or :variable
  info[:γ]          = γ0
  info[:solution]   = solution    # :hunt or :exact
  info[:ue]         = ue
  info[:B]           = D==3 ? B : nothing
  info[:uexBxB]      = D==3 ? uexBxB : nothing
  info[:is_periodic]         = is_periodic
  info[:is_stretched]        = is_stretched
  info[:domain]              = domain
  info[:scaling]             = scaling       # :small or :large
  info[:FESpace]             = fe_space      # :RT or :Qk
  info[:FE_order]            = fe_order
  info[:qdegree]             = qdegree 
  info[:cycle_type]          = cycle_type    # :v_cycle
  info[:prolongation]        = prolongation  # :patch_block_jacobi
  info[:restriction]         = restriction   # :projection
  info[:smoother]            = smoother # :patch
  info[:solver_atol]         = tolerances[1]
  info[:solver_rtol]         = tolerances[2]
  info[:projection_solver]   = projection_solver
  info[:projection_solver_atol] = tolerances[3]  # for :CG_Jacobi
  info[:projection_solver_rtol] = tolerances[4]  # for :CG_Jacobi
  info[:solver] = solver
  info[:vtk_out] = vtk_out

  niter = Matrix{Int}(undef,npar,nlev)
  error = Matrix{Float64}(undef,npar,nlev)
  time_setup = Matrix{Float64}(undef,npar,nlev)
  time_solve = Matrix{Float64}(undef,npar,nlev)
  time_model = Vector{Float64}(undef,nlev)

  # distribute = DebugArray
  with_mpi() do distribute
    parts = distribute(LinearIndices((prod(np),)))
    t = PTimer(parts,verbose=true)
    np_per_level = [np]
    for l=1:nlev
      push!(np_per_level,np)
      for nr=1:nruns
        tic!(t;barrier=true)
        global mh = create_model(parts,np_per_level,info)
        toc!(t,"model")
      end
      map_main(t.data) do data 
        time_model[l] = data["model"].max
      end
      for i=1:npar
        info[:ζ] = info[:ζs][i]
        if info[:γ_law] == :variable
           info[:γ] = γ0*info[:ζ]
        end
        i_am_main(parts) && println("Running for ζ=",info[:ζ])
        info[:vtk_name] = "solution_$i"
        for nr=1:nruns
          error[i,l], niter[i,l] = gmg_hdiv(parts,t,info,mh)
        end
        map_main(t.data) do data 
          time_setup[i,l] = data["setup"].max
          time_solve[i,l] = data["solve"].max
        end
      end
    end
    info[:niter] = niter
    info[:error] = error
    # Save results and print summary
    map_main(t.data) do data 
      info[:time_model] = time_model
      info[:time_setup] = time_setup
      info[:time_solve] = time_solve
      save(joinpath(path,"$title.bson"),info)
      # Print 
      for i=1:npar
        println("For ζ=",ζs[i]," the number of iterations per level: ",niter[i,:])
      end
    end
  end
end
