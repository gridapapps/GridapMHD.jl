# Testing GMG for the velocity block
using FileIO
using BSON
using PartitionedArrays
using Gridap
using Gridap.Helpers
using Gridap.Geometry
using GridapDistributed
using GridapSolvers
using GridapSolvers.MultilevelTools
using GridapSolvers.PatchBasedSmoothers

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
function get_scaling(info,β)
  if info[:scaling] == :small
    α = 1.0
    ν = 1.0/β
    s = 1.0/β
  elseif info[:scaling] == :large
    α = β
    ν = 1.0
    s = 1.0
  else
    error("Unknown scaling type")
  end
  return α,ν,s
end

function get_bilinear_form(model,info,β)
  order = info[:FE_order]
  Ω = Triangulation(model)
  if info[:FESpace] == :RT
    Λ = Skeleton(model)
    Γ = Boundary(Ω)  # Γ = Boundary(model) gives an error with #polytope branches
    return (u,v) -> biform_dg(u,v,Ω,Λ,Γ,order,get_scaling(info,β))
  elseif info[:FESpace] == :Qk
    return (u,v) -> biform(u,v,Ω,order,get_scaling(info,β))
  else
    error("Only RaviartThomas and Lagrangian FE spaces are implemented")
  end
  return a
end
# DG formulation
function biform_dg(u,v,Ω,Λ,Γ,order,p) 
  α,ν,s = p
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

  a = ∫( ν*(∇(v)⊙∇(u)) + s*(v⋅u) + α*(∇⋅v)*(∇⋅u))*dΩ +
      ∫((μ/h_Γ)*(v⋅u) - ν*(v⋅(n_Γ⋅∇(u))+(n_Γ⋅∇(v))⋅u) )*dΓ +
      ∫(
        (μ/h_Λ)*(jump(v⊗n_Λ)⊙jump(u⊗n_Λ)) -
          ν*(jump(v⊗n_Λ)⊙mean(∇(u)) + mean(∇(v))⊙jump(u⊗n_Λ))
      )*dΛ
  return a
end
# H1 formulation
function biform(u,v,Ω,order,p)
  α,ν,s = p
  dΩ = Measure(Ω,2*order)
  return ∫( ν*(∇(v)⊙∇(u)) + s*(v⋅u) + α*(∇⋅v)*(∇⋅u))*dΩ 
end

function liform_dg(v,Ω,Γ,order,ue,p)
  _,ν,s = p
  μ = order*(order+1)*ν
  f(x) = s*ue(x) - ν*Δ(ue)(x)
  # f(x) = s*VectorValue(0.,0.,1.) 

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

function liform(v,Ω,order,ue,p)
  _,ν,s = p
  f(x) = s*ue(x) - ν*Δ(ue)(x)
  dΩ = Measure(Ω,2*order)
  return ∫(v⋅f)dΩ
end

# Geometry
function create_model(parts,np_per_level,is_periodic)
    if is_periodic  # 3D periodic
      num_levels = length(np_per_level)
      nx = 4*2^(num_levels-1)
      Lz = 3.0/nx
      domain = (0,1,0,1,0,Lz)
      cells = (4,4,3)
      _np_per_level = map(x->(x[1],x[2],1),np_per_level)
      mh = CartesianModelHierarchy(parts,_np_per_level,domain,cells;nrefs = (2,2,1),isperiodic = (false,false,true))
    else
      D = length(np_per_level[1])
      domain = (D == 2) ? (0,1,0,1) : (0,1,0,1,0,1)
      cells = (D == 2) ? (4,4) : (4,4,4)
      mh = CartesianModelHierarchy(parts,np_per_level,domain,cells)
      # println("Created $(D)D model hierarchy")
    end
    return mh
end

function gmg_hdiv(parts,t,info,β,mh,ue)

    tic!(t;barrier=true)
    k = info[:FE_order]
    if info[:FESpace] == :RT
      RTk_FE = ReferenceFE(raviart_thomas,Float64,k-1)
      tests  = TestFESpace(mh,RTk_FE,dirichlet_tags="boundary")
    elseif info[:FESpace] == :Qk
      Qk_FE = ReferenceFE(lagrangian,VectorValue{info[:D],Float64},k)
      tests  = TestFESpace(mh,Qk_FE,dirichlet_tags="boundary")
      Pkm1_FE = ReferenceFE(lagrangian,Float64,k-1,space=:P)
      Πp = MultilevelTools.LocalProjectionMap(divergence,Pkm1_FE,2*(k-1)) # Is this the correct order?
      info[:Πp] = Πp
    else
      error("Only RaviartThomas and Lagrangian FE spaces are implemented")
    end
    trials = TrialFESpace(tests,ue)
    
    # Bilinear forms at all levels
    biforms = map(mhl -> get_bilinear_form(get_model(mhl),info,β),mh)

    # FE operator (finest level)
    model = get_model(mh,1)
    Ω = Triangulation(model)
    dΩ = Measure(Ω,2*k)
    a = get_bilinear_form(model,info,β)
    U = get_fe_space(trials,1)
    V = get_fe_space(tests,1)
    if info[:FESpace] == :RT
      Γ = Boundary(Ω)  # Γ = Boundary(model) gives an error with #polytope branches
      op = AffineFEOperator(a,v->liform_dg(v,Ω,Γ,k,ue,get_scaling(info,β)),U,V)
      # Nonlinear operator
      # res(u,v) = a(u,v) - liform_dg(v,Ω,Γ,k,β,ue)
      # jac(u0,u,v) = a(u,v)
      # op = FEOperator(res,jac,U,V)
    elseif info[:FESpace] == :Qk
      op = AffineFEOperator(a,v->liform(v,Ω,k,ue,get_scaling(info,β)),U,V)
    else
      error("Only RaviartThomas and Lagrangian FE spaces are implemented")
    end

    # Multigrid components
    if info[:projection_solver] == :CG_Jacobi
      projection_solver = CGSolver(
        JacobiLinearSolver();
        maxiter=20,
        atol=info[:projection_solver_atol],
        rtol=info[:projection_solver_rtol],
        verbose=i_am_main(parts),
        name = "Projection solver (CG_Jacobi)"
      )
      projection_solver.log.depth = 8
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
      smoothers = get_block_jacobi_smoothers(tests) 
    end

    # GMG solver
    gmg = GMGLinearSolver(
        trials,tests,biforms,
        prolongations,restrictions,
        pre_smoothers=smoothers,
        post_smoothers=smoothers,        
        maxiter=1,
        cycle_type=info[:cycle_type],
        verbose=i_am_main(parts)
        # ,is_nonlinear=true # Nonlinear solver
    )
    gmg.log.depth = 4

    atol = info[:solver_atol]
    rtol = info[:solver_rtol]
    if info[:solver] == :CG
      solver = CGSolver(gmg;maxiter=30,atol=atol,rtol=rtol,verbose=i_am_main(parts),name = "System (CG_GMG)")
    elseif info[:solver] == :FGMRES
      solver = FGMRESSolver(2,gmg;m_add=1,maxiter=30,atol=atol,rtol=rtol,verbose=i_am_main(parts),name="System (FGMRES_GMG)")
    end
    toc!(t,"setup")
    
    tic!(t;barrier=true)
    # Nonlinear solver
    # nl_solver = GridapSolvers.NewtonSolver(solver,maxiter=1,atol=atol,rtol=rtol,verbose=i_am_main(parts))
    # uh = zero(U)
    # uh, cache = solve!(uh, nl_solver,op)

    uh = solve(solver,op)
    toc!(t,"solve")

    eh = ue-uh
    e_l2  = sum(∫(eh⋅eh)dΩ)

    if i_am_main(parts)
      println("L2 error = ", e_l2)
    end

    e_l2, solver.log.num_iters
end

const βs = [1.0,1.0e2,1.0e4,1.0e6,1.0e8,1.0e10,1.0e12]
const nlev = 5
const npar = 7 # length(βs)
const nruns = 2 # To measure CPU times
const path = @__DIR__
# const ue(x) = VectorValue(x[1],-x[2]) 

function gmg_par_dep(;D=2,
  is_periodic=false,
  scaling=:small,
  fe_space=:RT,
  fe_order=1,
  qdegree=2*fe_order,
  prolongation=:default,
  restriction=:default,
  projection_solver=:LU,
  smoother=:block_jacobi,
  cycle_type=:v_cycle,
  solver=:CG,
  tolerances=[1e-11,1e-5,1e-12,1e-6],
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
    ue = x -> VectorValue(x[1],-x[2],0.0)
  end
  title = "$(name)_$(fe_space)$(fe_order)_S_$(scaling)_qdeg_$(qdegree)_$(cycle_type)_S_$(smoother)_P_$(prolongation)_R_$(restriction)_Ps_$(projection_solver)_nonlinear"

  info = Dict{Symbol,Any}()
  info[:title]      = title
  info[:max_levels] = nlev
  info[:max_params] = npar
  info[:β]          = βs
  info[:D]          = D
  info[:is_periodic]         = is_periodic
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
  # info[:projection_solver] = :CG_Jacobi
  info[:projection_solver_atol] = tolerances[3]  # for :CG_Jacobi
  info[:projection_solver_rtol] = tolerances[4]  # for :CG_Jacobi
  info[:solver] = solver

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
        global mh = create_model(parts,np_per_level,is_periodic)
        toc!(t,"model")
      end
      map_main(t.data) do data 
        time_model[l] = data["model"].max
      end
      for i=1:npar
        β=βs[i]
        i_am_main(parts) && println("Running for β=",β)
        for nr=1:nruns
          error[i,l], niter[i,l] = gmg_hdiv(parts,t,info,β,mh,ue)
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
        println("For β=",βs[i]," the number of iterations per level: ",niter[i,:])
      end
    end
  end
end
