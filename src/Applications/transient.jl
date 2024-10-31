
function transient(;
  backend = nothing,
  np      = 1,
  title   = "Transient",
  path    = datadir(),
  kwargs...)

  if isa(backend,Nothing)
    @assert np == 1
    info, t = _transient(;title=title,path=path,kwargs...)
  else
    @assert backend ∈ [:sequential,:mpi]
    if backend === :sequential
      info,t = with_debug() do distribute
        _transient(;distribute=distribute,np=np,title=title,path=path,kwargs...)
      end
    else
      info,t = with_mpi() do distribute
        _transient(;distribute=distribute,np=np,title=title,path=path,kwargs...)
      end
    end
  end
  info[:np] = np
  info[:backend] = backend
  info[:title] = title
  map_main(t.data) do data
    for (k,v) in data
      info[Symbol("time_$k")] = v.max
    end
    save(joinpath(path,"$title.bson"),info)
  end

  nothing
end

function _transient(;
  distribute = nothing,
  np = nothing,
  n=2,
  nc=(n,n,n),
  k = 2,
  t0 = 0.0,
  Δt = 1,
  tf = 1,
  ν=1.0,
  ρ=1.0,
  σ=1.0,
  B=(0.0,0.0,1.0),
  f=(0.0,0.0,1.0),
  ζ=0.0,
  L=1.0,
  u0=1.0,
  B0=norm(VectorValue(B)),
  nsums = 10,
  vtk=true,
  title = "transient",
  path  = datadir(),
  debug = false,
  res_assemble = false,
  jac_assemble = false,
  solve = true,
  solver = :julia,
  verbose = true,
  man_solution = nothing,
  max_error = 0.0,
  time_solver = :theta,
  θ = 0.5,
  μ = 0.0,
  )

  info = Dict{Symbol,Any}()
  params = Dict{Symbol,Any}(
    :debug=>debug,
    :solve=>solve,
    :res_assemble=>res_assemble,
    :jac_assemble=>jac_assemble,
  )

  # Communicator
  if isa(distribute,Nothing)
    @assert isa(np,Nothing)
    np = Tuple(fill(1,length(nc)))
    distribute = DebugArray
  end
  @assert length(np) == length(nc)
  parts = distribute(LinearIndices((prod(np),)))

  # Timer
  t = PTimer(parts,verbose=verbose)
  params[:ptimer] = t
  tic!(t,barrier=true)

  # Solver
  if isa(solver,Symbol)
    solver = default_solver_params(Val(solver))
  end
  params[:solver] = solver

  # Reduced quantities
  Re = u0*L/ν
  Ha = B0*L*sqrt(σ/(ρ*ν))
  N = Ha^2/Re
  α = 1.0
  β = 1.0/Re
  γ = N
  B̄ = (1/B0)*VectorValue(B)
  # f̄ = TimeSpaceFunction(t->x-> (L/(ρ*u0^2)) * f(t,x))

  is_manufactured = !isnothing(man_solution)
  if is_manufactured
    u,j,p,φ,f̄,g = _transient_solution(man_solution;B=B̄,α,β,γ)
  end

  domain = ( 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 )
  params[:model] = _box_model(parts,np,domain,nc)


  if debug && vtk
    writevtk(params[:model],joinpath(path,"$(title)_model"))
  end

  # FE Space parameters
  params[:fespaces] = Dict(
  :k => k,
  :p_constraint => :zeromean)

  # Fluid parameters
  params[:fluid] = Dict(
    :domain=>nothing,
    :α=>α,
    :β=>β,
    :γ=>γ,
    :f=>f̄,
    :B=>B̄,
    :ζ=>ζ,
    :g=>g,
  )

  # Boundary conditions
  params[:bcs] = Dict(
    :u=>Dict(:tags=>"boundary",:values=>u),
    :j=>Dict(:tags=>"top_bottom",:values=>j),
    :φ=>Dict(:domain=>"walls",:value=>φ),
  )

  if μ > 0
    params[:bcs][:stabilization] = Dict(:μ=>μ)
  end

  # Setup ODE solver
  params[:transient] = Dict{Symbol,Any}(
    :solver => ode_solver,
    :t0 => t0,
    :tf => tf,
    :Δt => Δt,
  )
  params[:transient][:solver] = Dict{Symbol,Any}(
    :solver => time_solver,
    :θ => θ
  )
  if is_manufactured
    params[:x0] = Dict{Symbol,Any}(
      :u => u(t0),
      :p => p(t0),
      :j => j(t0),
      :φ => φ(t0),
    )
  end

  toc!(t,"pre_process")

  # Lazy ODE solver
  xh,fullparams,info = main(params;output=info)

  Ω = Interior(params[:model])
  Ω_phys = L == 1.0 ? Ω : _warp(model,Ω,L)


  degree = k*2
  dΩ = Measure(Ω_phys,degree)
  l2(u,dΩ) = √(∑( ∫(u ⊙ u)dΩ ))
  h1(u,dΩ) = √(∑( ∫( ∇(u) ⊙ ∇(u) + u ⊙ u)dΩ ))

  nt = ceil(Int, (tf-t0)/Δt )
  results = Dict(
    :t => Δt:Δt:tf,
    :uh_el2 => zeros(nt),
    :uh_eh1 => zeros(nt),
    :jh_el2 => zeros(nt),
    :ph_el2 => zeros(nt),
    :φh_el2 => zeros(nt),
  )

  tic!(t,barrier=true)
  if vtk
    pvd = createpvd(parts,joinpath(path,title))
    pvd[t0] = createvtk(Ω_phys,joinpath(path,title,"$(title)_0"),
      order=2,
      cellfields=[
        "uh"=>u(t0),"ph"=>p(t0),"jh"=>j(t0),"phi"=>φ(t0),
        "u"=>u(t0),"j"=>j(t0)])
  end

  for (i,(t,xht)) in enumerate(xh)
    ūh,p̄h,j̄h,φ̄h = xht

    # Rescale quantities
    uh = u0*ūh
    ph = (ρ*u0^2)*p̄h
    jh = (σ*u0*B0)*j̄h
    φh = (u0*B0*L)*φ̄h

    if is_manufactured
      e_uh = uh - u(t)
      e_jh = jh - j(t)
      e_ph = (ph) - (p(t))
      e_φh = (φh) - (φ(t))
      uh_el2 = l2(e_uh,dΩ)
      jh_el2 = l2(e_jh,dΩ)
      ph_el2 = l2(e_ph,dΩ)
      φh_el2 = l2(e_φh,dΩ)
      uh_eh1 = h1(e_uh,dΩ)
      results[:uh_el2][i] = uh_el2
      results[:jh_el2][i] = jh_el2
      results[:ph_el2][i] = ph_el2
      results[:φh_el2][i] = φh_el2
      results[:uh_eh1][i] = uh_eh1
      if max_error > 0
        @assert uh_el2 < max_error
        @assert jh_el2 < max_error
        @assert ph_el2 < max_error
        @assert φh_el2 < max_error
        @assert uh_eh1 < max_error
      end
      _print_timestep_errors(results,i;verbose=verbose)
    end
    if vtk
      pvd[t] = createvtk(Ω_phys,joinpath(path,title,"$(title)_$i"),
        order=2,
        cellfields=[
        "uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh,
        "u"=>u(t),"j"=>j(t)])
    end
  end
  if vtk
    savepvd(pvd)
  end
  toc!(t,"time_stepping")

  info[:ncells] = num_cells(params[:model])
  info[:Re] = Re
  info[:Ha] = Ha
  info[:Δt] = Δt
  info[:t0] = t0
  info[:tf] = tf
  info[:θ] = θ

  merge(info, results), t
end

_transient_solution(sol;kwargs...) = _trasient_solution(Val(sol);kwargs...)

function _trasient_solution(::Val{:stationary_fespace};kwargs...)
  u = TimeSpaceFunction(t->x-> VectorValue( x[2], x[1], 0 ) )
  j = TimeSpaceFunction(t->x-> VectorValue( x[2], 1, 0 ) )
  p = TimeSpaceFunction(t->x-> 0 )
  φ = TimeSpaceFunction(t->x-> 1 )
  f = _transient_solution_f(;u,j,p,φ,kwargs...)
  g = _transient_solution_fj(;u,j,φ,kwargs...)
  u,j,p,φ,f,g
end

function _trasient_solution(::Val{:lineartime_fespace};kwargs...)
  u = TimeSpaceFunction(t->x-> VectorValue( x[2]*(1-t/2), x[1]*(1-t/2), 0 ) )
  j = TimeSpaceFunction(t->x-> VectorValue( x[2]*t, (1-t/2), 0 ) )
  p = TimeSpaceFunction(t->x-> 0 )
  φ = TimeSpaceFunction(t->x-> t*one(x[1]) )
  f = _transient_solution_f(;u,j,p,φ,kwargs...)
  g = _transient_solution_fj(;u,j,φ,kwargs...)
  u,j,p,φ,f,g
end

function _trasient_solution(::Val{:nonlineartime_fespace};kwargs...)
  u = TimeSpaceFunction(t->x-> VectorValue( x[2]*exp(-t), x[1]*cos(t),0) )
  j = TimeSpaceFunction(t->x-> VectorValue( sin(t), cos(t),0) )
  p = TimeSpaceFunction(t->x-> 0 )
  φ = TimeSpaceFunction(t->x-> cos(t)*one(x[1]) )
  f = _transient_solution_f(;u,j,p,φ,kwargs...)
  g = _transient_solution_fj(;u,j,φ,kwargs...)
  u,j,p,φ,f,g
end

function _trasient_solution(::Val{:stationary_nonfespace};kwargs...)
  u = TimeSpaceFunction(t->x-> VectorValue( cos(x[2]), cos(x[1]), 0 ) )
  j = TimeSpaceFunction(t->x-> VectorValue( cos(x[2]), 1, 0 ) )
  p = TimeSpaceFunction(t->x-> 0 )
  φ = TimeSpaceFunction(t->x-> cos(x[1]) )
  f = _transient_solution_f(;u,j,p,φ,kwargs...)
  g = _transient_solution_fj(;u,j,φ,kwargs...)
  u,j,p,φ,f,g
end

function _trasient_solution(::Val{:lineartime_nonfespace};kwargs...)
  u = TimeSpaceFunction(
    t->x-> VectorValue( cos(x[2])*(1-t/2), cos(x[1])*(1-t/2), 0 ) )
  j = TimeSpaceFunction(
    t->x-> VectorValue( cos(x[2])*t, (1-t/2), 0 ) )
  p = TimeSpaceFunction(t->x-> 0 )
  φ = TimeSpaceFunction(t->x-> cos(x[1]) )
  f = _transient_solution_f(;u,j,p,φ,kwargs...)
  g = _transient_solution_fj(;u,j,φ,kwargs...)
  u,j,p,φ,f,g
end

function _trasient_solution(::Val{:nonlineartime_nonfespace};kwargs...)
  u = TimeSpaceFunction(
    t->x-> VectorValue( cos(x[2])*exp(-t), cos(x[1])*cos(t), 0 ) )
  j = TimeSpaceFunction(
    t->x-> VectorValue( cos(x[2])*sin(t), cos(t), 0 ) )
  p = TimeSpaceFunction(t->x-> 0 )
  φ = TimeSpaceFunction(t->x-> cos(t)*cos(x[1]) )
  f = _transient_solution_f(;u,j,p,φ,kwargs...)
  g = _transient_solution_fj(;u,j,φ,kwargs...)
  u,j,p,φ,f,g
end

function _transient_solution_j(;u,φ,B,σ=1,kwargs...)
  jt = t-> x-> σ * (u(t,x)×B) - σ * ∇(φ)(t,x)
  TimeSpaceFunction(jt)
end

function _transient_solution_fj(;u,j,φ,B,σ=1,kwargs...)
  f = t->x-> j(t,x) + σ*∇(φ)(t,x) - σ*(u(t,x)×B)
  TimeSpaceFunction(f)
end

function _transient_solution_f(;u,j,p,B,α,β,γ,kwargs...)
  ft = t-> x->
    ∂t(u)(t,x) +
    α * u(t,x) ⋅ ∇(u)(t,x) +
    - β * Δ(u)(t,x) +
    ∇(p)(t,x) +
    - γ*( j(t,x) × B )
  TimeSpaceFunction(ft)
end

function _box_model(parts,np,domain,nc)
  model = CartesianDiscreteModel(parts,np,domain,nc)
  D = length(nc)
  d = D-1
  p = Polytope(Gridap.Helpers.tfill(HEX_AXIS,Val{D}()))
  facets_w = collect(2*(D-d)+1:2*(D))
  faces_w = facets_w .+ get_offset(p,D-1)
  facets_tb = setdiff(1:2*D,facets_w)
  faces_tb = facets_tb .+ get_offset(p,D-1)
  map(local_views(model)) do model
    labels = get_face_labeling(model)
    add_tag!(labels,"walls",faces_w)
    add_tag!(labels,"top_bottom",faces_tb)
  end
  model
end

function _print_timestep_errors(results,i;verbose=true)
  if verbose
    println("Time: $(results[:t][i])")
    println("  L²(uₕ-u): $(results[:uh_el2][i])")
    println("  L²(jₕ-j): $(results[:jh_el2][i])")
    println("  L²(pₕ-p): $(results[:ph_el2][i])")
    println("  L²(φₕ-φ): $(results[:φh_el2][i])")
    println("  H¹(uₕ-u): $(results[:uh_eh1][i])")
  end
end
