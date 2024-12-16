function toy(;
  backend = nothing,
  np      = nothing,
  title   = "toy",
  nruns   = 1,
  path    = datadir(),
  kwargs...
)

  for ir in 1:nruns
    _title = title*"_r$ir"
    if isa(backend,Nothing)
      @assert isa(np,Nothing)
      info, t = _toy(;title=_title,path=path,kwargs...)
    else
      @assert backend ∈ [:sequential,:mpi]
      @assert !isa(np,Nothing)
      if backend === :sequential
        info, t = with_debug() do distribute
          _toy(;distribute=distribute,rank_partition=np,title=_title,path=path,kwargs...)
        end
      else
        info,t = with_mpi() do distribute
          _toy(;distribute=distribute,rank_partition=np,title=_title,path=path,kwargs...)
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
      save(joinpath(path,"$_title.bson"),info)
    end
  end

  nothing
end

function _toy(;
  distribute=nothing,
  rank_partition=nothing,
  nc = (4,4),
  L  = 1.0,
  u0 = 1.0,
  ν  = 1.0,
  ρ  = 1.0,
  σ  = 1.0,
  B0 = 10.0,
  fz = 1.0,
  ζ  = 0.0, # Augmented Lagrangian weight  
  μ  = 0,   # Stabilization weight
  order = 2,
  order_j = order,
  nsums = 100,
  vtk=true,
  title = "test",
  path  = datadir(),
  debug = false,
  res_assemble = false,
  jac_assemble = false,
  solve = true,
  solver = :julia,
  formulation = :mhd,
  initial_value = :zero,
  rt_scaling = false,
  verbose = true,
  ranks_per_level = nothing,
  adaptivity_method = 0,
  fluid_disc = ifelse(iszero(adaptivity_method),:Qk_dPkm1,:SV),
  current_disc = :RT,
)
  @assert formulation ∈ [:cfd,:mhd]
  @assert initial_value ∈ [:zero,:solve]

  info = Dict{Symbol,Any}()
  params = Dict{Symbol,Any}(
    :debug=>debug,
    :solve=>solve,
    :res_assemble=>res_assemble,
    :jac_assemble=>jac_assemble,
  )

  # Communicator
  if isa(distribute,Nothing)
    @assert isa(rank_partition,Nothing)
    rank_partition = Tuple(fill(1,length(nc)))
    distribute = DebugArray
  end
  @assert length(rank_partition) == length(nc)
  parts = distribute(LinearIndices((prod(rank_partition),)))

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
  N  = Ha^2/Re
  f  = (L/(ρ*u0^2))*VectorValue(0.0,0.0,fz)
  B  = VectorValue(0.0,1.0,0.0)

  if formulation == :cfd # Option 1 (CFD)
    α = 1.0
    β = 1.0/Re
    γ = N
    ξ = 1.0
  elseif formulation == :mhd # Option 2 (MHD) is chosen in the experimental article
    α = (1.0/N)
    β = (1.0/Ha^2)
    γ = 1.0
    f = f/N
    ξ = 1.0/B0
  else
    error("Unknown formulation")
  end

  u(x) = analytical_hunt_u(L,L,ρ*ν,-fz/ρ,Ha,nsums,x)
  j(x) = ξ * analytical_hunt_j(L,L,σ,ρ*ν,-fz/ρ,Ha,nsums,x)

  # DiscreteModel in terms of reduced quantities

  model = toy_mesh(
    parts,params,nc,rank_partition;ranks_per_level,adaptivity_method
  )
  Ω = Interior(model)
  if debug && vtk
    writevtk(model,"hunt_model")
  end

  params[:fluid] = Dict{Symbol,Any}(
    :domain=>nothing,
    :α=>α,
    :β=>β,
    :γ=>γ,
    :f=>f,
    :B=>B,
    :ζ=>ζ,
  )

  params[:fespaces] = Dict{Symbol,Any}(
    :order_u => order,
    :order_j => order_j,
    :rt_scaling => rt_scaling ? 1.0/get_mesh_size(model) : nothing,
    :fluid_disc => fluid_disc,
    :current_disc => current_disc,
  )

  # Boundary conditions

  params[:bcs] = Dict(
    :u => Dict(:tags=>"boundary", :values => u),
    :j => Dict(:tags=>"boundary", :values => j),
  )
  params[:x0] = initial_value
  if μ > 0
    params[:bcs][:stabilization] = Dict(:μ=>μ)
  end

  toc!(t,"pre_process")

  # Solve it
  if !uses_petsc(params[:solver])
    xh,fullparams,info = main(params;output=info)
  else
    petsc_options = params[:solver][:petsc_options]
    xh,fullparams,info = GridapPETSc.with(args=split(petsc_options)) do
      xh,fullparams,info = main(params;output=info)
      GridapPETSc.gridap_petsc_gc() # Destroy all PETSc objects
      return xh,fullparams,info
    end
  end
  t = fullparams[:ptimer]
  
  # Rescale quantities

  tic!(t,barrier=true)
  ūh,p̄h,j̄h,φ̄h = xh
  uh = u0*ūh
  ph = (ρ*u0^2)*p̄h
  jh = (σ*u0*B0)*j̄h
  φh = (L*u0*B0)*φ̄h

  div_jh = ∇·jh
  div_uh = ∇·uh

  # Post process

  dΩ = Measure(Ω,4*(order+1))
  eu = u - uh
  ej = j - jh
  eu_h1 = sqrt(sum(∫( ∇(eu)⊙∇(eu) + eu⋅eu  )dΩ))
  eu_l2 = sqrt(sum(∫( eu⋅eu )dΩ))
  ej_l2 = sqrt(sum(∫( ej⋅ej )dΩ))
  uh_l2 = sqrt(sum(∫( uh⋅uh )dΩ))
  uh_h1 = sqrt(sum(∫( ∇(uh)⊙∇(uh) + uh⋅uh )dΩ))
  jh_l2 = sqrt(sum(∫( jh⋅jh )dΩ))
  toc!(t,"post_process")

  if vtk
    cellfields = [
      "uh" => uh,"ph" => ph,"jh" => jh,"phi" => φh,
      "div_jh" => div_jh,"div_uh" => div_uh,
      "u" => u, "j" => j
    ]
    writevtk(Ω,joinpath(path,title),order=max(order,order_j),cellfields=cellfields,append=false)
    toc!(t,"vtk")
  end
  if verbose 
    if i_am_main(parts)
      println(" >> H1 error for u = ", eu_h1)
      println(" >> L2 error for j = ", ej_l2)
    end
    display(t)
  end

  info[:nc] = nc
  info[:nsums] = nsums
  info[:ncells] = num_cells(model)
  info[:ndofs_u] = length(get_free_dof_values(ūh))
  info[:ndofs_p] = length(get_free_dof_values(p̄h))
  info[:ndofs_j] = length(get_free_dof_values(j̄h))
  info[:ndofs_φ] = length(get_free_dof_values(φ̄h))
  info[:ndofs] = length(get_free_dof_values(xh))
  info[:Re] = Re
  info[:Ha] = Ha
  info[:eu_h1] = eu_h1
  info[:eu_l2] = eu_l2
  info[:ej_l2] = ej_l2
  info[:uh_h1] = uh_h1
  info[:uh_l2] = uh_l2
  info[:jh_l2] = jh_l2

  info, t
end

function toy_mesh(
  parts,params,nc::NTuple{2,Int},np::NTuple{2,Int};
  ranks_per_level = nothing, adaptivity_method = 0
)
  Lz = 5.0*2.0/max(nc[1],nc[2]) # Preserve good aspect ratio
  _nc = (nc[1],nc[2],5)
  _np = (np[1],np[2],1)
  domain = (-1.0,1.0,-1.0,1.0,0.0,Lz)
  if isnothing(ranks_per_level) # Single grid
    model = CartesianDiscreteModel(parts,_np,domain,_nc)
    add_cavity_tags!(model)
  else # Multigrid
    mh = CartesianModelHierarchy(
      parts,np_per_level,domain,_nc;
      add_labels! = add_cavity_tags!,
    )
    params[:multigrid] = Dict{Symbol,Any}(
      :mh => mh,
      :num_refs_coarse => 0,
      :ranks_per_level => ranks_per_level,
    )
    model = get_model(mh,1)
  end
  model = Meshers.adapt_mesh(model,adaptivity_method)
  params[:model] = model
  return model
end
