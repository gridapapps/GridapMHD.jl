
function cavity(;
  backend = nothing,
  np      = 1,
  title   = "Cavity",
  path    = datadir(),
  kwargs...)

  if isa(backend,Nothing)
    @assert np == 1
    info, t = _cavity(;title=title,path=path,kwargs...)
  else
    @assert backend ∈ [:sequential,:mpi]
    if backend === :sequential
      info,t = with_debug() do distribute
        _cavity(;distribute=distribute,np=np,title=title,path=path,kwargs...)
      end
    else
      info,t = with_mpi() do distribute
        _cavity(;distribute=distribute,np=np,title=title,path=path,kwargs...)
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

function _cavity(;
  distribute=nothing,
  np=1,
  nc=(4,4,4),
  ν = 1.0,
  ρ = 1.0,
  σ = 1.0,
  ζ = 0.0,   # Augmented Lagrangian weight  
  μ = 0,     # Stabilization weight
  B = VectorValue(0.0, 0.0, 10.0),
  f = VectorValue(0.0, 0.0, 0.0),
  L = 1.0,
  u0 = 1.0,
  B0 = norm(B),
  order = 2,
  order_j = order,
  formulation = :mhd,
  initial_value = :zero,
  rt_scaling = false,
  title="Cavity",
  path=datadir(),
  solver=:julia,
  ranks_per_level=nothing,
  verbose=true,
  vtk=true,
  convection=:newton,
  closed_cavity=true
)
  @assert formulation ∈ [:cfd,:mhd]
  @assert initial_value ∈ [:zero,:solve]
  @assert convection ∈ [:newton,:picard,:none]

  info = Dict{Symbol,Any}()
  params = Dict{Symbol,Any}(
    :debug => false,
    :solve => true,
    :res_assemble => false,
    :jac_assemble => false,
    :check_valid => false,
  )

  if isa(distribute,Nothing)
    distribute = DebugArray
  end
  parts = distribute(LinearIndices((prod(np),)))
  
  t = PTimer(parts,verbose=verbose)
  params[:ptimer] = t
  tic!(t, barrier=true)

  # Solver
  if isa(solver,Symbol)
    solver = default_solver_params(Val(solver))
  end
  params[:solver] = solver

  # Model
  model = cavity_mesh(parts,params,nc,np,L,ranks_per_level)

  # Reduced quantities
  Re = u0 * L / ν
  Ha = B0 * L * sqrt(σ / (ρ * ν))
  N = Ha^2 / Re
  f̄ = (L / (ρ * u0^2)) * f
  B̄ = (1 / B0) * B

  if formulation == :cfd # Option 1 (CFD)
    α = 1.0
    β = 1.0/Re
    γ = N
  elseif formulation == :mhd # Option 2 (MHD) is chosen in the experimental article
    α = (1.0/N)
    β = (1.0/Ha^2)
    γ = 1.0
    f̄ = f̄ / N
  else
    error("Unknown formulation")
  end

  params[:fluid] = Dict(
    :domain => nothing,
    :α => α,
    :β => β,
    :γ => γ,
    :f => f̄,
    :B => B̄,
    :ζ => ζ,
    :convection => convection,
  )

  # FESpaces and Boundary conditions
  uw = VectorValue(0.0, 0.0, 0.0)
  ul = VectorValue(1.0, 0.0, 0.0)
  ji = VectorValue(0.0, 0.0, 0.0)
  params[:fespaces] = Dict{Symbol,Any}(
    :order_u => order,
    :order_j => order_j,
    :rt_scaling => rt_scaling ? 1.0/get_mesh_size(model) : nothing
  )

  params[:x0] = initial_value

  if closed_cavity
    params[:bcs] = Dict{Symbol,Any}(
      :u => Dict(:tags => ["cavity", "lid"], :values => [uw, ul]),
      :j => Dict(:tags => "insulating", :values => ji),
    )
    params[:fespaces][:p_constraint] = :zeromean
  else
    params[:bcs] = Dict{Symbol,Any}(
      :u => Dict(:tags => ["wall", "lid"], :values => [uw, ul]), # Bottom wall is Newman
      :j => Dict(:tags => "insulating", :values => ji),
    )
  end

  if μ > 0
    params[:bcs][:stabilization] = Dict(:μ=>μ)
  end

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

  if vtk
    tic!(t, barrier=true)
    Ω = Interior(model)
    ūh, p̄h, j̄h, φ̄h = xh
    uh = u0 * ūh
    ph = (ρ * u0^2) * p̄h
    jh = (σ * u0 * B0) * j̄h
    φh = (u0 * B0 * L) * φ̄h
    div_jh = ∇·jh
    div_uh = ∇·uh
    writevtk(
      Ω, joinpath(path,title), order=max(order,order_j), 
      cellfields=["uh" => uh, "ph" => ph, "jh" => jh, "phi" => φh, "div_jh" => div_jh, "div_uh" => div_uh],
      append=false
    )
    toc!(t, "vtk")
  end

  info[:ncells]  = num_cells(model)
  info[:ndofs_u] = length(get_free_dof_values(ūh))
  info[:ndofs_p] = length(get_free_dof_values(p̄h))
  info[:ndofs_j] = length(get_free_dof_values(j̄h))
  info[:ndofs_φ] = length(get_free_dof_values(φ̄h))
  info[:ndofs]   = sum([info[:ndofs_u], info[:ndofs_p], info[:ndofs_j], info[:ndofs_φ]])
  info[:Re]      = Re
  info[:Ha]      = Ha

  return info, t
end

function add_cavity_tags!(model::CartesianDiscreteModel)
  labels = get_face_labeling(model)
  add_cavity_tags!(labels)
end

function add_cavity_tags!(model::GridapDistributed.DistributedDiscreteModel)
  map(add_cavity_tags!, local_views(model))
end

function add_cavity_tags!(labels)
  Γl = [22]
  Γb = [21]
  Γw = vcat(collect(1:20),collect(23:26))
  add_tag_from_tags!(labels, "lid", Γl)
  add_tag_from_tags!(labels, "bottom", Γb)
  add_tag_from_tags!(labels, "wall", Γw)
  add_tag_from_tags!(labels, "cavity", ["wall","bottom"])
  add_tag_from_tags!(labels, "insulating", "boundary")
end

function cavity_mesh(parts,params,nc::Int,np,L,ranks_per_level)
  return cavity_mesh(parts,params,(nc,nc,nc),np,L,ranks_per_level)
end

function cavity_mesh(parts,params,nc::Tuple,np::Int,L,ranks_per_level)
  return cavity_mesh(parts,params,nc,(np,1,1),L,ranks_per_level)
end

function cavity_mesh(parts,params,nc::Tuple,np::Tuple,L,ranks_per_level)
  domain = (0.0,L,0.0,L,0.0,L)
  if isnothing(ranks_per_level) # Single grid
    model = CartesianDiscreteModel(parts,np,domain,nc)
    add_cavity_tags!(model)
    params[:model] = model
  else # Multigrid
    mh = CartesianModelHierarchy(
      parts,np_per_level,domain,nc;
      add_labels! = add_cavity_tags!,
    )
    params[:multigrid] = Dict{Symbol,Any}(
      :mh => mh,
      :num_refs_coarse => 0,
      :ranks_per_level => ranks_per_level,
    )
    model = get_model(mh,1)
    params[:model] = model
  end
  return model
end
