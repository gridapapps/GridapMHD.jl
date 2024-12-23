function expansion(;
  backend = nothing,
  np      = nothing,
  title   = "Expansion",
  mesh    = "710",
  path    = datadir(),
  kwargs...)

  if isa(backend,Nothing)
    @assert isa(np,Nothing)
    info, t = _expansion(;title=title,path=path,mesh=mesh,kwargs...)
  else
    @assert backend ∈ [:sequential,:mpi]
    @assert !isnothing(np)
    np = isa(np,Int) ? (np,) : np
    if backend === :sequential
      info,t = with_debug() do distribute
        _expansion(;distribute=distribute,rank_partition=np,title=title,path=path,mesh=mesh,kwargs...)
      end
    else
      info,t = with_mpi() do distribute
        _expansion(;distribute=distribute,rank_partition=np,title=title,path=path,mesh=mesh,kwargs...)
      end
    end
  end
  info[:np] = np
  info[:backend] = backend
  info[:title] = title
  info[:mesh] = mesh
  map_main(t.data) do data
    for (k,v) in data
      info[Symbol("time_$k")] = v.max
    end
    save(joinpath(path,"$title.bson"),info)
  end

  nothing
end

function _expansion(;
  distribute=nothing,
  rank_partition=nothing,
  title = "Expansion",
  mesh  = "710",
  vtk   = true,
  path  = datadir(),
  debug = false,
  verbose = true,
  solver  = :julia,
  formulation = :mhd,
  Z = 4.0,    # Expansion Ratio, it has to be consistent with the mesh
  b = 0.2,    # Outlet channel aspect ratio, it has to be consistent with the mesh
  N  = 1.0,
  Ha = 1.0,
  cw = 0.028,
  τ  = 100,
  ζ  = 0.0,   # Augmented Lagrangian weight  
  μ  = 0,     # Stabilization weight
  order = 2,
  order_j = order,
  inlet = :parabolic,
  initial_value = :zero,
  solid_coupling = :none,
  niter = nothing,
  convection = :newton,
  fluid_disc = :Qk_dPkm1,
  current_disc = :RT,
  rt_scaling = false,
  savelines = false,
  petsc_options = "",
)
  @assert solid_coupling ∈ [:none,:thin_wall,:solid]
  @assert inlet ∈ [:parabolic,:shercliff,:constant]
  @assert initial_value ∈ [:zero,:inlet,:solve]
  @assert formulation ∈ [:cfd,:mhd]
  @assert convection ∈ [:newton,:picard,:none]

  info   = Dict{Symbol,Any}()
  params = Dict{Symbol,Any}(
    :debug=>debug,
    :solve=>true,
    :res_assemble=>false,
    :jac_assemble=>false,
    :solver=> isa(solver,Symbol) ? default_solver_params(Val(solver)) : solver
  )

  if isa(distribute,Nothing)
    @assert isa(rank_partition,Nothing)
    rank_partition = (1,)
    distribute = DebugArray
  end
  parts = distribute(LinearIndices((prod(rank_partition),)))

  t = PTimer(parts,verbose=verbose)
  params[:ptimer] = t
  tic!(t,barrier=true)

  # Mesh
  if isa(mesh,String)
    mesh = Dict(:mesher => :gmsh, :base_mesh => mesh)
  end
  model = expansion_mesh(mesh,parts,params)
  if debug && vtk
    writevtk(model,path*"/expansion_model")
  end
  toc!(t,"model")

  # Parameters and bounday conditions
  # Fluid parameters

  Re = Ha^2/N
  if formulation == :cfd # Option 1 (CFD)
    α = 1.0
    β = 1.0/Re
    γ = N
  elseif formulation == :mhd # Option 2 (MHD) is chosen in the experimental article
    α = (1.0/N)
    β = (1.0/Ha^2)
    γ = 1.0
  else
    error("Unknown formulation")
  end

  params[:fespaces] = Dict{Symbol,Any}(
    :order_u => order,
    :order_j => order_j,
    :rt_scaling => rt_scaling ? 1.0/get_mesh_size(model) : nothing,
    :fluid_disc => fluid_disc,
    :current_disc => current_disc,
  )

  params[:fluid] = Dict{Symbol,Any}(
    :domain => "fluid",
    :α => α,
    :β => β,
    :γ => γ,
    :f => VectorValue(0.0,0.0,0.0),
    :B => VectorValue(0.0,1.0,0.0),
    :ζ => ζ,
    :convection => convection,
  )

  # Boundary conditions

  u_in = u_inlet(inlet,Ha,Z,b)

  params[:bcs] = Dict{Symbol,Any}()
  params[:bcs][:u] = Dict(
    :tags => ["inlet", "wall"],
    :values => [u_in, VectorValue(0.0, 0.0, 0.0)]
  )

  if solid_coupling == :none
    params[:bcs][:j] = Dict{Symbol,Any}(
      :tags => ["wall", "inlet", "outlet"],
      :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)],
    )
  elseif solid_coupling == :thin_wall
    params[:bcs][:j] = Dict{Symbol,Any}(
      :tags => ["inlet", "outlet"],
      :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)]
    )
    params[:bcs][:thin_wall] = Dict{Symbol,Any}(
      :τ=>τ,
      :cw=>cw,
      :domain => ["wall"],
    )
  else
    @assert solid_coupling == :solid
    params[:bcs][:j] = Dict(
      :tags => ["wall_exterior", "inlet", "outlet"],
      :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)]
    )
    params[:solid] = Dict(:domain => "wall",:σ => 1.0)
  end

  if μ > 0
    params[:bcs][:stabilization] = Dict(:μ=>μ)
  end

  # Initial conditions

  j_zero = VectorValue(0.0,0.0,0.0)
  if initial_value == :inlet
    params[:x0] = Dict(
      :u=>u_in,:j=>j_zero,:p=>0.0,:φ=>0.0
    )
  else
    params[:x0] = initial_value
  end

  # Solver options

  if params[:solver][:solver] == :petsc && petsc_options != ""
    params[:solver][:petsc_options] = petsc_options
  end
  if !isnothing(niter)
    params[:solver][:niter] = niter
    if params[:solver][:solver] == :petsc
      params[:solver][:petsc_options] *= " -snes_max_funcs $(niter+1)"
    end
  end

  toc!(t,"pre_process")

  # Main solve
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
  tic!(t,barrier=true)

  # Post-process
  uh,ph,jh,φh = xh
  div_jh = ∇·jh
  div_uh = ∇·uh
  Grad_p = ∇·ph

  if vtk
    Ω = Interior(model,tags="fluid")
    iorder = min(3,max(order,order_j))
    writevtk(
      Ω,joinpath(path,title),
      order=iorder,
      cellfields=[
        "uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh,"div_uh"=>div_uh,"div_jh"=>div_jh,"kp"=>Grad_p
      ],
      append=false
    )
    toc!(t,"vtk")
  end
  if savelines
    line = top_ha_line(model,Z)
    info[:line] = line
    info[:p_line] = evaluate_line(ph,line)
    toc!(t,"p-lines")
  end
  if verbose
    display(t)
  end

  info[:ncells] = num_cells(model)
  info[:ndofs_u] = length(get_free_dof_values(uh))
  info[:ndofs_p] = length(get_free_dof_values(ph))
  info[:ndofs_j] = length(get_free_dof_values(jh))
  info[:ndofs_φ] = length(get_free_dof_values(φh))
  info[:ndofs] = length(get_free_dof_values(xh))
  info[:Ha] = Ha
  info[:N]  = N
  info[:Re] = Ha^2/N
  return info, t
end

# Mesh

const meshes_dir = joinpath(@__DIR__,"..","..","meshes")

function expansion_mesh(mesh::Dict,ranks,params)
  expansion_mesh(Val(mesh[:mesher]),mesh,ranks,params)
end

function expansion_mesh(::Val{:gmsh},mesh::Dict,ranks,params)
  # The domain is of size L_out x 2 x 2/β and L_in x 2/Z x 2/β
  # after and before the expansion respectively.
  msh_file = mesh[:base_mesh]
  if !ispath(msh_file)
    msh_file = joinpath(meshes_dir,"Expansion_"*msh_file*".msh") |> normpath
  end
  model = UnstructuredDiscreteModel(GmshDiscreteModel(ranks,msh_file))
  setup_expansion_mesh_tags!(model)
  params[:model] = model
  return model
end

function expansion_mesh(::Val{:gridap_SG},mesh::Dict,ranks,params)
  @assert haskey(mesh,:num_refs)
  num_refs = mesh[:num_refs]
  if haskey(mesh,:base_mesh)
    msh_file = mesh[:base_mesh]
    if !ispath(msh_file)
      msh_file = joinpath(meshes_dir,"Expansion_"*msh_file*".msh") |> normpath
    end
    base_model = UnstructuredDiscreteModel(GmshDiscreteModel(ranks,msh_file))
  else
    base_model = expansion_generate_base_mesh()
  end
  setup_expansion_mesh_tags!(base_model)
  model = Meshers.generate_refined_mesh(ranks,base_model,num_refs)
  if haskey(mesh,:adaptivity_method)
    model = Meshers.adapt_mesh(model,mesh[:adaptivity_method])
  end
  params[:model] = model
  return model
end

function expansion_mesh(::Val{:p4est_SG},mesh::Dict,ranks,params)
  @assert haskey(mesh,:num_refs)
  num_refs = mesh[:num_refs]
  if haskey(mesh,:base_mesh)
    msh_file = mesh[:base_mesh]
    if !ispath(msh_file)
      msh_file = joinpath(meshes_dir,"Expansion_"*msh_file*".msh") |> normpath
    end
    base_model = UnstructuredDiscreteModel(GmshDiscreteModel(msh_file))
  else
    base_model = expansion_generate_base_mesh()
  end
  setup_expansion_mesh_tags!(base_model)
  model = Meshers.generate_p4est_refined_mesh(ranks,base_model,num_refs)
  if haskey(mesh,:adaptivity_method)
    model = Meshers.adapt_mesh(model,mesh[:adaptivity_method])
  end
  params[:model] = model
  return model
end

function expansion_mesh(::Val{:p4est_MG},mesh::Dict,ranks,params)
  @assert haskey(mesh,:num_refs_coarse) && haskey(mesh,:ranks_per_level)
  num_refs_coarse = mesh[:num_refs_coarse]
  ranks_per_level = mesh[:ranks_per_level]
  if haskey(mesh,:base_mesh)
    msh_file = mesh[:base_mesh]
    if !ispath(msh_file)
      msh_file = joinpath(meshes_dir,"Expansion_"*msh_file*".msh") |> normpath
    end
    base_model = UnstructuredDiscreteModel(GmshDiscreteModel(msh_file))
  else
    base_model = expansion_generate_base_mesh()
  end
  setup_expansion_mesh_tags!(base_model)

  mh = Meshers.generate_p4est_mesh_hierarchy(ranks,base_model,num_refs_coarse,ranks_per_level)
  params[:multigrid] = Dict{Symbol,Any}(
    :mh => mh,
    :num_refs_coarse => num_refs_coarse,
    :ranks_per_level => ranks_per_level,
  )

  model = get_model(mh,1)
  params[:model] = model
  return model
end

function setup_expansion_mesh_tags!(model::GridapDistributed.DistributedDiscreteModel)
  map(setup_expansion_mesh_tags!,local_views(model))
end

function setup_expansion_mesh_tags!(model::DiscreteModel)
  labels = get_face_labeling(model)

  tags = labels.tag_to_name
  solid = issubset(["wall_interior","wall_exterior","wall_volume"], tags)

  if "wall" ∉ tags
    @assert solid
    add_tag_from_tags!(labels,"wall",["wall_interior","wall_exterior","wall_volume"])
  end
  if "boundary" ∉ tags
    if solid
      add_tag_from_tags!(labels,"boundary",["inlet","outlet","wall_exterior"])
    else
      add_tag_from_tags!(labels,"boundary",["inlet","outlet","wall"])
    end
  end
  if "interior" ∉ tags
    if solid
      add_tag_from_tags!(labels,"interior",["fluid","wall_volume","wall_interior"])
    else
      add_tag_from_tags!(labels,"interior",["fluid"])
    end
  end
end

# It ensures avg(u) = 1 in the outlet channel in every case
function u_inlet(inlet,Ha,Z,β)
  u_inlet_parabolic((x,y,z)) = VectorValue(36.0*Z*(y-1/Z)*(y+1/Z)*(z-β*Z)*(z+β*Z),0,0)

  kp_inlet = GridapMHD.kp_shercliff_cartesian(β*Z,Ha/Z)
  function u_inlet_shercliff((x,y,z))
    x,y,z = (z*Z,y*Z,x)
    u_x = GridapMHD.analytical_GeneralHunt_u(β*Z, 0.0, -Z*kp_inlet, Ha/Z,200,(x,y,z))[3]
    if abs(x) > β*Z || abs(y) > 1.0
      u_x = 0.0
    end
    VectorValue(u_x,0.0,0.0)
  end

  u_inlet_cte = VectorValue(Z,0.0,0.0)

  if inlet == :parabolic
	  U = u_inlet_parabolic
  elseif inlet == :shercliff
    Ha > 10 || @warn "Shercliff inlet is not accurate for Ha=$Ha"
    U = u_inlet_shercliff
  else
    U = u_inlet_cte
  end

  return U
end

function evaluate_line(uh,line)
  model = get_background_model(get_triangulation(uh))
  # Simplexified model
  smodel = Gridap.Adaptivity.refine(model;refinement_method="simplexify")
  strian = Triangulation(smodel)
  # Distributed change of domain
  uhi = GridapDistributed.DistributedCellField(
    map((uhi,ti) -> change_domain(uhi,ti,ReferenceDomain()),local_views(uh),local_views(strian)),
    strian
  )
  return uhi(line)
end

function top_line(model,n=100)
  pmin,pmax = _get_bounding_box(model)
  xmin,xmax = pmin[1],pmax[1]
  zmax = pmax[3]
  line = map(x->Point(x,0.0,zmax), range(xmin,xmax,n+1))
  return line
end

function top_ha_line(model,Z,n=100)
  δ = 1.e6*eps(Float64)
  pmin,pmax = _get_bounding_box(model)
  xmin,xmax = pmin[1]+δ,pmax[1]-δ
  y1,y2 = pmax[2]-δ, pmax[2]/Z
  line = map( x -> x>0 ? Point(x,y1,δ) : Point(x,y2,δ), range(xmin,xmax,n+1))
  return line
end

function evaluation_lines(model,Z,n=100)
  pmin,pmax = _get_bounding_box(model)
  xmin,xmax = pmin[1],pmax[1]
  ymin,ymax = pmin[2]/Z,pmax[2]/Z
  zmin,zmax = pmin[3],pmax[3]
  pmin,pmax = Point(xmin,ymin,zmin),Point(xmax,ymax,zmax)

  xline = map(x->Point(x,0.0,0.0), range(xmin,xmax,n+1))
  yline = map(y->Point(xmin,y,0.0), range(ymin,ymax,n+1))
  zline = map(z->Point(xmin,0.0,z), range(zmin,zmax,n+1))
  xline,yline,zline
end

function _get_bounding_box(x::AbstractVector{<:Point{D}}) where D
  pmin = fill(Inf,D)
  pmax = fill(-Inf,D)
  for p in x
    for d in 1:D
      pmin[d] = min(pmin[d],p[d])
      pmax[d] = max(pmax[d],p[d])
    end
  end
  pmin, pmax
end

function _get_bounding_box(model::DiscreteModel)
  _get_bounding_box(get_node_coordinates(model))
end

function _get_bounding_box(
  model::GridapDistributed.DistributedDiscreteModel{Dc}
) where Dc
  _pmin, _pmax = map(_get_bounding_box,local_views(model)) |> tuple_of_arrays
  
  pmin = fill(Inf,Dc)
  pmax = fill(-Inf,Dc)
  for d in 1:Dc
    d_min = map(x -> x[d], _pmin)
    d_max = map(x -> x[d], _pmax)
    pmin[d] = PartitionedArrays.getany(reduction(min,d_min;destination=:all))
    pmax[d] = PartitionedArrays.getany(reduction(max,d_max;destination=:all))
  end
  
  return pmin, pmax
end
