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
  Z = 4.0,                  #Expansion Ratio, it has to be consistent with the mesh
  b = 0.2,                  #Outlet channel aspect ratio, it has to be consistent with the mesh
  N  = 1.0,
  Ha = 1.0,
  cw = 0.028,
  τ  = 100,
  ζ  = 0.0,
  order = 2,
  inlet = :parabolic,
  μ=0,
  initial_value=:zero,
  niter=nothing,
  convection=:true,
  savelines=false,
  petsc_options="",
  )

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
    writevtk(model,"data/expansion_model")
  end
  Ω = Interior(model,tags="PbLi")
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

  params[:fespaces] = Dict(
    :k => order
  )

  params[:fluid] = Dict(
    :domain => nothing, # whole domain
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

  if abs(cw) < 1.e-5
    params[:bcs][:j] = Dict(
      :tags => ["wall", "inlet", "outlet"],
      :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)],
    )

  else
    params[:bcs][:j] = Dict(
      :tags => ["inlet", "outlet"],
      :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)]
    )
    params[:bcs][:thin_wall] = Dict(
      :τ=>τ,
      :cw=>cw,
      :domain => ["wall"],
    )
  end

  if μ > 0
    params[:bcs][:stabilization] = Dict(:μ=>μ)
  end

  j_zero = VectorValue(0.0,0.0,0.0)
  if initial_value == :inlet
    params[:solver][:initial_values] = Dict(
      :u=>u_in,:j=>j_zero,:p=>0.0,:φ=>0.0)
  end
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
    writevtk(
      Ω,joinpath(path,title),
      order=order,
      cellfields=[
        "uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh,"div_uh"=>div_uh,"div_jh"=>div_jh,"kp"=>Grad_p],
      append=false
    )
    toc!(t,"vtk")
  end
  if savelines
    line = top_ha_line(model,Z)
    info[:line] = line
    info[:p_on_top] = ph(line)
    # xline,yline,zline = evaluation_lines(model,Z)
    # info[:xline] = xline
    # info[:yline] = yline
    # info[:zline] = zline
    # info[:uh_xline] = vector_field_eval(uh,xline)
    # info[:uh_yline] = vector_field_eval(uh,yline)
    # info[:uh_zline] = vector_field_eval(uh,zline)
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
  msh_name = mesh[:base_mesh]
  msh_file = joinpath(meshes_dir,"Expansion_"*msh_name*".msh") |> normpath
  model = GmshDiscreteModel(ranks,msh_file;has_affine_map=true)
  params[:model] = model
  return model
end

function epansion_mesh(::Val{:p4est_SG},mesh::Dict,ranks,params)
  @assert haskey(mesh,:num_refs)
  num_refs = mesh[:num_refs]
  if haskey(mesh,:base_mesh)
    msh_file = joinpath(meshes_dir,"Expansion_"*mesh[:base_mesh]*".msh") |> normpath
    base_model = GmshDiscreteModel(msh_file;has_affine_map=true)
    add_tag_from_tags!(get_face_labeling(base_model),"interior",["PbLi"])
    add_tag_from_tags!(get_face_labeling(base_model),"boundary",["inlet","outlet","wall"])
  else
    base_model = expansion_generate_base_mesh()
  end
  model = Meshers.generate_refined_mesh(ranks,base_model,num_refs)
  params[:model] = model
  return model
end

function expansion_mesh(::Val{:p4est_MG},mesh::Dict,ranks,params)
  @assert haskey(mesh,:num_refs_coarse) && haskey(mesh,:ranks_per_level)
  num_refs_coarse = mesh[:num_refs_coarse]
  ranks_per_level = mesh[:ranks_per_level]
  if haskey(mesh,:base_mesh)
    msh_file = joinpath(meshes_dir,"Expansion_"*mesh[:base_mesh]*".msh") |> normpath
    base_model = GmshDiscreteModel(msh_file;has_affine_map=true)
    add_tag_from_tags!(get_face_labeling(base_model),"interior",["PbLi"])
    add_tag_from_tags!(get_face_labeling(base_model),"boundary",["inlet","outlet","wall"])
  else
    base_model = expansion_generate_base_mesh()
  end

  mh = Meshers.generate_mesh_hierarchy(ranks,base_model,num_refs_coarse,ranks_per_level)
  params[:multigrid] = Dict{Symbol,Any}(
    :mh => mh,
    :num_refs_coarse => num_refs_coarse,
    :ranks_per_level => ranks_per_level,
  )

  model = get_model(mh,1)
  params[:model] = model
  return model
end

function u_inlet(inlet,Ha,Z,β) # It ensures avg(u) = 1 in the outlet channel in every case

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
 U
end

function top_line(model,n=100)
  pmin,pmax = _get_bounding_box(model)
  xmin,xmax = pmin[1],pmax[1]
  zmax = pmax[3]
  line = map(x->Point(x,0.0,zmax), range(xmin,xmax,n+1))
  return line
end

function top_ha_line(model,Z,n=100)
  pmin,pmax = _get_bounding_box(model)
  xmin,xmax = pmin[1],pmax[1]
  ymax = pmax[2]
  # line = map( x -> x>0 ? Point(x,ymax,0.0) : Point(x,ymax/Z,0.0), range(xmin,xmax,n+1))
  # line = map( x -> x>0 ? Point(x,1.0,0.0) : Point(x,0.25,0.0), range(xmin,xmax,n+1))
  line = map( x -> x>0 ? Point(x,0.99,0.0) : Point(x,0.25,0.0), range(xmin,xmax,n+1))
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

function _get_bounding_box(x::AbstractVector{<:Point})
  pmin,pmax = x[1],x[1]
  for p in x
    pmin = min.(pmin,p)
    pmax = max.(pmax,p)
  end
  pmin,pmax
end

function _get_bounding_box(model::DiscreteModel)
  _get_bounding_box(get_node_coordinates(model))
end

function _get_bounding_box(model::GridapDistributed.DistributedDiscreteModel)
  boxes = map(_get_bounding_box,local_views(model))
  boxes = gather(boxes,destination=:all)
  boxes = map(boxes) do boxes
    if ! isempty(boxes)
      pmin,pmax = boxes[1]
      for (bmin,bmax) in boxes
        pmin = min.(pmin,bmin)
        pmax = max.(pmax,bmax)
      end
      pmin,pmax
    else
      @unreachable
    end
  end
  pmin,pmax = PartitionedArrays.getany(boxes)
  pmin,pmax
end

# CellField evaluations until comments in
# https://github.com/gridap/GridapDistributed.jl/pull/146 are fixed

function vector_field_eval(f,x)
  fx = scalar_field_eval(f,x,1)
  fy = scalar_field_eval(f,x,2)
  fz = scalar_field_eval(f,x,3)
  map(VectorValue,fx,fy,fz)
end

function scalar_field_eval(f,x,comp)
  fx = (x->x[comp])∘f
  fx = emit(fx(x))
  PartitionedArrays.getany(fx)
end

function scalar_field_eval(f,x)
  fx = f(x)
  fx = emit(fx(x))
  PartitionedArrays.getany(fx)
end
