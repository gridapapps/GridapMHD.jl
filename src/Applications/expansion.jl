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
  mesh  = "720",            #Either a string (for pure gmsh) or, a dictionary for p4est?
  vtk   = true,
  path  = datadir(),
  debug = false,
  verbose = true,
  solver  = :julia,
  petsc_options = nothing,  #Added to have access to non-default petsc options
  formulation = :cfd,
  Z = 4.0,                  #Expansion Ratio, it has to be consistent with the mesh
  β = 0.2,                  #Outlet channel aspect ratio, , it has to be consistent with the mesh
  N  = 1.0,
  Ha = 1.0,
  cw = 0.028,
  τ  = 100,
  ζ  = 0.0, #What is this?
  order = 2,
  inlet = :parabolic
  )
  
  info   = Dict{Symbol,Any}()
  params = Dict{Symbol,Any}(
       :debug=>debug,
       :solve=>true,
       :res_assemble=>false,
       :jac_assemble=>false,
       :solver=> isa(solver,Symbol) ? default_solver_params(Val(solver)) : solver
    )
    
  if isa(petsc_options,Nothing)
    params[:solver][:petsc_options] = petsc_options
  end

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
  )

  # Boundary conditions
  # u_inlet((x,y,z)) = VectorValue(36.0*(y-1/4)*(y+1/4)*(z-1)*(z+1),0,0) # This gives mean(u_inlet)=1
  if abs(cw) < 1.e-5
    params[:bcs] = Dict( 
      :u => Dict(
        :tags => ["inlet", "wall"],
        :values => [u_inlet(inlet,Ha,Z,β), VectorValue(0.0, 0.0, 0.0)]
      ),
      :j => Dict(
		    :tags => ["wall", "inlet", "outlet"], 
        :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)],
      ) #Pending: Add J analyical formula in the inlet BC when Shercliff BC is selected at the inlet
    )
  else 
    params[:bcs] = Dict(
      :u => Dict(
        :tags => ["inlet", "wall"],
        :values => [u_inlet(inlet,Ha,Z,β), VectorValue(0.0, 0.0, 0.0)]
      ),
      :j => Dict(
        :tags => ["inlet", "outlet"], 
        :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)]
      ),
      :thin_wall => [Dict(
        :τ=>τ,
        :cw=>cw,
        :domain => ["wall"], #No necessary to specify that "wall" is a gridap boundary
      )]
    )
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
    writevtk(Ω,joinpath(path,title),
      order=2,
      cellfields=[
        "uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh,"div_uh"=>div_uh,"div_jh"=>div_jh,"kp"=>Grad_p])
    toc!(t,"vtk")
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

function expansion_mesh(mesh::Dict,ranks,params)
  expansion_mesh(Val(mesh[:mesher]),mesh,ranks,params)
end

function expansion_mesh(::Val{:gmsh},mesh::Dict,ranks,params)
  # The domain is of size L_out x 2 x 2/β and L_in x 2/Z x 2/β
  # after and before the expansion respectively.
  msh_name = mesh[:base_mesh]
  msh_file = joinpath(projectdir(),"meshes","Expansion_"*msh_name*".msh") |> normpath
  model = GmshDiscreteModel(ranks,msh_file)
  params[:model] = model
  return model
end

function epansion_mesh(::Val{:p4est_SG},mesh::Dict,ranks,params)
  @assert haskey(mesh,:num_refs)
  num_refs = mesh[:num_refs]
  if haskey(mesh,:base_mesh)
    msh_file = joinpath(projectdir(),"meshes","Expansion_"*mesh[:base_mesh]*".msh") |> normpath
    base_model = GmshDiscreteModel(msh_file)
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
    msh_file = joinpath(projectdir(),"meshes","Expansion_"*mesh[:base_mesh]*".msh") |> normpath
    base_model = GmshDiscreteModel(msh_file)
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
  u_inlet_shercliff((x,y,z)) = VectorValue(GridapMHD.analytical_GeneralHunt_u(β*Z, 0.0, -Z*kp_inlet, Ha/Z,200,(z*Z,y*Z,x))[3],0.0,0.0) 

  u_inlet_cte = VectorValue(Z,0.0,0.0)  

  if inlet == :parabolic
	U = u_inlet_parabolic
  elseif inlet == :shercliff
    U = u_inlet_shercliff
  else
    U = u_inlet_cte
  end
 U
end
