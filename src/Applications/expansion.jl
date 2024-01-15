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
    @assert !isa(np,Nothing)
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
  mesh  = "720",
  vtk   = true,
  path  = datadir(),
  debug = false,
  verbose = true,
  solver  = :julia,
  formulation = :cfd,
  N  = 1.0,
  Ha = 1.0,
  cw = 0.028,
  τ  = 100,
  ζ  = 0.0,
  )
  
  info   = Dict{Symbol,Any}()
  params = Dict{Symbol,Any}(
    :debug=>debug,
    :solve=>true,
    :res_assemble=>false,
    :jac_assemble=>false,
    :solver=>solver,
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
  model = expansion_mesh(mesh,parts,params)
  if debug && vtk
    writevtk(model,"expansion_model")
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
  u_inlet((x,y,z)) = VectorValue(36.0*(y-1/4)*(y+1/4)*(z-1)*(z+1),0,0) # This gives mean(u_inlet)=1
  if cw == 0.0
    params[:bcs] = Dict( 
      :u => Dict(
        :tags => ["inlet", "wall"],
        :values => [u_inlet, VectorValue(0.0, 0.0, 0.0)]
      ),
      :j => Dict(
		    :tags => ["wall", "inlet", "outlet"], 
        :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)],
      )
    )
  else 
    params[:bcs] = Dict(
      :u => Dict(
        :tags => ["inlet", "wall"],
        :values => [u_inlet, VectorValue(0.0, 0.0, 0.0)]
      ),
      :j => Dict(
        :tags => ["inlet", "outlet"], 
        :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)]
      ),
      :thin_wall => [Dict(
        :τ=>τ,
        :cw=>cw,
        :domain => ["wall"],
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
  if vtk
    writevtk(Ω,joinpath(path,title),
      order=2,
      cellfields=[
        "uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh,])
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

function expansion_mesh(mesh::String,ranks,params)
  # The domain is of size 8L x 2L x 2L and 8L x 2L/Z x 2L
  # after and before the expansion respectively (L=1).
  msh_file = joinpath(projectdir(),"meshes","Expansion_"*mesh*".msh") |> normpath
  model = GmshDiscreteModel(ranks,msh_file)
  params[:model] = model
  return model
end

function expansion_mesh(mesh::Dict,ranks,params)
  @assert haskey(mesh,:num_refs_coarse) && haskey(mesh,:ranks_per_level)
  num_refs_coarse = mesh[:num_refs_coarse]
  ranks_per_level = mesh[:ranks_per_level]
  mh = Meshers.expansion_generate_mesh_hierarchy(ranks,num_refs_coarse,ranks_per_level)
  params[:multigrid] = Dict{Symbol,Any}(
    :mh => mh,
    :num_refs_coarse => num_refs_coarse,
    :ranks_per_level => ranks_per_level,
  )

  model = get_model(mh,1)
  params[:model] = model
  return model
end
