function FullyDeveloped(;
  backend=nothing,
  np=nothing,
  parts=nothing,
  title = "FullyDeveloped",
  nruns=1,
  path=".",
  kwargs...)

  for ir in 1:nruns
    _title = title*"_r$ir"
    @assert parts === nothing
    if backend === nothing
      @assert np === nothing
      info, t = _FullyDeveloped(;title=_title,path=path,kwargs...)
    else
      @assert backend !== nothing
      info, t = with_backend(_find_backend(backend),(np...,1)) do _parts
        _FullyDeveloped(;parts=_parts,title=_title,path=path,kwargs...)
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

function _FullyDeveloped(;
  parts=nothing,
  nc=(3,3),
  Ha = 10.0,
  dir_B=(0.0,1.0,0.0),
  b=1.0,
  vtk=true,
  title="test",
  path=".",
  debug=false,
  res_assemble=false,
  jac_assemble=false,
  solve=true,
  solver=:julia,
  verbose=true,
  mesh = false,
  cw_Ha = 0.0,
  cw_s = 0.0,
  τ_Ha = 100.0,
  τ_s = 100.0,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
  )

  info = Dict{Symbol,Any}()

  if parts === nothing
    t_parts = get_part_ids(SequentialBackend(),1)
  else
    t_parts = parts
  end
  t = PTimer(t_parts,verbose=verbose)
  tic!(t,barrier=true)

  domain = (-b,b,-1.0,1.0,0.0,0.1)

  # Reduced quantities
  Re = 1.0
  N = Ha^2/Re
  dirB = (1/norm(VectorValue(dir_B)))*VectorValue(dir_B)
  α = 1.0/N
  β = 1.0/Ha^2
  γ = 1.0

  # Prepare problem in terms of reduced quantities

   strech_Ha = sqrt(Ha/(Ha-1))
   strech_side = sqrt(sqrt(Ha)/(sqrt(Ha)-1))

  function map1(coord)
     ncoord = GridapMHD.strechMHD(coord,domain=(0,-b,0,-1.0),factor=(strech_side,strech_Ha),dirs=(1,2))
     ncoord = GridapMHD.strechMHD(ncoord,domain=(0,b,0,1.0),factor=(strech_side,strech_Ha),dirs=(1,2))
     ncoord  
   end

 partition=(nc[1],nc[2],3)
 model = CartesianDiscreteModel(
     parts,domain,partition;isperiodic=(false,false,true),map=map1)
  Ω = Interior(model)
  
  labels = get_face_labeling(model)
  tags_u = append!(collect(1:20),[23,24,25,26])
  tags_j_Ha = append!(collect(1:20),[25,26])
  tags_j_side = append!(collect(1:20),[23,24])

  add_tag_from_tags!(labels,"noslip",tags_u)
  add_tag_from_tags!(labels,"Ha_walls",tags_j_Ha)
  add_tag_from_tags!(labels,"side_walls",tags_j_side)
  
  if mesh
    writevtk(model,"Mesh")
  end

  params = Dict(
    :ptimer=>t,
    :debug=>debug,
    :model=>model,
    :res_assemble=>res_assemble,
    :jac_assemble=>jac_assemble,
    :solve=>solve,
    :fluid=>Dict(
      :domain=>model,
      :α=>α,
      :β=>β,
      :γ=>γ,
      :f=>VectorValue(0.0,0.0,1.0),
      :B=>dirB,
    ),
    :bcs => Dict(
      :u=>Dict(:tags=>["noslip"]),
      :j=>Dict{Symbol,Vector{String}}(),
      :thin_wall=>Vector{Dict{Symbol,Any}}() 
    )
  )

  insulated_tags = Vector{String}()
  thinWall_options = Vector{Dict{Symbol,Any}}() 

  if cw_Ha == 0.0 
    push!(insulated_tags,"Ha_walls")
  else
    push!(thinWall_options, Dict(
        :cw => cw_Ha,
        :τ => τ_Ha, 
        :domain =>Boundary(model, tags="Ha_walls")
        )
      ) 
  end  

  if cw_s == 0.0
    push!(insulated_tags,"side_walls")
  else
    push!(thinWall_options, Dict(
        :cw => cw_s,
        :τ => τ_s,
        :domain =>Boundary(model, tags="side_walls")
        )
      )
  end
  
#  if insulated_tags != [] 
    params[:bcs][:j] = Dict(:tags=>insulated_tags)  
#  end

  if thinWall_options != []
    params[:bcs][:thin_wall] = thinWall_options
  end
 
  toc!(t,"pre_process")

  # Solve it
  if solver == :julia
    params[:solver] = NLSolver(show_trace=true,method=:newton)
    xh,fullparams = main(params)
  elseif solver == :petsc
    xh,fullparams = GridapPETSc.with(args=split(petsc_options)) do
      params[:matrix_type] = SparseMatrixCSR{0,PetscScalar,PetscInt}
      params[:vector_type] = Vector{PetscScalar}
      params[:solver] = PETScNonlinearSolver()
      params[:solver_postpro] = cache -> snes_postpro(cache,info)
      main(params)
    end
  else
    error()
  end
  t = fullparams[:ptimer]

  # Rescale quantities

  tic!(t,barrier=true)
  uh,ph,jh,φh = xh
  
  div_jh = ∇·jh 
  div_uh = ∇·uh  

  #Post process

  if vtk
    writevtk(Ω,joinpath(path,title),
      order=2,
      cellfields=[
        "uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh,"div_jh"=>div_jh,"div_uh"=>div_uh])
    toc!(t,"vtk")
  end
  if verbose
    display(t)
  end

  info[:nc] = nc
  info[:ncells] = num_cells(model)
  info[:ndofs_u] = length(get_free_dof_values(uh))
  info[:ndofs_p] = length(get_free_dof_values(ph))
  info[:ndofs_j] = length(get_free_dof_values(jh))
  info[:ndofs_φ] = length(get_free_dof_values(φh))
  info[:ndofs] = length(get_free_dof_values(xh))
  info[:Re] = Re
  info[:Ha] = Ha
  info[:N] = N

  info, t
end


