function hunt(;
  backend = nothing,
  np      = nothing,
  title   = "hunt",
  nruns   = 1,
  path    = datadir(),
  kwargs...)

  for ir in 1:nruns
    _title = title*"_r$ir"
    if isa(backend,Nothing)
      @assert isa(np,Nothing)
      info, t = _hunt(;title=_title,path=path,kwargs...)
    else
      @assert backend ∈ [:sequential,:mpi]
      @assert !isa(np,Nothing)
      if backend === :sequential
        info, t = with_debug() do distribute
          _hunt(;distribute=distribute,rank_partition=np,title=_title,path=path,kwargs...)
        end
      else
        info,t = with_mpi() do distribute
          _hunt(;distribute=distribute,rank_partition=np,title=_title,path=path,kwargs...)
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

function _hunt(;
  distribute=nothing,
  rank_partition=nothing,
  nc=(4,4),
  ν=1.0,
  ρ=1.0,
  σ=1.0,
  B=(0.0,10.0,0.0),
  f=(0.0,0.0,1.0),
  ζ=0.0,
  L=1.0,
  u0=1.0,
  B0=norm(VectorValue(B)),
  σw1=0.1,
  σw2=10.0,
  tw=0.0,
  nsums = 10,
  vtk=true,
  title = "test",
  path  = datadir(),
  debug = false,
  res_assemble = false,
  jac_assemble = false,
  solve = true,
  solver = :julia,
  verbose = true,
  BL_adapted = true,
  kmap_x = 1,
  kmap_y = 1,
  ranks_per_level = nothing
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
  N = Ha^2/Re
  f̄ = (L/(ρ*u0^2))*VectorValue(f)
  B̄ = (1/B0)*VectorValue(B)
  α = 1.0
  β = 1.0/Re
  γ = N
  σ̄1 = σw1/σ
  σ̄2 = σw2/σ

  # DiscreteModel in terms of reduced quantities

  model = hunt_mesh(parts,params,nc,rank_partition,L,tw,Ha,kmap_x,kmap_y,BL_adapted,ranks_per_level)
  Ω = Interior(model)
  if debug && vtk
    writevtk(model,"data/hunt_model")
  end

  params[:fluid] = Dict(
    :domain=>nothing,
    :α=>α,
    :β=>β,
    :γ=>γ,
    :f=>f̄,
    :B=>B̄,
    :ζ=>ζ,
  )

  if tw > 0.0
    σ_Ω = solid_conductivity(Ω,get_cell_gids(model),get_face_labeling(model))
    params[:solid] = Dict(:domain=>"solid",:σ=>σ_Ω)
  end

  # Boundary conditions

  params[:bcs] = Dict(
    :u=>Dict(:tags=>"noslip"),
    :j=>Dict(:tags=>"insulating"),
  )

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
  φh = (u0*B0*L)*φ̄h

  if L == 1.0
    Ω_phys = Ω
  else
    Ω_phys = _warp(model,Ω,L)
  end

  # Post process

  μ = ρ*ν
  grad_pz = -f[3]/ρ
  u(x) = analytical_hunt_u(L,L,μ,grad_pz,Ha,nsums,x)
  j(x) = analytical_hunt_j(L,L,σ,μ,grad_pz,Ha,nsums,x)
  u_ref(x) = analytical_hunt_u(L,L,μ,grad_pz,Ha,2*nsums,x)
  j_ref(x) = analytical_hunt_j(L,L,σ,μ,grad_pz,Ha,2*nsums,x)

  k = 2
  dΩ_phys = Measure(Ω_phys,2*(k+1))
  eu = u - uh
  ej = j - jh
  eu_h1 = sqrt(sum(∫( ∇(eu)⊙∇(eu) + eu⋅eu  )dΩ_phys))
  eu_l2 = sqrt(sum(∫( eu⋅eu )dΩ_phys))
  ej_l2 = sqrt(sum(∫( ej⋅ej )dΩ_phys))
  eu_ref = u_ref - uh
  ej_ref = j_ref - jh
  eu_ref_h1 = sqrt(sum(∫( ∇(eu_ref)⊙∇(eu_ref) + eu_ref⋅eu_ref  )dΩ_phys))
  eu_ref_l2 = sqrt(sum(∫( eu_ref⋅eu_ref )dΩ_phys))
  ej_ref_l2 = sqrt(sum(∫( ej_ref⋅ej_ref )dΩ_phys))
  uh_l2 = sqrt(sum(∫( uh⋅uh )dΩ_phys))
  uh_h1 = sqrt(sum(∫( ∇(uh)⊙∇(uh) + uh⋅uh )dΩ_phys))
  jh_l2 = sqrt(sum(∫( jh⋅jh )dΩ_phys))
  toc!(t,"post_process")

  if vtk
    if tw > 0.0
      writevtk(Ω_phys,joinpath(path,title),
        order=2,
        cellfields=[
          "uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh,
          "u"=>u,"j"=>j,"u_ref"=>u_ref,"j_ref"=>j_ref,"σ"=>σ_Ω])
    else
      writevtk(Ω_phys,joinpath(path,title),
        order=2,
        cellfields=[
          "uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh,
          "u"=>u,"j"=>j,"u_ref"=>u_ref,"j_ref"=>j_ref])
    end
    toc!(t,"vtk")
  end
  if verbose
    println(" >> H1 error for u = ", eu_ref_h1)
    println(" >> L2 error for j = ", ej_ref_l2)
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
  info[:eu_ref_h1] = eu_ref_h1
  info[:eu_ref_l2] = eu_ref_l2
  info[:ej_ref_l2] = ej_ref_l2
  info[:uh_h1] = uh_h1
  info[:uh_l2] = uh_l2
  info[:jh_l2] = jh_l2
  info[:kmap] = [kmap_x,kmap_y]

  info, t
end

function hunt_mesh(
  parts,params,
  nc::Tuple,np::Tuple,L::Real,tw::Real,Ha::Real,kmap_x::Number,kmap_y::Number,BL_adapted::Bool,
  ranks_per_level)
  if isnothing(ranks_per_level) # Single grid
    model = Meshers.hunt_generate_base_mesh(parts,np,nc,L,tw,Ha,kmap_x,kmap_y,BL_adapted)
    params[:model] = model
  else # Multigrid
    base_model = Meshers.hunt_generate_base_mesh(nc,L,tw,Ha,kmap_x,kmap_y,BL_adapted)
    mh = Meshers.generate_mesh_hierarchy(parts,base_model,0,ranks_per_level)
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

function solid_conductivity(Ω::GridapDistributed.DistributedTriangulation,cells,labels)
  D = num_cell_dims(Ω)
  fields =  map(local_views(labels),local_views(cells),local_views(Ω)) do labels,partition,trian
    cell_entity = labels.d_to_dface_to_entity[end]
    σ_field(labels,trian,cell_entity[partition.own_to_local])
  end
  GridapDistributed.DistributedCellField(fields)
end

function σ_field(labels,Ω,cell_entity)
  solid_1 = get_tag_entities(labels,"solid_1")
  solid_2 = get_tag_entities(labels,"solid_2")
  function entity_to_σ(entity)
    if entity == solid_1
      σ̄1
    elseif entity == solid_2
      σ̄2
    else
      1.0
    end
  end
  σ_field = CellField(map(entity_to_σ,cell_entity),Ω)
end

# This is not very elegant. This needs to be solved by Gridap and GridapDistributed
function _warp(model::DiscreteModel,Ω::Triangulation,L)
  grid_phys = UnstructuredGrid(get_grid(model))
  node_coords = get_node_coordinates(grid_phys)
  node_coords .= L .* node_coords
  Ω_phys = Gridap.Geometry.BodyFittedTriangulation(Ω.model,grid_phys,Ω.tface_to_mface)
end

function _warp(
  model::GridapDistributed.DistributedDiscreteModel,
  Ω::GridapDistributed.DistributedTriangulation,L)
  trians = map_parts(model.models,Ω.trians) do model,Ω
    grid_phys = UnstructuredGrid(get_grid(model))
    node_coords = get_node_coordinates(grid_phys)
    node_coords .= L .* node_coords
    gp = GridPortion(grid_phys,Ω.tface_to_mface)
    Ω_phys = Gridap.Geometry.BodyFittedTriangulation(model,gp,Ω.tface_to_mface)
  end
  GridapDistributed.DistributedTriangulation(trians,model)
end

function analytical_hunt_u(
  a::Float64,       # semi-length of side walls
  b::Float64,       # semi-length of Hartmann walls
  μ::Float64,       # fluid viscosity
  grad_pz::Float64, # presure gradient
  Ha::Float64,      # Hartmann number
  n::Int,           # number of sumands included in Fourier series
  x)                # evaluation point

  l = b/a
  ξ = x[1]/a
  η = x[2]/a

  V = 0.0; V0=0.0;
  for k in 0:n
    α_k = (k + 0.5)*π/l
    N = (Ha^2 + 4*α_k^2)^(0.5)
    r1_k = 0.5*( Ha + N)
    r2_k = 0.5*(-Ha + N)

    num = exp(-r1_k*(1-η))+exp(-r1_k*(1+η))
    den = 1+exp(-2*r1_k)
    V2 = (r2_k/N)*(num/den)

    num = exp(-r2_k*(1-η))+exp(-r2_k*(1+η))
    den = 1+exp(-2*r2_k)
    V3 = (r1_k/N)*(num/den)


    V += 2*(-1)^k*cos(α_k * ξ)/(l*α_k^3) * (1-V2-V3)
  end
  u_z = V/μ * (-grad_pz) * a^2

  VectorValue(0.0*u_z,0.0*u_z,u_z)
end


function analytical_hunt_j(
  a::Float64,       # semi-length of side walls
  b::Float64,       # semi-length of Hartmann walls
  σ::Float64,       # fluid conductivity
  μ::Float64,       # fluid viscosity
  grad_pz::Float64, # presure gradient
  Ha::Float64,      # Hartmann number
  n::Int,           # number of sumands included in Fourier series
  x)                # evaluation point

  l = b/a
  ξ = x[1]/a
  η = x[2]/a

  H_dx = 0.0; H_dy = 0.0
  for k in 0:n
    α_k = (k + 0.5)*π/l
    N = sqrt(Ha^2 + 4*α_k^2)
    r1_k = 0.5*( Ha + N)
    r2_k = 0.5*(-Ha + N)

    num = exp(-r1_k*(1-η)) - exp(-r1_k*(1+η))
    num_dy = exp(-r1_k*(1-η))*(r1_k/a) + exp(-r1_k*(1+η))*(r1_k/a)
    den = 1+exp(-2*r1_k)
    H2 = (r2_k/N)*(num/den)
    H2_dy = (r2_k/N)*(num_dy/den)

    num = exp(-r2_k*(1-η)) - exp(-r2_k*(1+η))
    num_dy = exp(-r2_k*(1-η))*(r2_k/a) + exp(-r2_k*(1+η))*(r2_k/a)
    den = 1+exp(-2*r2_k)
    H3 = (r1_k/N)*(num/den)
    H3_dy = (r1_k/N)*(num_dy/den)

    H_dx += -2*(-1)^k * sin(α_k * ξ)/(a*l*α_k^2) * (H2 - H3)
    H_dy += 2*(-1)^k * cos(α_k * ξ)/(l*α_k^3) * (H2_dy - H3_dy)
  end
  j_x = a^2*σ^0.5 / μ^0.5 * (-grad_pz) * H_dy
  j_y = a^2*σ^0.5 / μ^0.5 * (-grad_pz) * (-H_dx)

  VectorValue(j_x,j_y,0.0)
end
