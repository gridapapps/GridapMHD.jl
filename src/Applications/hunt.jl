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
        info,t = with_debug() do distribute
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
  L=1.0,
  u0=1.0,
  B0=norm(VectorValue(B)),
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
  mesh = false,
  BL_adapted = true,
  kmap_x = 1,
  kmap_y = 1,
  )

  info = Dict{Symbol,Any}()

  if isa(distribute,Nothing)
    @assert isa(rank_partition,Nothing)
    rank_partition = Tuple(fill(1,length(nc)))
    distribute = DebugArray
  end
  @assert length(rank_partition) == length(nc)
  parts = distribute(LinearIndices((prod(rank_partition),)))
  
  t = PTimer(parts,verbose=verbose)
  tic!(t,barrier=true)

  if isa(solver,Symbol)
    solver = default_solver_params(Val(solver))
  end

  domain_phys = (-L,L,-L,L,0.0*L,0.1*L)

  # Reduced quantities
  Re = u0*L/ν
  Ha = B0*L*sqrt(σ/(ρ*ν))
  N = Ha^2/Re
  f̄ = (L/(ρ*u0^2))*VectorValue(f)
  B̄ = (1/B0)*VectorValue(B)
  α = 1.0
  β = 1.0/Re
  γ = N
  domain = domain_phys ./ L

  # Prepare problem in terms of reduced quantities

  strech_Ha = sqrt(Ha/(Ha-1))
  strech_side = sqrt(sqrt(Ha)/(sqrt(Ha)-1))

  function map1(coord)
    ncoord = GridapMHD.strechMHD(coord,domain=(0,-L,0,-L),factor=(strech_side,strech_Ha),dirs=(1,2))
    ncoord = GridapMHD.strechMHD(ncoord,domain=(0,L,0,L),factor=(strech_side,strech_Ha),dirs=(1,2))
    ncoord  
  end
  layer(x,a) = sign(x)*abs(x)^(1/a)
  map2((x,y,z)) = VectorValue(layer(x,kmap_x),layer(y,kmap_y),z)

  mesh_partition = (nc[1],nc[2],3)
  mesh_rank_partition = (rank_partition[1],rank_partition[2],1)
  if BL_adapted
    model = CartesianDiscreteModel(
      parts,mesh_rank_partition,domain,mesh_partition;isperiodic=(false,false,true),map=map1)
  else
    model = CartesianDiscreteModel(
      parts,mesh_rank_partition,domain,mesh_partition;isperiodic=(false,false,true),map=map2)
  end
  Ω = Interior(model)
  labels = get_face_labeling(model)
  tags_u = append!(collect(1:20),[23,24,25,26])
  tags_j = append!(collect(1:20),[25,26])
  add_tag_from_tags!(labels,"noslip",tags_u)
  add_tag_from_tags!(labels,"insulating",tags_j)
  
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
      :f=>f̄,
      :B=>B̄,
    ),
    :bcs => Dict(
      :u=>Dict(:tags=>"noslip"),
      :j=>Dict(:tags=>"insulating"),
    ),
    :solver=>solver,
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
    writevtk(Ω_phys,joinpath(path,title),
      order=2,
      cellfields=[
        "uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh,
        "u"=>u,"j"=>j,"u_ref"=>u_ref,"j_ref"=>j_ref])
    toc!(t,"vtk")
  end
  if verbose
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

function snes_postpro(cache,info)
  snes = cache.snes[]
  i_petsc = Ref{PetscInt}()
  @check_error_code PETSC.SNESGetIterationNumber(snes,i_petsc)
  info[:nls_iters] = Int(i_petsc[])
  @check_error_code PETSC.SNESGetLinearSolveIterations(snes,i_petsc)
  info[:ls_iters] = Int(i_petsc[])
  nothing
end