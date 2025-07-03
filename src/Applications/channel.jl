
function channel(;
  backend = nothing,
  np      = 1,
  title   = "Channel",
  path    = datadir(),
  kwargs...)

  if isa(backend,Nothing)
    @assert np == 1
    info, t = _channel(;title=title,path=path,kwargs...)
  else
    @assert backend ∈ [:sequential,:mpi]
    if backend === :sequential
      info,t = with_debug() do distribute
        _channel(;distribute=distribute,np=np,title=title,path=path,kwargs...)
      end
    else
      info,t = with_mpi() do distribute
        _channel(;distribute=distribute,np=np,title=title,path=path,kwargs...)
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

function _channel(;
  distribute = nothing,
  np = nothing,
  nc=(4,2,2),
  L = 1.0, # Channel length (x-axis)
  w = 1.0, # Channel width (y/z-axis)
  sizes = (L,w,w),
  order = 2,
  order_j = 2,
  ν = 1.0,
  ρ = 1.0,
  σ = 1.0,
  B = (0.0,1.0,0.0),
  f = (0.0,0.0,0.0),
  ζ = 0.0,
  μ = 0.0,
  u0 = 1.0,
  B0 = norm(VectorValue(B)),
  vtk = true,
  title = "channel",
  path  = datadir(),
  debug = false,
  res_assemble = false,
  jac_assemble = false,
  solve = true,
  solver = :julia,
  verbose = true,
  nonuniform_B = false,
  inlet = :parabolic,
  γB = 0.45,
  bl_orders=(1,1,1),
  initial_value = :zero,
  convection = :newton,
  np_per_level = nothing,
  rt_scaling = false,
  formulation = :mhd,
  adaptivity_method = 0,
  fluid_disc = ifelse(iszero(adaptivity_method),:Qk_dPkm1,:SV),
  current_disc = :RT,
)
  @assert inlet ∈ [:parabolic,:shercliff,:constant]
  @assert initial_value ∈ [:zero,:inlet,:solve]
  @assert formulation ∈ [:cfd,:mhd]
  @assert convection ∈ [:newton,:picard,:none]

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
  Lx,Ly,Lz = sizes
  Re = u0*Ly/ν
  Ha = B0*Ly*sqrt(σ/(ρ*ν))
  N = Ha^2/Re

  if formulation == :cfd # Option 1 (CFD)
    α = 1.0
    β = 1.0/Re
    γ = N
    f̄ = VectorValue(f)
  elseif formulation == :mhd # Option 2 (MHD) is chosen in the experimental article
    α = (1.0/N)
    β = (1.0/Ha^2)
    γ = 1.0
    f̄ = VectorValue(f)/N
  else
    error("Unknown formulation")
  end

  γB = γB*(30/Lx)*(Ly/1)

  if nonuniform_B
    B̄ = channel_magnetic_field(;γ=γB,Ha,Re,d=Ly)
    Bx = x-> B0 * B̄(x)
  else
    B̄ = (1/norm(B))*VectorValue(B)
    Bx = x-> VectorValue(B)
  end

  Z_u = 2/Ly
  β_u = Lz/2
  ū = u_inlet(inlet,Ha,Z_u,β_u)

  disc_dirs = nonuniform_B ? 1 : []
  disc_factor = nonuniform_B ? γB : nothing
  model = channel_mesh(parts,np,sizes,nc,params;bl_orders,disc_dirs,disc_factor,np_per_level,adaptivity_method)

  if debug && vtk
    writevtk(params[:model],joinpath(path,"$(title)_model"))
  end

  # FE Space parameters
  params[:fespaces] = Dict(
    :order_u => order,
    :order_j => order_j,
    :rt_scaling => rt_scaling ? 1.0/get_mesh_size(model) : nothing,
    :fluid_disc => fluid_disc,
    :current_disc => current_disc,
  )

  # Fluid parameters
  params[:fluid] = Dict(
    :domain => nothing,
    :α=>α,
    :β=>β,
    :γ=>γ,
    :f=>f̄,
    :B=>B̄,
    :ζ=>ζ,
    :convection => convection,
  )

  # Boundary conditions
  u_zero = VectorValue(0.0,0.0,0.0)
  j_zero = VectorValue(0.0,0.0,0.0)
  params[:bcs] = Dict{Symbol,Any}()
  params[:bcs][:u] = Dict(:tags=>["inlet","walls"],:values=>[ ū, u_zero ] )
  params[:bcs][:j] = Dict(:tags=>["inlet","outlet","walls"])

  if μ > 0
    params[:bcs][:stabilization] = Dict(:μ=>μ)
  end

  toc!(t,"pre_process")

  if initial_value == :inlet
    params[:x0] = Dict(
      :u=>ū,:j=>j_zero,:p=>0.0,:φ=>0.0
    )
  end

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

  Ω = Interior(model)
  k = max(order,order_j)
  degree = k*2
  dΩ = Measure(Ω,degree)
  l2(u,dΩ) = sqrt(sum( ∫(u ⊙ u)dΩ ))
  h1(u,dΩ) = sqrt(sum( ∫( ∇(u) ⊙ ∇(u) + u ⊙ u)dΩ ))

  ūh,p̄h,j̄h,φ̄h = xh
  uh = u0*ūh
  ph = (ρ*u0^2)*p̄h
  jh = (σ*u0*B0)*j̄h
  φh = (u0*B0*Ly)*φ̄h
  div_jh = ∇·jh
  div_uh = ∇·uh
  Grad_p = ∇·ph

  info[:ncells] = num_cells(params[:model])
  info[:Re] = Re
  info[:Ha] = Ha

  if vtk
    vtk_order = min(k,3) 
    writevtk(
      Ω, joinpath(path,title),
      order = vtk_order,
      cellfields = [
        "uh"     => uh,
        "ph"     => ph,
        "jh"     => jh,
        "phi"    => φh,
        "B"      => Bx,
        "div_jh" => div_jh,
        "div_uh" => div_uh,
        "Grad_p" => Grad_p
      ],
      append = false
    )
    toc!(t,"vtk")
  end

  info, t
end

function channel_add_tags!(labels::GridapDistributed.DistributedFaceLabeling)
  map(channel_add_tags!,local_views(labels))
end

function channel_add_tags!(labels::FaceLabeling)
  D = length(labels.d_to_dface_to_entity) - 1
  p = (D == 2) ? QUAD : HEX

  facets_i = [2*D-1]
  facets_o = [2*D]
  facets_w = [1:2*(D-1);]
  facets_ins = [1:2;]

  faces_i = facets_i .+ get_offset(p,D-1)
  faces_o = facets_o .+ get_offset(p,D-1)
  faces_w = facets_w .+ get_offset(p,D-1)
  faces_ins = facets_ins .+ get_offset(p,D-1)
  faces_w = reduce(union, get_faces(p)[faces_w] ) |> sort
  faces_ins = reduce(union, get_faces(p)[faces_ins] ) |> sort
  
  add_tag!(labels,"inlet",faces_i)
  add_tag!(labels,"outlet",faces_o)
  add_tag!(labels,"walls",faces_w)
  add_tag!(labels,"insulating",faces_ins)
end

function channel_stretch_map(domain,bl_orders=(1,2,2),disc_dirs=1,disc_factor=1)
  map1 = boundary_layer_map(domain,bl_orders)
  map2 = discontinuity_map(domain,disc_dirs,disc_factor)
  return map1∘map2
end

function channel_mesh(
  parts,np,sizes,nc,params;
  bl_orders=(1,2,2),disc_dirs=1,disc_factor=1,
  np_per_level=nothing,adaptivity_method = 0
)
  D = length(nc)
  domain = ntuple(Val{D*2}()) do i
    isodd(i) ? -0.5*sizes[(i+1)÷2] : 0.5*sizes[(i+1)÷2]
  end
  coord_map = channel_stretch_map(domain,bl_orders,disc_dirs,disc_factor)
  if isnothing(np_per_level)
    model = CartesianDiscreteModel(parts,np,domain,nc,map=coord_map)
    channel_add_tags!(get_face_labeling(model))
  else
    mh = GridapSolvers.CartesianModelHierarchy(
      parts,np_per_level,domain,nc;
      nrefs = (1,2,2), # Only refine y-z plane
      map = coord_map, add_labels! = channel_add_tags!
    )
    params[:multigrid] = Dict{Symbol,Any}(
      :mh => mh,
      :num_refs_coarse => 0,
      :ranks_per_level => np_per_level,
    )
    model = get_model(mh,1)
  end
  model = Meshers.adapt_mesh(model,adaptivity_method)
  params[:model] = model
  return model
end

function channel_magnetic_field(;kwargs...)
  B = channel_magnetic_field_1d(;kwargs...)
  dB(x) = ForwardDiff.derivative(B,x)
  d²B(x) = ForwardDiff.derivative(dB,x)
  d³B(x) = ForwardDiff.derivative(d²B,x)
  d⁴B(x) = ForwardDiff.derivative(d³B,x)

  function magnetic_field_fun(X)
    x,y = X[1],X[2]
    Bx = dB(x)*y - d³B(x)*y^3/6
    By = B(x) - d²B(x)*y^2/2 + d⁴B(x)*y^4/24
    Bz = 0.0
    VectorValue(Bx,By,Bz)
  end
end

function channel_magnetic_field_1d(;γ=0.45,Ha=1,Re=0,d=1,c=0)
  Ha_Re_crit = 1/200
  bmin = Ha_Re_crit * Re / Ha
  b0 = 1-bmin
  function magnetic_field_1d_fun(x)
    b0*(( 1 - tanh(γ*(x/d -c)) ) / 2) + bmin
  end
end

function boundary_layer_map(domain,orders=(1,2,2))
  D = length(domain) ÷ 2
  function ref(x,d)
    xmin = domain[2*d-1]
    xmax = domain[2*d]
    (x-xmin) * 2 / (xmax-xmin) - 1
  end
  function phys(x,d)
    xmin = domain[2*d-1]
    xmax = domain[2*d]
    (x + 1)*(xmax-xmin) / 2 + xmin
  end
  ref_layer(x,a) = sign(x)*abs(x)^(1/a)
  function layer(x,d)
    x = ref(x,d)
    x = ref_layer(x,orders[d])
    phys(x,d)
  end
  function _boundary_layer_map(x)
    data = ntuple(d->layer(x[d],d),Val{D}())
    return VectorValue(data...)
  end
end

function discontinuity_map(domain,dirs=1,f=1)
  # TODO: Add support to change discontinuity offset
  D = length(domain) ÷ 2
  function ref(x,d)
    xmin = domain[2*d-1]
    xmax = domain[2*d]
    (x-xmin) * 2 / (xmax-xmin) - 1
  end
  function phys(x,d)
    xmin = domain[2*d-1]
    xmax = domain[2*d]
    (x + 1)*(xmax-xmin) / 2 + xmin
  end
  function ref_tanh_map(x,d)
    xmin = domain[2*d-1]
    xmax = domain[2*d]
    L = xmax - xmin
    fr = f*L
    tanh(x/fr)*abs(x) / tanh(1/fr)
  end
  function tanh_map(x,d)
    x = ref(x,d)
    x = ref_tanh_map(x,d)
    phys(x,d)
  end
  function _discontinuity_map(x)
    for d in dirs
      xd = tanh_map(x[d],d)
      data = Base.setindex(x.data,xd,d)
      x = VectorValue(data...)
    end
    x
  end
end
