
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
  n=2,
  nc=(n,n,n),
  a = 1.0,
  L = 1.0,
  sizes = (L,a,a),
  k = 2,
  ν=1.0,
  ρ=1.0,
  σ=1.0,
  B=(0.0,1.0,0.0),
  f=(0.0,0.0,0.0),
  ζ=0.0,
  u0=1.0,
  B0=norm(VectorValue(B)),
  nsums = 10,
  vtk=true,
  title = "channel",
  path  = datadir(),
  debug = false,
  res_assemble = false,
  jac_assemble = false,
  solve = true,
  solver = :julia,
  verbose = true,
  man_solution = nothing,
  nonuniform_B = false,
  inlet=:parabolic,
  γB = 0.45,
  μ = 0.0,
  niter=10,
  bl_orders=(1,1,1),
  initial_value=:zero,
  convection=true,
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
  ax,ay,az = sizes
  L = ay
  Re = u0*L/ν
  Ha = B0*L*sqrt(σ/(ρ*ν))
  N = Ha^2/Re
  α = 1.0
  β = 1.0/Re
  γ = N
  f̄ = VectorValue(f)
  γB = γB*(30/ax)*(ay/1)

  if nonuniform_B
    B̄ = magnetic_field(;γ=γB,Ha,Re,d=ay)
    Bx = x-> B0 * B̄(x)
  else
    B̄ = (1/norm(B))*VectorValue(B)
    Bx = x-> VectorValue(B)
  end
  h = ax/nc[1]

  Reh = Re*h/L
  Hah = Ha*h/L
  @show Re
  @show Ha
  @show Reh
  @show Hah

  Z_u = 2/ay
  β_u = az/2
  ū = u_inlet(inlet,Ha,Z_u,β_u)

  # ū = inlet_profile(sizes)


  if nonuniform_B
    disc_dirs = 1
    disc_factor = γB
  else
    disc_dirs = []
    disc_factor = nothing
  end

  params[:model] = _channel_model(parts,np,sizes,nc;bl_orders,disc_dirs,disc_factor)


  if debug && vtk
    writevtk(params[:model],joinpath(path,"$(title)_model"))
  end

  # FE Space parameters
  params[:fespaces] = Dict(
  :k => k)

  # Fluid parameters
  @show convection
  params[:fluid] = Dict(
    :domain=>nothing,
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


  params[:solver][:niter] = niter
  if initial_value == :inlet
    params[:solver][:initial_values] = Dict(
      :u=>ū,:j=>j_zero,:p=>0.0,:φ=>0.0)
  end
  if params[:solver][:solver] == :petsc
    params[:solver][:petsc_options] *= " -snes_max_funcs $(niter+1)"
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

  Ω = Interior(params[:model])
  # Ω_phys = L == 1.0 ? Ω : _warp(model,Ω,L)
  Ω_phys = Ω

  degree = k*2
  dΩ = Measure(Ω_phys,degree)
  l2(u,dΩ) = √(∑( ∫(u ⊙ u)dΩ ))
  h1(u,dΩ) = √(∑( ∫( ∇(u) ⊙ ∇(u) + u ⊙ u)dΩ ))

  ūh,p̄h,j̄h,φ̄h = xh
  uh = u0*ūh
  ph = (ρ*u0^2)*p̄h
  jh = (σ*u0*B0)*j̄h
  φh = (u0*B0*L)*φ̄h


  info[:ncells] = num_cells(params[:model])
  info[:Re] = Re
  info[:Ha] = Ha

  if vtk
    writevtk(Ω_phys,joinpath(path,title),
      order=2,
      cellfields=[
        "uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh,"B"=>Bx])
    toc!(t,"vtk")
  end

  info, t
end

function inlet_profile(L)
  u_2 = inlet_profile_1d(L[2])
  u_3 = inlet_profile_1d(L[3])
  function inlet_profile_fun(x)
    VectorValue( u_2(x[2])*u_3(x[3]), 0.0, 0.0  )
  end
end

function inlet_profile_1d(L)
  function inlet_profile_1d_fun(x)
    -((2/L)^2)*x*(x-L)
  end
end

function _channel_model(parts,np,sizes,nc;
  bl_orders=(1,2,2),disc_dirs=1,disc_factor=1)
  D = length(nc)
  domain = ntuple(Val{D*2}()) do i
    isodd(i) ? -0.5*sizes[(i+1)÷2] : 0.5*sizes[(i+1)÷2]
  end
  map1 = boundary_layer_map(domain,bl_orders)
  map2 = discontinuity_map(domain,disc_dirs,disc_factor)
  model = CartesianDiscreteModel(parts,np,domain,nc,map=map1∘map2)
  p = Polytope(Gridap.Helpers.tfill(HEX_AXIS,Val{D}()))
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
  map(local_views(model)) do model
    labels = get_face_labeling(model)
    add_tag!(labels,"inlet",faces_i)
    add_tag!(labels,"outlet",faces_o)
    add_tag!(labels,"walls",faces_w)
    add_tag!(labels,"insulating",faces_ins)
  end
  model
end


function magnetic_field(;kwargs...)
  B = magnetic_field_1d(;kwargs...)
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


function magnetic_field_1d(;γ=0.45,Ha=1,Re=0,d=1,c=0)
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
