
# MHD equations
# ∇⋅u = 0 in fluid
# u⋅∇(u) -ν*Δ(u) + (1/ρ)*∇(p) - (1/ρ)*(j×B) = (1/ρ)*f in fluid
# ∇⋅j = 0 in fluid and solid
# j + σ*∇(φ) - σ*(u×B) = 0 in fluid
# j + σ*∇(φ) = 0 in solid
# solving for u,p,j,φ for a given B,ν,ρ,σ
#
# One can provide characteristic quantities
# u0,B0,L
#
# To introduce these change of variables
# u = u0*ū, B = B0*B̄, j = σ*u0*B0*j̄, φ = u0*B0*L*φ̄ (I am not sure about L in φ)
# p = ρ*u0^2*p̄, f = (ρ*u0^2/L)*f̄ (option 1)
# p = σ*u0*B0*L*p̄, f = (σ*u0*B0^2)*f̄ (option 2)
#
# To solve the following equations in a scaled domain according to L
# ∇⋅ū = 0 in fluid
# ∇⋅j̄ = 0 in fluid and solid
# j̄ + ∇(φ̄) - ū×B̄ = 0 in fluid
# j̄ + ∇(φ̄) = 0 in solid
# ū⋅∇(ū) -(1/Re)*Δ(ū) + ∇(p̄) - N*(j̄×B̄) = f̄ in fluid (Option 1, CFD)
# (1/N)*ū⋅∇(ū) - (1/Ha^2)*Δ(ū) + ∇(p̄) - (j̄×B̄) = f̄ in fluid (Option 2,MHD)
#
# with
#  Re = u0*L/ν
#  Ha = B0*L*sqrt(σ/(ρ*ν))
#  N = Ha^2/Re

# In order to account for both options, the code solves these equations
#
# ∇⋅ū = 0 in fluid
# ∇⋅j̄ = 0 in fluid and solid
# j̄ + σ̄*∇(φ̄) - σ̄*(ū×B̄) = 0 in fluid
# j̄ + σ̄*∇(φ̄) = 0 in solid
# α*ū⋅∇(ū) - β*Δ(ū) + ∇(p̄) - γ*(j̄×B̄) = f̄ in fluid
#
# α = 1, β = (1/Re), γ = N (option 1,CFD)
# α = (1/N), β = (1/Ha^2), γ = 1 (option 2,MHD)
# σ̄ is 1, but is included just in case you want to solve
# the un-scaled equations directly
#
# Some boundary conditions
# (for simplicity we drop the bars, i.e. ū is simply u)
# 
# Velocity bc
# u = u (imposed strongly)
#
# Charge bc
# j = j (imposed strongly in the normal direction to the boundary)
#
# Traction bC
# # n⋅∇(u) - p*n = t
#
# Potential bc
# φ = φ (Imposed weakly)
#
# Thin wall bc
# j⋅n + cw*n⋅∇(j)⋅n = jw (imposed weakly via a penalty of value τ)

"""
    add_default_params(params::Dict{Symbol}) -> full_params

Return a new `Dict` with the contents of `params` plus
default values for the optional keys not provided in `params`.
It also checks the validity of the main parameter dictionary `params`.
"""
function add_default_params(_params)
  mandatory = Dict(
    :ptimer=>false,
    :debug=>false,
    :solver=>false,
    :matrix_type=>false,
    :vector_type=>false,
    :model=>true,
    :k=>false,
    :fluid=>true,
    :solid=>false,
    :bcs=>true,
    :check_valid=>false,
    :solver_postpro=>false,
  )
  _check_mandatory(_params,mandatory,"")
  optional = Dict(
    :solid=>nothing,
    :ptimer=>default_ptimer(_params[:model]),
    :debug=>false,
    :solver=>NLSolver(show_trace=true,method=:newton),
    :solver_postpro => (x->nothing),
    :matrix_type=>SparseMatrixCSC{Float64,Int},
    :vector_type=>Vector{Float64},
    :k=>2,
    :check_valid=>true,
  )
  params = _add_optional(_params,mandatory,optional,_params,"")
  _check_unused(params,mandatory,params,"")
  # Process sub-params
  params[:fluid] = params_fluid(params)
  if params[:solid] !== nothing
    params[:solid] = params_solid(params)
  end
  params[:bcs] = params_bcs(params)
  params
end

default_ptimer(model) = PTimer(get_part_ids(sequential,1))
default_ptimer(model::GridapDistributed.DistributedDiscreteModel) = PTimer(get_part_ids(model.models))

"""
Valid keys for `params[:fluid]` are the following.

# Mandatory keys
- `domain`: Domain where to solve the MHD problem.
  The domain is represented with a `Triangulation` or with a `Integer`/`String`
  tag in the underlying discrete model.
  It can also be a `Gridap.DiscreteModel` or `GridapDistributed.DistributedDiscreteModel`
  if `params[:solid]` is `nothing`. In this case,
  one needs to guarantee `params[:model]===params[:fluid][:domain]`.
-  `:α`: Value of the parameter `α`.
-  `:β`: Value of the parameter `β`.
-  `:γ`: Value of the parameter `γ`.
-  `:B`: Value of the parameter `B`.

# Optional keys
-  `:f=>VectorValue(0,0,0)`: Value of the parameter `f`.
-  `:σ=>1`: Value of the parameter `σ`.
"""
function params_fluid(params::Dict{Symbol,Any})
  mandatory = Dict(
   :domain=>true,
   :α=>true,
   :β=>true,
   :γ=>true,
   :B=>true,
   :f=>false,
   :σ=>false,
  )
  optional = Dict(:σ=>1,:f=>VectorValue(0,0,0))
  fluid = _check_mandatory_and_add_optional(params[:fluid],mandatory,optional,params,"[:fluid]")
  fluid
end

"""
Valid keys for `params[:solid]` are the following

# Mandatory keys
- `domain`: Domain occupied by the solid.
  The domain is represented either with a `Triangulation` or with a `Integer`/`String`
  tag in the underlying discrete model.
# Optional keys
-  `:σ=>1`: Value of the parameter `σ`.
"""
function params_solid(params::Dict{Symbol,Any})
  mandatory = Dict(
   :domain=>true,
   :σ=>false,
  )
  optional = Dict(:σ=>1)
  solid = _check_mandatory_and_add_optional(params[:solid],mandatory,optional,params,"[:solid]")
  solid
end

"""
Valid keys for `params[:bcs]` are the following

# Mandatory keys
- `:u`: A `Dict` defining strong Dirichlet conditions on the fluid velocity.
   See  [`params_bcs_u`](@ref) for further details.
- `:j`: A `Dict` defining strong Dirichlet conditions on the charge current.
   See  [`params_bcs_j`](@ref) for further details.

# Optional keys
- `:φ=>[]`: A `Dict` or a vector of `Dict`s.
   Each `Dict` defines a weak boundary condition for the potential.
   See  [`params_bcs_φ`](@ref) for further details.
- `:t=>[]`: A `Dict` or a vector of `Dict`s.
  Each `Dict` defines a boundary traction on the fluid.
   See  [`params_bcs_t`](@ref) for further details.
- `:thin_wall=>[]`: A `Dict` or a vector of `Dict`s.
  Each Dict defines a thin wall law.
   See  [`params_bcs_thin_wall`](@ref) for further details.
"""
function params_bcs(params)
  mandatory = Dict(
   :u=>true,
   :j=>true,
   :φ=>false,
   :t=>false,
   :thin_wall=>false,
   :f => false,
   :B => false,
  )
  optional = Dict(
   :φ=>[],
   :t=>[],
   :thin_wall=>[],
   :f =>[],
   :B =>[],
  )
  bcs = _check_mandatory_and_add_optional(params[:bcs],mandatory,optional,params,"[:bcs]")
  # Sub params
  bcs[:u] = params_bcs_u(params)
  bcs[:j] = params_bcs_j(params)
  if bcs[:φ] !== optional[:φ]
    bcs[:φ] = params_bcs_φ(params)
  end
  if bcs[:t] !== optional[:t]
    bcs[:t] = params_bcs_t(params)
  end
  if bcs[:thin_wall] !== optional[:thin_wall]
    bcs[:thin_wall] = params_bcs_thin_wall(params)
  end
  bcs
end

"""
Valid keys for `params[:bcs][:u]` are the following

# Mandatory keys
- `:tags`: Dirichlet tags where to impose strong boundary conditions for the fluid velocity.

# Optional keys
- `:values => zero_values(prams[:bcs][:u][:tags])`: The fluid velocity value or function
   to be imposed at each of the given tags.
"""
function params_bcs_u(params::Dict{Symbol,Any})
  mandatory = Dict(
   :tags=>true,
   :values=>false,
  )
  _check_mandatory(params[:bcs][:u],mandatory,"[:bcs][:u]")
  optional = Dict(
    :values=>zero_values(params[:bcs][:u][:tags]),
  )
  u = _add_optional(params[:bcs][:u],mandatory,optional,params,"[:bcs][:u]")
  _check_unused(u,mandatory,params,"[:bcs][:u]")
  u
end

zero_values(tags) = VectorValue(0,0,0)
zero_values(tags::Vector) = map(zero_values,tags)

"""
Valid keys for `params[:bcs][:j]` are the following

# Mandatory keys
- `:tags`: Dirichlet tags where to impose strong boundary conditions for the charge current in normal direction.

# Optional keys
- `:values => zero_values(params[:bcs][:j][:tags])`: The charge current value or function
   to be imposed at each of the given tags.
"""
function params_bcs_j(params::Dict{Symbol,Any})
  mandatory = Dict(
   :tags=>true,
   :values=>false,
  )
  _check_mandatory(params[:bcs][:j],mandatory,"[:bcs][:j]")
  optional = Dict(
    :values=>zero_values(params[:bcs][:j][:tags]),
  )
  j = _add_optional(params[:bcs][:j],mandatory,optional,params,"[:bcs][:j]")
  _check_unused(j,mandatory,params,"[:bcs][:j]")
  j
end

"""
Valid keys for the dictionaries in `params[:bcs][:φ]` are the following.

# Mandatory keys
- `:domain`: Domain where to impose the potential weakly.
  The domain is represented either with a `Triangulation` or with a `Integer`/`String`
  tag in the underlying discrete model.
- `:value`: Value of the electric potential to be imposed weakly.
"""
function params_bcs_φ(params::Dict{Symbol,Any})
  mandatory = Dict(
   :domain=>true,
   :value=>true,
  )
  optional = Dict()
  _check_mandatory_and_add_optional_weak(params[:bcs][:φ],mandatory,optional,params,"[:bcs][:φ]")
end

"""
Valid keys for the dictionaries in `params[:bcs][:t]` are the following.

# Mandatory keys
- `:domain`: Domain where to impose the fluid boundary traction weakly.
  The domain is represented either with a `Triangulation` or with a `Integer`/`String`
  tag in the underlying discrete model.
- `:value`: Value of the fluid traction to be imposed weakly.
"""
function params_bcs_t(params::Dict{Symbol,Any})
  mandatory = Dict(
   :domain=>true,
   :value=>true,
  )
  optional = Dict()
  _check_mandatory_and_add_optional_weak(params[:bcs][:t],mandatory,optional,params,"[:bcs][:t]")
end

"""
Valid keys for the dictionaries in `params[:bcs][:thin_wall]` are the following.

The thin wall law is

    j⋅n + cw*n⋅∇(j)⋅n = jw

where `j` is unknown, `n`  is the boundary outward normal, and `cw,jw` are parameters.
The thin wall law is imposed weakly via a penalty parameter `τ`.

# Mandatory keys
- `:domain`: Domain where to impose the thin wall law.
  The domain is represented either with a `Triangulation` or with a `Integer`/`String`
  tag in the underlying discrete model.
- `:cw`: Value of the parameter `cw`.
- `:τ`: Value of the parameter `τ`.

# Optional keys
- `:jw=>0`: Value of the parameter `jw`.
"""
function params_bcs_thin_wall(params::Dict{Symbol,Any})
  mandatory = Dict(
   :domain=>true,
   :cw=>true,
   :τ=>true,
   :jw=>false,
  )
  optional = Dict(:jw=>0)
  _check_mandatory_and_add_optional_weak([:bcs][:thin_wall],mandatory,optional,params,"[:bcs][:thin_wall]")
end

function _check_mandatory(params,mandatory,p)
  if !isa(params,Dict{Symbol})
    error("The params$p has to be a Dict{Symbol}")
  end
  for key in keys(mandatory)
    if mandatory[key] && !haskey(params,key)
      error("Key :$key is a mandatory key in params$p, but it is not provided.")
    end
  end
end

function _add_optional(_subparams,mandatory,optional,params,p)
  # New dict
  subparams = Dict{Symbol,Any}()
  merge!(subparams,_subparams)
  # Check that we have computed defaults for all optionals
  for key in keys(mandatory)
    if !mandatory[key] && !haskey(optional,key)
      error("Internal error")
    end
  end
  # Set default args
  for key in keys(optional)
    if !haskey(subparams,key)
      subparams[key] = optional[key]
    end
  end
  subparams
end

function _check_unused(subparams,mandatory,params,p)
  if params[:check_valid]
    for key in keys(subparams)
      if !haskey(mandatory,key)
        error("Key :$key is not a valid key in params$p. Set params[:check_valid] = false to ignore invalid keys.")
      end
    end
  end
end

function _check_mandatory_and_add_optional(_subparams,mandatory,optional,params,p)
  _check_mandatory(_subparams,mandatory,p)
  subparams = _add_optional(_subparams,mandatory,optional,params,p)
  _check_unused(subparams,mandatory,params,p)
  subparams
end

function _check_mandatory_and_add_optional_weak(t,mandatory,optional,params,p)
  function _params_bcs(_t)
    error("The value params$p has to be a Dict{Symbol} or a vector of Dict{Symbol}.")
  end
  function _params_bcs(_t::AbstractVector)
    map(i->_params_bcs(_t[i],"[$i]"),1:length(_t))
  end
  function _params_bcs(_t::Dict{Symbol},i="")
    _check_mandatory_and_add_optional(_t,mandatory,optional,params,"$p$i")
  end
  r = _params_bcs(t)
  isa(r,Dict) ? [r] : r
end

"""
    main(params::Dict{Symbol}) -> xh, full_params

Solve the MHD equations. Function `main` takes a dictionary `params`
that contains several parameters defining a MHD problem (see below)
and returns `xh` and `full_params`.
Object `xh` represents the solution of the MHD problem and it is of type `Gridap.MultiField.MultiFieldFEFunction`,
which can be unpacked to get access to the different fields of the MHD solution
(fluid velocity, fluid pressure, charge current, and electric potential respectively).
One can further post process these quantities using the tools provided by Gridap.
On the other hand `full_params` is a copy of `params` augmented with the default values
used in the computation.

In the fluid domain, the equations being solved are

     ∇⋅j = 0
     j + σ*∇(φ) - σ*(u×B) = 0
     ∇⋅u = 0
     α*u⋅∇(u) - β*Δ(u) + ∇(p) - γ*(j×B) = f

where `u,p,j,φ` are the unknowns and `α,β,γ,σ,B,f` are (possibly spatially-dependent) parameters.  In the solid domain, we solve

     ∇⋅j = 0
     j + σ*∇(φ) = 0

where `j,φ` are the unknowns and `σ` is a (possibly spatially-dependent) parameter.
These equations are augmented with suitable boundary conditions (see below).

The MHD problem is customized by setting the following keys in `params`.

# Mandatory keys
- `:model`: The finite element mesh representing the union of the fluid and solid domains. It should be either a `Gridap.DiscreteModel`
or a `GridapDistributed.DistributedDiscreteModel`
- `:fluid`: A `Dict` defining the fluid domain and fluid parameters.
   See [`params_fluid`](@ref) for further details.
- `:bcs`: A `Dict` defining the boundary conditions and other external loads.
  See [`params_bcs`](@ref) for further details.

# Optional keys
- `:solid => nothing`:
  A `Dict` defining the solid domain and solid parameters.
  If not provided or set to `nothing` the solid domain is not taken into account.
  See [`params_solid`](@ref) for further details.
- `:k => 2`:
  Maximum interpolation order (i.e., the order used for the fluid velocity).
- `:solver => default_solver()`:
  Nonlinear solver to compute the solution.
    It should be an instance of some type implementing the `NonlinearSolver` interface of Gridap.
- `:debug => false`:
  If true, setup the problem, but do not solve it. Otherwise, solve it.
- `:check_valid => true`: If `true`, check that all given keys are valid. Otherwise, silently ignore invalid keys.
- `:matrix_type => SparseMatrixCSC{Float64,Int}`:
   Matrix type to assemble the problem.
- `:vector_type => Vector{Float64}`:
  Vector type to assemble the problem.
- `:ptimer => default_ptimer(params[:model])`:
  Instance of `PTimer` used to monitor times. New time measurements are added to the given timer.
"""
function main(_params::Dict)

  params = add_default_params(_params)

  t = params[:ptimer]
  tic!(t;barrier=true)

  # ReferenceFEs
  k::Int = params[:k]
  T = Float64
  model = params[:model]
  D = num_cell_dims(model)
  reffe_u = ReferenceFE(lagrangian,VectorValue{D,T},k)
  reffe_p = ReferenceFE(lagrangian,T,k-1;space=:P)
  reffe_j = ReferenceFE(raviart_thomas,T,k-1)
  reffe_φ = ReferenceFE(lagrangian,T,k-1)

  # Test spaces
  Ωf = _fluid_mesh(model,params[:fluid][:domain])
  V_u = TestFESpace(Ωf,reffe_u;dirichlet_tags=params[:bcs][:u][:tags])
  V_p = TestFESpace(Ωf,reffe_p;conformity=p_conformity(Ωf))
  V_j = TestFESpace(model,reffe_j;dirichlet_tags=params[:bcs][:j][:tags])
  V_φ = TestFESpace(model,reffe_φ;conformity=:L2)
  V = MultiFieldFESpace([V_u,V_p,V_j,V_φ])

  # Trial spaces
  # TODO improve for parallel computations
  tic!(t;barrier=true)
  z = zero(VectorValue{D,Float64})
  u_bc = params[:bcs][:u][:values]
  j_bc = params[:bcs][:j][:values]
  U_u = u_bc == z ? V_u : TrialFESpace(V_u,u_bc)
  U_j = j_bc == z ? V_j : TrialFESpace(V_j,j_bc)
  U_p = TrialFESpace(V_p)
  U_φ = TrialFESpace(V_φ)
  U = MultiFieldFESpace([U_u,U_p,U_j,U_φ])
  toc!(t,"fe_spaces")

  tic!(t;barrier=true)
  if params[:debug]
    Random.seed!(1234)
    vt = get_vector_type(U)
    free_ids = get_free_dof_ids(U)
    free_vals = _rand(vt,free_ids)
    xh = FEFunction(U,free_vals)
  else
    res, jac = weak_form(params,k)
    Tm = params[:matrix_type]
    Tv = params[:vector_type]
    assem = SparseMatrixAssembler(Tm,Tv,U,V)
    op = FEOperator(res,jac,U,V,assem)
    solver = params[:solver]
    xh = zero(U)
    xh,cache = solve!(xh,solver,op)
    solver_postpro = params[:solver_postpro]
    solver_postpro(cache)
  end
  toc!(t,"solve")

  xh, params
end

function _fluid_mesh(
  model,domain::Union{Gridap.DiscreteModel,GridapDistributed.DistributedDiscreteModel})
  msg = "params[:fluid][:domain] is a discrete model, but params[:fluid][:domain]===params[:model] is not true."
  @assert model === domain msg
  domain
end

function _fluid_mesh(
  model,
  domain::Union{Gridap.Triangulation,GridapDistributed.DistributedTriangulation})
  domain
end

function _fluid_mesh( model, domain)
  Interior(model,tags=domain)
end

function _interior(
  model,domain::Union{Gridap.DiscreteModel,GridapDistributed.DistributedDiscreteModel})
  Interior(domain)
end

function _interior(
  model,
  domain::Union{Gridap.Triangulation,GridapDistributed.DistributedTriangulation})
  domain
end

function _interior( model, domain)
  Interior(model,tags=domain)
end

function _boundary(
  model,
  domain::Union{Gridap.Triangulation,GridapDistributed.DistributedTriangulation})
  domain
end

function _boundary( model, domain)
  Boundary(model,tags=domain)
end

function _rand(vt::Type{<:Vector{T}},r::AbstractUnitRange) where T
  rand(T,length(r))
end

function _rand(vt::Type{<:PVector{T,A}},ids::PRange) where {T,A}
  values = map_parts(ids.partition) do partition
    Tv = eltype(A)
    _rand(Tv,1:num_lids(partition))
  end
  PVector(values,ids)
end

function p_conformity(poly::Polytope)
  if is_simplex(poly)
    conf = :H1
  elseif is_n_cube(poly)
    conf = :L2
  else
    @unreachable "unsupported cell topology"
  end
  conf
end

function p_conformity(Ω::Triangulation)
  reffes = get_reffes(Ω)
  @assert length(reffes) == 1
  reffe = first(reffes)
  poly = get_polytope(reffe)
  p_conformity(poly)
end

function p_conformity(Ω::GridapDistributed.DistributedTriangulation)
  p = map_parts(Ω.trians) do Ω
    reffes = get_reffes(Ω)
    @assert length(reffes) == 1
    reffe = first(reffes)
    poly = get_polytope(reffe)
    poly
  end
  poly = get_part(p) # We assume same polytope in all parts
  p_conformity(poly)
end

function p_conformity(model::DiscreteModel)
  Ω = Interior(model)
  p_conformity(Ω)
end

function p_conformity(model::GridapDistributed.DistributedDiscreteModel)
  Ω = Interior(model)
  p_conformity(Ω)
end

function add_defaults!(params,defaults)
  for (key,val) in defaults
    if !haskey(params,key)
      params[key] = val
    elseif isa(val,AbstractDict)
      @assert isa(params[key],AbstractDict)
      add_defaults!(params[key],val)
    end
  end
end

function weak_form(params,k)

  fluid = params[:fluid]

  Ωf = _interior(params[:model],fluid[:domain])
  dΩf = Measure(Ωf,2*k)

  solid = params[:solid]
  if solid !== nothing
    Ωs = _interior(params[:model],solid[:domain])
    dΩs = Measure(Ωs,2*k)
    σs = solid[:σ]
  end

  α = fluid[:α]
  β = fluid[:β]
  γ = fluid[:γ]
  f = fluid[:f]
  B = fluid[:B]
  σf = fluid[:σ]

  bcs = params[:bcs]

  params_φ = []
  for i in 1:length(bcs[:φ])
    φ_i = bcs[:φ][i][:value]
    Γ = _boundary(params[:model],bcs[:φ][i][:domain])
    dΓ = Measure(Γ,2*k)
    n_Γ = get_normal_vector(Γ)
    push!(params_φ,(φ_i,n_Γ,dΓ))
  end

  params_thin_wall = []
  for i in 1:length(bcs[:thin_wall])
    τ_i = bcs[:thin_wall][i][:τ]
    cw_i = bcs[:thin_wall][i][:cw]
    jw_i = bcs[:thin_wall][i][:jw]
    Γ = _boundary(params[:model],bcs[:thin_wall][i][:domain])
    dΓ = Measure(Γ,2*k)
    n_Γ = get_normal_vector(Γ)
    push!(params_thin_wall,(τ_i,cw_i,jw_i,n_Γ,dΓ))
  end

  if length(bcs[:t]) != 0
    error("Boundary tranction not yet implemented")
  end

  params_f = []
  for i in 1:length(bcs[:f])
    f_i = bcs[:f][i][:value]
    Ω_i = _interior(params[:model],bcs[:f][i][:domain])
    dΩ_i = Measure(Ω_i,2*k)
    push!(params_f,(f_i,dΩ_i))
  end

  params_B = []
  for i in 1:length(bcs[:B])
    B_i = bcs[:B][i][:value]
    Ω_i = _interior(params[:model],bcs[:B][i][:domain])
    dΩ_i = Measure(Ω_i,2*k)
    push!(params_f,(γ,B_i,dΩ_i))
  end

  function a(x,dy)
    r = a_mhd(x,dy,β,γ,B,σf,dΩf)
    for p in params_thin_wall
      r = r + a_thin_wall(x,dy,p...)
    end
    for p in params_B
      r = r + a_B(x,dy,p...)
    end
    if solid !== nothing
      r = r + a_solid(x,dy,σs,dΩs)
    end
    r
  end

  function ℓ(dy)
    r = ℓ_mhd(dy,f,dΩf)
    for p in params_φ
      r = r + ℓ_φ(dy,p...)
    end
    for p in params_thin_wall
      r = r + ℓ_thin_wall(dy,p...)
    end
    for p in params_f
      r = r + ℓ_f(dy,p...)
    end
    r
  end

  function c(x,dy)
    r = c_mhd(x,dy,α,dΩf)
    r
  end

  function dc(x,dx,dy)
    r = dc_mhd(x,dx,dy,α,dΩf)
    r
  end

  res(x,dy) = c(x,dy) + a(x,dy) - ℓ(dy)
  jac(x,dx,dy) = dc(x,dx,dy) + a(dx,dy)

  res, jac
end

conv(u,∇u) = (∇u')⋅u

function a_solid(x,dy,σ,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( j⋅v_j - σ*φ*(∇⋅v_j) + (∇⋅j)*v_φ)dΩ
end

function a_mhd(x,dy,β,γ,B,σ,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫(
    β*(∇(u)⊙∇(v_u)) - p*(∇⋅v_u) -(γ*(j×B)⋅v_u) +
    (∇⋅u)*v_p +
    j⋅v_j - σ*φ*(∇⋅v_j) - σ*(u×B)⋅v_j +
    (∇⋅j)*v_φ ) * dΩ
end

function ℓ_mhd(dy,f,dΩ)
  v_u, v_p, v_j, v_φ = dy
  ∫( v_u⋅f )*dΩ
end

function c_mhd(x,dy,α,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( α*v_u⋅(conv∘(u,∇(u))) ) * dΩ
end

function dc_mhd(x,dx,dy,α,dΩ)
  u, p, j, φ = x
  du , dp , dj , dφ  = dx
  v_u, v_p, v_j, v_φ = dy
  ∫( α*v_u⋅( (conv∘(u,∇(du))) + (conv∘(du,∇(u))) ) ) * dΩ
end

function ℓ_φ(dy,φ,n_Γ,dΓ)
  v_u, v_p, v_j, v_φ = dy
  ∫( -(v_j⋅n_Γ)*φ )*dΓ
end

function ℓ_f(dy,f,dΩ)
  v_u, v_p, v_j, v_φ = dy
  ∫( v_u⋅f )*dΩ
end

function ℓ_thin_wall(dy,τ,cw,jw,n_Γ,dΓ)
  v_u, v_p, v_j, v_φ = dy
  ∫( τ*(v_j⋅n_Γ)*jw ) * dΓ
end

function a_thin_wall(x,dy,τ,cw,jw,n_Γ,dΓ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( τ*((v_j⋅n_Γ)*(j⋅n_Γ) + cw*(v_j⋅n_Γ)*(n_Γ⋅(∇(j)⋅n_Γ))) )*dΓ
end

function a_B(x,dy,γ,B,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( -(γ*(j×B)⋅v_u) - (u×B)⋅v_j )*dΩ
end

