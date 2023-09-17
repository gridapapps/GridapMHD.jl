
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
function main(_params::Dict;output::Dict=Dict{Symbol,Any}())

  params = add_default_params(_params)

  t = params[:ptimer]
  tic!(t;barrier=true)

  # ReferenceFEs
  k::Int = params[:fespaces][:k]
  T = Float64
  model = params[:model]
  D = num_cell_dims(model)
  reffe_u = ReferenceFE(lagrangian,VectorValue{D,T},k)
  reffe_p = ReferenceFE(lagrangian,T,k-1;space=params[:fespaces][:p_space])
  reffe_j = ReferenceFE(raviart_thomas,T,k-1)
  reffe_φ = ReferenceFE(lagrangian,T,k-1)

  # Test spaces
  mfs = _multi_field_style(params)
  Ωf  = _fluid_mesh(model,params[:fluid][:domain])
  V_u = TestFESpace(Ωf,reffe_u;dirichlet_tags=params[:bcs][:u][:tags])
  V_p = TestFESpace(Ωf,reffe_p)
  V_j = TestFESpace(model,reffe_j;dirichlet_tags=params[:bcs][:j][:tags])
  V_φ = TestFESpace(model,reffe_φ;conformity=:L2)
  V = MultiFieldFESpace([V_u,V_p,V_j,V_φ];style=mfs)

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
  U = MultiFieldFESpace([U_u,U_p,U_j,U_φ];style=mfs)
  toc!(t,"fe_spaces")

  tic!(t;barrier=true)
  if params[:debug]
    Random.seed!(1234)
    vt = get_vector_type(U)
    free_ids = get_free_dof_ids(U)
    free_vals = _rand(vt,free_ids)
    xh = FEFunction(U,free_vals)
    toc!(t,"solve")
  else
    op = _fe_operator(mfs,U,V,params)
    xh = zero(get_trial(op))
    if params[:solve]
      solver = _solver(op,params)
      xh,cache = solve!(xh,solver,op)
      solver_postpro = params[:solver][:solver_postpro]
      solver_postpro(cache,output)
      toc!(t,"solve")
      Base.finalize(solver)
    end
    if params[:res_assemble]
      tic!(t;barrier=true)
      r = residual(op,xh)
      toc!(t,"residual")
    end
    if params[:jac_assemble]
      tic!(t;barrier=true)
      j = jacobian(op,xh)
      toc!(t,"jacobian")
    end
  end

  xh, params, output
end

@inline function _solver(op,params)
  s = params[:solver][:solver]
  _solver(Val(s),op,params)
end
_solver(::Val{:julia},op,params) = NLSolver(show_trace=true,method=:newton)
_solver(::Val{:petsc},op,params) = PETScNonlinearSolver()
_solver(::Val{:block_gmres_li2019},op,params) = Li2019.Li2019Solver(op,params)

_multi_field_style(params) = _multi_field_style(Val(params[:solver][:solver]))
_multi_field_style(::Val{:julia}) = ConsecutiveMultiFieldStyle()
_multi_field_style(::Val{:petsc}) = ConsecutiveMultiFieldStyle()
_multi_field_style(::Val{:block_gmres_li2019}) = BlockMultiFieldStyle()

function _fe_operator(::ConsecutiveMultiFieldStyle,U,V,params)
  k = params[:fespaces][:k]
  res, jac = weak_form(params,k)
  Tm = params[:solver][:matrix_type]
  Tv = params[:solver][:vector_type]
  assem = SparseMatrixAssembler(Tm,Tv,U,V)
  return FEOperator(res,jac,U,V,assem)
end

function _fe_operator(::BlockMultiFieldStyle,U,V,params)
  k  = params[:fespaces][:k]
  Tm = params[:solver][:matrix_type]
  Tv = params[:solver][:vector_type]

  # Global operator
  res, jac = weak_form(params,k)
  assem = SparseMatrixAssembler(Tm,Tv,U,V)
  op_global = FEOperator(res,jac,U,V,assem)

  # u-u operator
  U_u, _, _, _ = U
  V_u, _, _, _ = V
  res_uu, jac_uu = weakform_uu(params,k)
  assem_uu = SparseMatrixAssembler(Tm,Tv,U_u,V_u)
  op_uu = FEOperator(res_uu,jac_uu,U_u,V_u,assem_uu)

  return FEOperatorMHD(op_global,op_uu)
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
  values = map(partition(ids)) do indices
    Tv = eltype(A)
    _rand(Tv,1:local_length(indices))
  end
  return PVector(values,partition(ids))
end

function weak_form(params,k)

  fluid = params[:fluid]
  Ωf  = _interior(params[:model],fluid[:domain])
  dΩf = Measure(Ωf,2*k)

  solid = params[:solid]
  if solid !== nothing
    Ωs  = _interior(params[:model],solid[:domain])
    dΩs = Measure(Ωs,2*k)
    σs  = solid[:σ]
  end

  α  = fluid[:α]
  β  = fluid[:β]
  γ  = fluid[:γ]
  f  = fluid[:f]
  B  = fluid[:B]
  σf = fluid[:σ]
  ζ  = params[:ζ]

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
    τ_i  = bcs[:thin_wall][i][:τ]
    cw_i = bcs[:thin_wall][i][:cw]
    jw_i = bcs[:thin_wall][i][:jw]
    Γ    = _boundary(params[:model],bcs[:thin_wall][i][:domain])
    dΓ   = Measure(Γ,2*k)
    n_Γ  = get_normal_vector(Γ)
    push!(params_thin_wall,(τ_i,cw_i,jw_i,n_Γ,dΓ))
  end

  if length(bcs[:t]) != 0
    error("Boundary tranction not yet implemented")
  end

  params_f = []
  for i in 1:length(bcs[:f])
    f_i  = bcs[:f][i][:value]
    Ω_i  = _interior(params[:model],bcs[:f][i][:domain])
    dΩ_i = Measure(Ω_i,2*k)
    push!(params_f,(f_i,dΩ_i))
  end

  params_B = []
  for i in 1:length(bcs[:B])
    B_i  = bcs[:B][i][:value]
    Ω_i  = _interior(params[:model],bcs[:B][i][:domain])
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
    if ζ !== nothing
      r = r + a_augmented_lagragian(x,dy,ζ,dΩf)
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

  return res, jac
end

conv(u,∇u) = (∇u')⋅u

function a_solid(x,dy,σ,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( j⋅v_j - σ*φ*(∇⋅v_j) + (∇⋅j)*v_φ)dΩ
end

function a_augmented_lagragian(x,dy,ζ,dΩf)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  return ∫( ζ*(∇⋅u)*(∇⋅v_u) ) * dΩf
end

function a_mhd(x,dy,β,γ,B,σ,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫(
    β*(∇(u)⊙∇(v_u)) - p*(∇⋅v_u) -(γ*(j×B)⋅v_u) -
    (∇⋅u)*v_p +
    γ*j⋅v_j - γ*σ*φ*(∇⋅v_j) - γ*σ*(u×B)⋅v_j +
    - γ*(∇⋅j)*v_φ ) * dΩ
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



