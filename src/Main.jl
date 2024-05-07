
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

  # FESpaces
  tic!(t;barrier=true)
  U, V = _fe_spaces(params)
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
    op = _fe_operator(U,V,params)
    xh = _allocate_solution(op)
    if params[:solve]
      solver = _solver(op,params)
      toc!(t,"solver_setup")
      tic!(t;barrier=true)
      xh,cache = _solve(xh,solver,op,params)
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

# Solver

_solver(op,params) = _solver(Val(params[:solver][:solver]),Val(params[:transient]),op,params)
_solver(val::Val,::Val{false},op,params) = _solver(val,op,params)

#_solver(::Val{:julia},op,params) = NLSolver(show_trace=true,method=:newton)
_solver(::Val{:julia},op,params) = GridapSolvers.NewtonSolver(LUSolver(),maxiter=10,rtol=1.e-6,verbose=true)
_solver(::Val{:petsc},op,params) = PETScNonlinearSolver()
_solver(::Val{:li2019},op,params) = Li2019Solver(op,params)

function _solver(val::Val,::Val{true},op,params)
  solver = _solver(val,op,params)
  _ode_solver(solver,params)
end

_ode_solver(solver,params) = _ode_solver(Val(params[:ode][:solver]),solver,params)

function _ode_solver(::Val{:theta},solver,params)
  Δt = params[:ode][:Δt]
  θ = params[:ode][:solver_params][:θ]
  ThetaMethod(solver,Δt,θ)
end

function _ode_solver(::Val{:forward},solver,params)
  Δt = params[:ode][:Δt]
  ForwardEuler(solver,Δt)
end

# MultiFieldStyle

_multi_field_style(params) = _multi_field_style(Val(params[:solver][:solver]))
_multi_field_style(::Val{:julia}) = ConsecutiveMultiFieldStyle()
_multi_field_style(::Val{:petsc}) = ConsecutiveMultiFieldStyle()
_multi_field_style(::Val{:li2019}) = BlockMultiFieldStyle(4,(1,1,1,1),(3,1,2,4)) # (j,u,p,φ)

# FESpaces

function _fe_spaces(params)
  k = params[:fespaces][:k]
  T = Float64
  model = params[:model]

  # ReferenceFEs
  D = num_cell_dims(model)
  reffe_u = ReferenceFE(lagrangian,VectorValue{D,T},k)
  reffe_p = ReferenceFE(lagrangian,T,k-1;space=params[:fespaces][:p_space])
  reffe_j = ReferenceFE(raviart_thomas,T,k-1)
  reffe_φ = ReferenceFE(lagrangian,T,k-1)

  # Test spaces
  mfs = _multi_field_style(params)
  Ωf  = _fluid_mesh(model,params[:fluid][:domain])
  V_u = TestFESpace(Ωf,reffe_u;dirichlet_tags=params[:bcs][:u][:tags])
  V_p = TestFESpace(Ωf,reffe_p;
    conformity=p_conformity(model,params[:fespaces]),
    constraint=params[:fespaces][:p_constraint])
  V_j = TestFESpace(model,reffe_j;dirichlet_tags=params[:bcs][:j][:tags])
  V_φ = TestFESpace(model,reffe_φ;conformity=:L2)
  V = MultiFieldFESpace([V_u,V_p,V_j,V_φ];style=mfs)

  # Trial spaces
  z = zero(VectorValue{D,Float64})
  u_bc = params[:bcs][:u][:values]
  j_bc = params[:bcs][:j][:values]
  if !params[:transient]
    U_u = u_bc == z ? V_u : TrialFESpace(V_u,u_bc)
    U_j = j_bc == z ? V_j : TrialFESpace(V_j,j_bc)
    U_p = TrialFESpace(V_p)
    U_φ = TrialFESpace(V_φ)
    U = MultiFieldFESpace([U_u,U_p,U_j,U_φ];style=mfs)
  else
    U_u = u_bc == z ? V_u : TransientTrialFESpace(V_u,u_bc)
    U_j = j_bc == z ? V_j : TransientTrialFESpace(V_j,j_bc)
    U_p = TransientTrialFESpace(V_p)
    U_φ = TransientTrialFESpace(V_φ)

    U = TransientMultiFieldFESpace([U_u,U_p,U_j,U_φ];style=mfs)
  end

  return U, V
end

# FEOperator

_fe_operator(U,V,params) = __fe_operator(_multi_field_style(params),U,V,params)

function __fe_operator(T,U,V,params)
  if ! params[:transient]
    _fe_operator(T,U,V,params)
  else
   _ode_fe_operator(T,U,V,params)
  end
end

function _fe_operator(::ConsecutiveMultiFieldStyle,U,V,params)
  k = params[:fespaces][:k]
  res, jac = weak_form(params,k)
  Tm = params[:solver][:matrix_type]
  Tv = params[:solver][:vector_type]
  assem = SparseMatrixAssembler(Tm,Tv,U,V)
  return FEOperator(res,jac,U,V,assem)
end

function _fe_operator(::BlockMultiFieldStyle,U,V,params)
  # TODO: BlockFEOperator, which only updates nonlinear blocks (only important for high Re)
  k = params[:fespaces][:k]
  res, jac = weak_form(params,k)
  Tm = params[:solver][:matrix_type]
  Tv = params[:solver][:vector_type]
  assem = SparseMatrixAssembler(Tm,Tv,U,V)
  return FEOperator(res,jac,U,V,assem)
end

function _ode_fe_operator(::ConsecutiveMultiFieldStyle,U,V,params)
  k = params[:fespaces][:k]
  res, jac, jac_t = weak_form(params,k)
  Tm = params[:solver][:matrix_type]
  Tv = params[:solver][:vector_type]
  assem = SparseMatrixAssembler(Tm,Tv,U(0),V(0))
  return TransientFEOperator(res,(jac,jac_t),U,V,assembler=assem)
end

# Sub-triangulations

const DiscreteModelTypes = Union{Gridap.DiscreteModel,GridapDistributed.DistributedDiscreteModel}
const TriangulationTypes = Union{Gridap.Triangulation,GridapDistributed.DistributedTriangulation}

function _fluid_mesh(model,domain::DiscreteModelTypes)
  msg = "params[:fluid][:domain] is a discrete model, but params[:fluid][:domain]===params[:model] is not true."
  @assert model === domain msg
  return domain
end
_fluid_mesh(model,domain::TriangulationTypes) = domain
_fluid_mesh(model,domain::Nothing) = model # This should be removed, but Gridap needs fixes
_fluid_mesh(model,domain) = Interior(model,tags=domain)

_interior(model,domain::DiscreteModelTypes) = Interior(domain)
_interior(model,domain::TriangulationTypes) = domain
_interior(model,domain::Nothing) = Triangulation(model) # This should be removed, but Gridap needs fixes
_interior(model,domain) = Interior(model,tags=domain)

_boundary(model,domain::TriangulationTypes) = domain
_boundary(model,domain) = Boundary(model,tags=domain)

_skeleton(model,domain::TriangulationTypes) = SkeletonTriangulation(domain)
_skeleton(model,domain::Nothing) = SkeletonTriangulation(model)
_skeleton(model,domain) = _skeleton(model,_interior(model,domain))

# Random vector generation

function _rand(vt::Type{<:Vector{T}},r::AbstractUnitRange) where T
  rand(T,length(r))
end

function _rand(vt::Type{<:PVector{VT,A}},ids::PRange) where {VT,A}
  T = eltype(VT)
  prand(T,partition(ids))
end

# Solve

_solve(xh,solver,op,params) = _solve(Val(params[:transient]),xh,solver,op,params)

_solve(::Val{false},xh,solver,op,params) = solve!(xh,solver,op)

function _solve(::Val{true},xh,solver,op,params)
  xh0 = initial_value(op,params)
  t0,tf = time_interval(params)
  cache = nothing
  solve(solver,op,t0,tf,xh0), cache
end

initial_value(op,params) = initial_value(Val(params[:ode][:U0]),op,params)

initial_value(::Val{:zero},op,params) = zero(get_trial(op))
initial_value(::Val{:solve},op,params) = @notimplemented

function initial_value(::Val{:value},op,params)
  t0 = params[:ode][:t0]
  U0 = get_trial(op)(t0)
  v = params[:ode][:initial_values]
  interpolate([v[:u],v[:p],v[:j],v[:φ]],U0)
end

time_interval(params) = ( params[:ode][:t0], params[:ode][:tf] )

function _allocate_solution(op::FEOperator)
  zero(get_trial(op))
end

function _allocate_solution(op::TransientFEOperator)
  nothing
end

function _get_cell_size(t::Triangulation)
  meas = get_cell_measure(t)
  d = num_dims(t)
  map(m->m^(1/d),meas)
end

function _get_cell_size(t::GridapDistributed.DistributedTriangulation)
  map(_get_cell_size,local_views(t))
end
