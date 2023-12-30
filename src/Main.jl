
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
    xh = zero(get_trial(op))
    if params[:solve]
      solver = _solver(op,params)
      toc!(t,"solver_setup")
      tic!(t;barrier=true)
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

# Solver

_solver(op,params) = _solver(Val(params[:solver][:solver]),op,params)
_solver(::Val{:julia},op,params) = NLSolver(show_trace=true,method=:newton)
_solver(::Val{:petsc},op,params) = PETScNonlinearSolver()
_solver(::Val{:li2019},op,params) = Li2019Solver(op,params)
_solver(::Val{:badia2024},op,params) = Badia2024Solver(op,params)

# MultiFieldStyle

_multi_field_style(params) = _multi_field_style(Val(params[:solver][:solver]))
_multi_field_style(::Val{:julia}) = ConsecutiveMultiFieldStyle()
_multi_field_style(::Val{:petsc}) = ConsecutiveMultiFieldStyle()
_multi_field_style(::Val{:li2019}) = BlockMultiFieldStyle(4,(1,1,1,1),(3,1,2,4)) # (j,u,p,φ)
_multi_field_style(::Val{:badia2024}) = BlockMultiFieldStyle(3,(2,1,1),(1,3,2,4)) # ([u,j],p,φ)

# FESpaces

_fe_spaces(params) = _fe_spaces(Val(uses_multigrid(params[:solver])),params)

function _fe_spaces(::Val{false},params)
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
  V_p = TestFESpace(Ωf,reffe_p)
  V_j = TestFESpace(model,reffe_j;dirichlet_tags=params[:bcs][:j][:tags])
  V_φ = TestFESpace(model,reffe_φ;conformity=:L2)
  V = MultiFieldFESpace([V_u,V_p,V_j,V_φ];style=mfs)

  # Trial spaces
  tic!(t;barrier=true)
  z = zero(VectorValue{D,Float64})
  u_bc = params[:bcs][:u][:values]
  j_bc = params[:bcs][:j][:values]
  U_u = u_bc == z ? V_u : TrialFESpace(V_u,u_bc)
  U_j = j_bc == z ? V_j : TrialFESpace(V_j,j_bc)
  U_p = TrialFESpace(V_p)
  U_φ = TrialFESpace(V_φ)
  U = MultiFieldFESpace([U_u,U_p,U_j,U_φ];style=mfs)

  return U, V
end

function _fe_spaces(::Val{true},params)
  # TODO: Add fluid/solid mesh support
  k = params[:fespaces][:k]
  T = Float64
  model = params[:model]
  mh    = params[:multigrid][:mh]
  Ωf    = _fluid_mesh(model,params[:fluid][:domain])
  @assert get_model(mh,1) == model

  uses_mg = space_uses_multigrid(params[:solver])
  trians  = map((m,a,b) -> m ? a : b,uses_mg,[mh,mh,mh,mh],[Ωf,model,Ωf,model])

  # ReferenceFEs
  D = num_cell_dims(model)
  reffe_u = ReferenceFE(lagrangian,VectorValue{D,T},k)
  reffe_p = ReferenceFE(lagrangian,T,k-1;space=params[:fespaces][:p_space])
  reffe_j = ReferenceFE(raviart_thomas,T,k-1)
  reffe_φ = ReferenceFE(lagrangian,T,k-1)

  # Test spaces
  mfs = _multi_field_style(params)
  V_u = TestFESpace(trians[1],reffe_u;dirichlet_tags=params[:bcs][:u][:tags])
  V_j = TestFESpace(trians[2],reffe_j;dirichlet_tags=params[:bcs][:j][:tags])
  V_p = TestFESpace(trians[3],reffe_p)
  V_φ = TestFESpace(trians[4],reffe_φ;conformity=:L2)
  
  # Trial spaces
  z = zero(VectorValue{D,Float64})
  u_bc = params[:bcs][:u][:values]
  j_bc = params[:bcs][:j][:values]
  U_u = (u_bc == z) ? V_u : TrialFESpace(V_u,u_bc)
  U_j = (j_bc == z) ? V_j : TrialFESpace(V_j,j_bc)
  U_p = TrialFESpace(V_p)
  U_φ = TrialFESpace(V_φ)

  # Sort spaces
  trials, tests, sh_trials, sh_tests = map(uses_mg,[U_u,U_p,U_j,U_φ],[V_u,V_p,V_j,V_φ]) do m,trial,test
    if m
      GridapSolvers.get_fe_space(trial,1), GridapSolvers.get_fe_space(test,1), trial, test
    else
      trial, test, nothing, nothing
    end
  end

  params[:multigrid][:trials] = sh_trials
  params[:multigrid][:tests]  = sh_tests
  V = MultiFieldFESpace(tests;style=mfs)
  U = MultiFieldFESpace(trials;style=mfs)
  return U, V
end

# FEOperator

_fe_operator(U,V,params) = _fe_operator(_multi_field_style(params),U,V,params)

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

# Sub-triangulations

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

# Random vector generation

function _rand(vt::Type{<:Vector{T}},r::AbstractUnitRange) where T
  rand(T,length(r))
end

function _rand(vt::Type{<:PVector{VT,A}},ids::PRange) where {VT,A}
  T = eltype(VT)
  prand(T,partition(ids))
end
