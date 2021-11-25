
# MHD equations
# ∇⋅u = 0
# u⋅∇(u) -ν*Δ(u) + (1/ρ)*∇(p) - (1/ρ)*(j×B) = (1/ρ)*f
# ∇⋅j = 0
# j + σ*∇(φ) - σ*(u×B) = 0
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
# ∇⋅ū = 0
# ∇⋅j̄ = 0
# j̄ + ∇(φ̄) - ū×B̄ = 0
# ū⋅∇(ū) -(1/Re)*Δ(ū) + ∇(p̄) - N*(j̄×B̄) = f̄ (Option 1, CFD)
# (1/N)*ū⋅∇(ū) - (1/Ha^2)*Δ(ū) + ∇(p̄) - (j̄×B̄) = f̄ (Option 2,MHD)
#
# In order to account for both options, the code solves these equations
#
# ∇⋅ū = 0
# ∇⋅j̄ = 0
# j̄ + ∇(φ̄) - ū×B̄ = 0
# α*ū⋅∇(ū) - β*Δ(ū) + ∇(p̄) - γ*(j̄×B̄) = f̄
#
# α = 1, β = (1/Re), γ = N (option 1,CFD)
# α = (1/N), β = (1/Ha^2), γ = 1 (option 2,MHD)
#

abstract type Action end
abstract type BoundaryAction <: Action end
abstract type BodyAction <: Action end

@with_kw struct ConductingFluid{A,B,C} <: BodyAction
  domain::A
  α::B
  β::C
end

@with_kw struct MagneticField{A,X,D} <: BodyAction
  domain::A
  B::X
  γ::D
end

@with_kw struct FluidForce{A,B} <: BodyAction
  domain::A
  f::B
end

# u = u
@with_kw struct VelocityBc{A,B} <: BoundaryAction
  domain::A
  u::B = VectorValue(0.,0.,0.)
end

# n⋅∇(u) - p*n = t
@with_kw struct TractionBc{A,B} <: BoundaryAction
  domain::A
  t::B
end

# j = j (imposed in normal direction only)
@with_kw struct CurrentBc{A,B} <: BoundaryAction
  domain::A
  j::B = VectorValue(0.,0.,0.)
end

# φ = φ
@with_kw struct PotentialBc{A,B} <: BoundaryAction
  domain::A
  φ::B = 0.0
end

# j⋅n + cw*n⋅∇(j)⋅n = jw
# imposed via a penalty of value τ
@with_kw struct ConductingThinWall{A,B,C,D} <: BoundaryAction
  domain::A
  cw::B
  jw::C = 0.0
  τ::D = 1.0
end

function main(
  model::DiscreteModel,
  actions::Vector{<:Action};
  debug=false,
  title="test",
  solver = NLSolver(show_trace=true,method=:newton))

  info = Dict{String,Float64}()

  @check count(i->isa(i,ConductingFluid),actions) == 1 "Only one instance of ConductingFluid allowed"
  ifluid = findall(i->isa(i,ConductingFluid),actions) |> first
  fluid = actions[ifluid]
  Ω = get_domain(model,fluid)
  u_tags, u_vals = find_strong_bcs_u(actions)
  j_tags, j_vals = find_strong_bcs_j(actions)
  info["cells fluid"] = num_cells(Ω)

  # Reference FES
  k = 2
  T = Float64
  reffe_u = ReferenceFE(lagrangian,VectorValue{3,T},k)
  reffe_p = ReferenceFE(lagrangian,T,k-1;space=:P)
  reffe_j = ReferenceFE(raviart_thomas,T,k-1)
  reffe_φ = ReferenceFE(lagrangian,T,k-1)

  # Test spaces
  V_u = TestFESpace(Ω,reffe_u;dirichlet_tags=u_tags)
  V_p = TestFESpace(Ω,reffe_p)
  V_j = TestFESpace(Ω,reffe_j;dirichlet_tags=j_tags)
  V_φ = TestFESpace(Ω,reffe_φ;conformity=:L2)
  V = MultiFieldFESpace([V_u,V_p,V_j,V_φ])
  info["dofs u"] = num_free_dofs(V_u)
  info["dofs p"] = num_free_dofs(V_p)
  info["dofs j"] = num_free_dofs(V_j)
  info["dofs φ"] = num_free_dofs(V_φ)
  info["dofs total"] = num_free_dofs(V)

  # Trial Spaces TODO improve for parallel computations
  U_u = TrialFESpace(V_u,u_vals)
  U_p = TrialFESpace(V_p)
  U_j = TrialFESpace(V_j,j_vals)
  U_φ = TrialFESpace(V_φ)
  U = MultiFieldFESpace([U_u,U_p,U_j,U_φ])

  if debug
    # Create some plot data
    for (i,action) in enumerate(actions)
      trian = get_domain(model,action)
      pn = propertynames(action)
      cellfields = [ string(p)=>getproperty(action,p) for p in pn if p != :domain ]
      writevtk(trian,"$(title)_action_$i",order=k,cellfields=cellfields)
    end
    # Create a random solution (useful to debug the FESpaces)
    Random.seed!(1234)
    xh = FEFunction(U,rand(num_free_dofs(U)))
  else
    # Build and solve the problem
    trian_and_meas = []
    for action in actions
      trian = get_domain(model,action)
      m = Measure(trian,2*k)
      push!(trian_and_meas,(trian,m))
    end
    context = nothing
    res(x,dy) =
      c(actions,x,dy,context,trian_and_meas) +
      a(actions,x,dy,context,trian_and_meas) -
      ℓ(actions,dy,context,trian_and_meas)
    jac(x,dx,dy) =
      dc(actions,x,dx,dy,context,trian_and_meas) +
      a(actions,dx,dy,context,trian_and_meas)
    op=FEOperator(res,jac,U,V)
    xh = solve(solver,op)
  end

  out = (solution=xh,info=info)
  out
end

function a(actions::Vector{<:Action},dx,dy,context,trian_and_meas)
  cont = DomainContribution()
  for (i,action) in enumerate(actions)
    dc = a(action,dx,dy,context,trian_and_meas[i])
    if dc !== nothing
      cont = cont + dc
    end
  end
  cont
end

function a(action::Action,dx,dy,context,trian_and_meas)
  nothing
end

function a(action::ConductingFluid,dx,dy,context,trian_and_meas)
  u, p, j, φ = dx
  v_u, v_p, v_j, v_φ = dy
  Ω, dΩ = trian_and_meas
  β = action.β
  ∫(
    β*(∇(u)⊙∇(v_u)) - p*(∇⋅v_u)  +
    (∇⋅u)*v_p +
    j⋅v_j - φ*(∇⋅v_j) +
    (∇⋅j)*v_φ ) * dΩ
end

function a(action::MagneticField,dx,dy,context,trian_and_meas)
  u, p, j, φ = dx
  v_u, v_p, v_j, v_φ = dy
  Ω, dΩ = trian_and_meas
  B = action.B
  γ = action.γ
  ∫( -(γ*(j×B)⋅v_u) - (u×B)⋅v_j )*dΩ
end

function a(action::ConductingThinWall,dx,dy,context,trian_and_meas)
  u, p, j, φ = dx
  v_u, v_p, v_j, v_φ = dy
  Γ, dΓ = trian_and_meas
  n_Γ = get_normal_vector(Γ)
  cw = action.cw
  τ = action.τ
  ∫( τ*((v_j⋅n_Γ)*(j⋅n_Γ) + cw*(v_j⋅n_Γ)*(n_Γ⋅(∇(j)⋅n_Γ))) )*dΓ
end

function ℓ(actions::Vector{<:Action},dy,context,trian_and_meas)
  cont = DomainContribution()
  for (i,action) in enumerate(actions)
    dc = ℓ(action,dy,context,trian_and_meas[i])
    if dc !== nothing
      cont = cont + dc
    end
  end
  cont
end

function ℓ(action::Action,dy,context,trian_and_meas)
  nothing
end

function ℓ(action::ConductingThinWall,dy,context,trian_and_meas)
  u, p, j, φ = dx
  v_u, v_p, v_j, v_φ = dy
  Γ, dΓ = trian_and_meas
  n_Γ = get_normal_vector(Γ)
  jw = action.jw
  τ = action.τ
  ∫( τ*(v_j⋅n_Γ)*jw ) * dΓ
end

function ℓ(action::FluidForce,dy,context,trian_and_meas)
  v_u, v_p, v_j, v_φ = dy
  Ω, dΩ = trian_and_meas
  f = action.f
  ∫( v_u⋅f )*dΩ
end

function ℓ(action::PotentialBc,dy,context,trian_and_meas)
  v_u, v_p, v_j, v_φ = dy
  Γ, dΓ = trian_and_meas
  φ = action.φ
  ∫( -(v_j⋅n_Γ)*φ )*dΓ
end

function c(actions::Vector{<:Action},x,dy,context,trian_and_meas)
  cont = DomainContribution()
  for (i,action) in enumerate(actions)
    d = c(action,x,dy,context,trian_and_meas[i])
    if d !== nothing
      cont = cont + d
    end
  end
  cont
end

function c(action::Action,x,dy,context,trian_and_meas)
  nothing
end

conv(u,∇u) = (∇u')⋅u

function c(action::ConductingFluid,x,dy,context,trian_and_meas)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  Ω, dΩ = trian_and_meas
  α = action.α
  ∫( α*v_u⋅(conv∘(u,∇(u))) ) * dΩ
end

function dc(actions::Vector{<:Action},x,dx,dy,context,trian_and_meas)
  cont = DomainContribution()
  for (i,action) in enumerate(actions)
    d = dc(action,x,dx,dy,context,trian_and_meas[i])
    if d !== nothing
      cont = cont + d
    end
  end
  cont
end

function dc(action::Action,x,dx,dy,context,trian_and_meas)
  nothing
end

function dc(action::ConductingFluid,x,dx,dy,context,trian_and_meas)
  u, p, j, φ = x
  du , dp , dj , dφ  = dx
  v_u, v_p, v_j, v_φ = dy
  Ω, dΩ = trian_and_meas
  α = action.α
  ∫( α*v_u⋅( (conv∘(u,∇(du))) + (conv∘(du,∇(u))) ) ) * dΩ
end

function find_strong_bcs_u(bcs)
  tags = String[]
  vals = []
  for bc in bcs
    if isa(bc,VelocityBc)
      push!(tags,bc.domain)
      push!(vals,bc.u)
    end
  end
  tags, vals
end

function find_strong_bcs_j(bcs)
  tags = String[]
  vals = []
  for bc in bcs
    if isa(bc,CurrentBc)
      push!(tags,bc.domain)
      push!(vals,bc.j)
    end
  end
  tags, vals
end

function get_domain(model,action)
  get_domain(model,action,action.domain)
end

function get_domain(model,action,domain::Triangulation)
  domain
end

function get_domain(model,action::BodyAction,tag::String)
  Triangulation(model,tags=tag)
end

function get_domain(model,action::BoundaryAction,tag::String)
  Boundary(model,tags=tag)
end

