
# MHD equations
# ∇⋅u = 0
# u⋅∇(u) -ν*Δ(u) + (1/ρ)*∇(p) - (1/ρ)*(j×B) = (1/ρ)*f
# ∇⋅j = 0
# j + σ*∇(φ) - σ*(u×B) = 0
#
# One can provide scaling factors
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

@with_kw struct ConductingFluid{A,B,C,D}
  domain::A
  α::B
  β::C
  γ::D
end

@with_kw struct MagneticField{A,B}
  domain::A
  value::B
end

@with_kw struct FluidForce{A,B}
  domain::A
  value::B
end

# u = value
@with_kw struct VelocityBc{A,B}
  domain::A
  value::B = VectorValue(0.,0.,0.)
end

# n⋅∇(u) - p*n = value
@with_kw struct TractionBc{A,B}
  domain::A
  value::B
end

# j = value (imposed in normal direction only)
@with_kw struct InsulatingBc{A,B}
  domain::A
  value::B = VectorValue(0.,0.,0.)
end

# φ = value
@with_kw struct ConductingBc{A,B}
  domain::A
  value::B = 0.0
end

# j⋅n + cw*n⋅∇(j)⋅n = value
@with_kw struct SemiConductingBc{A,B,C}
  domain::A
  cw::B
  value::C = 0.0
end

function find_strong_bcs_u(bcs)
  tags = String[]
  vals = []
  for bc in bcs
    if isa(bc,VelocityBc)
      push!(tags,bc.domain)
      push!(vals,bc.value)
    end
  end
  tags, vals
end

function find_strong_bcs_j(bcs)
  tags = String[]
  vals = []
  for bc in bcs
    if isa(bc,InsulatingBc)
      push!(tags,bc.domain)
      push!(vals,bc.value)
    end
  end
  tags, vals
end

function main(
  fluid::ConductingFluid,
  actions::Vector;
  debug_setup=false,
  title="test")

  Ω = fluid.domain
  model = get_background_model(Ω)
  u_tags, u_vals = find_strong_bcs_u(actions)
  j_tags, j_vals = find_strong_bcs_j(actions)

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

  # Trial Spaces TODO improve for parallel computations
  U_u = TrialFESpace(V_u,u_vals)
  U_p = TrialFESpace(V_p)
  U_j = TrialFESpace(V_j,j_vals)
  U_φ = TrialFESpace(V_φ)
  U = MultiFieldFESpace([U_u,U_p,U_j,U_φ])
  xh = zero(U)

  if debug_setup
    rh = FEFunction(U,rand(num_free_dofs(U)))
    uh,ph,jh,φh = rh
    writevtk(Ω,"$(title)_Ω_fluid",
      order=2,
      cellfields=["uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh])
    for (i,action) in enumerate(actions)
      debug_setup_action(title,model,i,action)
    end
  end

  out = (solution=xh,)

  out
end

function debug_setup_action(title,model,i,action)
end

function debug_setup_action(title,model,i,action::MagneticField)
  Ω = get_domain(model,action)
  writevtk(Ω,order=2,"$(title)_action_$(i)",cellfields=["B"=>action.value])
end

function debug_setup_action(title,model,i,action::ConductingBc)
  Ω = get_domain(model,action)
  writevtk(Ω,order=2,"$(title)_action_$(i)",cellfields=["φ_bc"=>action.value])
end

function debug_setup_action(title,model,i,action::VelocityBc)
  Ω = get_domain(model,action)
  writevtk(Ω,order=2,"$(title)_action_$(i)",cellfields=["u_bc"=>action.value])
end

function debug_setup_action(title,model,i,action::InsulatingBc)
  Ω = get_domain(model,action)
  writevtk(Ω,order=2,"$(title)_action_$(i)",cellfields=["j_bc"=>action.value])
end

function get_domain(model,action)
  get_domain(model,action,action.domain)
end

function get_domain(model,action,domain::Triangulation)
  domain
end

function get_domain(model,action::MagneticField,tag::String)
  Triangulation(model,tags=tag)
end

function get_domain(model,action::ConductingBc,tag::String)
  Boundary(model,tags=tag)
end

function get_domain(model,action::InsulatingBc,tag::String)
  Boundary(model,tags=tag)
end

function get_domain(model,action::VelocityBc,tag::String)
  Boundary(model,tags=tag)
end


