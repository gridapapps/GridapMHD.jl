
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

u_physical(ū,u0) = u0*ū
p_physical_cfd(p̄,ρ,u0) = ρ*u0^2*p̄
p_physical_mhd(p̄,σ,u0,B0,L) = σ*u0*B0*L*p̄
p_physical(p̄,σ,u0,B0,L) = p_physical_mhd(p̄,σ,u0,B0,L)
j_physical(j̄,σ,u0,B0) = σ*u0*B0*j̄
φ_physical(φ̄,u0,B0,L) = u0*B0*L*φ̄
f_solver_cfd(f,ρ,u0,L) = (L/ρ*u0^2)*f
f_solver_mhd(f,σ,u0,B0) = f/(σ*u0*B0^2)
f_solver(f,σ,u0,B0) = f_solver_mhd(f,σ,u0,B0)

re_number(ν,u0,L) = (u0*L)/ν
n_number(σ,ρ,u0,B0,L) = (σ*L*B0^2)/(ρ*u0)
ha2_number(σ,ρ,ν,B0,L) = ((L*B0)^2*σ)/(ρ*ν)
ha_number(σ,ρ,ν,B0,L) = (L*B0)*sqrt(σ/(ρ*ν))

@with_kw struct PhysicalParams{A,B,C}
  ρ::A
  ν::B
  σ::C
end

@with_kw struct ScalingFactors
  u0::Float64
  B0::Float64
  L::Float64
end

@with_kw struct ReducedParams{A,B,C}
  Re::A
  Ha::B
  N::C
end

@with_kw struct SolverParams{A,B,C}
  α::A
  β::B
  γ::C
end

function ReducedParams(p::PhysicalParams,f::ScalingFactors)
  Re = re_number(p.ν,f.u0,f.L)
  Ha = ha_number(p.σ,p.ρ,p.ν,f.B0,f.L)
  N = n_number(p.σ,p.ρ,f.u0,f.B0,f.L)
  ReducedParams(Re=Re,Ha=Ha,N=N)
end

function EquationParams(p::PhysicalParams,f::ScalingFactors)
  Re = re_number(p.ν,f.u0,f.L)
  Ha2 = ha2_number(p.σ,p.ρ,p.ν,f.B0,f.L)
  N = n_number(p.σ,p.ρ,f.u0,f.B0,f.L)
  α = (1/N)
  β = (1/Ha2)
  γ = 1
  SolverParams(α=α,β=β,γ=γ)
end

function u_physical(ū,p::PhysicalParams,f::ScalingFactors)
  u_physical(ū,f.u0)
end

function p_physical(p̄,a::PhysicalParams,f::ScalingFactors)
  p_physical(p̄,a.σ,f.u0,f.B0,f.L)
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
@with_kw struct WallLaw{A,B,C}
  domain::A
  cw::B
  value::C = 0.0
end

# source term for the momentum equation
@with_kw struct FluidForce{A,B}
  domain::A
  value::B
end

@with_kw struct Params{A,B,C,D}
  domain::A
  material::B
  actions::C = nothing[]
  B0::D = VectorValue(0.,0.,0.)
end

@with_kw struct Output{A}
  solution::A
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

function main(params::Params)

  Ω = params.domain
  actions = params.actions
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

  out = Output(solution=xh)

  out
end

