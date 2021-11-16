module GridapMHD

using Parameters
using Gridap

export Material
export VelocityBc
export TractionBc
export InsulatingBc
export ConductingBc
export WallLaw
export Params

@with_kw struct Material{A,B}
  N::A
  Re::B
end

# u = value
@with_kw struct VelocityBc{A,B}
  domain::A
  value::B
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

@with_kw struct Params{A,B,C}
  domain::A
  material::B
  bcs::C = nothing[]
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
  bcs = params.bcs
  u_tags, u_vals = find_strong_bcs_u(bcs)
  j_tags, j_vals = find_strong_bcs_j(bcs)

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


end # module
