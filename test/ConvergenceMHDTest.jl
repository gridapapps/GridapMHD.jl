module ConvergenceMHDTest

using Gridap
using Test
using Polynomials: fit

u(x) = VectorValue(x[1]^4+x[2]^4,-x[1]^2-x[2]^3,one(x[1])+sin(x[2]^4))
p(x) = sin(x[1]+x[2]+x[3])
j(x) = VectorValue(x[1]^4+x[2],-x[1]-cos(x[2]),one(x[1])+x[3]^4)
φ(x) = sin(x[1]+x[2]+x[3])
B = VectorValue(1.0,1.0,1.0)

∇u(x) = ∇(u)(x)
Δu(x) = Δ(u)(x)
∇p(x) = ∇(p)(x)
∇φ(x) = ∇(φ)(x)


f_u(x) = (∇u(x)')⋅u(x) - Δu(x) + ∇p(x) - j(x)×B
f_p(x) = (∇⋅u)(x)
f_j(x) = j(x) + ∇φ(x) - u(x)×B
f_φ(x) = (∇⋅j)(x)

g_u(x) = u(x)
g_j(x) = j(x)

eu_l2 = Vector{Float64}()
eu_h1 = Vector{Float64}()
ep_l2 = Vector{Float64}()
ej_l2 = Vector{Float64}()
ej_hdiv = Vector{Float64}()
eφ_l2 = Vector{Float64}()

order = 2
nxs = 2:5
domain = (-0.5,0.5,-0.5,0.5,-0.5,0.5)

for n=nxs
  partition = (n,n,n)
  model = CartesianDiscreteModel(domain,partition)

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels,"dirichlet_u",collect(1:24))
  add_tag_from_tags!(labels,"dirichlet_j",collect(1:24))
  add_tag_from_tags!(labels,"neumann_u",collect(25:26))
  add_tag_from_tags!(labels,"neumann_j",collect(25:26))


  Vu = FESpace(model, ReferenceFE(:Lagrangian,VectorValue{3,Float64},order);
      conformity=:H1, dirichlet_tags="dirichlet_u")

  Vp = FESpace(model, ReferenceFE(:Lagrangian,Float64,order-1,space=:P);
      conformity=:L2)

  Vj = FESpace(model, ReferenceFE(:RaviartThomas,Float64,order-1);
      conformity=:Hdiv, dirichlet_tags="dirichlet_j")

  Vφ = FESpace(model, ReferenceFE(:Lagrangian,Float64,order-1,space=:Q);
      conformity=:L2)

  U = TrialFESpace(Vu,g_u)
  P = TrialFESpace(Vp)
  J = TrialFESpace(Vj,g_j)
  Φ = TrialFESpace(Vφ)

  Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
  X = MultiFieldFESpace([U, P, J, Φ])

  trian = Triangulation(model)
  degree = 2*(order)
  dΩ = Measure(trian,degree)

  uk = interpolate(u,U)
  function a(X,Y)
    u  , p  , j  , φ   = X
    v_u, v_p, v_j, v_φ = Y

    ∫((∇(u)'⋅uk)⋅v_u + inner(∇(u),∇(v_u)) - p*(∇⋅v_u) - (j×B)⋅v_u +
      (∇⋅u)*v_p +
      j⋅v_j - φ*(∇⋅v_j) - (u×B)⋅v_j +
      (∇⋅j)*v_φ )*dΩ
  end

  function l(Y)
    v_u, v_p, v_j, v_φ = Y

    ∫( v_u⋅f_u + v_p*f_p + v_j⋅f_j + v_φ*f_φ )*dΩ
  end

  btrian_u = BoundaryTriangulation(model,tags=["neumann_u"])
  dΓ_u = Measure(btrian_u,degree)
  nb_u = get_normal_vector(btrian_u)

  function l_Γ_u(Y)
    v_u, v_p, v_j, v_φ = Y

    ∫( v_u⋅(nb_u⋅∇u) - (nb_u⋅v_u)*p)*dΓ_u
  end

  btrian_j = BoundaryTriangulation(model,tags=["neumann_j"])
  dΓ_j = Measure(btrian_j,degree)
  nb_j = get_normal_vector(btrian_j)

  function l_Γ_j(Y)
    v_u, v_p, v_j, v_φ = Y

    ∫( -(v_j⋅nb_j)*φ)*dΓ_j
  end

  b(Y) = l(Y) + l_Γ_u(Y) + l_Γ_j(Y)

  op  = AffineFEOperator(a,b,X,Y)

  xh = solve(op)
  uh, ph, jh, φh = xh

  eu = uh - u
  ep = ph - p
  ej = jh - j
  eφ = φh - φ

  l2(v) = sqrt(sum(∫(v⋅v)*dΩ ))
  h1(v) = sqrt(sum(∫(v⋅v + inner(∇(v),∇(v)) )*dΩ ))
  hdiv(v) = sqrt(sum(∫(v⋅v + inner((∇⋅v),(∇⋅v)) )*dΩ ))

  append!(eu_l2, l2(eu) )
  append!(eu_h1, h1(eu) )
  append!(ep_l2, l2(ep) )
  append!(ej_l2, l2(ej) )
  append!(ej_hdiv, hdiv(ej) )
  append!(eφ_l2, l2(eφ) )

end


slope_eu_l2 = fit([log(x) for x in nxs], [log(y) for y in eu_l2], 1)[end]
slope_eu_h1 = fit([log(x) for x in nxs], [log(y) for y in eu_h1], 1)[end]
slope_ep_l2 = fit([log(x) for x in nxs], [log(y) for y in ep_l2], 1)[end]
slope_ej_l2 = fit([log(x) for x in nxs], [log(y) for y in ej_l2], 1)[end]
slope_ej_hdiv = fit([log(x) for x in nxs], [log(y) for y in ej_hdiv], 1)[end]
slope_eφ_l2 = fit([log(x) for x in nxs], [log(y) for y in eφ_l2], 1)[end]

@test abs(slope_eu_l2 + order + 1.0) < 2e-1
@test abs(slope_eu_h1 + order) < 2e-1
@test abs(slope_ep_l2 + order) < 2e-1
@test abs(slope_ej_l2 + order + 1.0) < 2e-1
@test abs(slope_ej_hdiv + order) < 2e-1
@test abs(slope_eφ_l2 + order) < 2e-1

end #module
