module ConvergenceWithPeriodicBCMHDTest

using GridapMHD

using Gridap
using Test
using Polynomials: fit


u(x) = VectorValue(one(x[1]),x[3]^4+x[2]^4,-x[3]^4-x[2]^4)
p(x) = sin(x[3]+x[2])
j(x) = VectorValue(one(x[1]),(x[3])+x[2]^4,-x[3]^4-x[2])
φ(x) = sin(x[3]+x[2])
B = VectorValue(1.0,1.0,1.0)

∇u(x) = ∇(u)(x)
Δu(x) = Δ(u)(x)
∇p(x) = ∇(p)(x)
∇φ(x) = ∇(φ)(x)

function ∇u_n(x)
  if (x[2] < -0.49999)
    n = VectorValue(0.0,-1.0,0.0)
  elseif (x[2] > 0.49999)
    n = VectorValue(0.0,1.0,0.0)
  elseif (x[3] < -0.49999)
    n = VectorValue(0.0,0.0,-1.0)
  elseif (x[3] > 0.49999)
    n = VectorValue(0.0,0.0,1.0)
  elseif (x[1] < -0.49999)
    n = VectorValue(-1.0,0.0,0.0)
  else
    n = VectorValue(1.0,0.0,0.0)
  end
  n*∇u(x)
end

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
  partition = (3,n,n)
  model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false,false))

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels,"dirichlet_u",collect(1:22))
  add_tag_from_tags!(labels,"dirichlet_j",collect(1:22))
  add_tag_from_tags!(labels,"neumann_u",collect(23:24))
  add_tag_from_tags!(labels,"neumann_j",collect(23:24))


  Vu = FESpace(
      reffe=:Lagrangian, order=order, valuetype=VectorValue{3,Float64},
      conformity=:H1, model=model, dirichlet_tags="dirichlet_u")

  Vp = FESpace(
      reffe=:PLagrangian, order=order-1, valuetype=Float64,
      conformity=:L2, model=model)

  Vj = FESpace(
      reffe=:RaviartThomas, order=order-1, valuetype=VectorValue{3,Float64},
      conformity=:Hdiv, model=model, dirichlet_tags="dirichlet_j")

  Vφ = FESpace(
      reffe=:QLagrangian, order=order-1, valuetype=Float64,
      conformity=:L2, model=model)

  U = TrialFESpace(Vu,g_u)
  P = TrialFESpace(Vp)
  J = TrialFESpace(Vj,g_j)
  Φ = TrialFESpace(Vφ)

  Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
  X = MultiFieldFESpace([U, P, J, Φ])

  trian = Triangulation(model)
  degree = 2*(order)
  quad = CellQuadrature(trian,degree)

  uk = interpolate(U,u)
  function a(X,Y)
    u  , p  , j  , φ   = X
    v_u, v_p, v_j, v_φ = Y

    (∇(u)'⋅uk)⋅v_u + inner(∇(u),∇(v_u)) - p*(∇⋅v_u) - (j×B)⋅v_u +
    (∇⋅u)*v_p +
    j⋅v_j - φ*(∇⋅v_j) - (u×B)⋅v_j +
    (∇⋅j)*v_φ
  end

  function l(Y)
    v_u, v_p, v_j, v_φ = Y

    v_u⋅f_u + v_p*f_p + v_j⋅f_j + v_φ*f_φ
  end

  btrian_u = BoundaryTriangulation(model,"neumann_u")
  bquad_u = CellQuadrature(btrian_u,degree)
  nb_u = get_normal_vector(btrian_u)

  function l_Γ_u(Y)
    v_u, v_p, v_j, v_φ = Y

    v_u⋅(nb_u⋅∇u) - (nb_u⋅v_u)*p
  end

  btrian_j = BoundaryTriangulation(model,"neumann_j")
  bquad_j = CellQuadrature(btrian_j,degree)
  nb_j = get_normal_vector(btrian_j)

  function l_Γ_j(Y)
    v_u, v_p, v_j, v_φ = Y

    -(v_j⋅nb_j)*φ
  end

  t_Ω = AffineFETerm(a,l,trian,quad)
  t_Γ_u = FESource(l_Γ_u,btrian_u,bquad_u)
  t_Γ_j = FESource(l_Γ_j,btrian_j,bquad_j)
  op  = AffineFEOperator(X,Y,t_Ω, t_Γ_u, t_Γ_j )

  xh = solve(op)
  uh, ph, jh, φh = xh

  eu = uh - u
  ep = ph - p
  ej = jh - j
  eφ = φh - φ

  l2(v) = v⋅v
  h1(v) = v⋅v + inner(∇(v),∇(v))
  hdiv(v) = v⋅v + inner((∇⋅v),(∇⋅v))

  append!(eu_l2, sqrt(sum(integrate(l2(eu),trian,quad))))
  append!(eu_h1, sqrt(sum(integrate(h1(eu),trian,quad))))
  append!(ep_l2, sqrt(sum(integrate(l2(ep),trian,quad))))
  append!(ej_l2, sqrt(sum(integrate(l2(ej),trian,quad))))
  append!(ej_hdiv, sqrt(sum(integrate(hdiv(ej),trian,quad))))
  append!(eφ_l2, sqrt(sum(integrate(l2(eφ),trian,quad))))

  # writevtk(trian,"results", nsubcells=10,cellfields=["uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh,
  #                                                    "eu"=>eu,"ep"=>ep,"ej"=>ej,"eφ"=>eφ])

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
