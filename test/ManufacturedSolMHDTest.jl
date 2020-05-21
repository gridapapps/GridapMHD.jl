module ManufacturedSolMHDTest

using Gridap
using Gridap.FESpaces: residual
using Gridap.MultiField: MultiFieldFEFunction
using Gridap: ∇, Δ
using LinearAlgebra
using Test

u(x) = VectorValue(x[1]^2+x[2]^2,-x[1]-x[2]^2,one(x[1])+x[2]^2)
p(x) = x[1]+x[2]+x[3]
j(x) = VectorValue(x[1]^2+x[2],-x[1]-x[2]^2,one(x[1])+x[3]^2)
φ(x) = x[1]+x[2]+x[3]
B = VectorValue(1.0,1.0,1.0)

∇u(x) = ∇(u)(x)
Δu(x) = Δ(u)(x)
∇p(x) = ∇(p)(x)
∇φ(x) = ∇(φ)(x)

@law vprod(a,b) = VectorValue(a[2]b[3]-a[3]b[2], a[1]b[3]-a[3]b[1], a[1]b[2]-a[2]b[1])


f_u(x) = (∇u(x)')*u(x) - Δu(x) + ∇p(x) - vprod(j(x),B)
f_p(x) = (∇*u)(x)
f_j(x) = j(x) + ∇φ(x) - vprod(u(x),B)
f_φ(x) = (∇*j)(x)

g_u(x) = u(x)
g_j(x) = j(x)

n = 4
domain = (-0.5,0.5,-0.5,0.5,-0.5,0.5)
partition = (n,n,n)
model = CartesianDiscreteModel(domain,partition)

dirichlet_u = collect(1:24)
dirichlet_j = collect(1:24)
neumann_u = collect(25:26)
neumann_j = collect(25:26)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet_u",dirichlet_u)
add_tag_from_tags!(labels,"dirichlet_j",dirichlet_j)
add_tag_from_tags!(labels,"neumann_u",neumann_u)
add_tag_from_tags!(labels,"neumann_j",neumann_j)

order = 2

Vu = FESpace(
    reffe=:Lagrangian, order=order, valuetype=VectorValue{3,Float64},
    conformity=:H1, model=model, dirichlet_tags="dirichlet_u")

if length(neumann_u) == 0
  Vp = FESpace(
      reffe=:PLagrangian, order=order-1, valuetype=Float64,
      conformity=:L2, model=model, constraint=:zeromean)
else
  Vp = FESpace(
      reffe=:PLagrangian, order=order-1, valuetype=Float64,
      conformity=:L2, model=model)
end

Vj = FESpace(
    reffe=:RaviartThomas, order=order-1, valuetype=VectorValue{3,Float64},
    conformity=:Hdiv, model=model, dirichlet_tags="dirichlet_j")

if length(neumann_j) == 0
  Vφ = FESpace(
      reffe=:QLagrangian, order=order-1, valuetype=Float64,
      conformity=:L2, model=model, constraint=:zeromean)
else
  Vφ = FESpace(
      reffe=:QLagrangian, order=order-1, valuetype=Float64,
      conformity=:L2, model=model)
end

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

  (∇(u)'*uk)*v_u + inner(∇(u),∇(v_u)) - p*(∇*v_u) - vprod(j,B)*v_u +
  (∇*u)*v_p +
  j*v_j - φ*(∇*v_j) - vprod(u,B)*v_j +
  (∇*j)*v_φ
end

function l(Y)
  v_u, v_p, v_j, v_φ = Y

  v_u*f_u + v_p*f_p + v_j*f_j + v_φ*f_φ
end
t_Ω = AffineFETerm(a,l,trian,quad)

if length(neumann_u) > 0
  btrian_u = BoundaryTriangulation(model,"neumann_u")
  bquad_u = CellQuadrature(btrian_u,degree)
  nb_u = get_normal_vector(btrian_u)

  function l_Γ_u(Y)
    v_u, v_p, v_j, v_φ = Y

    v_u*(nb_u*∇u) - (nb_u*v_u)*p
  end
  t_Γ_u = FESource(l_Γ_u,btrian_u,bquad_u)
end

if length(neumann_j) > 0
  btrian_j = BoundaryTriangulation(model,"neumann_j")
  bquad_j = CellQuadrature(btrian_j,degree)
  nb_j = get_normal_vector(btrian_j)

  function l_Γ_j(Y)
    v_u, v_p, v_j, v_φ = Y

    -(v_j*nb_j)*φ
  end
  t_Γ_j = FESource(l_Γ_j,btrian_j,bquad_j)
end

if (length(neumann_u) > 0) & (length(neumann_j) > 0)
  op  = AffineFEOperator(X,Y,t_Ω, t_Γ_u, t_Γ_j )
elseif length(neumann_u) > 0
  op  = AffineFEOperator(X,Y,t_Ω, t_Γ_u )
elseif length(neumann_j) > 0
  op  = AffineFEOperator(X,Y,t_Ω, t_Γ_j )
else
  op  = AffineFEOperator(X,Y,t_Ω)
end


xh = solve(op)
uh, ph, jh, φh = xh

r = residual(op,xh)

rh = FEFunction(Y,r)
ruh, rph, rjh, rφh = rh


eu = uh - u
ep = ph - p
ej = jh - j
eφ = φh - φ

l2(v) = v*v
h1(v) = v*v + inner(∇(v),∇(v))
hdiv(v) = v*v + inner((∇*v),(∇*v))

eu_l2 = sqrt(sum(integrate(l2(eu),trian,quad)))
eu_h1 = sqrt(sum(integrate(h1(eu),trian,quad)))
ep_l2 = sqrt(sum(integrate(l2(ep),trian,quad)))
ej_l2 = sqrt(sum(integrate(l2(ej),trian,quad)))
ej_hdiv = sqrt(sum(integrate(hdiv(ej),trian,quad)))
eφ_l2 = sqrt(sum(integrate(l2(eφ),trian,quad)))

@test eu_l2 < 1e-12
@test eu_h1 < 1e-12
@test ep_l2 < 1e-12
@test ej_l2 < 1e-12
@test ej_hdiv < 1e-12
@test eφ_l2 < 1e-12

end
