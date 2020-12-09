module ManufacturedSolMHDTest

using Gridap
using Gridap.FESpaces: residual
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


f_u(x) = (∇u(x)')⋅u(x) - Δu(x) + ∇p(x) - j(x)×B
f_p(x) = (∇⋅u)(x)
f_j(x) = j(x) + ∇φ(x) - u(x)×B
f_φ(x) = (∇⋅j)(x)

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

Vu = FESpace(model, ReferenceFE(:Lagrangian,VectorValue{3,Float64},order);
        conformity=:H1, dirichlet_tags="dirichlet_u")

if length(neumann_u) == 0
  Vp = FESpace(model, ReferenceFE(:Lagrangian,Float64,order-1,space=:P);
          conformity=:L2, constraint=:zeromean)
else
  Vp = FESpace(model, ReferenceFE(:Lagrangian,Float64,order-1,space=:P);
          conformity=:L2)
end

Vj = FESpace(model, ReferenceFE(:RaviartThomas,Float64,order-1);
    conformity=:Hdiv, dirichlet_tags="dirichlet_j")

if length(neumann_j) == 0
  Vφ = FESpace(model, ReferenceFE(:Lagrangian,Float64,order-1,space=:Q);
      conformity=:L2, constraint=:zeromean)
else
  Vφ = FESpace(model, ReferenceFE(:Lagrangian,Float64,order-1,space=:Q);
      conformity=:L2)
end

U = TrialFESpace(Vu,g_u)
P = TrialFESpace(Vp)
J = TrialFESpace(Vj,g_j)
Φ = TrialFESpace(Vφ)

Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
X = MultiFieldFESpace([U, P, J, Φ])

trian = Triangulation(model)
degree = 2*(order)
dΩ = LebesgueMeasure(trian,degree)

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

  ∫(v_u⋅f_u + v_p*f_p + v_j⋅f_j + v_φ*f_φ)*dΩ
end

if length(neumann_u) > 0
  btrian_u = BoundaryTriangulation(model,tags=["neumann_u"])
  dΓ_u = LebesgueMeasure(btrian_u,degree)
  nb_u = get_normal_vector(btrian_u)

  function l_Γ_u(Y)
    v_u, v_p, v_j, v_φ = Y

    ∫(v_u⋅(nb_u⋅∇u) - (nb_u⋅v_u)*p)*dΓ_u
  end
end

if length(neumann_j) > 0
  btrian_j = BoundaryTriangulation(model,tags=["neumann_j"])
  dΓ_j = LebesgueMeasure(btrian_j,degree)
  nb_j = get_normal_vector(btrian_j)

  function l_Γ_j(Y)
    v_u, v_p, v_j, v_φ = Y

    ∫(-(v_j⋅nb_j)*φ)*dΓ_j
  end
end

if (length(neumann_u) > 0) & (length(neumann_j) > 0)
  b(Y) = l(Y) + l_Γ_u(Y) + l_Γ_j(Y)
elseif length(neumann_u) > 0
  b(Y) = l(Y) + l_Γ_u(Y)
elseif length(neumann_j) > 0
  b(Y) = l(Y) + l_Γ_j(Y)
else
  b(Y) = l(Y)
end

op  = AffineFEOperator(a,b,X,Y)

xh = solve(op)
uh, ph, jh, φh = xh

r = residual(op,xh)

rh = FEFunction(Y,r)
ruh, rph, rjh, rφh = rh


eu = uh - u
ep = ph - p
ej = jh - j
eφ = φh - φ

l2(v) = sqrt(sum(∫(v⋅v)*dΩ ))
h1(v) = sqrt(sum(∫(v⋅v + inner(∇(v),∇(v)) )*dΩ ))
hdiv(v) = sqrt(sum(∫(v⋅v + inner((∇⋅v),(∇⋅v)) )*dΩ ))

eu_l2 = l2(eu)
eu_h1 = h1(eu)
ep_l2 = l2(ep)
ej_l2 = l2(ej)
ej_hdiv = hdiv(ej)
eφ_l2 = l2(eφ)

@test eu_l2 < 1e-12
@test eu_h1 < 1e-12
@test ep_l2 < 1e-12
@test ej_l2 < 1e-12
@test ej_hdiv < 1e-12
@test eφ_l2 < 1e-12

end #module
