module TMP

using Gridap
using Gridap: ∇, Δ
using LinearAlgebra
using Test

u(x) = VectorValue(x[1]+x[2],-x[1]-x[2],1.0)
p(x) = x[1]+x[2]+x[3]-1.5
j(x) = VectorValue(x[2],-x[1],0.0)
φ(x) = x[1]+x[2]+x[3]-1.5

B = VectorValue(0.0,1.0,0.0)

# f_u(x) = ((∇(u)(x))')*u(x) - Δ(u)(x) + ∇(p)(x)
# f_p(x) = (∇*u)(x)
# f_j(x) = j(x) + ∇(φ)(x)
# f_φ(x) = (∇*j)(x)

@law vprod(a,b) = VectorValue(a[2]b[3]-a[3]b[2], a[1]b[3]-a[3]b[1], a[1]b[2]-a[2]b[1])

f_u(x) = 2.0*VectorValue(x[1]+x[2],x[1]+x[2],0.0) + VectorValue(1.0,1.0,1.0)
f_p(x) = 0.0
f_j(x) = j(x) + VectorValue(1.0,1.0,1.0) #- vprod(u(x),B)
f_φ(x) = 0.0

g_u(x) = u(x)
g_j(x) = j(x)

n = 4
domain = (0,1,0,1,0,1)
partition = (n,n,n)
model = CartesianDiscreteModel(domain,partition)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet_u","boundary")
add_tag_from_tags!(labels,"dirichlet_j","boundary")

order = 2

Vu = FESpace(
    reffe=:Lagrangian, order=order, valuetype=VectorValue{3,Float64},
    conformity=:H1, model=model, dirichlet_tags="dirichlet_u")

Vp = FESpace(
    reffe=:PLagrangian, order=order-1, valuetype=Float64,
    conformity=:L2, model=model, constraint=:zeromean)

Vj = FESpace(
    reffe=:RaviartThomas, order=order-1, valuetype=VectorValue{3,Float64},
    conformity=:Hdiv, model=model, dirichlet_tags="dirichlet_j")

Vφ = FESpace(
    reffe=:QLagrangian, order=order-1, valuetype=Float64,
    conformity=:L2, model=model, constraint=:zeromean)

U = TrialFESpace(Vu,g_u)
P = TrialFESpace(Vp)
J = TrialFESpace(Vj,g_j)
Φ = TrialFESpace(Vφ)

Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
X = MultiFieldFESpace([U, P, J, Φ])

trian = Triangulation(model)
degree = 2*(order)
quad = CellQuadrature(trian,degree)

uk = interpolate_everywhere(Vu,u)
function a(X,Y)
  u  , p  , j  , φ   = X
  v_u, v_p, v_j, v_φ = Y

  uk*(∇(u)*v_u) + inner(∇(u),∇(v_u)) - p*(∇*v_u) +
  (∇*u)*v_p +
  j*v_j - φ*(∇*v_j) + # vprod(u,B)*v_j
  (∇*j)*(v_φ)
end

function l(Y)
  v_u, v_p, v_j, v_φ = Y

  v_u*f_u + v_p*f_p + v_j*f_j + v_φ*f_φ
end

btrian = BoundaryTriangulation(model,"boundary")
bquad = CellQuadrature(btrian,degree)

nb = get_normal_vector(btrian)
# function a_Γ(X,Y)
#   u  , p  , j  , φ   = X
#   v_u, v_p, v_j, v_φ = Y
#
#   v_φ*(j*nb)
# end

function l_Γ(Y)
  v_u, v_p, v_j, v_φ = Y

  -(v_j*nb)*φ
end

t_Ω = AffineFETerm(a,l,trian,quad)
t_Γ = FESource(l_Γ,btrian,bquad)
op  = AffineFEOperator(X,Y,t_Ω)

uh, ph, jh, φh = solve(op)

eu = uh - u
ep = ph - p
ej = jh - j
eφ = φh - φ

writevtk(trian,"results",cellfields=["uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh,
                                     "eu"=>eu,"ep"=>ep,"ej"=>ej,"eφ"=>eφ])


writevtk(btrian,"btrian",cellfields=["uh"=>restrict(uh,btrian),"jhn"=>nb*restrict(jh,btrian)])

end
