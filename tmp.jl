module TMP

using Gridap
using LinearAlgebra
using Test

u(x) = VectorValue(x[1]+x[2],-x[1]-x[2],1.0)
p(x) = x[1]+x[2]+x[3]
j(x) = VectorValue(x[1]+x[2]+x[3],x[1]+x[2]+x[3],x[1]+x[2]+x[3])
φ(x) = x[1]+x[2]+x[3]

f_u(x) = - Δ(u)(x) + ∇(p)(x)
f_p(x) = (∇*u)(x)
f_j(x) = j(x) + ∇(φ)(x)
f_φ(x) = (∇*j)(x)

g_u(x) = u(x)
g_j(x) = j(x)

n = 3
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

uh, ph, jh, φh = FEFunction(X,rand(num_free_dofs(X)))

trian = Triangulation(model)

writevtk(trian,"trian",cellfields=["uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh])

btrian = BoundaryTriangulation(model,"boundary")
nb = get_normal_vector(btrian)

writevtk(btrian,"btrian",cellfields=["uh"=>restrict(uh,btrian),"jhn"=>nb*restrict(jh,btrian)])

end
