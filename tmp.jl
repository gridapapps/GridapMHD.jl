module TMP

using Gridap
using Gridap.FESpaces: residual
using Gridap.MultiField: MultiFieldFEFunction
using Gridap: ∇, Δ
using LinearAlgebra
using Test

u(x) = VectorValue(x[1]+x[2],-x[1]-x[2],1.0)
p(x) = x[1]+x[2]+x[3]
j(x) = VectorValue(x[1]+x[2],-x[1]-x[2],1.0)
φ(x) = x[1]+x[2]+x[3]

∇u(x) = TensorValue(1.0,1.0,0.0, -1.0,-1.0,0.0, 0.0,0.0,0.0)

B = VectorValue(1.0,1.0,1.0)

# f_u(x) = ((∇(u)(x))')*u(x) - Δ(u)(x) + ∇(p)(x) - cross(j,B)
# f_p(x) = (∇*u)(x)
# f_j(x) = j(x) + ∇(φ)(x) - cross(u,B)
# f_φ(x) = (∇*j)(x)

@law vprod(a,b) = VectorValue(a[2]b[3]-a[3]b[2], a[1]b[3]-a[3]b[1], a[1]b[2]-a[2]b[1])

f_u(x) = 2.0*VectorValue(x[1]+x[2],x[1]+x[2],0.0) + VectorValue(1.0,1.0,1.0) - vprod(j(x),B)
f_p(x) = 0.0
f_j(x) = j(x) + VectorValue(1.0,1.0,1.0) - vprod(u(x),B)
f_φ(x) = 0.0

g_u(x) = u(x)
g_j(x) = j(x)

n = 4
domain = (-0.5,0.5,-0.5,0.5,-0.5,0.5)
partition = (n,n,n)
model = CartesianDiscreteModel(domain,partition)

none = Array{Int,1}(undef,0)
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet_u",collect(1:24))
add_tag_from_tags!(labels,"dirichlet_j",collect(1:24))

order = 2

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

  uk*(∇(u)*v_u) + inner(∇(u),∇(v_u)) - p*(∇*v_u) - vprod(j,B)*v_u +
  (∇*u)*v_p +
  j*v_j - φ*(∇*v_j) - vprod(u,B)*v_j +
  (∇*j)*v_φ
end

function l(Y)
  v_u, v_p, v_j, v_φ = Y

  v_u*f_u + v_p*f_p + v_j*f_j + v_φ*f_φ
end

btrian = BoundaryTriangulation(model,collect(25:26))
bquad = CellQuadrature(btrian,degree)
nb = get_normal_vector(btrian)

function l_Γ(Y)
  v_u, v_p, v_j, v_φ = Y

  -(v_j*nb)*φ +
  v_u*(nb*∇u) - (nb*v_u)*p
end

t_Ω = AffineFETerm(a,l,trian,quad)
t_Γ = FESource(l_Γ,btrian,bquad)
op  = AffineFEOperator(X,Y,t_Ω, t_Γ )

uh = interpolate(U,u)
ph = interpolate(P,p)
jh = interpolate(J,j)
φh = interpolate(Φ,φ)

fv = Vector{Float64}(undef,0)
append!(fv, uh.free_values)
append!(fv, ph.free_values)
append!(fv, jh.free_values)
append!(fv, φh.free_values)
xh = FEFunction(X,fv)
r = residual(op,xh)

rh = FEFunction(Y,r)
ruh, rph, rjh, rφh = rh

uh, ph, jh, φh = solve(op)

eu = uh - u
ep = ph - p
ej = jh - j
eφ = φh - φ

writevtk(trian,"results", nsubcells=10,cellfields=["uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh,
                                     "ruh"=>ruh,"rph"=>rph,"rjh"=>rjh,"rφh"=>rφh,
                                     "eu"=>eu,"ep"=>ep,"ej"=>ej,"eφ"=>eφ])


# writevtk(btrian,"btrian",cellfields=["uh"=>restrict(uh,btrian),"jhn"=>nb*restrict(jh,btrian)])

end
