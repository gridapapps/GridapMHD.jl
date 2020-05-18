module DBG

using Gridap
using Gridap: ∇, Δ
using LinearAlgebra
using Test

j(x) = VectorValue(x[1]+x[2],-x[1]-x[2],1.0)
f_φ(x) = 0.0
g_j(x) = j(x)

n = 4
domain = (-0.5,0.5,-0.5,0.5,-0.5,0.5)
partition = (n,n,n)
model = CartesianDiscreteModel(domain,partition)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet_j",collect(1:24))
add_tag_from_tags!(labels,"neumann_j",collect(25:26))

order = 2

Vj = FESpace(
    reffe=:RaviartThomas, order=order-1, valuetype=VectorValue{3,Float64},
    conformity=:Hdiv, model=model, dirichlet_tags="dirichlet_j")

Vφ = FESpace(
    reffe=:QLagrangian, order=order-1, valuetype=Float64,
    conformity=:H1, model=model)

Uj = TrialFESpace(Vj,g_j)
Uϕ = TrialFESpace(Vφ)

jh = interpolate(Uj,j)

v_φ = FEFunction(Vφ,rand(num_free_dofs(Vφ)))


trian = Triangulation(model)

trian_Γn = BoundaryTriangulation(model,"neumann_j")
n_Γn = get_normal_vector(trian_Γn)

trian_Γd = BoundaryTriangulation(model,"dirichlet_j")
n_Γd = get_normal_vector(trian_Γd)

jh_Γn = restrict(jh,trian_Γn)
jh_Γd = restrict(jh,trian_Γd)

v_φ_Γn = restrict(v_φ,trian_Γn)
v_φ_Γd = restrict(v_φ,trian_Γd)

quad = CellQuadrature(trian,2*order)
quad_Γn = CellQuadrature(trian_Γn,2*order)
quad_Γd = CellQuadrature(trian_Γd,2*order)

a_Ω(j,v_φ) = v_φ*(∇*j)
b_Ω(j,v_φ) = - ∇(v_φ)*j
c_Γn(j,v_φ) = n_Γn*(v_φ*j)
d_Γd(v_φ) = - n_Γd*(v_φ*g_j)

ia = sum( integrate( a_Ω(jh,v_φ), trian, quad) )
ib = sum( integrate( b_Ω(jh,v_φ), trian, quad) )
ic = sum( integrate( c_Γn(jh_Γn,v_φ_Γn), trian_Γn, quad_Γn) )
id = sum( integrate( d_Γd(v_φ_Γd), trian_Γd, quad_Γd) )

@show ia
@show ib + ic - id



end # mdoule
