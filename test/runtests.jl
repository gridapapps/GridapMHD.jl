# using GridapMHD
include("../src/GridapMHD.jl")
using .GridapMHD

# include("../../../.julia/dev/Gridap/src/FESpaces/FESpaces.jl")
# using .FESpaces
using Gridap
import .Gridap: ∇

using Test
using Gridap.Geometry: CartesianDiscreteModel
# using Gridap.Geometry: get_cell_vertices, get_vertex_coordinates, get_cell_nodes
# using Gridap.Geometry: get_node_coordinates


model = CartesianDiscreteModel((0,1,0,1),(32,32))

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",[1,2,3,4,5,6])

trian = get_triangulation(model)
degree = 1
quad = CellQuadrature(trian,degree)

Vu = FESpace(
      reffe=:Lagrangian, order=1, valuetype=Float64,
      conformity=:H1, model=model, dirichlet_tags=[1,2,3,4,5,6])

g(x) = 0
f(x) = x[1]*sin(5.0π*x[2])+exp(-((x[1]-0.5)^2+(x[2]-0.5)^2)/0.02)
U = TrialFESpace(Vu,g)

aa(u,v) = inner(∇(v),∇(u))
l(v) = v*f
t_Ω = AffineFETerm(aa,l,trian,quad)

op = AffineFEOperator(U,Vu,t_Ω)

uh = solve(op)

writevtk(trian,"results",cellfields=["u"=>uh])
