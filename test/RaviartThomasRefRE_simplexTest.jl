#module RaviartThomasRefRE_simplexTest

using Gridap
using GridapMHD
using GridapGmsh
using GridapMHD.RaviartThomasRefFE_simplex: RaviartThomasRefFE_t
using Test

# using Gridap.Integrationi
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.FESpaces

xi = Point(4,2)
np = 1
x = fill(xi,np)

order = 2
D = 2
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = PCurlGradMonomialBasis{D}(T,order)

v = V[
  (1.0, 0.0), (0.0, 1.0), (4.0, 0.0), (0.0, 2.0), (16.0, 0.0), (0.0, 4.0),
  (2.0, 0.0), (0.0, 4.0), (8.0, 0.0), (0.0, 8.0), (4.0, 0.0), (0.0, 16.0),
  (64.0, 32.0), (32.0, 16.0), (16.0, 8.0)]

g = G[
  (0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0),
  (0.0, 0.0, 0.0, 1.0), (8.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 4.0),
  (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0), (2.0, 4.0, 0.0, 0.0),
  (0.0, 0.0, 2.0, 4.0), (0.0, 4.0, 0.0, 0.0), (0.0, 0.0, 8.0, 0.0),
  (48.0, 0.0, 16.0, 16.0), (16.0, 16.0, 4.0, 16.0), (4.0, 16.0, 0.0, 12.0)]

vb = evaluate_field(b,x)

for (vi,vbi) in zip(v,vb)
  @test vi == vbi
end

∇b = field_gradient(b)
gvb = evaluate_field(∇b,x)
for (vi,vbi) in zip(g,gvb)
  @test vi == vbi
end

@test num_terms(b) == 15
@test get_order(b) == 2

xi = Point(2,3,5)
np = 3
x = fill(xi,np)

order = 1
D = 3
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = PCurlGradMonomialBasis{D}(T,order)

@test num_terms(b) == 15
@test get_order(b) == 1


p = TRI
D = num_dims(TRI)
et = Float64
order = 0

reffe = RaviartThomasRefFE_t(et,p,order)
@show test_reference_fe(reffe)
@show num_terms(get_prebasis(reffe))
@show num_dofs(reffe)
@show get_order(get_prebasis(reffe))

model = GmshDiscreteModel("./test_2d.msh")
labels = get_face_labeling(model)
dir_tags = Array{Integer}(undef,0)
trian = Triangulation(model)
strian = SkeletonTriangulation(model)
btrian = BoundaryTriangulation(model)


quad = CellQuadrature(trian,2*order+1)
V = ConformingFESpace([reffe],DivConformity(),model,labels,dir_tags)
free_values = ones(num_free_dofs(V))
uh = FEFunction(V,free_values)

#writevtk(strian,"test",cellfields=["nv"=>nv])

I = integrate(uh,trian,quad)
@show I
@show sum(I)
#end #module
