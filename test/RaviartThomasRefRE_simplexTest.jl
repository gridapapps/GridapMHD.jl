module RaviartThomasRefRE_simplexTest

using Gridap
using GridapMHD
using Test

using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials

xi = Point(2,3)
np = 5
x = fill(xi,np)

order = 0
D = 2
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = PCurlGradMonomialBasis{D}(T,order)

@test num_terms(b) == 3
@test get_order(b) == 0

xi = Point(2,3,5)
np = 5
x = fill(xi,np)

order = 1
D = 3
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = PCurlGradMonomialBasis{D}(T,order)

@test num_terms(b) == 15
@test get_order(b) == 1

end
