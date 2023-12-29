using PartitionedArrays
using Gridap, GridapSolvers, GridapDistributed

using GridapSolvers.LinearSolvers
using GridapSolvers.MultilevelTools
using GridapSolvers.BlockSolvers 
using GridapSolvers.NonlinearSolvers

using Gridap.MultiField
using Gridap.Arrays, Gridap.Fields

np = (2,2,1)
ranks = with_debug() do distribute
  distribute(LinearIndices((prod(np),)))
end

n = 16
model = CartesianDiscreteModel(ranks,np,(0,1,0,1,0,1),(n,n,n))

k = 1
T = Float64
D = num_cell_dims(model)
reffe_u = ReferenceFE(lagrangian,VectorValue{D,T},k)
reffe_p = ReferenceFE(lagrangian,T,k-1;space=:P)
reffe_j = ReferenceFE(raviart_thomas,T,k-1)
reffe_φ = ReferenceFE(lagrangian,T,k-1)

# Test spaces
mfs = BlockMultiFieldStyle(4,(1,1,1,1),(3,1,2,4))
Ωf  = Triangulation(model)
V_u = TestFESpace(model,reffe_u;dirichlet_tags="boundary")
V_p = TestFESpace(model,reffe_p)
V_j = TestFESpace(model,reffe_j;dirichlet_tags="boundary")
V_φ = TestFESpace(model,reffe_φ;conformity=:L2)
V = MultiFieldFESpace([V_u,V_p,V_j,V_φ];style=mfs)

dΩ = Measure(Ωf,3*k)
f = VectorValue(1.0,1.0,1.0)
a((u,p,j,φ),(v,q,ψ,χ)) = ∫(u⋅v + p*q + j⋅ψ + φ*χ)dΩ
l((v,q,ψ,χ)) = ∫(v⋅f + ψ⋅f)dΩ

assem = SparseMatrixAssembler(V,V)
op = AffineFEOperator(a,l,V,V)
A, b = get_matrix(op), get_vector(op);
