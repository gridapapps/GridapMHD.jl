using PartitionedArrays
using Gridap, GridapSolvers, GridapDistributed
using BlockArrays

using GridapSolvers.LinearSolvers
using GridapSolvers.MultilevelTools
using GridapSolvers.BlockSolvers 
using GridapSolvers.NonlinearSolvers

using Gridap.MultiField
using Gridap.Arrays, Gridap.Fields

using GridapMHD

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
mfs = BlockMultiFieldStyle(3,(2,1,1),(1,3,2,4))
Ωf  = Triangulation(model)
V_u = TestFESpace(model,reffe_u;dirichlet_tags="boundary")
V_p = TestFESpace(model,reffe_p)
V_j = TestFESpace(model,reffe_j;dirichlet_tags="boundary")
V_φ = TestFESpace(model,reffe_φ;conformity=:L2)
V = MultiFieldFESpace([V_u,V_p,V_j,V_φ];style=mfs)

params = Dict(
  :model=>model,
  :fluid=>Dict(
    :domain=>model,
    :α=>1.0,
    :β=>1.0,
    :γ=>1.0,
    :f=>VectorValue(0.0,0.0,0.0),
    :B=>VectorValue(0.0,0.0,1.0),
    :σ=>1.0,
    :ζ=>1.0,
  ),
  :solid=>nothing,
  :bcs => Dict(
    :u=>Dict(:tags=>"noslip"),
    :j=>Dict(:tags=>"insulating"),
    :φ=>[],
    :t=>[],
    :thin_wall=>[],
    :f =>[],
    :B =>[],
  ),
)

Ω = Triangulation(model)
dΩ = Measure(Ω,2*k)
x0 = zero(V)
res, jac = GridapMHD.weak_form(params,2)
biform(x,y) = jac(x0,x,y)

assem = SparseMatrixAssembler(V,V,FullyAssembledRows())
A = assemble_matrix(biform,assem,V,V)

a_Ip(p,v_p) = ∫(p*v_p)*dΩ
a_Iφ(φ,v_φ) = ∫(φ*v_φ)*dΩ
Ip = assemble_matrix(a_Ip,V_p,V_p)
Iφ = assemble_matrix(a_Iφ,V_φ,V_φ)

A_uj_p = A[Block(1,2)]
A_uj_φ = A[Block(1,3)]

using PartitionedArrays: matching_local_indices, matching_own_indices, matching_ghost_indices
using Gridap.Algebra

x = allocate_in_domain(A)
y = allocate_in_domain(A)
p = x[Block(2)]
φ = x[Block(3)]

mul!(y[Block(2)],Ip,p)
mul!(y[Block(3)],Ip,φ)

ghost_to_local(axes(Ip,2))
ghost_to_local(axes(p,1))

# NOTE: 
# If we assemble the global matrix using SubAssembledRows(), the ghost indices do not match.
# This is because column ghost indices include dofs which are needed to assemble ghost rows. 
# This is further proof that we need to do something about this....
