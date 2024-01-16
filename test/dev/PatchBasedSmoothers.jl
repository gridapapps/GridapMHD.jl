using PartitionedArrays
using Gridap, GridapSolvers, GridapDistributed

using GridapSolvers.LinearSolvers
using GridapSolvers.MultilevelTools
using GridapSolvers.BlockSolvers 
using GridapSolvers.NonlinearSolvers
import GridapSolvers.PatchBasedSmoothers as PBS

using Gridap.MultiField, Gridap.Arrays, Gridap.Fields
using Gridap.ReferenceFEs, Gridap.FESpaces, Gridap.Algebra

using PartitionedArrays: getany

using GridapMHD

np = (1,1,1)
ranks = with_mpi() do distribute
  distribute(LinearIndices((prod(np),)))
end

num_refs = 1
#model = GridapMHD.Meshers.expansion_generate_mesh(ranks,num_refs)
model = GridapMHD.Meshers.expansion_generate_mesh_hierarchy(ranks,0,[1,1]);
PD = PBS.PatchDecomposition(model)

k = 1
T = Float64
D = 3
reffe_u = ReferenceFE(lagrangian,VectorValue{D,T},k)
reffe_j = ReferenceFE(raviart_thomas,T,k-1)

#conf_u = GridapSolvers.MultilevelTools._cell_conformity(model,reffe_u;conformity=H1Conformity())
#conf_j = GridapSolvers.MultilevelTools._cell_conformity(model,reffe_j;conformity=DivConformity())

V_u = TestFESpace(model,reffe_u;dirichlet_tags=["inlet","wall"])
V_j = TestFESpace(model,reffe_j;dirichlet_tags=["inlet","outlet"])
V = MultiFieldFESpace([V_u,V_j])

P = PBS.PatchFESpace(V,PD)#,[conf_u,conf_j])

PD = PD[1]
P = GridapSolvers.MultilevelTools.get_fe_space(P,1)
V = GridapSolvers.MultilevelTools.get_fe_space(V,1)

Ωp = Interior(PD)
Γp = Boundary(PD;tags="wall")

dΩp = Measure(Ωp,2*k)
dΓp = Measure(Γp,2*k-1)

n_Γ = get_normal_vector(Γp)
x0 = zero(V)

conv(u,∇u) = (∇u')⋅u
a((u0,j0),(u,j),(v,q),dΩp,dΓp,n_Γ) = ∫(u⋅v + j⋅q)*dΩp + ∫(v⋅((conv∘(u0,∇(u))) + (conv∘(u,∇(u0)))))*dΩp + ∫((q⋅n_Γ)*(j⋅n_Γ) + (q⋅n_Γ)*(n_Γ⋅(∇(j)⋅n_Γ)))*dΓp

a_glob(x,y) = a(x0,x,y,dΩp,dΓp,n_Γ)
Ap = assemble_matrix(a_glob,P,P)

Api = map(local_views(P),local_views(dΩp),local_views(dΓp),local_views(n_Γ),local_views(x0)) do P, dΩp, dΓp, n_Γ, x0
  a_loc(x,y) = a(x0,x,y,dΩp,dΓp,n_Γ)
  assem = SparseMatrixAssembler(P,P)
  assemble_matrix(a_loc,assem,P,P)
end


Pi = getany(local_views(P))
dΩpi = getany(local_views(dΩp))
dΓpi = getany(local_views(dΓp))
nΓpi = getany(local_views(n_Γ))
x0i = getany(local_views(x0))

ai(x,y) = a(x0i,x,y,dΩpi,dΓpi,nΓpi)
assem_i = SparseMatrixAssembler(Pi,Pi)

ui = get_trial_fe_basis(Pi)
vi = get_fe_basis(Pi)
contr = ai(ui,vi)
matdata = collect_cell_matrix(Pi,Pi,contr)

assemble_matrix(assem_i,matdata)

