
using PartitionedArrays
using Gridap, GridapPETSc, GridapSolvers, GridapDistributed

using GridapSolvers.LinearSolvers
import GridapSolvers.PatchBasedSmoothers as PBS
using Gridap.ReferenceFEs

np = (2,2,1)
ranks = with_debug() do distribute
  distribute(LinearIndices((prod(np),)))
end

# Geometry
n = 4
model = CartesianDiscreteModel(ranks,np,(0,1,0,1,0,1),(n,n,n))

labels = get_face_labeling(model);
tags_j = append!(collect(1:20),[25,26])
add_tag_from_tags!(labels,"insulating",tags_j)
j_bc(x) = VectorValue(0.0,0.0,0.0)

# FESpaces
k = 2
reffe_j = ReferenceFE(raviart_thomas,Float64,k-1)
V_j = TestFESpace(model,reffe_j;dirichlet_tags="insulating")
U_j = TrialFESpace(V_j,j_bc)

γ = 1.0
Ω = Interior(model)
dΩ = Measure(Ω,2*k)
a_j(j,v_j) = ∫(γ*j⋅v_j + γ*(∇⋅j)⋅(∇⋅v_j))*dΩ
D_j = assemble_matrix(a_j,U_j,V_j)

# Patch solver
PD  = PBS.PatchDecomposition(model)
P_j = PBS.PatchFESpace(model,reffe_j,DivConformity(),PD,V_j)

Ωp = Triangulation(PD)
dΩp = Measure(Ωp,2*k)
ap_j(j,v_j) = ∫(γ*j⋅v_j + γ*(∇⋅j)⋅(∇⋅v_j))*dΩp

local_solver = LUSolver()
patch_solver = PatchBasedLinearSolver(ap_j,P_j,U_j,local_solver)
smoother = RichardsonSmoother(patch_solver,100,1.0)

ns_smoother = numerical_setup(symbolic_setup(smoother,D_j),D_j)

b = GridapSolvers.allocate_col_vector(D_j)
x = GridapSolvers.allocate_col_vector(D_j)
y = GridapSolvers.allocate_row_vector(D_j)

fill!(b,1.0)
y = b - D_j*x
solve!(x,ns_smoother,y)
err = norm(b - D_j*x)
