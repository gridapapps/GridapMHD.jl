
using PartitionedArrays
using Gridap, GridapPETSc, GridapSolvers, GridapDistributed

using GridapSolvers.LinearSolvers
import GridapSolvers.PatchBasedSmoothers as PBS
using Gridap.ReferenceFEs

function test_smoother(s,D_j)
  ns = numerical_setup(symbolic_setup(s,D_j),D_j)

  b = GridapSolvers.allocate_in_domain(D_j)
  x = GridapSolvers.allocate_in_domain(D_j)
  y = GridapSolvers.allocate_in_range(D_j)

  fill!(b,1.0)
  y = b - D_j*x
  solve!(x,ns,y)
  err = norm(b - D_j*x)

  return err
end

function test_solver(s,D_j)
  ns = numerical_setup(symbolic_setup(s,D_j),D_j)

  b = GridapSolvers.allocate_in_domain(D_j)
  x = GridapSolvers.allocate_in_domain(D_j)

  fill!(b,1.0)
  solve!(x,ns,b)
  err = norm(b - D_j*x)

  return err
end

np = (2,2,1)
ranks = with_debug() do distribute
  distribute(LinearIndices((prod(np),)))
end

# Geometry
n = 16
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

#smoother = RichardsonSmoother(patch_solver,10,1.0)
#test_smoother(smoother,D_j)

solver = GMRESSolver(100,patch_solver;tol=1e-6,verbose=true)
test_solver(solver,D_j)
