using Gridap

using Test
using LinearAlgebra
using FillArrays, BlockArrays

using Gridap
using Gridap.ReferenceFEs, Gridap.Algebra, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.MultiField, Gridap.Algebra
using PartitionedArrays
using GridapDistributed

using GridapSolvers
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools, GridapSolvers.PatchBasedSmoothers

function get_patch_smoothers(mh,tests,biform,patch_decompositions,qdegree)
  patch_spaces = PatchFESpace(tests,patch_decompositions)
  nlevs = num_levels(mh)
  smoothers = map(view(tests,1:nlevs-1),patch_decompositions,patch_spaces) do tests, PD, Ph
    Vh = get_fe_space(tests)
    Ω  = Triangulation(PD)
    dΩ = Measure(Ω,qdegree)
    ap = (u,v) -> biform(u,v,dΩ)
    patch_smoother = PatchBasedLinearSolver(ap,Ph,Vh)
    return RichardsonSmoother(patch_smoother,10,0.2)
  end
  return smoothers
end

function get_bilinear_form(mh_lev,biform,qdegree)
  model = get_model(mh_lev)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,qdegree)
  return (u,v) -> biform(u,v,dΩ)
end

u(x) = VectorValue(x[1],-x[2])
nc = (4,4)

np = (1,1)
np_per_level = fill(np,2)
parts = with_mpi() do distribute
  distribute(LinearIndices((prod(np),)))
end

Dc = length(nc)
domain = (Dc == 2) ? (0,1,0,1) : (0,1,0,1,0,1)
mh = CartesianModelHierarchy(parts,np_per_level,domain,nc)
model = get_model(mh,1)

order = 2
qdegree = 2*(order+1)
reffe_u = ReferenceFE(lagrangian,VectorValue{Dc,Float64},order)
reffe_p = ReferenceFE(lagrangian,Float64,order-1;space=:P)

tests_u  = TestFESpace(mh,reffe_u,dirichlet_tags="boundary");
trials_u = TrialFESpace(tests_u,u);
U, V = get_fe_space(trials_u,1), get_fe_space(tests_u,1)
Q = TestFESpace(model,reffe_p;conformity=:L2,constraint=:zeromean) 

α = 1.0
β = 1.0
Π_Qh = LocalProjectionMap(divergence,reffe_p,qdegree)

biform(u,v,dΩ) = ∫(∇(v)⊙∇(u) + α*(∇⋅v)⋅Π_Qh(u) + β*(u×B)⋅(v×B))dΩ 

Ω = Triangulation(model)
dΩ = Measure(Ω,qdegree)

a(u,v) = biform(u,v,dΩ)
l(v) = liform(v,dΩ)
op = AffineFEOperator(a,l,X,Y)
A, b = get_matrix(op), get_vector(op);

biforms = map(mhl -> get_bilinear_form(mhl,biform_u,qdegree),mh)
patch_decompositions = PatchDecomposition(mh)
smoothers = get_patch_smoothers(
  mh,tests_u,biform_u,patch_decompositions,qdegree
)
prolongations = setup_patch_prolongation_operators(
  tests_u,biform_u,graddiv,qdegree
)
restrictions = setup_patch_restriction_operators(
  tests_u,prolongations,graddiv,qdegree;solver=CGSolver(JacobiLinearSolver())
)
gmg = GMGLinearSolver(
  trials_u,tests_u,biforms,
  prolongations,restrictions,
  pre_smoothers=smoothers,
  post_smoothers=smoothers,
  coarsest_solver=LUSolver(),
  maxiter=4,mode=:preconditioner,verbose=i_am_main(parts)
)

solver_u = gmg
solver_p = CGSolver(JacobiLinearSolver();maxiter=20,atol=1e-14,rtol=1.e-6,verbose=i_am_main(parts))
solver_u.log.depth = 2
solver_p.log.depth = 2

diag_blocks  = [LinearSystemBlock(),BiformBlock((p,q) -> ∫(-1.0/α*p*q)dΩ,Q,Q)]
bblocks = map(CartesianIndices((2,2))) do I
  (I[1] == I[2]) ? diag_blocks[I[1]] : LinearSystemBlock()
end
coeffs = [1.0 1.0;
          0.0 1.0]  
P = BlockTriangularSolver(bblocks,[solver_u,solver_p],coeffs,:upper)
solver = FGMRESSolver(20,P;atol=1e-10,rtol=1.e-12,verbose=i_am_main(parts))
ns = numerical_setup(symbolic_setup(solver,A),A)

x = allocate_in_domain(A); fill!(x,0.0)
solve!(x,ns,b)
xh = FEFunction(X,x)
