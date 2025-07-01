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

Dc = 3
B0 = VectorValue(0.,0.,1.)
u(x) = VectorValue(x[1],-x[2],0.)
nc = Tuple(fill(4,Dc))
np = Tuple(fill(1,Dc))
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

α = 1.e4
β = 1.e4
Π_Qh = LocalProjectionMap(divergence,reffe_p,qdegree)
utb(x) = u(x)×B0
biform(u,v,dΩ) = ∫(∇(v)⊙∇(u) + α*(∇⋅v)⋅Π_Qh(u) + β*(u×B0)⋅(v×B0))dΩ 
liform(v,dΩ) = ∫(∇(v)⊙∇(u) + β*(utb⋅(v×B0)))dΩ 

Ω = Triangulation(model)
dΩ = Measure(Ω,qdegree)

a(u,v) = biform(u,v,dΩ)
l(v) = liform(v,dΩ)
op = AffineFEOperator(a,l,U,V)
A, b = get_matrix(op), get_vector(op);

biforms = map(mhl -> get_bilinear_form(mhl,biform,qdegree),mh)
patch_decompositions = PatchDecomposition(mh)
smoothers = get_patch_smoothers(
  mh,tests_u,biform,patch_decompositions,qdegree
)
prolongations = setup_patch_prolongation_operators(
  tests_u,biform,biform,qdegree
)
restrictions = setup_patch_restriction_operators(
  tests_u,prolongations,biform,qdegree
)
gmg = GMGLinearSolver(
  trials_u,tests_u,biforms,
  prolongations,restrictions,
  pre_smoothers=smoothers,
  post_smoothers=smoothers,
  coarsest_solver=LUSolver(),
  maxiter=4,mode=:preconditioner,verbose=i_am_main(parts)
)
gmg.log.depth = 4

solver = FGMRESSolver(20,gmg;atol=1e-10,rtol=1.e-12,verbose=i_am_main(parts))
ns = numerical_setup(symbolic_setup(solver,A),A)

x = allocate_in_domain(A)
fill!(x,0.0)
solve!(x,ns,b)
xh = FEFunction(U,x)

eh = xh - u
err_l2 = sqrt(sum(∫(eh⋅eh)dΩ))


