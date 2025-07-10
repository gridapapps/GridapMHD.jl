
using LinearAlgebra
using FillArrays

using Gridap
using Gridap.ReferenceFEs, Gridap.Algebra, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.MultiField, Gridap.Algebra, Gridap.Adaptivity
using PartitionedArrays
using GridapDistributed

using GridapSolvers
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools, GridapSolvers.PatchBasedSmoothers

using GridapMHD: PatchModel, get_cell_size

function get_patch_smoothers(mh,tests,biform;w=0.2)
  nlevs = num_levels(mh)
  smoothers = map(view(mh,1:nlevs-1),view(tests,1:nlevs-1)) do mhl,tests
    model = get_model(mhl)
    if isa(model,Union{Gridap.Adaptivity.AdaptedDiscreteModel,GridapDistributed.DistributedAdaptedDiscreteModel})
      model = get_model(model)
    end
    Vh = get_fe_space(tests)
    ptopo = Geometry.PatchTopology(ReferenceFE{0},model)
    ap = biform(PatchModel(model,ptopo))
    solver = PatchBasedSmoothers.PatchSolver(ptopo,Vh,Vh,ap;assembly=:star,collect_factorizations=true)
    if w > 0.0
      return RichardsonSmoother(solver,10,w)
    else
      return FGMRESSolver(10,solver;maxiter=10)
    end
  end
  return smoothers
end

function get_bilinear_form(mh_lev,biform)
  model = get_model(mh_lev)
  return biform(model)
end

np = (2,2)
np_per_level = [np,np]
parts = with_debug() do distribute
  distribute(LinearIndices((prod(np),)))
end

nc = (4,4)
Dc = length(nc)
domain = (Dc == 2) ? (0,1,0,1) : (0,1,0,1,0,1)
# mh = CartesianModelHierarchy(parts,np_per_level,domain,nc)

cmodel = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
fmodel = refine(cmodel)
mh = ModelHierarchy([fmodel,cmodel])

model = get_model(mh,1)

order = 2
qdegree = 2*order
reffe = ReferenceFE(raviart_thomas,Float64,order-1)

μ = 10.0
β = 1.0
γ = 1000.0
ζᵤ = γ
ζⱼ = γ
B0 = VectorValue(0.0,0.0,1.0)

u_exact(x) = VectorValue(x[1],-2*x[2],x[3])
j_exact(x) = VectorValue(-2*x[1]+x[2],x[2],-x[3])

f(x) = -Δ(u_exact)(x) - γ*j_exact(x)×B0
g(x) = j_exact(x) - u_exact(x)×B0

function biform(model)

  Ω = Triangulation(model)
  Γ = Boundary(model)
  Λ = Skeleton(model)

  h_Γ = get_cell_size(Γ)
  h_Λ = get_cell_size(Λ)
  n_Γ = get_normal_vector(Γ)
  n_Λ = get_normal_vector(Λ)

  dΩ = Measure(Ω,qdegree)
  dΓ = Measure(Γ,qdegree)
  dΛ = Measure(Λ,qdegree)

  function jac(dx,dy)
    du, dj = dx
    v, s = dy

    ∇du, ∇v = ∇(du), ∇(v)
    div_du, div_dj = divergence(du), divergence(dj)
    div_v, div_s = divergence(v), divergence(s)

    uᵗ, vᵗ = jump(du⊗n_Λ), jump(v⊗n_Λ)
    αΛ, αΓ = μ/h_Λ, μ/h_Γ

    u_block = β*(∇du⊙∇v) + ζᵤ*(div_du*div_v)
    j_block = dj⋅s + ζⱼ*(div_dj*div_s)
    c = ∫(u_block - γ*(dj×B0)⋅v + j_block - (du×B0)⋅s )dΩ
    c += ∫(αΛ*vᵗ⊙uᵗ - vᵗ⊙mean(∇du) - mean(∇v)⊙uᵗ)dΛ
    c += ∫(αΓ*v⋅du - v⋅(∇du⋅n_Γ) - (∇v⋅n_Γ)⋅du)dΓ
    return c
  end

  return jac
end

liform((v,s),dΩ) = ∫(v⋅f + s⋅g)dΩ

tests_u  = TestFESpace(mh,reffe,dirichlet_tags=["boundary"]);
trials_u = TrialFESpace(tests_u,[u_exact]);
tests_j = TestFESpace(mh,reffe,dirichlet_tags=["boundary"]);
trials_j = TrialFESpace(tests_j,[j_exact]);
tests = MultiFieldFESpace([tests_u, tests_j])
trials = MultiFieldFESpace([trials_u, trials_j])
U, V = get_fe_space(trials,1), get_fe_space(tests,1)

Ω = Triangulation(model)
dΩ = Measure(Ω,qdegree)
a = biform(model)
l(v) = liform(v,dΩ)
op = AffineFEOperator(a,l,U,V)
A, b = get_matrix(op), get_vector(op);

uh, jh = solve(op)

eh_u = u_exact - uh
eh_j = j_exact - jh
error_u = sum(∫(eh_u ⋅ eh_u)dΩ)
error_j = sum(∫(eh_j ⋅ eh_j)dΩ)

biforms = map(mhl -> get_bilinear_form(mhl,biform), mh)

w = 0.2
smoothers = get_patch_smoothers(
  mh,tests,biform;w
)
prolongations = setup_prolongation_operators(
  tests,qdegree;mode=:residual
)
restrictions = setup_restriction_operators(
  tests,qdegree;mode=:residual
)

gmg = GMGLinearSolver(
  trials,tests,biforms,
  prolongations,restrictions,
  pre_smoothers=smoothers,
  post_smoothers=smoothers,
  coarsest_solver=LUSolver(),
  maxiter=1,mode=:preconditioner,
  verbose=i_am_main(parts),
);

solver = FGMRESSolver(10,gmg;verbose=i_am_main(parts),atol=1e-6,rtol=1e-10)
ns = numerical_setup(symbolic_setup(solver,A),A)

x = Algebra.allocate_in_domain(A)
fill!(x,0.0)
solve!(x,ns,b)
