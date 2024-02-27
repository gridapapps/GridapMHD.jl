
using Gridap, GridapSolvers, GridapDistributed, PartitionedArrays, GridapPETSc
using BlockArrays, SparseArrays

using Gridap.Algebra, Gridap.FESpaces, Gridap.MultiField
using GridapDistributed: i_am_in
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools, GridapSolvers.PatchBasedSmoothers

using GridapP4est
using GridapMHD

function get_hierarchy_matrices(
  trials::FESpaceHierarchy,
  tests::FESpaceHierarchy,
  a::Function,
  qdegree::Integer;
  is_nonlinear::Bool=false
)
  nlevs = num_levels(trials)
  mh    = trials.mh

  mats = Vector{PSparseMatrix}(undef,nlevs)
  for lev in 1:nlevs
    parts = get_level_parts(mh,lev)
    if i_am_in(parts)
      model = get_model(mh,lev)
      U = GridapSolvers.get_fe_space(trials,lev)
      V = GridapSolvers.get_fe_space(tests,lev)
      Ω = Triangulation(model)
      dΩ = Measure(Ω,qdegree)
      ai(u,v) = is_nonlinear ? a(zero(U),u,v,dΩ) : a(u,v,dΩ)
      mats[lev] = assemble_matrix(ai,U,V)
    end
  end
  return mats
end

function get_patch_smoothers(tests,patch_spaces,patch_decompositions,biform,qdegree)
  mh = tests.mh
  nlevs = num_levels(mh)
  smoothers = Vector{RichardsonSmoother}(undef,nlevs-1)
  for lev in 1:nlevs-1
    parts = get_level_parts(mh,lev)
    if i_am_in(parts)
      PD = patch_decompositions[lev]
      Ph = GridapSolvers.get_fe_space(patch_spaces,lev)
      Vh = GridapSolvers.get_fe_space(tests,lev)
      Ω  = Triangulation(PD)
      dΩ = Measure(Ω,qdegree)
      ap(u,du,v) = biform(u,du,v,dΩ)
      patch_smoother = PatchBasedLinearSolver(ap,Ph,Vh;is_nonlinear=true)
      smoothers[lev] = RichardsonSmoother(patch_smoother,10,0.2)
    end
  end
  return smoothers
end

function u_exact_3d(x)
  VectorValue(x[1]^2,-2.0*x[1]*x[2],1.0)
end
function u_exact_2d(x)
  VectorValue(x[1]^2,-2.0*x[1]*x[2])
end
function j_exact_3d(x)
  VectorValue(-x[2]*x[1]^2,x[1]*x[2]^2,0.0)
end
function j_exact_2d(x)
  VectorValue(-x[2]*x[1]^2,x[1]*x[2]^2)
end

np = 1
Dc = 3
nc = 2
Re = 0.1
Ha = 10.0
η_u, η_j = 10.0,10.0
B = VectorValue(0.0,0.0,1.0)

u_exact(x) = (Dc==2) ? u_exact_2d(x) : u_exact_3d(x)
j_exact(x) = (Dc==2) ? j_exact_2d(x) : j_exact_3d(x)

α = 1.0
β = 1.0/Re
γ = Ha^2/Re
f(x) = α*u_exact(x)⋅∇(u_exact)(x) - β*Δ(u_exact)(x) - γ*(j_exact(x)×B)

ranks = with_mpi() do distribute
  distribute(LinearIndices((np,)))
end

domain = Tuple(vcat(repeat([0.0,1.0],Dc)))
mesh_partition = Tuple(fill(nc,Dc))
base_model = CartesianDiscreteModel(domain,mesh_partition)
mh = GridapMHD.Meshers.generate_mesh_hierarchy(ranks,base_model,0,[np,np])
model = get_model(mh,1)


order = 2
reffe_u  = ReferenceFE(lagrangian,VectorValue{Dc,Float64},order)
tests_u  = FESpace(mh,reffe_u;dirichlet_tags="boundary");
trials_u = TrialFESpace(tests_u);

reffe_j  = ReferenceFE(raviart_thomas,Float64,order-1)
tests_j  = FESpace(mh,reffe_j;dirichlet_tags="boundary");
trials_j = TrialFESpace(tests_j);

trials = MultiFieldFESpace([trials_u,trials_j]);
tests  = MultiFieldFESpace([tests_u,tests_j]);
spaces = tests, trials


conv(u,∇u) = (∇u')⋅u
a_al((u,j),(v_u,v_j),dΩ) = ∫(η_u*(∇⋅u)⋅(∇⋅v_u))*dΩ + ∫(η_j*(∇⋅j)⋅(∇⋅v_j))*dΩ
a_mhd((u,j),(v_u,v_j),dΩ) = ∫(β*∇(u)⊙∇(v_u) -γ*(j×B)⋅v_u + j⋅v_j - (u×B)⋅v_j)dΩ
c_mhd((u,j),(v_u,v_j),dΩ) = ∫( α*v_u⋅(conv∘(u,∇(u))) ) * dΩ
dc_mhd((u,j),(du,dj),(v_u,v_j),dΩ) = ∫(α*v_u⋅( (conv∘(u,∇(du))) + (conv∘(du,∇(u)))))dΩ
rhs((u,j),(v_u,v_j),dΩ) = ∫(f⋅v_u)dΩ

jac(x0,x,y,dΩ) = a_mhd(x,y,dΩ) + a_al(x,y,dΩ) + dc_mhd(x0,x,y,dΩ)
res(x0,y,dΩ) = a_mhd(x0,y,dΩ) + a_al(x0,y,dΩ) + c_mhd(x0,y,dΩ) - rhs(x0,y,dΩ)

qdegree = 2*(order+1)
patch_decompositions = PatchDecomposition(mh)
patch_spaces = PatchFESpace(tests,patch_decompositions);
smoothers = get_patch_smoothers(tests,patch_spaces,patch_decompositions,jac,qdegree)

smatrices = get_hierarchy_matrices(trials,tests,jac,qdegree;is_nonlinear=true);
A = smatrices[1]

dΩ = Measure(Triangulation(get_model(mh,1)),qdegree)
x0 = zero(GridapSolvers.get_fe_space(trials,1))
b = assemble_vector(v -> res(x0,v,dΩ), GridapSolvers.get_fe_space(tests,1))

coarse_solver = LUSolver()
restrictions, prolongations = setup_transfer_operators(trials,
                                                        qdegree;
                                                        mode=:residual,
                                                        solver=LUSolver());

#######################################
# A) GMG as solver 

gmg_solver = GMGLinearSolver(mh,
                      smatrices,
                      prolongations,
                      restrictions,
                      pre_smoothers=smoothers,
                      post_smoothers=smoothers,
                      coarsest_solver=LUSolver(),
                      maxiter=20,
                      rtol=1.0e-8,
                      verbose=true,
                      mode=:preconditioner)
gmg_solver.log.depth += 1
gmg_ns = numerical_setup(symbolic_setup(gmg_solver,A),A)

x = pfill(0.0,partition(axes(A,2)))
r = b - A*x
solve!(x,gmg_ns,r)

############################################
# B) GMG as preconditioner for GMRES

gmg = GMGLinearSolver(mh,
                      smatrices,
                      prolongations,
                      restrictions,
                      pre_smoothers=smoothers,
                      post_smoothers=smoothers,
                      coarsest_solver=LUSolver(),
                      maxiter=3,
                      rtol=1.0e-8,
                      verbose=true,
                      mode=:preconditioner)
gmg.log.depth += 1

gmres_solver = FGMRESSolver(10,gmg;m_add=5,maxiter=30,rtol=1.0e-6,verbose=i_am_main(ranks))
gmres_ns = numerical_setup(symbolic_setup(gmres_solver,A),A)

x = pfill(0.0,partition(axes(A,2)))
solve!(x,gmres_ns,b)
