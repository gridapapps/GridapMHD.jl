
using Gridap, GridapSolvers, GridapDistributed, PartitionedArrays, GridapPETSc
using BlockArrays, SparseArrays

using Gridap.Algebra, Gridap.FESpaces, Gridap.MultiField
using GridapDistributed: i_am_in
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools, GridapSolvers.PatchBasedSmoothers

using GridapP4est

using GridapMHD

np = 1
ranks = with_mpi() do distribute
  distribute(LinearIndices((np,)))
end

params = Dict{Symbol,Any}(
  :debug=>false,
  :solve=>true,
  :res_assemble=>false,
  :jac_assemble=>false,
)

t = PTimer(ranks,verbose=true);
params[:ptimer] = t;

L = 1.0
u0 = 1.0
f = VectorValue(0.0,0.0,1.0)
B = VectorValue(0.0,10.0,0.0)
B0 = norm(B)

Re = u0*L
Ha = B0*L
N = Ha^2/Re

params[:solver] = Dict(
  :solver         => :badia2024,
  :matrix_type    => SparseMatrixCSC{Float64,Int64},
  :vector_type    => Vector{Float64},
  :petsc_options  => "-ksp_error_if_not_converged true -ksp_converged_reason",
  :block_solvers  => [:gmg,:cg_jacobi,:cg_jacobi],
)
params[:fluid] = Dict(
  :domain=>nothing,
  :α=>1.0,
  :β=>1.0/Re,
  :γ=>N,
  :f=>(L/(u0^2))*f,
  :B=>(1/B0)*B,
  :ζ=>100.0,
)
params[:bcs] = Dict(
  :u=>Dict(:tags=>"noslip"),
  :j=>Dict(:tags=>"insulating"),
)

model = GridapMHD.hunt_mesh(ranks,params,(2,2),(1,1),L,Ha,1.0,1.0,false,[1,1])

params = GridapMHD.add_default_params(params);

U, V = GridapMHD._fe_spaces(params)

op = GridapMHD._fe_operator(U,V,params)
xh = zero(get_trial(op))
nlsolver = GridapMHD._solver(op,params)
solve!(xh,nlsolver,op);

A = jacobian(op,xh)
b = residual(op,xh)
global_ls = nlsolver.ls
x = get_free_dof_values(xh)
global_ns = numerical_setup(symbolic_setup(global_ls,A,x),A,x)
solve!(x,global_ns,b)


P_ns = global_ns.Pr_ns
gmg_ns = P_ns.block_ns[1]

b_uj = copy(blocks(b)[1])
x_uj = copy(blocks(x)[1])
solve!(x_uj,gmg_ns,b_uj)



###########################################################################################

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
      U = GridapSolvers.MultilevelTools.get_fe_space(trials,lev)
      V = GridapSolvers.MultilevelTools.get_fe_space(tests,lev)
      Ω = Triangulation(model)
      dΩ = Measure(Ω,qdegree)
      ai(u,v) = is_nonlinear ? a(zero(U),u,v,dΩ) : a(u,v,dΩ)
      mats[lev] = assemble_matrix(ai,U,V)
    end
  end
  return mats
end

function get_patch_smoothers(tests,patch_decompositions,biform,qdegree)
  mh = tests.mh
  patch_spaces = PatchFESpace(tests,patch_decompositions);
  nlevs = num_levels(mh)
  smoothers = Vector{RichardsonSmoother}(undef,nlevs-1)
  for lev in 1:nlevs-1
    parts = get_level_parts(mh,lev)
    if i_am_in(parts)
      PD = patch_decompositions[lev]
      Ph = GridapSolvers.MultilevelTools.get_fe_space(patch_spaces,lev)
      Vh = GridapSolvers.MultilevelTools.get_fe_space(tests,lev)
      Ω  = Triangulation(PD)
      dΩ = Measure(Ω,qdegree)
      ap(u,du,v) = biform(u,du,v,dΩ)
      patch_smoother = PatchBasedLinearSolver(ap,Ph,Vh;is_nonlinear=true)
      smoothers[lev] = RichardsonSmoother(patch_smoother,10,0.2)
    end
  end
  return smoothers
end

function get_mesh_hierarchy(parts,cmodel,num_refs_coarse,np_per_level)
  num_levels   = length(np_per_level)
  cparts       = generate_subparts(parts,np_per_level[num_levels])
  coarse_model = OctreeDistributedDiscreteModel(cparts,cmodel,num_refs_coarse)
  mh = ModelHierarchy(parts,coarse_model,np_per_level)
  return mh
end

function add_hunt_tags!(model)
  labels = get_face_labeling(model)
  tags_u = append!(collect(1:20),[23,24,25,26])
  tags_j = append!(collect(1:20),[25,26])
  add_tag_from_tags!(labels,"noslip",tags_u)
  add_tag_from_tags!(labels,"insulating",tags_j)
end

α = 1.0
β = 1.0/Re
γ = N
ζ = params[:fluid][:ζ]
_B = (1/B0)*B
η_u = η_j = params[:fluid][:ζ]
biform((u,j),(du,dj),(v_u,v_j),dΩ) = GridapMHD.a_mhd_u_u(du,v_u,β,dΩ) + GridapMHD.a_mhd_u_j(dj,v_u,γ,_B,dΩ) + 
                                  GridapMHD.a_mhd_j_u(du,v_j,1.0,B,dΩ) + GridapMHD.a_mhd_j_j(dj,v_j,dΩ) + GridapMHD.dc_mhd_u_u(u,du,v_u,α,dΩ) +
                                  GridapMHD.a_al_u_u(du,v_u,ζ,dΩ) + GridapMHD.a_al_j_j(dj,v_j,ζ,dΩ)

conv(u,∇u) = (∇u')⋅u
a_al((u,j),(v_u,v_j),dΩ) = ∫(η_u*(∇⋅u)⋅(∇⋅v_u))*dΩ + ∫(η_j*(∇⋅j)⋅(∇⋅v_j))*dΩ
a_mhd((u,j),(v_u,v_j),dΩ) = ∫(β*∇(u)⊙∇(v_u) -γ*(j×_B)⋅v_u + j⋅v_j - (u×_B)⋅v_j)dΩ
c_mhd((u,j),(v_u,v_j),dΩ) = ∫( α*v_u⋅(conv∘(u,∇(u))) ) * dΩ
dc_mhd((u,j),(du,dj),(v_u,v_j),dΩ) = ∫(α*v_u⋅((conv∘(u,∇(du))) + (conv∘(du,∇(u)))))dΩ
rhs((u,j),(v_u,v_j),dΩ) = ∫(f⋅v_u)dΩ

jac(x0,x,y,dΩ) = a_mhd(x,y,dΩ) + a_al(x,y,dΩ) + dc_mhd(x0,x,y,dΩ)
res(x0,y,dΩ) = a_mhd(x0,y,dΩ) + a_al(x0,y,dΩ) + c_mhd(x0,y,dΩ) - rhs(x0,y,dΩ)

if false
  mh = params[:multigrid][:mh]
else
  domain = (-1.0,1.0,-1.0,1.0,0.0,1.0)
  cmodel = CartesianDiscreteModel(domain,(4,4,4);isperiodic=(false,false,true))
  add_hunt_tags!(cmodel)
  mh = get_mesh_hierarchy(ranks,cmodel,0,[1,1]);
end

if false
  tests = MultiFieldFESpace([params[:multigrid][:tests][[1,3]]...])
  trials = MultiFieldFESpace([params[:multigrid][:trials][[1,3]]...])
else
  Dc = 3
  order = 2
  reffe_u  = ReferenceFE(lagrangian,VectorValue{Dc,Float64},order)
  tests_u  = FESpace(mh,reffe_u;dirichlet_tags="noslip");
  trials_u = TrialFESpace(tests_u);

  reffe_j  = ReferenceFE(raviart_thomas,Float64,order-1)
  tests_j  = FESpace(mh,reffe_j;dirichlet_tags="insulating");
  trials_j = TrialFESpace(tests_j);

  trials = MultiFieldFESpace([trials_u,trials_j]);
  tests  = MultiFieldFESpace([tests_u,tests_j]);
end

patch_decompositions = PatchDecomposition(mh)
smoothers = get_patch_smoothers(trials,patch_decompositions,jac,2*order+1)

smatrices = get_hierarchy_matrices(trials,tests,jac,2*order+1;is_nonlinear=true);
A_uj = smatrices[1]

dΩ = Measure(Triangulation(get_model(mh,1)),2*order+1)
x0 = zero(GridapSolvers.get_fe_space(trials,1))
b_uj = assemble_vector(v -> res(x0,v,dΩ),GridapSolvers.get_fe_space(tests,1))

coarse_solver = LUSolver()
restrictions, prolongations = setup_transfer_operators(tests,
                                                        2*order+1;
                                                        mode=:residual,
                                                        solver=LUSolver()); 

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
gmg_ns = numerical_setup(symbolic_setup(gmg_solver,A_uj),A_uj)

x = pfill(0.0,partition(axes(A_uj,2)))
r = b_uj - A_uj*x
r = prandn(partition(axes(A_uj,2)))
solve!(x,gmg_ns,r)


using Gridap, Gridap.Geometry
using GridapMHD

base = CartesianDiscreteModel((-1.0,1.0,-1.0,1.0,0.0,0.1),(4,4,3);isperiodic=(false,false,true))
writevtk(base,"data/hunt_meshes/base")

Ha = 1.0
base_blt_1 = GridapMHD.Meshers.hunt_generate_base_mesh((4,4),1.0,Ha,1.0,1.0,true)
writevtk(base_blt_1,"data/hunt_meshes/base_blt_1")
base_noblt_1 = GridapMHD.Meshers.hunt_generate_base_mesh((4,4),1.0,Ha,1.0,1.0,false)
writevtk(base_noblt_1,"data/hunt_meshes/base_noblt_1")

Ha = 10.0
base_blt_10 = GridapMHD.Meshers.hunt_generate_base_mesh((4,4),1.0,Ha,1.0,1.0,true)
writevtk(base_blt_10,"data/hunt_meshes/base_blt_10")
base_noblt_10 = GridapMHD.Meshers.hunt_generate_base_mesh((4,4),1.0,Ha,1.0,1.0,false)
writevtk(base_noblt_10,"data/hunt_meshes/base_noblt_10")

Gridap.Geometry.get_node_coordinates(base_blt_10)
