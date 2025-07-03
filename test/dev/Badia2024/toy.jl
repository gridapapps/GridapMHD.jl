
using Gridap, GridapSolvers, GridapDistributed, PartitionedArrays, GridapPETSc
using BlockArrays, SparseArrays

using Gridap.Algebra, Gridap.FESpaces, Gridap.MultiField
using GridapDistributed: i_am_in
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools, GridapSolvers.PatchBasedSmoothers

using GridapP4est
using GridapMHD

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
ζ = 100.0
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

params = Dict{Symbol,Any}(
  :debug=>false,
  :solve=>true,
  :res_assemble=>false,
  :jac_assemble=>false,
)

t = PTimer(ranks,verbose=true);
params[:ptimer] = t;

params[:solver] = Dict(
  :solver         => :badia2024,
  :matrix_type    => SparseMatrixCSC{Float64,Int64},
  :vector_type    => Vector{Float64},
  :petsc_options  => "-ksp_error_if_not_converged true -ksp_converged_reason",
  :block_solvers  => [:gmg,:cg_jacobi,:cg_jacobi],
)
params[:fluid] = Dict(
  :domain=>nothing,
  :α=>α,
  :β=>β,
  :γ=>γ,
  :f=>f,
  :B=>B,
  :ζ=>ζ,
)
params[:bcs] = Dict(
  :u=>Dict(:tags=>"boundary",:values=>u_exact),
  :j=>Dict(:tags=>"boundary",:values=>j_exact),
)

domain = Tuple(vcat(repeat([0.0,1.0],Dc)))
mesh_partition = Tuple(fill(nc,Dc))
base_model = CartesianDiscreteModel(domain,mesh_partition)
mh = GridapMHD.Meshers.generate_mesh_hierarchy(ranks,base_model,0,[np,np]);
params[:multigrid] = Dict{Symbol,Any}(
  :mh => mh,
  :num_refs_coarse => 0,
  :ranks_per_level => [np,np],
)
model = get_model(mh,1)
params[:model] = model

params = GridapMHD.add_default_params(params);

U, V = GridapMHD._fe_spaces(params)

op = GridapMHD._fe_operator(U,V,params)
xh = zero(get_trial(op))
nlsolver = GridapMHD._solver(op,params)
solve!(xh,nlsolver,op);

A = jacobian(op,xh)
b = residual(op,xh)
global_ls = nlsolver.ls
x = mortar(map(ax -> pfill(0.0,partition(ax)),blocks(axes(A,2))))
#b = mortar(map(ax -> prandn(partition(ax)),blocks(axes(A,2))))
global_ns = numerical_setup(symbolic_setup(global_ls,A,x),A,x)
solve!(x,global_ns,b)

P_ns = global_ns.Pr_ns
gmg_ns = P_ns.block_ns[1]
b_uj = copy(blocks(b)[1])
b_uj = prand(partition(axes(blocks(A)[1,1],2)))
x_uj = pfill(0.0,partition(axes(blocks(A)[1,1],2)))
solve!(x_uj,gmg_ns,b_uj)
