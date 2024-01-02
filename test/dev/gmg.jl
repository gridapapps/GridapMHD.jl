
using Gridap, GridapSolvers, GridapDistributed, PartitionedArrays
using GridapMHD

using Gridap.Algebra

np = 1
ranks = with_mpi() do distribute
  distribute(LinearIndices((np,)))
end

model = GridapMHD.Meshers.expansion_generate_mesh_hierarchy(ranks,1,[1,1]);

T = Float64
k = 1
D = num_cell_dims(get_model(model,1))
reffe_u = ReferenceFE(lagrangian,VectorValue{D,T},k)
reffe_j = ReferenceFE(raviart_thomas,T,k-1)

V_u = TestFESpace(model,reffe_u;dirichlet_tags=["inlet", "wall"])
V_j = TestFESpace(model,reffe_j;dirichlet_tags=["inlet", "outlet"])
V = MultiFieldFESpace([V_u,V_j]);

u_inlet((x,y,z)) = VectorValue(36.0*(y-1/4)*(y+1/4)*(z-1)*(z+1),0,0)
u_bc = [u_inlet, VectorValue(0.0, 0.0, 0.0)]
j_bc = [VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)]
U_u = TrialFESpace(V_u,u_bc);
U_j = TrialFESpace(V_j,j_bc);
U = MultiFieldFESpace([U_u,U_j]);

using GridapMHD: a_mhd_u_u, a_mhd_u_j, a_mhd_j_u, a_mhd_j_j, a_al_u_u, a_al_j_j
B = VectorValue(0.,0.,1.); σf = 1.0; β = 1.0; γ = 1.e3; ζ = 1.e3
function jac_uj(x,y,dΩ)
  u, j = x
  v_u, v_j = y
  # TODO: Add thin_wall bcs... how we deal with more than one triangulation?
  # TODO: Add nonlinear terms to u-u
  r = a_mhd_u_u(u,v_u,β,dΩ) + a_mhd_u_j(j,v_u,γ,B,dΩ) + a_mhd_j_u(u,v_j,σf,B,dΩ) + a_mhd_j_j(j,v_j,dΩ)
  if abs(ζ) > eps(typeof(ζ))
    r = r + a_al_u_u(u,v_u,ζ,dΩ) + a_al_j_j(j,v_j,ζ,dΩ)
  end
  return r
end

Ω = Triangulation(get_model(model,1))
dΩ = Measure(Ω,2)

a(x,y) = jac_uj(x,y,dΩ)
l((v_u,v_j)) = GridapMHD.ℓ_mhd_u(v_u,VectorValue(0.,0.,0.),dΩ)
op = AffineFEOperator(a,l,get_fe_space(U,1),get_fe_space(V,1))
A, b = get_matrix(op), get_vector(op);

nlevs = num_levels(model)
qdegree = map(lev -> 2*k,1:nlevs)
solver = GridapMHD.gmg_solver(model,U,V,jac_uj,qdegree)

ns = numerical_setup(symbolic_setup(solver,A),A)

x = allocate_in_domain(A)
fill!(x,0.0)
solve!(x,ns,b)

norm(b-A*x)
