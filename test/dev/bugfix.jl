
using Gridap
using Gridap.Geometry
using Gridap.MultiField
using Gridap.FESpaces
using Gridap.Algebra

using GridapDistributed, PartitionedArrays
using GridapGmsh

using GridapMHD

np = (2,2)
ranks = with_debug() do distribute
  distribute(LinearIndices((prod(np),)))
end

mesh = Dict{Symbol,Any}(
  :mesher => :gmsh,
  :base_mesh => "solid",
)

mesh = Dict{Symbol,Any}(
  :mesher => :p4est_SG,
  :base_mesh => "solid",
  :num_refs => 0,
)

params = Dict{Symbol,Any}(
  :debug=>false,
  :solve=>true,
  :res_assemble=>false,
  :jac_assemble=>false,
  :solver=> GridapMHD.default_solver_params(Val(:badia2024))
)
params[:fluid] = Dict{Symbol,Any}(
  :domain => "fluid",
  :α => 1.0,
  :β => 1.0,
  :γ => 0.0,
  :f => VectorValue(0.0,0.0,0.0),
  :B => VectorValue(0.0,1.0,0.0),
  :ζ => 1.0,
  :convection => :newton,
)
u_in = GridapMHD.u_inlet(:parabolic,1.0,4.0,1.0)
params[:bcs] = Dict{Symbol,Any}(
  :u => Dict(
    :tags => ["inlet", "wall"],
    :values => [u_in, VectorValue(0.0, 0.0, 0.0)]
  ),
  :j => Dict(
    :tags => ["inlet", "outlet"],
    :values => [VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)],
  ),
)
params[:solid] = Dict(:domain => "wall",:σ => 1.0)
params[:fespaces] = Dict{Symbol,Any}(
  :order_u => 2,
  :order_j => 3,
)

model = GridapMHD.expansion_mesh(mesh,ranks,params)

params = GridapMHD.add_default_params(params);
GridapMHD._setup_trians!(params)

U_u, V_u = GridapMHD._fe_space(Val(:u),params)
U_p, V_p = GridapMHD._fe_space(Val(:p),params)
U_j, V_j = GridapMHD._fe_space(Val(:j),params)
U_φ, V_φ = GridapMHD._fe_space(Val(:φ),params)

mfs = GridapMHD._multi_field_style(params)
V = MultiFieldFESpace([V_u,V_p,V_j,V_φ];style=mfs)
U = MultiFieldFESpace([U_u,U_p,U_j,U_φ];style=mfs)

op = GridapMHD._fe_operator(U,V,params)
xh = GridapMHD._allocate_solution(op,params)
r = residual(op,xh)
j = jacobian(op,xh)

u = get_trial_fe_basis(U);
v = get_fe_basis(V);

contr = jac(xh,u,v)
map(local_views(contr)) do a 
  for strian in Gridap.CellData.get_domains(a)
    cell_mat = Gridap.CellData.get_contribution(a,strian)
    #for mat in cell_mat
    #  println(mat)
    #end
    println(first(cell_mat))
  end
end

res, jac = GridapMHD.weak_form(params)

# res, jac = GridapMHD.weak_form(params)
# xh = zero(U);
# u = get_trial_fe_basis(U);
# v = get_fe_basis(V);
# 
# contr = jac(xh,u,v)
# map(local_views(contr)) do a 
#   for strian in Gridap.CellData.get_domains(a)
#     scell_mat = Gridap.CellData.get_contribution(a,strian)
#     for mat in scell_mat
#       println(mat)
#     end
#   end
# end
#cell_mat = collect_cell_matrix(U, V, contr)
