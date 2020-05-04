module GridapMHD

using Gridap

# using Gridap: ∇, divergence

using Gridap.Arrays: Table, LocalToGlobalArray
using Gridap.Geometry: GridTopology
using Gridap.Geometry: CartesianDescriptor
using Gridap.Geometry: UnstructuredGrid
using Gridap.Geometry: FaceLabeling
using Gridap.Geometry: num_nodes, get_cell_nodes, get_cell_type,get_reffes
using Gridap.Geometry: num_faces, _generate_cell_to_vertices
using Gridap.Geometry: _generate_grid_topology_from_grid
using Gridap.Geometry: _generate_cell_to_vertices
using Gridap.Geometry: _generate_cell_to_vertices_count!
using Gridap.Geometry: _generate_cell_to_vertices_fill!
using Gridap.Geometry: _fill_cartesian_face_labeling!
import Gridap.Geometry: _generate_grid_topology_from_grid
import Gridap.Geometry: _generate_cell_to_vertices
import Gridap.Geometry: _generate_cell_to_vertices_count!
import Gridap.Geometry: _generate_cell_to_vertices_fill!
import Gridap.Geometry: _fill_cartesian_face_labeling!
import Gridap.Geometry: CartesianDescriptor
import Gridap.Geometry: CartesianDiscreteModel
import Gridap.Geometry: UnstructuredGridTopology

using Gridap.FESpaces: get_local_item, find_local_index, length_to_ptrs!
using Gridap.FESpaces: _generate_face_to_own_dofs_count_d!
using Gridap.FESpaces: _generate_face_to_own_dofs
import Gridap.FESpaces: _generate_face_to_own_dofs

export main
export analytical_solution

# Geometry/CartesianDiscreteModels.jl:1:31
#
# struct CartesianDiscreteModel2{D,T,F} <: CartesianDiscreteModel{D,T,F}
#   # grid::CartesianGrid{D,T,F}
#   # grid_topology::UnstructuredGridTopology{D,D,T,true}
#   # face_labeling::FaceLabeling
#   # @doc """
#   #     CartesianDiscreteModel(desc::CartesianDescriptor)
#   #
#   # Inner constructor
#   # """
#   function CartesianDiscreteModel2(grid::CartesianGrid{D,T,F},topo::UnstructuredGridTopology{D,D,T,true},labels::FaceLabeling) where {D,T,F}
#     new{D,T,F}(grid,topo,labels)
#   end
# end
#
#
# function CartesianDiscreteModel2(desc::CartesianDescriptor{D,T,F}) where {D,T,F}
#   grid = CartesianGrid(desc)
#   _grid = UnstructuredGrid(grid)
#   topo = UnstructuredGridTopology(_grid)
#   nfaces = [num_faces(topo,d) for d in 0:num_cell_dims(topo)]
#   labels = FaceLabeling(nfaces)
#   _fill_cartesian_face_labeling!(labels,topo)
#   CartesianDiscreteModel2(grid,topo,labels)
# end


function CartesianDescriptor(domain,partition,periodic::Array{Int,1},map::Function=identity)
  D = length(partition)
  limits = [(domain[2*d-1],domain[2*d]) for d in 1:D]
  sizes = [(limits[d][2]-limits[d][1])/partition[d] for d in 1:D]
  origin = [ limits[d][1] for d in 1:D]
  return CartesianDescriptor(origin,sizes,partition,map), periodic
end

function CartesianDiscreteModel(args::Tuple{CartesianDescriptor{D,T,F},Array{Int,1}}) where {D,T,F}
  desc, periodic = args
  grid = CartesianGrid(desc)
  _grid = UnstructuredGrid(grid)
  topo = UnstructuredGridTopology(_grid,desc.partition, periodic)
  nfaces = fill(0,num_cell_dims(topo)+1)
  for d in 0:num_cell_dims(topo)
    nfaces[d+1] = num_faces(topo,d)
  end
  labels = FaceLabeling(nfaces)
  _fill_cartesian_face_labeling!(labels,topo)
  CartesianDiscreteModel(grid,topo,labels)
end

function ijk_to_index(ijk::Array{Int,1},partition)
  index = ijk[1] -1
  for (i,ind) in enumerate(ijk[2:end])
    index += (ind-1) * prod(partition[1:i])
  end
  index + 1
end

function index_to_ijk(index::Int,partition)
  D = length(partition)
  ijk = Array{Int,1}(undef,D)
  ijk .= 0
  i = index -1
  for d in D:-1:1
    ijk[d] = i ÷ prod(partition[1:d-1])
    i = i % prod(partition[1:d-1])
  end
  ijk .+ 1
end

function UnstructuredGridTopology(grid::UnstructuredGrid,partition,periodic::Array{Int,1})
  cell_to_vertices, vertex_to_node, = _generate_cell_to_vertices_from_grid(grid,partition,periodic)
  _generate_grid_topology_from_grid(grid,cell_to_vertices,vertex_to_node)
end

function _generate_cell_to_vertices_from_grid(grid::UnstructuredGrid,partition,periodic::Array{Int,1})
  if is_first_order(grid)
    nodes = get_cell_nodes(grid)
    cell_to_vertices = copy(nodes)

    nnodes = num_nodes(grid)

    # compute periodic nodes
    D = length(partition)
    num_nodes_x_dir = Array{Int,1}(undef,D)
    num_nodes_x_dir = [partition[i]+1 for i in 1:D]

    point_to_isperiodic = fill(false,nnodes)
    slave_point_to_master_point = Array{Int,1}(undef,0)

    for periodic_dir in periodic
      free_dirs = collect(1:D)
      deleteat!(free_dirs, periodic_dir)

      slave_point_ijk = Array{Int,1}(undef,D)
      master_point_ijk = Array{Int,1}(undef,D)
      slave_point_ijk .= 1
      slave_point_ijk[free_dirs[1]] = 0
      slave_point_ijk[periodic_dir] = num_nodes_x_dir[periodic_dir]
      for i in collect(1:prod(num_nodes_x_dir[free_dirs]))
        slave_point_ijk[free_dirs[1]] += 1
        for j in collect(1:D-1)
          if slave_point_ijk[free_dirs[j]] > num_nodes_x_dir[free_dirs[j]]
            slave_point_ijk[free_dirs[j+1]] += 1
            slave_point_ijk[free_dirs[j]] = 1
          end
        end
        slave_point_index = ijk_to_index(slave_point_ijk,num_nodes_x_dir)
        point_to_isperiodic[slave_point_index] = true
        master_point_ijk .= slave_point_ijk
        master_point_ijk[periodic_dir] = 1
        append!(slave_point_to_master_point, ijk_to_index(master_point_ijk,num_nodes_x_dir))
      end
    end

    vertex_to_point = findall( .! point_to_isperiodic)
    slave_point_to_point = findall( point_to_isperiodic)

    point_to_vertex = fill(-1,length(point_to_isperiodic))
    point_to_vertex[vertex_to_point] = 1:length(vertex_to_point)
    point_to_vertex[slave_point_to_point] = point_to_vertex[slave_point_to_master_point]

    cell_to_vertices = Table(LocalToGlobalArray(nodes,point_to_vertex))

    vertex_to_node = vertex_to_point
    node_to_vertex = point_to_vertex
  else
    cell_to_nodes = get_cell_nodes(grid)
    cell_to_cell_type = get_cell_type(grid)
    reffes = get_reffes(grid)
    cell_type_to_lvertex_to_lnode = map(get_vertex_node, reffes)
    cell_to_vertices, vertex_to_node, node_to_vertex = _generate_cell_to_vertices(
      cell_to_nodes,
      cell_to_cell_type,
      cell_type_to_lvertex_to_lnode,
      num_nodes(grid))
  end
  (cell_to_vertices, vertex_to_node, node_to_vertex)
end



function _generate_cell_to_vertices(
  cell_to_nodes::Table,
  cell_to_cell_type::AbstractVector{<:Integer},
  cell_type_to_lvertex_to_lnode::Vector{Vector{Int}},
  nnodes::Int=maximum(cell_to_nodes.data))

  data, ptrs, vertex_to_node, node_to_vertex = _generate_cell_to_vertices(
    cell_to_nodes.data,
    cell_to_nodes.ptrs,
    cell_to_cell_type,
    cell_type_to_lvertex_to_lnode,
    nnodes)

  (Table(data,ptrs), vertex_to_node)
end

function _generate_cell_to_vertices(
  cell_to_nodes_data,
  cell_to_nodes_ptrs,
  cell_to_cell_type,
  cell_type_to_lvertex_to_lnode,
  nnodes)

  cell_to_vertices_ptrs = similar(cell_to_nodes_ptrs)

  cell_type_to_nlvertices = map(length,cell_type_to_lvertex_to_lnode)

  _generate_cell_to_vertices_count!(
    cell_to_vertices_ptrs,
    cell_to_cell_type,
    cell_type_to_nlvertices)

  T = eltype(cell_to_nodes_data)

  node_to_vertex = fill(T(UNSET),nnodes)

  length_to_ptrs!(cell_to_vertices_ptrs)

  ndata = cell_to_vertices_ptrs[end]-1
  cell_to_vertices_data = zeros(T,ndata)

  _generate_cell_to_vertices_fill!(
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    cell_to_nodes_data,
    cell_to_nodes_ptrs,
    node_to_vertex,
    cell_to_cell_type,
    cell_type_to_lvertex_to_lnode)

  vertex_to_node = find_inverse_index_map(node_to_vertex)

  (cell_to_vertices_data, cell_to_vertices_ptrs, vertex_to_node, node_to_vertex)

end

function  _generate_cell_to_vertices_count!(
    cell_to_vertices_ptrs,
    cell_to_cell_type,
    cell_type_to_nlvertices)

  cells = 1:length(cell_to_cell_type)
  for cell in cells
    cell_type = cell_to_cell_type[cell]
    nlvertices = cell_type_to_nlvertices[cell_type]
    cell_to_vertices_ptrs[1+cell] = nlvertices
  end
end

function  _generate_cell_to_vertices_fill!(
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    cell_to_nodes_data,
    cell_to_nodes_ptrs,
    node_to_vertex,
    cell_to_cell_type,
    cell_type_to_lvertex_to_lnode)

  cells = 1:length(cell_to_cell_type)

  vertex = 1

  for cell in cells
    cell_type = cell_to_cell_type[cell]
    a = cell_to_nodes_ptrs[cell]-1
    b = cell_to_vertices_ptrs[cell]-1

    lvertex_to_lnode = cell_type_to_lvertex_to_lnode[cell_type]
    for (lvertex, lnode) in enumerate(lvertex_to_lnode)
      node = cell_to_nodes_data[a+lnode]
      if node_to_vertex[node] == UNSET
        node_to_vertex[node] = vertex
        vertex += 1
      end
      cell_to_vertices_data[b+lvertex] = node_to_vertex[node]
    end

  end
end


function _generate_face_to_own_dofs(
  n_faces,
  cell_to_ctype,
  d_to_cell_to_dfaces::Vector{Table{T,P}},
  d_to_dface_to_cells::Vector{Table{T,P}},
  d_to_offset,
  d_to_ctype_to_ldface_to_own_ldofs) where {T,P}

  face_to_own_dofs_ptrs = zeros(P,n_faces+1)

  D = length(d_to_offset)-1

  icell = 1
  d_to_dface_to_cell = [ get_local_item(d_to_dface_to_cells[d+1],icell)  for d in 0:D ]

  d_to_dface_to_ldface = [
    find_local_index(d_to_dface_to_cell[d+1],d_to_cell_to_dfaces[d+1]) for d in 0:D ]

  for d in 0:D
    cell_to_dfaces = d_to_cell_to_dfaces[d+1]
    dface_to_cells = d_to_dface_to_cells[d+1]
    offset = d_to_offset[d+1]
    ctype_to_ldface_to_own_ldofs = d_to_ctype_to_ldface_to_own_ldofs[d+1]
    ctype_to_ldface_to_num_own_ldofs = map( (x) -> length.(x) ,ctype_to_ldface_to_own_ldofs)
    dface_to_cell_owner = d_to_dface_to_cell[d+1]
    dface_to_ldface = d_to_dface_to_ldface[d+1]

    if any( ctype_to_ldface_to_num_own_ldofs .!= 0)

      _generate_face_to_own_dofs_count_d!(
        face_to_own_dofs_ptrs,
        offset,
        cell_to_ctype,
        dface_to_cell_owner,
        dface_to_ldface,
        ctype_to_ldface_to_num_own_ldofs)
    end
  end

  length_to_ptrs!(face_to_own_dofs_ptrs)

  n_dofs = face_to_own_dofs_ptrs[end]-1
  face_to_own_dofs_data = collect(T(1):T(n_dofs))

  face_to_own_dofs = Table(face_to_own_dofs_data,face_to_own_dofs_ptrs)
  (face_to_own_dofs, n_dofs, d_to_dface_to_cell, d_to_dface_to_ldface)
end

function writePVD(filename,timeSteps)
  rm(filename,force=true,recursive=true)
  mkdir(filename)
  pvdcontent  = """<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>\n"""
  for t in timeSteps
    pvdcontent *= """    <DataSet timestep=""" * '"'
    pvdcontent *= string(t) * '"' * """ group="" part="0" file=""" * '"'
    pvdcontent *= filename*"""/time_"""*string(t)*""".vtu"/>\n"""
  end
  pvdcontent  *= "  </Collection>\n</VTKFile>"
  f = open(filename * ".pvd", "w")
  write(f,pvdcontent)
  close(f)
end

function analytical_solution(a::Float64,       # semi-length of side walls
                             b::Float64,       # semi-length of Hartmann walls
                             t_w::Float64,     # wall thickness
                             σ_w::Float64,     # wall conductivity
                             σ::Float64,       # fluid conductivity
                             μ::Float64,       # fluid viscosity
                             grad_pz::Float64, # presure gradient
                             Ha::Float64,      # Hartmann number
                             n::Int,           # number of sumands included in Fourier series
                             x)                # evaluation point
  l = b/a
  ξ = x[1]/a
  η = x[2]/a

  d_B = t_w*σ_w/(a*σ)

  V = 0.0
  dH_dx = 0.0; dH_dy = 0.0
  for k in 0:n
    α_k = (k + 0.5)*π/l
    r1_k = 0.5*( Ha + (Ha^2 + 4*α_k^2)^0.5)
    r2_k = 0.5*(-Ha + (Ha^2 + 4*α_k^2)^0.5)
    N = (Ha^2 + 4*α_k^2)^0.5

    V2 = ((d_B * r2_k + (1-exp(-2*r2_k))/(1+exp(-2*r2_k))) * 0.5 * (exp(-r1_k*(1-η))+exp(-r1_k*(1+η))))/
         (0.5*(1+exp(-2*r1_k))*d_B*N + (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r2_k)))

    V3 = ((d_B * r1_k + (1-exp(-2*r1_k))/(1+exp(-2*r1_k))) * 0.5 * (exp(-r2_k*(1-η))+exp(-r2_k*(1+η))))/
         (0.5*(1+exp(-2*r2_k))*d_B*N + (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r1_k)))

    V += 2*(-1)^k*cos(α_k * ξ)/(l*α_k^3)*(1-V2-V3)

    H2 = ((d_B * r2_k + (1-exp(-2*r2_k))/(1+exp(-2*r2_k))) * 0.5 * (exp(-r1_k*(1-η))-exp(-r1_k*(1+η))))/
         (0.5*(1+exp(-2*r1_k))*d_B*N + (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r2_k)))

    H3 = ((d_B * r1_k + (1-exp(-2*r1_k))/(1+exp(-2*r1_k))) * 0.5 * (exp(-r2_k*(1-η))-exp(-r2_k*(1+η))))/
         (0.5*(1+exp(-2*r2_k))*d_B*N + (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r1_k)))

    H2_dy = ((d_B * r2_k + (1-exp(-2*r2_k))/(1+exp(-2*r2_k))) * 0.5 * (exp(-r1_k*(1-η))*(r1_k/a)-exp(-r1_k*(1+η))*(-r1_k/a)))/
         (0.5*(1+exp(-2*r1_k))*d_B*N + (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r2_k)))

    H3_dy = ((d_B * r1_k + (1-exp(-2*r1_k))/(1+exp(-2*r1_k))) * 0.5 * (exp(-r2_k*(1-η))*(r2_k/a)-exp(-r2_k*(1+η))*(-r2_k/a)))/
         (0.5*(1+exp(-2*r2_k))*d_B*N + (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r1_k)))


    dH_dx += -2*(-1)^k*sin(α_k * ξ)/(a*l*α_k^3)*(H2-H3)
    dH_dy += 2*(-1)^k*cos(α_k * ξ)/(l*α_k^3)*(H2_dy-H3_dy)

  end
  u_z = V/μ * (-grad_pz) * a^2
  j_x = dH_dy / μ^0.5 * (-grad_pz) * a^2*σ^0.5
  j_y = -dH_dx / μ^0.5 * (-grad_pz) * a^2*σ^0.5

  u = VectorValue(0.0,0.0,u_z)
  j = VectorValue(j_x,j_y,0.0)
  return u,j
end

function main(partition=(4,4,3),Δt=1e-4)


t0 = 0.0
tf = Δt*4
domain = (-0.5,0.5,-0.5,0.5,0.0,0.01)
order = 2

map(x) = VectorValue(sign(x[1])*(abs(x[1])*0.5)^0.5,   sign(x[2])*(abs(x[2])*0.5)^0.5,  x[3])
model = CartesianDiscreteModel(domain,partition,[3],map)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",append!(collect(1:20),[23,24,25,26]))

Vu = FESpace(
    reffe=:Lagrangian, order=order, valuetype=VectorValue{3,Float64},
    conformity=:H1, model=model, dirichlet_tags="dirichlet")

Vp = FESpace(
    reffe=:PLagrangian, order=order-1, valuetype=Float64,
    conformity=:L2, model=model, constraint=:zeromean)

Vj = FESpace(
    reffe=:RaviartThomas, order=order-1, valuetype=VectorValue{3,Float64},
    conformity=:Hdiv, model=model, dirichlet_tags="dirichlet")

Vφ = FESpace(
    reffe=:QLagrangian, order=order-1, valuetype=Float64,
    conformity=:L2, model=model, constraint=:zeromean)

u0 = VectorValue(0.0,0.0,0.010)
gu = VectorValue(0.0,0.0,0.0)
gj = VectorValue(0.0,0.0,0.0)

U = TrialFESpace(Vu,gu)
P = TrialFESpace(Vp)
j = TrialFESpace(Vj,gj)
φ = TrialFESpace(Vφ)

Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
X = MultiFieldFESpace([U, P, j, φ])

trian = Triangulation(model)
degree = 2*(order+1)
quad = CellQuadrature(trian,degree)

# # xh = FEFunction(X,rand(num_free_dofs(X)))
# uh = FEFunction(U,rand(num_free_dofs(U)))
# ph = FEFunction(P,rand(num_free_dofs(P)))
# jh = FEFunction(j,rand(num_free_dofs(j)))
# φh = FEFunction(φ,rand(num_free_dofs(φ)))
#
# # uh, ph, jh, φh = xh
#
#
# # writevtk(trian,"results",cellfields=["u"=>uh, "p"=>ph, "j"=>jh, "phi"=>φh])
# writevtk(trian,"results",cellfields=["u"=>uh, "p"=>ph, "j"=>jh, "phi"=>φh])
#
# @assert false

neumanntags = [21,22]
btrian = BoundaryTriangulation(model,neumanntags)
degree = 2*order
bquad = CellQuadrature(btrian,degree)
nb = get_normal_vector(btrian)

# B? uk?
x = get_physical_coordinate(trian)
B(x) = VectorValue(0.0,10.0,0.0)
ρ = 1.0
ν = 1.0
σ = 1.0
Re = 10.0 # U = 10.0, L = 1.0/ ν = 1.0
Ha = 10.0
K = Ha / (1-0.825*Ha^(1/2)-Ha^(-1))
f(x) = VectorValue(0.0,0.0,K / Re)
analytical_u(x) = analytical_solution(0.5,  # semi-length of side walls
                                      0.5,  # semi-length of Hartmann walls
                                      0.0,  # wall conductivity
                                      1.0,  # wall thickness
                                      1.0,  # fluid conductivity
                                      1.0,  # fluid viscosity
                                      K/Re, # presure gradient
                                      Ha,   # Hartmann number
                                      10,   # number of sumands included in Fourier series
                                      x)[1]
un = interpolate(Vu, analytical_u)

@law vprod(a,b) = VectorValue(a[2]b[3]-a[3]b[2], a[1]b[3]-a[3]b[1], a[1]b[2]-a[2]b[1])

function a(X,Y)
  u  , p  , j  , φ   = X
  v_u, v_p, v_j, v_φ = Y

  # uk*(∇(u)*v_u) + ν*inner(∇(u),∇(v_u)) - p*(∇*v_u) + (∇*u)*v_p
  (1/Δt)*u*v_u + un*(∇(u)*v_u) + ν*inner(∇(u),∇(v_u)) - p*(∇*v_u) + (∇*u)*v_p - 1/ρ * vprod(j,B(x))*v_u +
  j*v_j + σ*(∇(φ)*v_j) - σ*vprod(u,B(x))*v_j - ∇(v_φ)*j
end

function l(y)
 v_u, v_p, v_j, v_φ = y
 (1/Δt)*un*v_u + v_u*f + v_p*0.0 + v_j*VectorValue(0.0,0.0,0.0) + v_φ*0.0
end

u_nbc = VectorValue(0.0,0.0,0.0)
p_nbc = 0.0
j_nbc = VectorValue(0.0,0.0,0.0)
φ_nbc = 0.0

function l_Γ(y)
 v_u, v_p, v_j, v_φ = y
 u_nbc * v_u + p_nbc * v_p + j_nbc * v_j + φ_nbc * v_φ
end

t_Ω = AffineFETerm(a,l,trian,quad)
t_Γ = FESource(l_Γ,btrian,bquad)
op  = AffineFEOperator(X,Y,t_Ω)

ls = LUSolver()
solver = LinearFESolver(ls)



timeSteps = collect(t0:Δt:tf)

writePVD("results",timeSteps[1:end])
writevtk(trian, "results/time_"*string(t0)*".vtu",cellfields=["u"=>un])

for t in timeSteps[2:end]

  xh = solve(solver,op)
  un, pn, jn, φn = xh

  writevtk(trian,"results/time_"*string(t)*".vtu",cellfields=["u"=>un, "p"=>pn, "j"=>jn, "phi"=>φn])

  @show t

  # Update operator
  op = AffineFEOperator(X,Y,t_Ω)
end


end
end # module
