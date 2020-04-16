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
    #
    # cell_index = Array{Int,1}(undef,D)
    # for (icell,cell_nodes) in enumerate(nodes)
    #   for (inode,node) in enumerate(cell_nodes)
    #     index = findfirst(x->x>node,list_deleted)
    #     if isnothing(index)
    #       index = 0
    #     else
    #       index -= 1
    #     end
    #     cell_to_vertices[icell][inode] = node - index
    #   end
    #   cell_ijk = index_to_ijk(icell,num_cells)
    #   source_ijk = Array{Int,1}(undef,D)
    #   for dir in periodic
    #     if cell_ijk[dir] == partition[dir]
    #       source_ijk .= cell_ijk
    #       source_ijk[dir] = 1
    #       source_cell = ijk_to_index(source_ijk, num_cells)
    #       vertex_map = Array{Int,2}(undef,(2^(D-1),2))
    #       l = 1
    #       for k in collect(1:2^(D-dir))
    #         for kj in collect(1:2)
    #           for ki in collect(1:2^(dir-1))
    #             vertex_map[ki+(k-1)*2^(dir-1),kj] = l
    #             l=l+1
    #           end
    #         end
    #       end
    #       cell_to_vertices[icell][vertex_map[:,2]] .= cell_to_vertices[source_cell][vertex_map[:,1]]
    #     end
    #   end
    # end
    # cell_to_vertices = Table(cell_to_vertices)
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
      # n = length(dface_to_ldface)
      # for i in 1:length(dface_to_ldface)
      #   if !(isassigned(dface_to_ldface,i))
      #     n = i - 1
      #     break
      #   end
      # end

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



function main()

domain = (0,1,0,1,0,0.01)
partition = (3,3,3)
order = 2
model = CartesianDiscreteModel(domain,partition)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",append!(collect(1:20),[23,24,25,26]))

Vu = FESpace(
    reffe=:Lagrangian, order=order, valuetype=VectorValue{3,Float64},
    conformity=:H1, model=model, dirichlet_tags="dirichlet")

Vp = FESpace(
    reffe=:PLagrangian, order=order-1, valuetype=Float64,
    conformity=:L2, model=model)

Vj = FESpace(
    reffe=:RaviartThomas, order=order-1, valuetype=VectorValue{3,Float64},
    conformity=:Hdiv, model=model, dirichlet_tags="dirichlet")

Vφ = FESpace(
    reffe=:QLagrangian, order=order-1, valuetype=Float64,
    conformity=:L2, model=model)

u0 = VectorValue(0.0,0.0,0.0)
gu = VectorValue(0.0,0.0,0.0)
gj = VectorValue(0.0,0.0,0.0)

U = TrialFESpace(Vu,gu)
P = TrialFESpace(Vp)
j = TrialFESpace(Vj,gj)
φ = TrialFESpace(Vφ)

Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
X = MultiFieldFESpace([U, P, j, φ])

trian = Triangulation(model)
degree = 2
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
uk = interpolate(Vu, u0)
B(x) = VectorValue(0.0,10.0,0.0)
ρ = 1.0
ν = 1.0
σ = 1.0
Re = 10.0 # U = 10.0, L = 1.0/ ν = 1.0
Ha = 10.0
K = Ha / (1-0.825Ha^(1/2)-Ha^(-1))
f(x) = VectorValue(0.0,0.0,K / Re)

@law vprod(a,b) = VectorValue(a[2]b[3]-a[3]b[2], a[1]b[3]-a[3]b[1], a[1]b[2]-a[2]b[1])

function a(X,Y)
  u  , p  , j  , φ   = X
  v_u, v_p, v_j, v_φ = Y

  # uk*(∇(u)*v_u) + ν*inner(∇(u),∇(v_u)) - p*(∇*v_u) + (∇*u)*v_p
  uk*(∇(u)*v_u) + ν*inner(∇(u),∇(v_u)) - p*(∇*v_u) + (∇*u)*v_p - 1/ρ * vprod(j,B(x))*v_u +
  j*v_j + σ*(∇(φ)*v_j) - σ*vprod(u,B(x))*v_j - ∇(v_φ)*j
end

function l(y)
 v_u, v_p, v_j, v_φ = y
 v_u*f
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

Δx = 1.0

while Δx > 1e-3

 xh = solve(solver,op)
 uk, pk, jk, φk = xh

 Δx = (xh - xk)*(xh - xk)/(xh*xh)
 xk = xh
 @show Δx
 # Update operator
 op = AffineFEOperator(X,Y,t_Ω,t_Γ)
 writevtk(trian,"results",cellfields=["u"=>uk, "p"=>pk, "j"=>jk, "phi"=>φk])
end



end
end # module
