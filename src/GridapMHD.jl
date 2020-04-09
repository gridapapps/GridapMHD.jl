module GridapMHD

using Gridap
import .Gridap: ∇, divergence

using Gridap.Arrays: Table
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
  nfaces = [num_faces(topo,d) for d in 0:num_cell_dims(topo)]
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

    # compute periodic nodes
    D = length(partition)
    num_cells = Array{Int,1}(undef,D)
    num_cells = [partition[i] for i in 1:D]
    list_deleted = Array{Int,1}(undef,0)
    for dir in periodic
      dir_x = collect(1:D)
      deleteat!(dir_x, dir)

      facet_index = Array{Int,1}(undef,D)
      facet_index .= 1
      facet_index[dir_x[1]] = 0
      facet_index[dir] = num_cells[dir]
      for i in collect(1:prod(num_cells[dir_x]))
        facet_index[dir_x[1]] += 1
        for j in collect(1:D-1)
          if facet_index[dir_x[j]] > num_cells[dir_x[j]]
            facet_index[dir_x[j+1]] += 1
            facet_index[dir_x[j]] = 0
          end
        end
        cell = ijk_to_index(facet_index,num_cells)

        vertex_map = Array{Int,2}(undef,(2^(D-1),2))
        l = 1
        for k in collect(1:2^(D-dir))
          for kj in collect(1:2)
            for ki in collect(1:2^(dir-1))
              vertex_map[ki+(k-1)*2^(dir-1),kj] = l
              l=l+1
            end
          end
        end

        vertices = vertex_map[:,2]
        append!(list_deleted,nodes[cell][vertices])
      end
    end
    sort!(list_deleted)
    unique!(list_deleted)

    vertex_to_node = collect(1:num_nodes(grid))
    node_to_vertex = vertex_to_node

    cell_index = Array{Int,1}(undef,D)
    for (icell,cell_nodes) in enumerate(nodes)
      for (inode,node) in enumerate(cell_nodes)
        index = findfirst(x->x>node,list_deleted)
        if isnothing(index)
          index = 0
        else
          index -= 1
        end
        cell_to_vertices[icell][inode] = node - index
      end
      cell_ijk = index_to_ijk(icell,num_cells)
      source_ijk = Array{Int,1}(undef,D)
      for dir in periodic
        if cell_ijk[dir] == partition[dir]
          source_ijk .= cell_ijk
          source_ijk[dir] = 1
          source_cell = ijk_to_index(source_ijk, num_cells)
          vertex_map = Array{Int,2}(undef,(2^(D-1),2))
          l = 1
          for k in collect(1:2^(D-dir))
            for kj in collect(1:2)
              for ki in collect(1:2^(dir-1))
                vertex_map[ki+(k-1)*2^(dir-1),kj] = l
                l=l+1
              end
            end
          end
          cell_to_vertices[icell][vertex_map[:,2]] .= cell_to_vertices[source_cell][vertex_map[:,1]]
        end
      end
    end
    cell_to_vertices = Table(cell_to_vertices)
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


#
# function main()
#
# u(x) = VectorValue(2*x[1],x[1]+x[2])
# divergence(::typeof(u)) = (x) -> 3
#
# p(x) = x[1]-x[2]
#
# ∇p(x) = VectorValue(1,-1)
#
# ∇(::typeof(p)) = ∇p
#
# f(x) = u(x) + ∇p(x)
#
# domain = (0,1,0,1)
# partition = (4,4)
# order = 2
# model = CartesianDiscreteModel(domain,partition)
#
# Vu = FESpace(
#      reffe=:Lagrangian, order=order, valuetype=Float64,
#      conformity=:H1, model=model)
#
# Vp = FESpace(
#      reffe=:QLagrangian, order=order-1, valuetype=Float64,
#      conformity=:L2, model=model)
#
# Vj = FESpace(
#      reffe=:RaviartThomas, order=order, valuetype=VectorValue{2,Float64},
#      conformity=:Hdiv, model=model, dirichlet_tags=[5,6])
#
# Vφ = FESpace(
#      reffe=:QLagrangian, order=order-1, valuetype=Float64,
#      conformity=:L2, model=model)
#
# U = TrialFESpace(Vu,u)
# P = TrialFESpace(Vp)
# j = TrialFESpace(Vj,u)
# φ = TrialFESpace(Vφ)
#
# Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
# X = MultiFieldFESpace([U, P, j, φ])
#
# trian = Triangulation(model)
# degree = 2
# quad = CellQuadrature(trian,degree)
#
# neumanntags = [7,8]
# btrian = BoundaryTriangulation(model,neumanntags)
# degree = 2*order
# bquad = CellQuadrature(btrian,degree)
# nb = get_normal_vector(btrian)
#
# # B? uk?
# x = get_physical_coordinate(trian)
# uk = interpolate(Vu, u0)
# B(x) = VectorValue(1,1)
#
# function a(x,y)
#   u  , p  , j  , φ   = x
#   v_u, v_p, v_j, v_φ = y
#   uk*(∇(u)*v_u) + ν*(∇(u)*∇(v_u)) - p*(∇*v_u) + v_p*(∇*u) - 1/ρ * (j×B(x))*v_u +
#   j*v_j + σ*(∇(φ)*v_j) - σ*(u×B(x))*v_j - ∇(v_φ)*j
# end
#
# function l(y)
#   v_u, v_p, v_j, v_φ = y
#   v_u*f
# end
#
# u_nbc =
# p_nbc =
# j_nbc =
# φ_nbc =
#
# function l_Γ(y)
#   v_u, v_p, v_j, v_φ = y
#   u_nbc * v_u + p_nbc * v_p + j_nbc * v_j + φ_nbc * v_φ
# end
#
# t_Ω = AffineFETerm(a,l,trian,quad)
# t_Γ = FESource(l_Γ,btrian,bquad)
# op  = AffineFEOperator(X,Y,t_Ω,t_Γ)
#
# ls = LUSolver()
# solver = LinearFESolver(ls)
#
# while Δx > 1e-3
#
#   xh = solve(solver,op)
#   uk, pk, jk, φk = xh
#
#   Δx = (xh - xk)*(xh - xk)/(xh*xh)
#   xk = xh
#
#   # Update operator
#   op = AffineFEOperator(X,Y,t_Ω,t_Γ)
#   writevtk(trian,"results",cellfields=["u"=>uk, "p"=>pk, "j"=>jk, "phi"=>φk])
# end
#
#
#
# end

end # module
