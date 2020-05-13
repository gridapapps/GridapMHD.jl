module PeriodicBC

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

export CartesianDiscreteModel

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
  sizes = Tuple([(limits[d][2]-limits[d][1])/partition[d] for d in 1:D])
  origin = Point([ limits[d][1] for d in 1:D]...)
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

end # module
