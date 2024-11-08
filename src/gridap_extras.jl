
# Get polytopes

function Geometry.get_polytopes(model::GridapDistributed.DistributedDiscreteModel)
  polys = map(get_polytopes,local_views(model))
  return getany(polys)
end

function Geometry.get_polytopes(trian::Triangulation)
  reffes = get_reffes(trian)
  unique(map(get_polytope,reffes))
end

function Geometry.get_polytopes(trian::GridapDistributed.DistributedTriangulation)
  polys = map(get_polytopes,local_views(trian))
  return getany(polys)
end

# Face/cell sizes

get_cell_size(t::TriangulationTypes) = CellField(_get_cell_size(t),t)

get_cell_size(m::DiscreteModelTypes) = get_cell_size(Triangulation(m))
_get_cell_size(m::DiscreteModelTypes) = _get_cell_size(Triangulation(m))

function _get_cell_size(t::Triangulation) :: Vector{Float64}
  if iszero(num_cells(t))
    return Float64[]
  end
  meas = get_cell_measure(t)
  d = num_dims(t)
  return collect(Float64, meas .^ (1/d))
end

function _get_cell_size(t::GridapDistributed.DistributedTriangulation)
  map(_get_cell_size,local_views(t))
end

get_mesh_size(m::DiscreteModel) = minimum(_get_cell_size(m))

function get_mesh_size(m::GridapDistributed.DistributedDiscreteModel)
  h = map(get_mesh_size,local_views(m))
  return reduce(min,h;init=one(eltype(h)))
end

