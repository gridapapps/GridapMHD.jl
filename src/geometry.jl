
const DiscreteModelTypes = Union{Gridap.DiscreteModel,GridapDistributed.DistributedDiscreteModel}
const TriangulationTypes = Union{Gridap.Triangulation,GridapDistributed.DistributedTriangulation}

function setup_geometry!(params)
  params[:geometry] = Dict{Symbol,Any}(
    :interiors  => Dict{UInt64,Any}(),
    :boundaries => Dict{UInt64,Any}(),
    :skeletons  => Dict{UInt64,Any}(),
    :normals    => Dict{UInt64,Any}(),
    :measures   => Dict{UInt64,Any}(),
    :other      => Dict{UInt64,Any}(),
  )

  Ω = interior(params,nothing)
  params[:Ω] = Ω

  fluid = params[:fluid]
  Ωf = interior(params,fluid[:domain])
  params[:Ωf] = Ωf

  solid = params[:solid]
  Ωs = !isnothing(solid) ? interior(params,solid[:domain]) : nothing
  params[:Ωs] = Ωs

  if uses_multigrid(params[:solver])
    params[:multigrid][:Ω]  = params[:multigrid][:mh]
    params[:multigrid][:Ωf] = params[:multigrid][:mh]
    params[:multigrid][:Ωs] = params[:multigrid][:mh]
  end
end

domain_hash(domain::Union{Nothing,DiscreteModelTypes,TriangulationTypes}) = objectid(domain)
domain_hash(domain::Union{String,Vector{String}}) = hash(join(domain))

domain_hash(model,domain::Union{String,Vector{String}}) = hash(join(domain),objectid(model))
domain_hash(model,domain) = domain_hash(domain)

interior(params::Dict,domain) = interior(params,params[:model],domain)

function interior(params::Dict,model,domain)
  interiors = params[:geometry][:interiors]
  key = domain_hash(model,domain)
  if haskey(interiors,key)
    trian = interiors[key]
  else
    trian = _interior(model,domain)
    interiors[domain_hash(trian)] = trian
    interiors[key] = trian
  end
  return trian
end

boundary(params::Dict,domain) = boundary(params,params[:model],domain)

function boundary(params::Dict,model,domain)
  boundaries = params[:geometry][:boundaries]
  key = domain_hash(model,domain)
  if haskey(boundaries,key)
    trian = boundaries[key]
  else
    trian = _boundary(model,domain)
    boundaries[domain_hash(trian)] = trian
    boundaries[key] = trian
  end
  return trian
end

skeleton(params::Dict,domain) = skeleton(params,params[:model],domain)

function skeleton(params::Dict,model,domain)
  skeletons = params[:geometry][:skeletons]
  key = domain_hash(model,domain)
  if haskey(skeletons,key)
    trian = skeletons[key]
  else
    trian = _skeleton(model,domain)
    skeletons[domain_hash(trian)] = trian
    skeletons[key] = trian
  end
  return trian
end

function normal_vector(params::Dict,trian::TriangulationTypes)
  normals = params[:geometry][:normals]
  key = domain_hash(trian)
  if haskey(normals,key)
    n = normals[key]
  else
    n = get_normal_vector(trian)
    normals[key] = n
  end
  return n
end

function measure(params::Dict,trian::TriangulationTypes)
  measures = params[:geometry][:measures]
  key = domain_hash(trian)
  if haskey(measures,key)
    m = measures[key]
  else
    d = num_cell_dims(trian)
    q = params[:fespaces][:quadratures][d]
    m = Measure(trian,q)
    measures[key] = m
  end
  return m
end

_interior(model,domain::DiscreteModelTypes) = Triangulation(domain)
_interior(model,domain::TriangulationTypes) = domain
_interior(model,domain::Nothing) = Triangulation(model)
_interior(model,domain) = Triangulation(model,tags=domain)

_boundary(model,domain::TriangulationTypes) = domain
_boundary(model,domain::Nothing) = Boundary(model)
_boundary(model,domain) = Boundary(model,tags=domain)

_skeleton(model,domain::TriangulationTypes) = Skeleton(domain)
_skeleton(model,domain::Nothing) = Skeleton(model)
_skeleton(model,domain) = Skeleton(_interior(model,domain))

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
