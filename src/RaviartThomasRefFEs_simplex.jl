module RaviartThomasRefFE_simplex

using Gridap
using Gridap.Helpers
using Gridap.Fields
using Gridap.Arrays
using Gridap.Integration
using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.Polynomials: _prepare_perms
using Gridap.Polynomials: _evaluate_1d!
using Gridap.Polynomials: _gradient_1d!
using Gridap.ReferenceFEs
using Gridap.ReferenceFEs: _RT_face_moments
# using Gridap.ReferenceFEs: _face_own_dofs_from_moments

import Gridap.Polynomials: get_order
import Gridap.Polynomials: get_orders
import Gridap.Polynomials: get_value_type
import Gridap.Polynomials: num_terms
import Gridap.Polynomials: field_cache
import Gridap.Polynomials: gradient_cache
import Gridap.Polynomials: evaluate_field!
import Gridap.Polynomials: evaluate_gradient!

# import Gridap.ReferenceFEs: _RT_face_values
import Gridap.ReferenceFEs: _RT_cell_values
# import Gridap.ReferenceFEs: _face_own_dofs_from_moments

export PGradMonomialBasis
export PCurlGradMonomialBasis

export get_order
export get_orders
export get_value_type
export num_terms

export field_cache
export gradient_cache
export evaluate_field!
export evaluate_gradient!

export RaviartThomasRefFE_t


# struct PGradMonomialBasis{D,T} <: Field
#   order::Int
#   terms::Array{CartesianIndex{D},1}
#   perms::Matrix{Int}
#   function PGradMonomialBasis(::Type{T},order::Int,terms::Array{CartesianIndex{D},1},perms::Matrix{Int}) where {D,T}
#     new{D,T}(order,terms,perms)
#   end
# end
#
# """
#     num_terms(f::PGradMonomialBasis{D,T}) where {D,T}
# """
# function num_terms(f::PGradMonomialBasis{D,T}) where {D,T}
#   # (_p_dim(f.order,D) + _p_dim(f.order+1,D-1))*D - _p_dim(f.order+2)
#   dim = f.order+1
#   for d in 2:D
#     dim *= (f.order+1+d)
#   end
#   Int(dim/factorial(dim-1))
# end


struct PCurlGradMonomialBasis{D,T} <: Field
  order::Int
  pterms::Array{CartesianIndex{D},1}
  sterms::Array{CartesianIndex{D},1}
  perms::Matrix{Int}
  function PCurlGradMonomialBasis(::Type{T},order::Int,
      pterms::Array{CartesianIndex{D},1},sterms::Array{CartesianIndex{D},1},
      perms::Matrix{Int}) where {D,T}
    new{D,T}(order,pterms,sterms,perms)
  end
end

function PCurlGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}
  @assert T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
  P_k = MonomialBasis{D}(T, order, _p_filter)
  S_k = MonomialBasis{D}(T, order, _s_filter)
  pterms = P_k.terms
  sterms = S_k.terms
  perms = _prepare_perms(D)
  PCurlGradMonomialBasis(T,order,pterms,sterms,perms)
end

function field_cache(f::PCurlGradMonomialBasis{D,T},x) where {D,T}
  @assert D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = _ndofs_pgrad(f)
  n = 1 + f.order+1
  V = VectorValue{D,T}
  r = CachedArray(zeros(V,(np,ndof)))
  v = CachedArray(zeros(V,(ndof,)))
  c = CachedArray(zeros(T,(D,n)))
  (r, v, c)
end

function evaluate_field!(cache,f::PCurlGradMonomialBasis{D,T},x) where {D,T}
  r, v, c = cache
  np = length(x)
  ndof = _ndofs_pgrad(f)
  n = 1 + f.order+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
      _evaluate_nd_pcurlgrad!(v,xi,f.order+1,f.pterms,f.sterms,f.perms,c)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r
end

function gradient_cache(f::PCurlGradMonomialBasis{D,T},x)  where {D,T}
  @assert D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = _ndofs_pgrad(f)
  n = 1 + f.order+1
  xi = testitem(x)
  V = VectorValue{D,T}
  G = gradient_type(V,xi)
  r = CachedArray(zeros(G,(np,ndof)))
  v = CachedArray(zeros(G,(ndof,)))
  c = CachedArray(zeros(T,(D,n)))
  g = CachedArray(zeros(T,(D,n)))
  (r, v, c, g)
end

function evaluate_gradient!(cache,f::PCurlGradMonomialBasis{D,T},x) where {D,T}
  r, v, c, g = cache
  np = length(x)
  ndof = _ndofs_pgrad(f)
  n = 1 + f.order+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  setsize!(g,(D,n))
  V = VectorValue{D,T}
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd_pcurlgrad!(v,xi,f.order+1,f.pterms,f.sterms,f.perms,c,g,V)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r
end

get_value_type(::PCurlGradMonomialBasis{D,T}) where {D,T} = T

"""
    num_terms(f::QCurlGradMonomialBasis{D,T}) where {D,T}
"""
function num_terms(f::PCurlGradMonomialBasis{D,T}) where {D,T}
  Int(_p_dim(f.order,D)*D + _p_dim(f.order,D-1))
end

get_order(f::PCurlGradMonomialBasis{D,T}) where {D,T} = f.order

# Helpers

_p_filter(e,order) = (sum(e) <= order)
_s_filter(e,order) = (sum(e) == order)

function _p_dim(order,D)
  dim = 1
  for d in 1:D
    dim *= order+d
  end
  dim/factorial(D)
end

_ndofs_pgrad(f::PCurlGradMonomialBasis{D}) where D = num_terms(f)



function _evaluate_nd_pcurlgrad!(
  v::AbstractVector{V},
  x,
  order,
  pterms::Array{CartesianIndex{D},1},
  sterms::Array{CartesianIndex{D},1},
  perms::Matrix{Int},
  c::AbstractMatrix{T}) where {V,T,D}

  dim = D
  for d in 1:dim
    _evaluate_1d!(c,x,order,d)
  end

  o = one(T)
  k = 1
  m = zero(mutable(V))
  js = eachindex(m)
  z = zero(T)

  for ci in pterms
    for j in js

      @inbounds for i in js
        m[i] = z
      end

      s = o
      @inbounds for d in 1:dim
        s *= c[d,ci[perms[d,j]]]
      end

      m[j] = s
      v[k] = m
      k += 1
    end
  end

  for ci in sterms
    @inbounds for i in js
      m[i] = z
    end
    for j in js

      s = c[j,2]
      @inbounds for d in 1:dim
        s *= c[d,ci[d]]
      end

      m[j] = s

    end
    v[k] = m
    k += 1
  end
end

function _gradient_nd_pcurlgrad!(
  v::AbstractVector{G},
  x,
  order,
  pterms::Array{CartesianIndex{D},1},
  sterms::Array{CartesianIndex{D},1},
  perms::Matrix{Int},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  ::Type{V}) where {G,T,D,V}

  dim = D
  for d in 1:dim
    _evaluate_1d!(c,x,order,d)
    _gradient_1d!(g,x,order,d)
  end

  z = zero(mutable(V))
  m = zero(mutable(G))
  js = eachindex(z)
  mjs = eachindex(m)
  o = one(T)
  zi = zero(T)
  k = 1

  for ci in pterms
    for j in js

      s = z
      for i in js
        s[i] = o
      end

      for q in 1:dim
        for d in 1:dim
          if d != q
            @inbounds s[q] *= c[d,ci[perms[d,j]]]
          else
            @inbounds s[q] *= g[d,ci[perms[d,j]]]
          end
        end
      end

      @inbounds for i in mjs
        m[i] = zi
      end

      for i in js
        @inbounds m[i,j] = s[i]
      end
      @inbounds v[k] = m
      k += 1
    end
  end

  for ci in sterms

    @inbounds for i in mjs
      m[i] = zi
    end

    for j in js

      s = z
      for i in js
        s[i] = c[j,2]
      end

      for q in 1:dim
        for d in 1:dim
          if d != q
            @inbounds s[q] *= c[d,ci[d]]
          else
            @inbounds s[q] *= g[d,ci[d]]
          end
        end
      end
      aux = o
      @inbounds for d in 1:dim
        aux *= c[d,ci[d]]
      end
      s[j] += aux

      for i in js
        @inbounds m[i,j] = s[i]
      end
    end
    @inbounds v[k] = m
    k += 1
  end
end




##

function RaviartThomasRefFE_t(::Type{et},p::Polytope,order::Integer) where et

  D = num_dims(p)

  prebasis = PCurlGradMonomialBasis{D}(et,order)

  # @show prebasis

  nf_nodes, nf_moments = _RT_nodes_and_moments(et,p,order)

  face_own_dofs = _face_own_dofs_from_moments(nf_moments)

  face_dofs = face_own_dofs

  dof_basis = MomentBasedDofBasis(nf_nodes, nf_moments)

  ndofs = num_dofs(dof_basis)

  metadata = nothing

  reffe = GenericRefFE(
    ndofs,
    p,
    prebasis,
    dof_basis,
    DivConformity(),
    metadata,
    face_dofs)

  reffe
end

function _RT_nodes_and_moments(::Type{et}, p::Polytope, order::Integer) where et

  D = num_dims(p)
  ft = VectorValue{D,et}
  pt = Point{D,et}

  nf_nodes = [ zeros(pt,0) for face in 1:num_faces(p)]
  nf_moments = [ zeros(ft,0,0) for face in 1:num_faces(p)]

  fcips, fmoments = _RT_face_values(p,et,order)

  frange = get_dimrange(p,D-1)
  nf_nodes[frange] = fcips
  nf_moments[frange] = fmoments

  if (order > 0)
    if is_n_cube(p)
      ccips, cmoments = _RT_cell_values(p,et,order)
    else
      ccips, cmoments = _RT_cell_values_t(p,et,order)
    end
    crange = get_dimrange(p,D)
    nf_nodes[crange] = ccips
    nf_moments[crange] = cmoments
  end

  nf_nodes, nf_moments
end



# Ref FE to faces geomaps
function _ref_face_to_faces_geomap(p,fp)
  cfvs = get_face_coordinates(p,num_dims(fp))
  nc = length(cfvs)
  freffe = LagrangianRefFE(Float64,fp,1)
  fshfs = get_shapefuns(freffe)
  cfshfs = fill(fshfs, nc)
  fgeomap = lincomb(cfshfs,cfvs)
end

function _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)
  nc = length(fgeomap)
  c_fips = fill(fips,nc)
  c_wips = fill(wips,nc)
  pquad = evaluate(fgeomap,c_fips)
  c_fips, pquad, c_wips
end

function _broadcast(::Type{T},n,b) where T
  c = Array{T}(undef,size(b))
  for (ii, i) in enumerate(b)
    c[ii] = i*n
  end
  return c
end

function _RT_face_moments_t(p, fshfs, c_fips, fcips, fwips)
  nc = length(c_fips)
  cfshfs = fill(fshfs, nc)
  cvals = evaluate(cfshfs,c_fips)
  cvals = [fwips[i].*cvals[i] for i in 1:nc]
  # fns, os = get_facet_normals(p)
  fns = get_facet_normals(p)
  os = get_facet_orientations(p)
  # @santiagobadia : Temporary hack for making it work for structured hex meshes
  cvals = [ _broadcast(typeof(n),n,b) for (n,b) in zip(fns,cvals)]
  #cvals = [ _broadcast(typeof(n),n,b) for (n,b) in zip(fns,cvals)]
  return cvals
end

# It provides for every face the nodes and the moments arrays
function _RT_face_values(p,et,order)

  # Reference facet
  @assert is_simplex(p) || is_n_cube(p) "We are assuming that all n-faces of the same n-dim are the same."
  fp = Polytope{num_dims(p)-1}(p,1)

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp)

  # Nodes are integration points (for exact integration)
  # Thus, we define the integration points in the reference
  # face polytope (fips and wips). Next, we consider the
  # n-face-wise arrays of nodes in fp (constant cell array c_fips)
  # the one of the points in the polytope after applying the geopmap
  # (fcips), and the weights for these nodes (fwips, a constant cell array)
  # Nodes (fcips)
  degree = (order+1)*2
  fquad = Quadrature(fp,degree)
  fips = get_coordinates(fquad)
  wips = get_weights(fquad)

  c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Moments (fmoments)
  # The RT prebasis is expressed in terms of shape function
  fshfs = MonomialBasis{num_dims(fp)}(et,order, _p_filter)

  # Face moments, i.e., M(Fi)_{ab} = q_RF^a(xgp_RFi^b) w_Fi^b n_Fi ⋅ ()
  fmoments = _RT_face_moments_t(p, fshfs, c_fips, fcips, fwips)

  return fcips, fmoments

end

function _RT_cell_moments(p, cbasis, ccips, cwips)
  # Interior DOFs-related basis evaluated at interior integration points
  ishfs_iips = evaluate(cbasis,ccips)
  return cwips.*ishfs_iips
end

# It provides for every cell the nodes and the moments arrays
function _RT_cell_values_t(p,et,order)
  # Compute integration points at interior
  degree = 2*(order+1)
  iquad = Quadrature(p,degree)
  ccips = get_coordinates(iquad)
  cwips = get_weights(iquad)

  # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
  T = VectorValue{num_dims(p),et}
  cbasis = MonomialBasis{num_dims(p)}(T,order-1, _p_filter)
  cmoments = _RT_cell_moments(p, cbasis, ccips, cwips )
  return [ccips], [cmoments]

end

function _face_own_dofs_from_moments(f_moments)
  face_dofs = Vector{Int}[]
  o = 1
  for moments in f_moments
    ndofs = size(moments,2)
    dofs = collect(o:(o+ndofs-1))
    push!(face_dofs,dofs)
    o += ndofs
  end
  face_dofs
end

function _trivial_face_own_dofs_permutations(face_own_dofs)
  [ [collect(Int,1:length(dofs)),]  for dofs in face_own_dofs  ]
end

end # module
