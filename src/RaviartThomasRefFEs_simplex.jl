module RaviartThomasRefFE_simplex

using Gridap
using Gridap.Helpers
using Gridap.Fields
using Gridap.Polynomials
using Gridap.Polynomials: _prepare_perms
using Gridap.ReferenceFEs
using Gridap.ReferenceFEs: _RT_face_values
using Gridap.ReferenceFEs: _face_own_dofs_from_moments

import Gridap.Polynomials: get_order
import Gridap.Polynomials: get_orders
import Gridap.Polynomials: get_value_type
import Gridap.Polynomials: num_terms

export PGradMonomialBasis
export PCurlGradMonomialBasis

export get_order
export get_orders
export get_value_type
export num_terms

export RaviartThomasRefFE_t


struct PGradMonomialBasis{D,T} <: Field
  order::Int
  terms::Array{CartesianIndex{D},1}
  perms::Matrix{Int}
  function PGradMonomialBasis(::Type{T},order::Int,terms::Array{CartesianIndex{D},1},perms::Matrix{Int}) where {D,T}
    new{D,T}(order,terms,perms)
  end
end

"""
    num_terms(f::PGradMonomialBasis{D,T}) where {D,T}
"""
num_terms(f::PGradMonomialBasis{D,T}) where {D,T} = length(f.terms)*D

get_order(f::PGradMonomialBasis) = f.order

function field_cache(f::PGradMonomialBasis{D,T},x) where {D,T}
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

function evaluate_field!(cache,f::PGradMonomialBasis{D,T},x) where {D,T}
  r, v, c = cache
  np = length(x)
  ndof = _ndofs_pgrad(f)
  n = 1 + f.order+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
      _evaluate_nd_pgrad!(v,xi,f.order+1,f.terms,f.perms,c)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r
end

function gradient_cache(f::PGradMonomialBasis{D,T},x) where {D,T}
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

function evaluate_gradient!(cache,f::PGradMonomialBasis{D,T},x) where {D,T}
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
    _gradient_nd_pgrad!(v,xi,f.order+1,f.terms,f.perms,c,g,V)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r
end

# Helpers

_ndofs_pgrad(f::PGradMonomialBasis{D}) where D = D*(length(f.terms))

function _evaluate_nd_pgrad!(
  v::AbstractVector{V},
  x,
  order,
  terms::Array{CartesianIndex{D},1},
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

  for ci in terms

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

end

function _gradient_nd_pgrad!(
  v::AbstractVector{G},
  x,
  order,
  terms::Array{CartesianIndex{D},1},
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

  for ci in terms

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

end


struct PCurlGradMonomialBasis{D,T} <: Field
  pgrad::PGradMonomialBasis{D,T}
  function PCurlGradMonomialBasis(::Type{T},order::Int,terms::Array{CartesianIndex{D},1},perms::Matrix{Int}) where {D,T}
    pgrad = PGradMonomialBasis(T,order,terms,perms)
    new{D,T}(pgrad)
  end
end

_p_filter(e,order) = (sum(e) <= order)
_t_filter(e,order) = (sum(e) == order)

function PCurlGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}
  @assert T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
  P_k = MonomialBasis{D}(T, order, _p_filter)
  T_k = MonomialBasis{D}(T, order, _t_filter)
  g = Vector(undef,D)
  g .= 0 ; g[1] = 1
  g = CartesianIndex(Tuple(g))
  for i in 1:length(T_k.terms)
    T_k.terms[i] += g
  end
  terms = [P_k.terms; T_k.terms ]
  perms = _prepare_perms(D)
  PCurlGradMonomialBasis(T,order,terms,perms)
end



function field_cache(f::PCurlGradMonomialBasis,x)
  field_cache(f.pgrad,x)
end

@inline function evaluate_field!(cache,f::PCurlGradMonomialBasis,x)
  evaluate_field!(cache,f.pgrad,x)
end

function gradient_cache(f::PCurlGradMonomialBasis,x)
  gradient_cache(f.pgrad,x)
end

@inline function evaluate_gradient!(cache,f::PCurlGradMonomialBasis,x)
  evaluate_gradient!(cache,f.pgrad,x)
end

get_value_type(::PCurlGradMonomialBasis{D,T}) where {D,T} = T

"""
    num_terms(f::QCurlGradMonomialBasis{D,T}) where {D,T}
"""
function num_terms(f::PCurlGradMonomialBasis{D,T}) where {D,T}
  if D == 2
    return (length(f.pgrad.terms)-f.pgrad.order-1)*D + f.pgrad.order+1
  elseif D == 3
    return (f.pgrad.order+4)*(f.pgrad.order+2)*(f.pgrad.order+1)/2
  else
    @notimplemented "Not implented for arbitrary dimensions yet"
  end
end

get_order(f::PCurlGradMonomialBasis{D,T}) where {D,T} = get_order(f.pgrad)


function RaviartThomasRefFE_t(::Type{et},p::Polytope,order::Integer) where et

  D = num_dims(p)

  prebasis = PCurlGradMonomialBasis{D}(et,order)

  @show prebasis

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

  if is_n_cube(p)
    fcips, fmoments = _RT_face_values(p,et,order)
  else
    fcips, fmoments = _RT_face_values(p,et,order)
  end
  frange = get_dimrange(p,D-1)
  nf_nodes[frange] = fcips
  nf_moments[frange] = fmoments

  if (order > 0)
    ccips, cmoments = _RT_cell_values(p,et,order)
    crange = get_dimrange(p,D)
    nf_nodes[crange] = ccips
    nf_moments[crange] = cmoments
  end

  nf_nodes, nf_moments
end

end # module
