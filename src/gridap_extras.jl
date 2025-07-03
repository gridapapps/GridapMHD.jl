
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

# MacroReferenceFEs

function conformity_from_symbol(sym::Symbol)
  if sym == :H1
    return H1Conformity()
  elseif sym == :L2
    return L2Conformity()
  else
    @assert false
  end
end

function Gridap.Adaptivity.MacroReferenceFE(
  rrule::Gridap.Adaptivity.RefinementRule,
  reffes::AbstractVector{<:Tuple};
  macro_kwargs...
)
  polys = Gridap.Adaptivity.get_cell_polytopes(rrule)
  _reffes = map(polys,reffes) do p,r
    basis, args, kwargs = r
    ReferenceFE(p,basis,args...;kwargs...)
  end
  return Gridap.Adaptivity.MacroReferenceFE(rrule,_reffes;macro_kwargs...)
end

function Gridap.Adaptivity.MacroReferenceFE(
  rrule::Gridap.Adaptivity.RefinementRule,
  reffe::Tuple;
  macro_kwargs...
)
  basis, args, kwargs = reffe
  reffes = ReferenceFE(rrule.ref_grid,basis,args...;kwargs...)
  return Gridap.Adaptivity.MacroReferenceFE(rrule,reffes;macro_kwargs...)
end

function Gridap.Adaptivity.MacroReferenceFE(
  rrule::Gridap.Adaptivity.RefinementRule,
  reffe::ReferenceFE;
  macro_kwargs...
)
  reffes = Fill(reffe,Gridap.Adaptivity.num_subcells(rrule))
  return Gridap.Adaptivity.MacroReferenceFE(rrule,reffes;macro_kwargs...)
end
