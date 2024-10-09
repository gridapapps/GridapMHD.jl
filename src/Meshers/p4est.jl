
"""
    expansion_generate_mesh(ranks,num_refs)

  ranks      :: Number of processors of the distributed mesh
  base_model :: Base serial model
  num_refs   :: Number of refinements to perform. 
  
  Final mesh will have nc = num_cells(base_model)*2^(3*num_refs) cells.
"""
function generate_refined_mesh(
  ranks,
  base_model::DiscreteModel,
  num_refs::Integer
)
  model = GridapP4est.with(ranks) do
    OctreeDistributedDiscreteModel(ranks,base_model,num_refs)
  end
  return model.dmodel
end

"""
    expansion_generate_mesh_hierarchy(level_ranks,num_refs)

  ranks           :: Number of processors
  base_model      :: Base serial model
  ranks_per_level :: Number of processors per level
  num_refs_coarse :: Number of refinements to perform for the coarsest model.

  Finest mesh will have 
     num_cells = num_cells(base_model)*2^(3*(num_refs_coarse+nlevels))
  where 
     nlevels = length(ranks_per_level)
"""
function generate_mesh_hierarchy(
  ranks,
  base_model::DiscreteModel,
  num_refs_coarse::Integer,
  ranks_per_level::Vector{<:Integer}
)
  mh = GridapP4est.with(ranks) do
    nlevels = length(ranks_per_level)
    cparts  = generate_subparts(ranks,ranks_per_level[nlevels])
    model   = OctreeDistributedDiscreteModel(cparts,base_model,num_refs_coarse)
    return ModelHierarchy(ranks,model,ranks_per_level)
  end
  return mh
end
