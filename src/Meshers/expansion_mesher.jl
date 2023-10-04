
# Coordinate transformations

rotate(x) = rotate(x,[2,3,1])
function rotate(x::VectorValue{3},p)
  return VectorValue{3}(x[p[1]],x[p[2]],x[p[3]])
end
function stretch(x::VectorValue{3})
  offset = VectorValue{3}(0.0,-1.0,-1.0)
  α = 8.0; β = 2.0/3.0; γ = 2.0;
  return VectorValue{3}(α⋅x[1],β⋅x[2],γ⋅x[3]) + offset
end
coordinate_transformation(x) = stretch(rotate(VectorValue{3}(x...)))

# Generation of the base mesh

function expansion_generate_face_labeling()
  num_faces = [20,36,21,4]
  PbLi_dfaces = [
    [],       # Nodes
    [],       # Edges
    [2,7,8],  # Faces
    [1,2,3,4] # Cells
  ]
  inlet_dfaces = [
    [17,18,19,20],
    [29,30,36,35],
    [19],
    []
  ]
  outlet_dfaces = [
    [3,4,7,8,11,12,15,16],
    [2,4,11,12,14,19,20,22,27,28],
    [4,9,14],
    []
  ]
  wall_dfaces = [
    [1,2,5,6,9,10,13,14],
    [1,3,5,6,7,8,9,10,13,15,16,17,18,21,23,24,25,26,31,32,33,34],
    [1,3,5,6,10,11,12,13,15,16,17,18,20,21],
    []
  ]
  d_to_dface_to_entity = map([0,1,2,3]) do d
    entity_to_dface = [PbLi_dfaces[d+1],inlet_dfaces[d+1],outlet_dfaces[d+1],wall_dfaces[d+1]]
    dface_to_entity = map(f -> findfirst(faces -> f ∈ faces, entity_to_dface),1:num_faces[d+1])
    return dface_to_entity
  end
  tag_to_name = ["PbLi","inlet","outlet","wall","interior","boundary"]
  tag_to_entities = [[1],[2],[3],[4],[1],[2,3,4]]
  return FaceLabeling(d_to_dface_to_entity,tag_to_entities,tag_to_name)
end

function expansion_asymmetric_refinement(model)
  vertices = [
    (0,0,0)  ,(1,0,0)  ,(0,1/3,0)  ,(1,1/3,0)  ,(0,2/3,0)  ,(1,2/3,0)  ,(0,1,0)  ,(1,1,0),
    (0,0,1)  ,(1,0,1)  ,(0,1/3,1)  ,(1,1/3,1)  ,(0,2/3,1)  ,(1,2/3,1)  ,(0,1,1)  ,(1,1,1),
    (0,0,2)  ,(1,0,2)  ,(0,1/3,2)  ,(1,1/3,2)  ,(0,2/3,2)  ,(1,2/3,2)  ,(0,1,2)  ,(1,1,2),
    (0,0,3)  ,(1,0,3)  ,(0,1/3,3)  ,(1,1/3,3)  ,(0,2/3,3)  ,(1,2/3,3)  ,(0,1,3)  ,(1,1,3),
    (0,-1,1) ,(1,-1,1) ,(0,-2/3,1) ,(1,-2/3,1) ,(0,-1/3,1) ,(1,-1/3,1),
    (0,-1,2) ,(1,-1,2) ,(0,-2/3,2) ,(1,-2/3,2) ,(0,-1/3,2) ,(1,-1/3,2),
  ]
  vertices = map(coordinate_transformation,vertices)

  c2n_map = Table([
    [ 1, 2, 3, 4, 9,10,11,12],[ 3, 4, 5, 6,11,12,13,14],[ 5, 6, 7, 8,13,14,15,16],
    [ 9,10,11,12,17,18,19,20],[11,12,13,14,19,20,21,22],[13,14,15,16,21,22,23,24],
    [17,18,19,20,25,26,27,28],[19,20,21,22,27,28,29,30],[21,22,23,24,29,30,31,32],
    [33,34,35,36,39,40,41,42],[35,36,37,38,41,42,43,44],[37,38, 9,10,43,44,17,18],
  ])

  n2o_cells_map = [1,1,1,2,2,2,3,3,3,4,4,4]
  n2o_cell_to_child_id = [1,2,3,1,2,3,1,2,3,1,2,3]
  refinement_rules = fill(RefinementRule(HEX,(1,3,1)),4)
  n2o_faces_map = [Int[],Int[],Int[],n2o_cells_map]
  glue = AdaptivityGlue(n2o_faces_map,n2o_cell_to_child_id,refinement_rules)

  cell_type = fill(1,length(c2n_map)); polys = [HEX]; orientation = Oriented()
  reffes    = map(p->LagrangianRefFE(Float64,p,1),polys)
  new_grid  = UnstructuredGrid(vertices,c2n_map,reffes,cell_type,orientation)
  new_topo  = UnstructuredGridTopology(vertices,c2n_map,cell_type,polys,orientation)
  new_labels = Gridap.Adaptivity._refine_face_labeling(
    get_face_labeling(model),
    glue,
    get_grid_topology(model),
    new_topo
  )
  new_model = UnstructuredDiscreteModel(new_grid,new_topo,new_labels)
  return new_model
end

function expansion_generate_base_mesh()
  vertices = [
    (0,0,0),(1,0,0),(0,1,0),(1,1,0), # Main block, 1st layer
    (0,0,1),(1,0,1),(0,1,1),(1,1,1), # Main block, 2nd layer
    (0,0,2),(1,0,2),(0,1,2),(1,1,2), # Main block, 3rd layer
    (0,0,3),(1,0,3),(0,1,3),(1,1,3), # Main block, 4th layer
    (0,-1,1),(1,-1,1),               # Extension,  2nd layer
    (0,-1,2),(1,-1,2),               # Extension,  3rd layer
  ]
  vertices = map(coordinate_transformation,vertices)

  c2n_map = Table([
    [ 1, 2, 3, 4, 5, 6, 7, 8],
    [ 5, 6, 7, 8, 9,10,11,12],
    [ 9,10,11,12,13,14,15,16],
    [17,18, 5, 6,19,20, 9,10]
  ])

  cell_type = [1,1,1,1]; polys = [HEX]; orientation = Oriented()
  reffes = map(p->LagrangianRefFE(Float64,p,1),polys)
  grid   = UnstructuredGrid(vertices,c2n_map,reffes,cell_type,orientation)
  topo   = UnstructuredGridTopology(vertices,c2n_map,cell_type,polys,orientation)
  labels = expansion_generate_face_labeling()
  model  = UnstructuredDiscreteModel(grid,topo,labels)
  asym_model = expansion_asymmetric_refinement(model)
  return asym_model
end

# Mesh refinement 

"""
  expansion_generate_mesh(ranks,num_refs)

  ranks    :: Number of processors of the distributed mesh
  num_refs :: Number of refinements to perform. Final mesh will have 
                num_cells = 12*2^(3*num_refs)
"""
function expansion_generate_mesh(ranks,num_refs)
  model = GridapP4est.with(ranks) do
    cmodel = expansion_generate_base_mesh()
    return OctreeDistributedDiscreteModel(ranks,cmodel,num_refs)
  end
  return model
end

"""
  expansion_generate_mesh_hierarchy(level_ranks,num_refs)

  ranks           :: Number of processors
  ranks_per_level :: Number of processors per level
  num_refs_coarse :: Number of refinements to perform for the coarsest model. 
                     Finest mesh will have 
                        num_cells = 12*2^(3*(num_refs+num_levels))
                     where 
                        num_levels = length(ranks_per_level)
"""
function expansion_generate_mesh_hierarchy(ranks,num_refs_coarse::Integer,ranks_per_level::Vector{<:Integer})
  mh = GridapP4est.with(ranks) do
    num_levels = length(ranks_per_level)
    cparts     = generate_subparts(ranks,ranks_per_level[num_levels])
    cmodel     = expansion_generate_base_mesh()
    model      = OctreeDistributedDiscreteModel(cparts,cmodel,num_refs_coarse)
    return ModelHierarchy(ranks,model,ranks_per_level)
  end
  return mh
end
