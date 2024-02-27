
using DrWatson
using PartitionedArrays
using Gridap, GridapP4est

using Gridap.Arrays, Gridap.Geometry, Gridap.ReferenceFEs, Gridap.Adaptivity

function asymmetric_refinement_bis(model)
  vertices = [
    (0,0,0)  ,(1,0,0)  ,(0,1/3,0)  ,(1,1/3,0)  ,(0,2/3,0)  ,(1,2/3,0)  ,(0,1,0)  ,(1,1,0),
    (0,0,1/2),(1,0,1/2),(0,1/3,1/2),(1,1/3,1/2),(0,2/3,1/2),(1,2/3,1/2),(0,1,1/2),(1,1,1/2),
    (0,0,1)  ,(1,0,1)  ,(0,1/3,1)  ,(1,1/3,1)  ,(0,2/3,1)  ,(1,2/3,1)  ,(0,1,1)  ,(1,1,1),
    (0,0,2)  ,(1,0,2)  ,(0,1/3,2)  ,(1,1/3,2)  ,(0,2/3,2)  ,(1,2/3,2)  ,(0,1,2)  ,(1,1,2),
    (0,0,5/2),(1,0,5/2),(0,1/3,5/2),(1,1/3,5/2),(0,2/3,5/2),(1,2/3,5/2),(0,1,5/2),(1,1,5/2),
    (0,0,3)  ,(1,0,3)  ,(0,1/3,3)  ,(1,1/3,3)  ,(0,2/3,3)  ,(1,2/3,3)  ,(0,1,3)  ,(1,1,3),
    (0,-1,1) ,(1,-1,1) ,(0,-2/3,1) ,(1,-2/3,1) ,(0,-1/3,1) ,(1,-1/3,1),
    (0,-1,2) ,(1,-1,2) ,(0,-2/3,2) ,(1,-2/3,2) ,(0,-1/3,2) ,(1,-1/3,2),
  ]
  vertices = map(x -> rotate(stretch(VectorValue{3,Float64}(x...))),vertices)

  c2n_map = Table([
    [ 1, 2, 3, 4, 9,10,11,12],[ 3, 4, 5, 6,11,12,13,14],[ 5, 6, 7, 8,13,14,15,16],
    [ 9,10,11,12,17,18,19,20],[11,12,13,14,19,20,21,22],[13,14,15,16,21,22,23,24],
    [17,18,19,20,25,26,27,28],[19,20,21,22,27,28,29,30],[21,22,23,24,29,30,31,32],
    [25,26,27,28,33,34,35,36],[27,28,29,30,35,36,37,38],[29,30,31,32,37,38,39,40],
    [33,34,35,36,41,42,43,44],[35,36,37,38,43,44,45,46],[37,38,39,40,45,46,47,48],
    [49,50,51,52,55,56,57,58],[51,52,53,54,57,58,59,60],[53,54,17,18,59,60,25,26],
  ])

  n2o_cells_map = [
    1,1,1,1,1,1,
    2,2,2,
    3,3,3,3,3,3,
    4,4,4
  ]
  n2o_cell_to_child_id = [
    1,2,3,4,5,6,
    1,2,3,
    1,2,3,4,5,6,
    1,2,3
  ]

  rr_1x3x1 = RefinementRule(HEX,(1,3,1))
  rr_1x3x2 = RefinementRule(HEX,(1,3,2))
  refinement_rules = [rr_1x3x2,rr_1x3x1,rr_1x3x2,rr_1x3x1]
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

function generate_face_labeling()
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
  tag_to_name = ["PbLi","inlet","outlet","wall"]
  tag_to_entities = [[1],[2],[3],[4]]
  return FaceLabeling(d_to_dface_to_entity,tag_to_entities,tag_to_name)
end

function generate_base_mesh()
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
  labels = generate_face_labeling()
  model  = UnstructuredDiscreteModel(grid,topo,labels)
  asym_model = asymmetric_refinement(model)
  return asym_model
end

function asymmetric_refinement(model)
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

function refine_mesh(ranks,model,num_refs)
  ref_model = GridapP4est.with(ranks) do
    ref_model = OctreeDistributedDiscreteModel(ranks,model,num_refs)
    return ref_model
  end
  return ref_model
end

function write_bcs(model)
  order = 2
  reffe = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(model,reffe,dirichlet_tags=["inlet","outlet","wall"])
  U = TrialFESpace(V,[1.0,2.0,3.0])
  uh = interpolate(-1.0,U)
  writevtk(Triangulation(model),datadir("P4est/expansion_bcs"),cellfields=["uh"=>uh],order=order)
end

function write_mesh(model)
  trian = Triangulation(model)
  writevtk(trian,datadir("P4est/expansion"))
end

function main(distribute,np)
  ranks  = distribute(LinearIndices((prod(np),)))
  cmodel = generate_base_mesh()
  model  = refine_mesh(ranks,cmodel,3)
  #write_mesh(model)
  write_bcs(model)
  return nothing
end

with_mpi() do distribute
  main(distribute,(1,1,1))
end
