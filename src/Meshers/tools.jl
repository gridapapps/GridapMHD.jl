function add_non_slip_at_solid_entity!(model,solid_entities,fluid_entity,name)
  D = num_cell_dims(model)
  labels = get_face_labeling(model)
  topo = get_grid_topology(model)
  cell_entity = get_cell_entity(labels)
  for d in 0:D-1  
    dface_entity = labels.d_to_dface_to_entity[d+1]
    dface_cells = get_faces(topo,d,D)
    cache = array_cache(dface_cells)
    for dface in 1:length(dface_cells)
      cells = getindex!(cache,dface_cells,dface)
      solid_found = false
      fluid_found = false
      for cell in cells
        solid_found = solid_found || cell_entity[cell] in solid_entities
        fluid_found = fluid_found || cell_entity[cell] == fluid_entity
      end
      if solid_found && fluid_found
        dface_entity[dface] = name
      end
    end
  end
end
