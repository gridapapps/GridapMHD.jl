
"""
  Mesh streching according to Smolentsev formula
"""
function strechMHD(coord;domain=(0.0,1.0,0.0,1.0,0.0,1.0),factor=(1.0,1.0,1.0),dirs=(1,2,3))
  ncoord = collect(coord.data)
  for (i,dir) in enumerate(dirs)
    ξ0 = domain[i*2-1]
    ξ1 = domain[i*2]
    l =  ξ1 - ξ0
    c = (factor[i] + 1)/(factor[i] - 1)

    if l > 0
      if ξ0 <= coord[dir] <= ξ1
        ξx = (coord[dir] - ξ0)/l                     # ξx from 0 to 1 uniformly distributed
        ξx_streched = factor[i]*(c^ξx-1)/(1+c^ξx)    # ξx streched from 0 to 1 towards 1
        ncoord[dir] =  ξx_streched*l + ξ0            # coords streched towards ξ1
      end
    else
      if ξ1 <= coord[dir] <= ξ0
        ξx = (coord[dir] - ξ0)/l                     # ξx from 0 to 1 uniformly distributed
        ξx_streched = factor[i]*(c^ξx-1)/(1+c^ξx)    # ξx streched from 0 to 1 towards 1
        ncoord[dir] =  ξx_streched*l + ξ0            # coords streched towards ξ1
      end
    end
  end
  return VectorValue(ncoord)
end

function hunt_stretch_map(
  L::Real,Ha::Real,kmap_x::Number,kmap_y::Number,BL_adapted::Bool
)
  strech_Ha = sqrt(Ha/(Ha-1))
  strech_side = sqrt(sqrt(Ha)/(sqrt(Ha)-1))
  function map1(x)
    y = strechMHD(x,domain=(0,-L,0,-L),factor=(strech_side,strech_Ha),dirs=(1,2))
    z = strechMHD(y,domain=(0,L,0,L),factor=(strech_side,strech_Ha),dirs=(1,2))
    return z  
  end
  function map2(x)
    layer(x,a) = sign(x)*abs(x)^(1/a)
    return VectorValue(layer(x[1],kmap_x),layer(x[2],kmap_y),x[3])
  end
  coord_map = BL_adapted ? map1 : map2
  return coord_map
end

function hunt_add_tags!(model)
  labels = get_face_labeling(model)
  tags_u = append!(collect(1:20),[23,24,25,26])
  tags_j = append!(collect(1:20),[25,26])
  add_tag_from_tags!(labels,"noslip",tags_u)
  add_tag_from_tags!(labels,"insulating",tags_j)
end

"""
    hunt_generate_base_mesh(ranks,nc::Tuple,np::Tuple,L::Real,Ha::Real,kmap_x::Number,kmap_y::Number,BL_adapted::Bool)

  Generate a mesh for the Hunt problem, distributed amongst the communicator `ranks`.
"""
function hunt_generate_base_mesh(
  ranks,np::Tuple,
  nc::Tuple,L::Real,Ha::Real,kmap_x::Number,kmap_y::Number,BL_adapted::Bool
)
  coord_map = hunt_stretch_map(L,Ha,kmap_x,kmap_y,BL_adapted)
  _nc = (nc[1],nc[2],3)
  _np = (np[1],np[2],1)
  domain = (-1.0,1.0,-1.0,1.0,0.0,0.1)
  model = CartesianDiscreteModel(ranks,_np,domain,_nc;isperiodic=(false,false,true),map=coord_map)
  hunt_add_tags!(model)
  return model
end

"""
    hunt_generate_base_mesh(nc::Tuple,np::Tuple,L::Real,Ha::Real,kmap_x::Number,kmap_y::Number,BL_adapted::Bool)

  Generate a serial mesh for the Hunt problem.
"""
function hunt_generate_base_mesh(
  nc::Tuple,L::Real,Ha::Real,kmap_x::Number,kmap_y::Number,BL_adapted::Bool
)
  coord_map = hunt_stretch_map(L,Ha,kmap_x,kmap_y,BL_adapted)
  _nc = (nc[1],nc[2],3)
  domain = (-1.0,1.0,-1.0,1.0,0.0,0.1)
  model = CartesianDiscreteModel(domain,_nc;isperiodic=(false,false,true),map=coord_map)
  hunt_add_tags!(model)
  return model
end
