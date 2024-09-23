##Functions for mesh manipulation

function ChangeDensity(coord;domain=(0.0,1.0,0.0,1.0,0.0,1.0),subDomain=(0.0,1.0,0.0,1.0,0.0,1.0),
			     nodesTot=(1.0,1.0,1.0), nodesSub=(1.0,1.0,1.0), dirs=(1,2,3))
  ncoord = collect(coord.data)
  for (i,dir) in enumerate(dirs)
    ξ0 = domain[i*2-1]
    ξ1 = domain[i*2]
    Ltot =  ξ1 - ξ0
    Lsub = subDomain[i*2] - subDomain[i*2-1]

    alpha = (Lsub/Ltot)*(nodesTot[i]/nodesSub[i])
    betta = ((Ltot-Lsub)/Ltot)*(nodesTot[i]/(nodesTot[i]-nodesSub[i]))

    if Ltot != Lsub
      if ξ0 <= coord[dir] <= ξ1
        if coord[dir] <= (Lsub/alpha + ξ0)
          ncoord[dir] = alpha*coord[dir] + ξ0*(1-alpha)
        else
          ncoord[dir] = betta*coord[dir] + ξ0*(1-betta) + Lsub*(1-betta/alpha)
        end
      end
    end
  end
  return VectorValue(ncoord)
end

## Functions for geometry generation

function cuboid(;dx=1,dy=1,dz=1,x0=Point(0,0,0),name="cuboid",faces=["face$i" for i in 1:6])
  e1 = VectorValue(1,0,0)
  e2 = VectorValue(0,1,0)
  e3 = VectorValue(0,0,1)

  plane1 = plane(x0=x0-0.5*dz*e3,v=-e3,name=faces[1])
  plane2 = plane(x0=x0+0.5*dz*e3,v=+e3,name=faces[2])
  plane3 = plane(x0=x0-0.5*dy*e2,v=-e2,name=faces[3])
  plane4 = plane(x0=x0+0.5*dy*e2,v=+e2,name=faces[4])
  plane5 = plane(x0=x0-0.5*dx*e1,v=-e1,name=faces[5])
  plane6 = plane(x0=x0+0.5*dx*e1,v=+e1,name=faces[6])

  geo12 = intersect(plane1,plane2)
  geo34 = intersect(plane3,plane4)
  geo56 = intersect(plane5,plane6)

  return intersect(intersect(geo12,geo34),geo56,name=name)
end

# Function that given a filter set a name to those elements that satisfy the
# conditions of the filter
function add_entity!(model,in,name)
  labels = get_face_labeling(model)
  node_coordinates = get_node_coordinates(model)
  entity = num_entities(labels) + 1
  for d in 0:num_dims(model)-1
    facets = get_face_nodes(model,d)
    for (i,facet) in enumerate(facets)
      coord = sum(node_coordinates[facet])/length(facet)
      if in(coord)
        labels.d_to_dface_to_entity[d+1][i] = entity
      end
    end
  end
  add_tag!(labels,name,[entity])
end

# Analytical formulas for pressure drop gradients

function kp_shercliff_cartesian(b,Ha)
  kp = 1/(Ha*(1-0.852*Ha^(-0.5)/b-1/Ha))
  return kp
end

function kp_shercliff_cylinder(Ha)
  kp = (3/8)*pi/(Ha-(3/2)*pi)
  return kp
end

function kp_hunt(b,Ha)
  kp = 1/(Ha*(1-0.956*Ha^(-0.5)/b-1/Ha))
  return kp
end

function kp_tillac(b,Ha,cw_s,cw_Ha)
  k_s = (1/(3*b))*(Ha^(0.5)/(1+cw_s*Ha^(0.5)))
  k_Ha = (1+cw_Ha)/(1/Ha + cw_Ha)
  kp = 1/(k_s+k_Ha)
  return kp
end

function kp_glukhih(Ha,cw)
  kp = (3/8)*pi*(1+0.833*cw*Ha-0.019*(cw*Ha)^2)/Ha
  return kp
end
