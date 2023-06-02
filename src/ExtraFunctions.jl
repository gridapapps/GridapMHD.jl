##Functions for mesh manipulation

function strechMHD(coord;domain=(0.0,1.0,0.0,1.0,0.0,1.0),factor=(1.0,1.0,1.0),dirs=(1,2,3))~
##Mesh streching according to Smolentsev formula
ncoord = collect(coord.data)
for (i,dir) in enumerate(dirs)
  ξ0 = domain[i*2-1]
  ξ1 = domain[i*2]
  l =  ξ1 - ξ0
  c = (factor[i] + 1)/(factor[i] - 1)

  if l > 0
        if ξ0 <= coord[dir] <= ξ1
           ξx = (coord[dir] - ξ0)/l                     #ξx from 0 to 1 uniformly distributed
           ξx_streched = factor[i]*(c^ξx-1)/(1+c^ξx)    #ξx streched from 0 to 1 towards 1
           ncoord[dir] =  ξx_streched*l + ξ0            #coords streched towards ξ1
        end
  else
        if ξ1 <= coord[dir] <= ξ0
           ξx = (coord[dir] - ξ0)/l                        #ξx from 0 to 1 uniformly distributed
           ξx_streched = factor[i]*(c^ξx-1)/(1+c^ξx)       #ξx streched from 0 to 1 towards 1
           ncoord[dir] =  ξx_streched*l + ξ0               #coords streched towards ξ1

        end
  end
end
VectorValue(ncoord)
end

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
VectorValue(ncoord)
end

##Functions for geometry generation

function cuboid(;dx=1,dy=1,dz=1,x0=Point(0,0,0),name="cuboid",
  faces=["face$i" for i in 1:6])
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

  intersect(intersect(geo12,geo34),geo56,name=name)
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
