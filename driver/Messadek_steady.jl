
using Gridap
using Gridap.Geometry:get_grid_topology, get_node_coordinates, get_face_nodes
using GridapEmbedded

using GridapMHD: driver_inductionless_MHD

# Problem setting
ρ = 1.0
ν = 1.0
σ = 1.0

U0 = 10.0
B0 = 10.0
L = 1.0
Re = U0 * L / ν
Ha = B0 * L * sqrt(σ/(ρ*ν))

function g_u(x)
  if abs(x[1] + 0.05) < 1e-8
    y = x[2]; z = x[3]
    return VectorValue(10.0/0.0125^4 * (y-0.0125) * (y+0.0125) *
                      (z-0.0125) * (z+0.0125), 0.0, 0.0)
  else
    return VectorValue(0.0, 0.0, 0.0)
  end
end

function g_j(x)
  return VectorValue(0.0,0.0,0.0)
end

function f_u(x)
  return VectorValue(0.0,0.0,0.0)
end



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

faces = ["face$i" for i in 1:6]
in_faces = copy(faces)
in_faces[3] = "inlet"
out_faces = copy(faces)
out_faces[4] = "outlet"

inlet = cuboid(dx=0.05-1e-6,dy=0.025-1e-6,dz=0.025,x0=Point(-0.025,0.0,0.0),faces=in_faces);
out = cuboid(dx=0.2-1e-6,dy=0.1-1e-6,dz=0.025,x0=Point(0.1,0.0,0.0),faces=out_faces);

out2 = cuboid(dx=0.15+1e-6,dy=0.0125+1e-6,dz=0.025,x0=Point(0.125,0.01875,0.0));
out3 = cuboid(dx=0.15+1e-6,dy=0.0125+1e-6,dz=0.025,x0=Point(0.125,-0.01875,0.0));
out1 = setdiff(out,out3);
outlet = setdiff(out1,out2);

geo = union(inlet,outlet,name="csg");

pmin=Point(-0.05,-0.05,-0.0125);
pmax=Point(0.2,0.05,0.0125);
partition = (20,8,2);
bgmodel=CartesianDiscreteModel(pmin,pmax,partition);

cutgeo = cut(bgmodel,geo);
model = DiscreteModel(cutgeo,"csg");

function walls(coord)
  x = coord[1]
  y = coord[2]
  z = coord[3]
  tol = 1e-6
  if (-0.05+tol < x < -tol) & ((0.0125-tol < y) | (y < -0.0125+tol ))
    return true
  elseif (z > 0.0125 - tol) | ( z < -0.0125 + tol)
    return true
  elseif (-tol < x < tol)
    if (-0.05-tol < y < -0.0125+tol) | (0.0125-tol < y < 0.5+tol)
      return true
    end
  elseif (0.05-tol < x < 0.2-tol) & ((0.0125-tol < y < 0.025+tol) | (-0.025-tol < y < -0.0125+tol ))
      return true
  elseif (y < -0.05 + tol ) | ( y > 0.05 - tol)
    return true
  end
  false
end

function inlet_boundary(coord)
  x = coord[1]
  y = coord[2]
  tol = 1e-6

  if (x<-0.05+tol)
    return true
  end
  false
end

function outlet_boundary(coord)
  x = coord[1]
  y = coord[2]
  tol = 1e-6

  if (x>0.2-tol)
    return true
  end
  false
end


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

add_entity!(model,walls,"walls")
add_entity!(model,inlet_boundary,"inlet")
add_entity!(model,outlet_boundary,"outlet")




xh, trian, quad = driver_inductionless_MHD(;
  Re=10.0,
  Ha=10.0,
  model=model,
  fluid_dirichlet_tags = ["inlet","walls"],
  fluid_neumann_tags = ["outlet"],
  magnetic_dirichlet_tags = ["walls"],
  magnetic_neumann_tags = ["inlet","outlet"],
  fluid_dirichlet_conditions = g_u ,
  magnetic_dirichlet_conditions = g_j,
  fluid_body_force = f_u,
  constraint_presures = (false,false),
  precond_tau = 1e-10
)

uh, ph, jh, φh = xh

uh_r = restrict(uh,trian)
ph_r = restrict(ph,trian)
jh_r = restrict(jh,trian)
φh_r = restrict(φh,trian)
writevtk(trian,"results",cellfields=["u"=>uh_r,"p"=>ph_r,"j"=>jh_r,"phi"=>φh_r])
