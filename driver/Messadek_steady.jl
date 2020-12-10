
using Gridap
using GridapEmbedded
using LineSearches: Static

using GridapMHD
using GridapMHD: driver_inductionless_MHD

# Problem setting
ρ = 1.0
ν = 1.0
σ = 1.0

U0 = 10.0
B0 = 20.0
L = 1.0
Re = U0 * L / ν
Ha = B0 * L * sqrt(σ/(ρ*ν))

function g_u(x)
  if abs(x[1] + 0.05) < 1e-8
    y = x[2]; z = x[3]
    return VectorValue(1.0/0.0125^4 * (y-0.0125) * (y+0.0125) *
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

function map1(coord)
  ncoord = streching(coord,  domain=(0.0,0.025,-0.05,-0.0375,-0.0125,0.0),factors=(3.0,3.0,3.0),dirs=(1,2,3))
  ncoord = streching(ncoord, domain=(0.025,0.05,-0.0375,-0.025,0.0,0.0125),factors=(1.25,1.25,1.25),dirs=(1,2,3))
  ncoord = streching(ncoord, domain=(0.175,0.2,-0.0125,0.0),factors=(1.25,3.0),dirs=(1,2))
  ncoord = streching(ncoord, domain=(0.0,0.0125),factors=(1.25,),dirs=(2,))
  ncoord = streching(ncoord, domain=(0.025,0.0375),factors=(3.0,),dirs=(2,))
  ncoord = streching(ncoord, domain=(0.0375,0.05),factors=(1.25,),dirs=(2,))
  ncoord
end

# Background mesh definition
partition = (20,16,5);
#partition = (60,24,6);
domain = (-0.05,0.2,-0.05,0.05,-0.0125,0.0125)
bgmodel=CartesianDiscreteModel(domain,partition,map=map1);

# Build the final domain using simple geomtries
inlet = cuboid(dx=0.05-1e-6,dy=0.025-1e-6,dz=0.025,x0=Point(-0.025,0.0,0.0));
out = cuboid(dx=0.2-1e-6,dy=0.1-1e-6,dz=0.025,x0=Point(0.1,0.0,0.0));

out2 = cuboid(dx=0.15+1e-6,dy=0.0125+1e-6,dz=0.025,x0=Point(0.125,0.01875,0.0));
out3 = cuboid(dx=0.15+1e-6,dy=0.0125+1e-6,dz=0.025,x0=Point(0.125,-0.01875,0.0));
out1 = setdiff(out,out3);
outlet = setdiff(out1,out2);

# Add the simple geomtries and "cut"
geo = union(inlet,outlet,name="csg");

# "Cut" the final geomtry over the Background mesh
cutgeo = cut(bgmodel,geo);
model = DiscreteModel(cutgeo,"csg");

# Filter to select all walls
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

# Filter to select inflow boundary
function inlet_boundary(coord)
  x = coord[1]
  y = coord[2]
  tol = 1e-6

  if (x<-0.05+tol)
    return true
  end
  false
end

# Filter to select outflow boundary
function outlet_boundary(coord)
  x = coord[1]
  y = coord[2]
  tol = 1e-6

  if (x>0.2-tol)
    return true
  end
  false
end

# Name inlet, outlet, and walls element sets.
add_entity!(model,walls,"walls")
add_entity!(model,inlet_boundary,"inlet")
add_entity!(model,outlet_boundary,"outlet")

# Call main MHD driver
xh, trian, dΩ = driver_inductionless_MHD(;
  Re=Re,
  Ha=Ha,
  model=model,
  fluid_dirichlet_tags = ["inlet","walls"],
  fluid_neumann_tags = ["outlet"],
  magnetic_dirichlet_tags = ["walls"],
  magnetic_neumann_tags = ["inlet","outlet"],
  fluid_dirichlet_conditions = g_u ,
  magnetic_dirichlet_conditions = g_j,
  fluid_body_force = f_u,
  linesearch=Static(),
  constraint_presures = (false,false),
  precond_tau = 1e-10
)

uh, ph, jh, φh = xh

# Scale unknowns
uh = uh * U0
jh = jh * σ * B0 * U0

uh_r = restrict(uh,trian)
ph_r = restrict(ph,trian)
jh_r = restrict(jh,trian)
φh_r = restrict(φh,trian)
writevtk(trian,"results",cellfields=["u"=>uh_r,"p"=>ph_r,"j"=>jh_r,"phi"=>φh_r])
