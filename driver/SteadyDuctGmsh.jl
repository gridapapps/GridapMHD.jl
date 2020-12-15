
using Gridap
using GridapGmsh

using GridapMHD
using GridapMHD: analytical_shercliff_u
using GridapMHD: analytical_shercliff_j
using GridapMHD: transient_driver_inductionless_MHD
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

# Mesh load
model = GmshDiscreteModel("driver/duct.msh")
partition = (12,12,2)
domain = (-1,1,-1,1,0,0.1)
#model0 = CartesianDiscreteModel(domain,partition)
#model = simplexify(model0)

# Boundary conditions location

function hartmann_walls(coord)
  y = coord[2]
  tol = 1e-6
  if (1.0-tol < y ) | ( y < -1.0+tol )
    return true
  end
  false
end

function side_walls(coord)
  x = coord[1]
  tol = 1e-6

  if (1.0-tol < x ) | ( x < -1.0+tol )
    return true
  end
  false
end

function inlet(coord)
  z = coord[3]
  tol = 1e-6
  if (z < tol )
    return true
  end
  false
end

function outlet(coord)
  z = coord[3]
  tol = 1e-6

  if (0.1 - tol < z )
    return true
  end
  false
end

# Name inlet, outlet, and walls element sets.
add_entity!(model,hartmann_walls,"hartmann_walls")
add_entity!(model,side_walls,"side_walls")
add_entity!(model,inlet,"inlet")
add_entity!(model,outlet,"outlet")

side_wall_semilength = 1.0
hartmann_wall_semilength = 1.0
number_fourier_sumands = 10
wall_thickness = 1.0
wall_conductivity = 0.0
K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))
∂p∂z = -Re * K / L^3
u0(x) = analytical_shercliff_u(side_wall_semilength, hartmann_wall_semilength,
  wall_thickness, wall_conductivity, σ, ν*ρ, ∂p∂z, Ha, number_fourier_sumands, x)/U0

j0(x) = analytical_shercliff_j(side_wall_semilength, hartmann_wall_semilength,
  wall_thickness, wall_conductivity, σ, ν*ρ, ∂p∂z, Ha, number_fourier_sumands, x)/(σ * B0 * U0)

# Boundary conditions values
function g_u(coord,t)
  if (abs(coord[3]) < 1e-6) | (abs(coord[3]-0.1) < 1e-6)
    return u0(coord)
  else
    return VectorValue(0.0, 0.0, 0.0)
  end
end
g_u(t::Real) = x -> g_u(x,t)
g_u(x::VectorValue) =  g_u(x,0.0)

function g_j(x,t)
  #return j0(x)
  return VectorValue(0.0,0.0,0.0)
end
g_j(t::Real) = x -> g_j(x,t)
g_j(x::VectorValue) = g_j(x,0.0)

# Body force
function f_u(x)
#  VectorValue(0.0,0.0, -∂p∂z) * L/U0^2
  VectorValue(0.0,0.0, 0.0)
end

# Initial conditions
#function u0(coord)
#  x = coord[1]; y = coord[2]
#  return VectorValue(0.0, 0.0, (x-1)*(x+1) * (y-1)*(y+1) )
#end

xh_t, trian, dΩ = driver_inductionless_MHD(;
  #t0=0.0,
  #tF=500.0,
  #Δt=500.0,
  #θ=1.0,
  Re=Re,
  Ha=Ha,
  model=model,
  fluid_dirichlet_tags = ["inlet","side_walls","hartmann_walls"],
  fluid_neumann_tags = ["outlet"],
  magnetic_dirichlet_tags = ["inlet","side_walls","hartmann_walls","outlet"],
  magnetic_neumann_tags = [],
  c_w = 0.0,
  α = 10.0,
  fluid_dirichlet_conditions = g_u ,
  magnetic_dirichlet_conditions = g_j,
  fluid_body_force = f_u,
  #u0 = u0,
  #j0 = j0,
  constraint_presures = (false,false),
  usegmres = false,
  resultsfile = nothing
)

uh, ph, jh, φh = xh_t
divj = ∇⋅jh
writevtk(trian, "TransientDuct", order=2,
  cellfields=["uh"=>uh, "ph"=>ph, "jh"=>jh, "φh"=>φh, "divj"=>divj])


# xh, tf = xh_t[end]
# uh, ph, jh, φh = xh
# divj = (∇⋅jh)
#
# # Scale unknowns
# uh = uh * U0
# jh = jh * σ * B0 * U0
#
# # Analyical solutions
# side_wall_semilength = 1.0
# hartmann_wall_semilength = 1.0
# number_fourier_sumands = 10
# wall_thickness = 1.0
# wall_conductivity = 0.0
# K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))
# ∂p∂z = -Re * K / L^3
