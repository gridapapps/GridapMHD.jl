
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
partition = (4,4,3)
domain = (-1,1,-1,1,0,0.1)
model = CartesianDiscreteModel(domain,partition,isperiodic=(false,false,true))

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

K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))
∂p∂z = -Re * K / L^3

# Boundary conditions values
g_u(coord,t) = VectorValue(0.0, 0.0, 0.0)
g_u(t::Real) = x -> g_u(x,t)
g_u(x::VectorValue) =  g_u(x,0.0)

g_j(x,t) = VectorValue(0.0,0.0,0.0)
g_j(t::Real) = x -> g_j(x,t)
g_j(x::VectorValue) = g_j(x,0.0)

# Body force
f_u(x) = VectorValue(0.0,0.0, -∂p∂z) * L/U0^2

xh_t, trian, dΩ = transient_driver_inductionless_MHD(;
  t0=0.0,
  tF=1.0,
  Δt=0.25,
  θ=0.5,
  Re=Re,
  Ha=Ha,
  model=model,
  fluid_dirichlet_tags = ["side_walls","hartmann_walls"],
  fluid_neumann_tags = [],
  magnetic_dirichlet_tags = ["side_walls"],
  magnetic_non_perfectly_conducting_walls_tag = ["hartmann_walls",],
  magnetic_neumann_tags = [],
  c_w = 0.0,
  α = 10.0,
  fluid_dirichlet_conditions = g_u ,
  magnetic_dirichlet_conditions = g_j,
  fluid_body_force = f_u,
  constraint_presures = (false,false),
  usegmres = true,
  verbosity = Verbosity(1),
  resultsfile = "output_test"
)
