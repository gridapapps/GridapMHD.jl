module DriverInductionlessMHDTest

using Gridap
using Test
using LineSearches: Static

using GridapMHD
using GridapMHD: driver_inductionless_MHD
using GridapMHD: transient_driver_inductionless_MHD
using GridapMHD: analytical_shercliff_u
using GridapMHD: analytical_shercliff_j
using GridapMHD: compute_u_j_errors

# Problem setting
ρ = 1.0
ν = 1.0
σ = 1.0

U0 = 10.0
B0 = 10.0
L = 1.0
Re = U0 * L / ν
Ha(x) = B0 * L * sqrt(σ/(ρ*ν))


K = Ha(0) / (1-0.825*Ha(0)^(-1/2)-Ha(0)^(-1))
∂p∂z = -Re * K / L^3

f_u(x) = VectorValue(0.0,0.0, -∂p∂z) * L/U0^2
g_u(x) = VectorValue(0.0,0.0,0.0)
g_j(x) = VectorValue(0.0,0.0,0.0)

function map1(coord)
  ncoord = streching(coord,  domain=(-1.0,0.0,-1.0,0.0,-1.0,0.0),factors=(3.0,3.0,3.0),dirs=(1,2,3))
  ncoord = streching(ncoord, domain=( 0.0,1.0, 0.0,1.0, 0.0,1.0),factors=(1.25,1.25,1.25),dirs=(1,2,3))
  ncoord
end

# Background mesh definition
partition = (5,5,3);
domain = (-1.0,1.0,-1.0,1.0,0.0,0.1)
model=CartesianDiscreteModel(domain,partition,map=map1,isperiodic=(false,false,true));


# Filter to select all walls
function hartmann_walls(coord)
  y = coord[2]
  tol = 1e-6
  if (1.0-tol < y ) | ( y < -1.0+tol )
    return true
  end
  false
end

# Filter to select inflow boundary
function side_walls(coord)
  x = coord[1]
  tol = 1e-6

  if (1.0-tol < x ) | ( x < -1.0+tol )
    return true
  end
  false
end

# Name inlet, outlet, and walls element sets.
add_entity!(model,hartmann_walls,"hartmann_walls")
add_entity!(model,side_walls,"side_walls")

# Call main MHD driver
xh, trian, dΩ = driver_inductionless_MHD(;
  Re=Re,
  Ha=Ha,
  model=model,
  fluid_dirichlet_tags = ["hartmann_walls","side_walls"],
  fluid_neumann_tags = [],
  magnetic_dirichlet_tags = ["hartmann_walls","side_walls"],
  magnetic_neumann_tags = [],
  fluid_dirichlet_conditions = g_u,
  magnetic_dirichlet_conditions = g_j,
  fluid_body_force = f_u,
  linesearch=Static(),
  constraint_presures = (true,true),
  precond_tau = 1e-10
)

uh, ph, jh, φh = xh

# Scale unknowns
uh = uh * U0
jh = jh * σ * B0 * U0

# Analyical solutions
side_wall_semilength = 1.0
hartmann_wall_semilength = 1.0
number_fourier_sumands = 10
wall_thickness = 1.0
wall_conductivity = 0.0
u0(x) = analytical_shercliff_u(side_wall_semilength, hartmann_wall_semilength,
  wall_thickness, wall_conductivity, σ, ν*ρ, ∂p∂z, Ha(0), number_fourier_sumands, x)

j0(x) = analytical_shercliff_j(side_wall_semilength, hartmann_wall_semilength,
  wall_thickness, wall_conductivity, σ, ν*ρ, ∂p∂z, Ha(0), number_fourier_sumands, x)

eu_l2, ej_l2 = compute_u_j_errors(uh, jh, u0, j0, dΩ)

@test eu_l2 < 0.3
@test ej_l2 < 5

Ha_ = Ha(0) # test for constant Hartmann

# Test transient driver
g_u(t::Real) = x -> g_u(x)
g_u(x,t) =  g_u(x)
g_j(t::Real) = x -> g_j(x)
g_j(x,t) = g_j(x)
xh_t, trian, dΩ = transient_driver_inductionless_MHD(;
  t0=0.0,
  tF=0.25,
  Δt=0.25,
  θ=0.5,
  Re=Re,
  Ha=Ha_,
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
  verbosity = Verbosity(3)
)

end #module
