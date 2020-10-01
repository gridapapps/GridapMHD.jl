
using Gridap

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

xh, trian, quad = driver_inductionless_MHD(;
  Re=10.0,
  Ha=10.0,
  meshfile="Messadek_coarse.msh",
  fluid_dirichlet_tags = ["inlet","wall_external","wall_internal"],
  fluid_neumann_tags = ["outlet_top","outlet_mid","outlet_bot"],
  magnetic_dirichlet_tags = ["wall_external","wall_internal"],
  magnetic_neumann_tags = ["inlet","outlet_top","outlet_mid","outlet_bot"],
  fluid_dirichlet_conditions = g_u ,
  magnetic_dirichlet_conditions = g_j,
  fluid_body_force = f_u,
  constraint_presures = (false,false),
  resultsfile = "Messadek"
)

uh, ph, jh, φh = xh
divj = (∇⋅jh)

# Scale unknowns
uh = uh * U0
jh = jh * σ * B0 * U0

# Analyical solutions
side_wall_semilength = 1.0
hartmann_wall_semilength = 1.0
number_fourier_sumands = 10
wall_thickness = 1.0
wall_conductivity = 0.0
K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))
∂p∂z = -Re * K / L^3
