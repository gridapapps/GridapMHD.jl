module ShercliffTest

using Gridap
using Test

using GridapMHD: shercliff
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
Ha = B0 * L * sqrt(σ/(ρ*ν))

xh, trian, quad = shercliff(nx=5,ny=5,Re=10.0,Ha=10.0)

uh, ph, jh, φh = xh
divj = (∇*jh)

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
∂p∂z = -L^3 * K / Re
u0(x) = analytical_shercliff_u(side_wall_semilength, hartmann_wall_semilength,
  wall_thickness, wall_conductivity, σ, ν*ρ, ∂p∂z, Ha, number_fourier_sumands, x)

j0(x) = analytical_shercliff_j(side_wall_semilength, hartmann_wall_semilength,
  wall_thickness, wall_conductivity, σ, ν*ρ, ∂p∂z, Ha, number_fourier_sumands, x)

eu_l2, ej_l2 = compute_u_j_errors(uh, jh, u0, j0, trian, quad)
e_divj = sqrt(sum(integrate(divj*divj,trian,quad)))

@test eu_l2 < 0.003
@test ej_l2 < 0.05
@test e_divj < 1e-12

end #module
