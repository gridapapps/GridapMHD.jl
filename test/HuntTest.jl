module HuntTest

using Gridap
using Test

using GridapMHD: hunt
using GridapMHD: analytical_hunt_u
using GridapMHD: analytical_hunt_j
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

xh, trian, quad = hunt(nx=5,ny=5,Re=10.0,Ha=10.0)

# Postprocess
uh, ph, jh, φh = xh
divj = (∇⋅jh)

uh = uh * U0
jh = jh * σ * B0 * U0

# Analyical solutions
side_wall_semilength = 1.0
hartmann_wall_semilength = 1.0
number_fourier_sumands = 10
K = Ha / (1-0.95598*Ha^(-1/2)-Ha^(-1))
∂p∂z = -Re * K / L^3
u0(x) = analytical_hunt_u(side_wall_semilength, hartmann_wall_semilength,
  ρ*ν, ∂p∂z, Ha, number_fourier_sumands, x)

j0(x) = analytical_hunt_j(side_wall_semilength, hartmann_wall_semilength,
  σ, ρ*ν, ∂p∂z, Ha, number_fourier_sumands, x)

eu_l2, ej_l2 = compute_u_j_errors(uh, jh, u0, j0, trian, quad)
e_divj = sqrt(sum(integrate(divj*divj,trian,quad)))

@test eu_l2 < 0.06
@test ej_l2 < 0.8
@test e_divj < 2e-10
end #module
