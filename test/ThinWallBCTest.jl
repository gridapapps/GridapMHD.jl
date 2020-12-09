module ThinWallBCTest

using Gridap
using Test
using LineSearches: BackTracking

using GridapMHD: conductive_thin_wall
using GridapMHD: analytical_shercliff_u
using GridapMHD: analytical_shercliff_j
using GridapMHD: compute_u_j_errors

# Problem setting
ρ = 1.0
ν = 1.0
σ = 1.0
wall_thickness = 0.10
wall_conductivity = 1.0
c_w = wall_conductivity * wall_thickness/σ

U0 = 10.0
B0 = 10.0
L = 1.0
Re = U0 * L / ν
Ha = B0 * L * sqrt(σ/(ρ*ν))

xh, trian, dΩ = conductive_thin_wall(nx=5,ny=5,c_w=c_w,Re=10.0,Ha=10.0,α=10.0)

uh, ph, jh, φh = xh
divj = (∇⋅jh)

# Scale unknowns
uh = uh * U0
jh = jh * σ * B0 * U0


# Analyical solutions
side_wall_semilength = 1.0
hartmann_wall_semilength = 1.0
number_fourier_sumands = 20
K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))
∂p∂z = -Re * K / L^3
u0(x) = analytical_shercliff_u(side_wall_semilength, hartmann_wall_semilength,
  wall_thickness, wall_conductivity, σ, ν*ρ, ∂p∂z, Ha, number_fourier_sumands, x)

j0(x) = analytical_shercliff_j(side_wall_semilength, hartmann_wall_semilength,
  wall_thickness, wall_conductivity, σ, ν*ρ, ∂p∂z, Ha, number_fourier_sumands, x)

writevtk(trian, "results.vtu",
  cellfields=["uh"=>uh, "ph"=>ph, "jh"=>jh, "φh"=>φh, "u"=>u0, "j"=>j0])

eu_l2, ej_l2 = compute_u_j_errors(uh, jh, u0, j0, dΩ)
e_divj = sqrt(sum(∫(divj*divj)*dΩ))

@test eu_l2 < 0.3
@test ej_l2 < 2
@test e_divj < 1e-12

end #module
